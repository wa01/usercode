#include "TF1.h"
#include "TH1F.h"
#include "TEfficiency.h"
#include "TRandom2.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"

#include "HtMetTreeReader.h"

#include "Minuit2/FCNBase.h"
#include "TFitterMinuit.h"

#include <vector>
#include <string>
using namespace std;

//
// response function: error function with shape parameter (skewness) +
//   rescaling to range [emin,emax]
//
double responseFunction (double met, double location, double scale, double shape, 
			 double emin, double emax) {
  double y = (met-location)/scale;
  if ( fabs(shape)>1.e-3 )  y = -1/shape*log(1-shape*y);
  double res = TMath::Freq(y);
  res = res*(emax-emin) + emin;
  return res;
}
//
// alternative interface to the response function
//
double responseFunction (double met, vector<double> pars) {
  return responseFunction(met,pars[0],pars[1],pars[2],pars[3],pars[4]);
}
//
// alternative interface to the response function
//
double responseFunction (double* xs, double* pars) {
  return responseFunction(*xs,pars[0],pars[1],pars[2],pars[3],pars[4]);
}

//
// -logL for the unbinned fit to the response function
//
class FcnLogL : public ROOT::Minuit2::FCNBase {
public:
  // constructor from event vector (HT,MET) and HT cut
  FcnLogL (const vector< pair<float,float> >& events, float htMin) : 
    events_(events), htMin_(htMin), errorDef_(0.5) {}
  // error definition
  double Up() const {return errorDef_;}
  // (re)set htMin
  void setHTmin (float htMin) {htMin_ = htMin;}
  // calculation of the LL
  double operator() (const vector<double>& pars) const
  {
    // initialization
    double result(0.);
    // parameters
    double location = pars[0];
    double scale = pars[1];
    double shape = pars[2];
    double emin = pars[3];
    double emax = pars[4];
    //
    // event loop
    //
    for ( size_t i=0; i<events_.size(); ++i ) {
      // MET value and efficiency at this point
      double met = events_[i].second;
      double eff = responseFunction(met,location,scale,shape,emin,emax);
      // avoid problems with log(eff)
      if ( eff>0.999999 )  eff = 0.999999;
      if ( eff<0.000001 )  eff = 0.000001;
      // probability "eff" to pass the HT cut at MET value "met"
      if ( events_[i].first>htMin_ )  result -= log(eff);
      else  result -= log(1.-eff);
    }
    
    return result;
  }
private:
  const vector< pair<float,float> >& events_;
  float htMin_;
  double errorDef_;
};


class LogLErfTest {
public:
  LogLErfTest (const string& path, float metMin = 150, float metMax = 1500) : 
    path_(path), metMin_(metMin), metMax_(metMax), minuitx_(0), fcn_(0) {
    if ( !path.empty() && path[path.size()-1]!='/' )   path_ += "/";
  }
  LogLErfTest (float metMin = 150, float metMax = 1500) : 
    metMin_(metMin), metMax_(metMax), minuitx_(0), fcn_(0) {}

  void fitSingleHT (TH1* hAll, TH1* hAllS1, float htMin);

  void fitMultiHT (const vector<string>& filenames, const vector<float>& fileweights,
		   vector<float> htMin = vector<float>());
  
  void fitMultiHT (const vector<string>& filenames, const vector<float>& fileweights, float htMin) {
    vector<float> htMins(1,htMin);
    fitMultiHT(filenames,fileweights,htMins);
  }
  
  void fitMultiHT (const string& filename, vector<float> htMin = vector<float>()) {
    vector<string> filenames(1,filename);
    vector<float> fileweights(1,1.);
    fitMultiHT (filenames,fileweights,htMin);
  }
  
  void fitMultiHT (const string& filename, float htMin) {
    vector<float> htMins(1,htMin);
    fitMultiHT(filename,htMins);
  }
  
  string fullFileName (const string& name) const {
    string fullName(name);
    if ( !path_.empty() ) fullName = path_ + fullName;
    return fullName;
  }

  string fileTTJets () const {
    return fullFileName("8TeV-TTJets.root");
  }
  vector<string> filesWJets () const {
    vector<string> filenames;
    filenames.push_back(fullFileName("8TeV-WJets-HT250to300.root"));
    filenames.push_back(fullFileName("8TeV-WJets-HT300to400.root"));
    filenames.push_back(fullFileName("8TeV-WJets-HT400.root"));
    return filenames;
  }
  vector<float> weightsWJets () const {
    vector<float> fileweights;
    fileweights.push_back(0.1496);
    fileweights.push_back(0.1127);
    fileweights.push_back(0.0765);
    return fileweights;
  }

  void fitTTJets (float htMin) {
    fitMultiHT(fileTTJets(),htMin);
  }
  
  void fitTTJets (vector<float> htMins = vector<float>()) {
    fitMultiHT(fileTTJets(),htMins);
  }
  
  void fitWJets (float htMin) {
    fitMultiHT(filesWJets(),weightsWJets(),htMin);
  }
  
  void fitWJets (vector<float> htMins = vector<float>()) {
    fitMultiHT(filesWJets(),weightsWJets(),htMins);
  }

  void drawGraphs () const;

  void writeGraphs (string name="logLErfTest_Graphs.root") const;

private:
  string path_;
  float metMin_;
  float metMax_;

  TFitterMinuit* minuitx_;
  FcnLogL* fcn_;

  typedef vector< pair<float,float> >  DataContainer;
  DataContainer data_;

  vector<TGraphErrors*> parGraphs_;
};
//
// Fit
//
void 
LogLErfTest::fitSingleHT (TH1* hAll, TH1* hAllS1, float htMin) {

  // TRandom2 rgen;
  // for ( size_t i=0; i<data_.size(); ++i ) {
  //   double met = rgen.Exp(250.);
  //   double eff = responseFunction(met,230.,130.,0.02);
  //   double ht = rgen.Uniform()<eff ? 1.1*htMin : 0.9*htMin;
  //   data_[i].first = ht;
  //   data_[i].second = met;
  // }

  //
  // histograms
  //
  // int nhBins = int((metMax_-metMin_)/50.+0.5);
  int ihtMin(htMin+0.5);
  char hName[128], hTitle[256];
  sprintf(hName,"hCut%4.4d",ihtMin); 
  sprintf(hTitle,"MET (HT > %d)",ihtMin);
  TH1* hCut = (TH1*)hAll->Clone(hName);
  hCut->SetTitle(hTitle);
  hCut->Reset();
  for ( size_t i=0; i<data_.size(); ++i ) {
    if ( data_[i].first>htMin ) hCut->Fill(data_[i].second);
  }
  // hAllS1->Divide(hAll); // mean MET / MET bin

  char name[128];
  sprintf(name,"eff%4.4d",ihtMin);
  TEfficiency* eff = new TEfficiency(*hCut,*hAll);
  eff->SetNameTitle(name,name);
  eff->SetMarkerStyle(24);
  eff->Draw();
  //
  // LL function
  //
  minuitx_->SetParameter(0,"loc",htMin,100.,-metMax_,metMax_);
  minuitx_->SetParameter(1,"scale",100.,20.,0.,500.);
  minuitx_->SetParameter(2,"shape",0.,0.05,-0.2,0.2);
  minuitx_->SetParameter(3,"min",0.,0.05,0,1);
  minuitx_->SetParameter(4,"max",1.,0.05,0,1);
  //
  // first fit: location and scale
  //
  minuitx_->FixParameter(2);
  minuitx_->FixParameter(3);
  minuitx_->FixParameter(4);
  // minuitx_->SetMaxIterations(10000);
  // minuitx_->SetPrintLevel(2);
  // minuitx_->CreateMinimizer();
  int ierr = minuitx_->Minimize();
  // if ( ierr==0 ) {
  //   cout << "Results = " << minuitx_->GetParameter(0)
  // 	 << " " << minuitx_->GetParameter(1) << endl;    
  // }
  //
  // second fit: skewness
  //
  minuitx_->FixParameter(0);
  minuitx_->FixParameter(1);
  minuitx_->ReleaseParameter(2);
  double shapeMax = minuitx_->GetParameter(1)/max(metMax_-minuitx_->GetParameter(0),minuitx_->GetParameter(0)-metMin_);
  minuitx_->SetParameter(2,"shape",minuitx_->GetParameter(2),0.05,-0.2,shapeMax);
  ierr = minuitx_->Minimize();
  // if ( ierr==0 ) {
  //   cout << "Results = " << minuitx_->GetParameter(0)
  // 	 << " " << minuitx_->GetParameter(1) << endl;    
  // }
  //
  // third fit: location, scale, shape
  minuitx_->ReleaseParameter(0);
  minuitx_->ReleaseParameter(1);
  ierr = minuitx_->Minimize();
  // if ( ierr==0 ) {
  //   cout << "Results = " << minuitx_->GetParameter(0)
  // 	 << " " << minuitx_->GetParameter(1) << endl;
  // }
  //
  // new Graph with updated mean METs / bin
  //
  gPad->Update();
  TGraphAsymmErrors* graph = eff->GetPaintedGraph();
  TGraphAsymmErrors* graph1 = (TGraphAsymmErrors*)graph->Clone();
  TAxis* axis = hAll->GetXaxis();
  int nGraph = graph1->GetN();
  double* xsGraph = graph1->GetX();
  double* ysGraph = graph1->GetY();
  for ( int i=0; i<nGraph; ++i ) {
    int ibin = axis->FindBin(xsGraph[i]);
    double xl = axis->GetBinLowEdge(ibin);
    double xh = axis->GetBinUpEdge(ibin);
    double xmean = hAllS1->GetBinContent(ibin);
    // cout << i << " " << xsGraph[i] << " " << xmean << endl;
    graph1->SetPoint(i,xmean,ysGraph[i]);
    graph1->SetPointEXlow(i,xmean-xl);
    graph1->SetPointEXhigh(i,xh-xmean);
  }
  graph1->SetMarkerStyle(20); graph1->SetLineColor(2); graph1->SetMarkerColor(2);
  // graph1->Draw("same P");
  graph1->Draw("AP");
  //
  // function corresponding to the result of the unbinned LL fit
  //
  sprintf(name,"fLL%4.4d",ihtMin);
  TF1* f = new TF1(name,responseFunction,hAll->GetXaxis()->GetXmin(),
		   hAll->GetXaxis()->GetXmax(),5);
  f->SetParameters(minuitx_->GetParameter(0),minuitx_->GetParameter(1),
		   minuitx_->GetParameter(2),minuitx_->GetParameter(3),
		   minuitx_->GetParameter(4));
  f->SetLineColor(2);
  f->SetLineWidth(2);
  f->Draw("same");
  //
  // fit to the graph
  //
  sprintf(name,"fG%4.4d",ihtMin);
  TF1* g = new TF1(name,responseFunction,hAll->GetXaxis()->GetXmin(),
		   hAll->GetXaxis()->GetXmax(),5);
  g->SetParameters(minuitx_->GetParameter(0),minuitx_->GetParameter(1),
		   minuitx_->GetParameter(2),0.,1.);
  g->SetParLimits(0,-metMax_,metMax_);
  g->SetParLimits(1,0.,500.);
  shapeMax = minuitx_->GetParameter(1)/max(metMax_-minuitx_->GetParameter(0),minuitx_->GetParameter(0)-metMin_);
  g->SetParLimits(2,-0.2,shapeMax);
  g->FixParameter(3,0.);
  g->FixParameter(4,1.);
  g->SetLineColor(4);
  g->SetLineWidth(2);
  g->Print();
  // g->Draw("same");
  graph1->Fit(g);

  // TCanvas* c = new TCanvas();
  // c->Divide(1,2);
  // c->cd(1);
  // hAllHT->DrawCopy();
  // c->cd(2);
  // hAll->DrawCopy();
}


//
// Fit
//
void 
LogLErfTest::fitMultiHT (const vector<string>& filenames, const vector<float>& fileweights,
		 vector<float> htMin) {
  //
  // (HT,MET) pairs from tree
  //
  assert (!filenames.empty() && (filenames.size()==fileweights.size()));

  double maxFileWgt(0.);
  for ( size_t i=0; i<fileweights.size(); ++i ) {
    if ( fileweights[i]>maxFileWgt )  maxFileWgt = fileweights[i];
  }

  TRandom2 rgen;
  data_.clear();
  for ( size_t i=0; i<filenames.size(); ++i ) {
    DataContainer dataTmp;
    HtMetTreeReader reader(filenames[i].c_str());
    reader.loop(dataTmp);
    for ( size_t j=0; j<dataTmp.size(); ++j ) {
      if ( dataTmp[j].second<metMin_ )  continue;
      if ( rgen.Uniform()<fileweights[i]/maxFileWgt )  data_.push_back(dataTmp[j]);
    }
  }
  cout << "Number of events = " << data_.size() << endl;
  // //
  // // impose lower limit on MET
  // //
  // DataContainer data;
  // data.reserve(data1.size());
  // for ( size_t i=0; i<data1.size(); ++i ) {
  //   if ( data1[i].second>=metMin_ )  data.push_back(data1[i]);
  // }
  // cout << "Number of events after cut = " << data.size() << endl;

  if ( htMin.empty() ) {
    htMin.push_back(400.);
    htMin.push_back(450.);
    htMin.push_back(500.);
    htMin.push_back(550.);
    htMin.push_back(600.);
    htMin.push_back(650);
    htMin.push_back(700.);
    htMin.push_back(750.);
    htMin.push_back(800.);
    htMin.push_back(1000.);
    htMin.push_back(1200.);
    htMin.push_back(1500.);
  }

  TCanvas* cnv = new TCanvas("c","c",800,600);
  int npad = int(sqrt(htMin.size()/12.)+0.5);
  cnv->Divide(4*npad,3*npad);
  //
  // LL function
  //
  delete minuitx_;
  delete fcn_;
  fcn_ = new FcnLogL(data_,0.);
  minuitx_ = new TFitterMinuit();
  minuitx_->SetMinuitFCN(fcn_);
  minuitx_->SetParameter(0,"loc",0.,100.,0.,metMax_);
  minuitx_->SetParameter(1,"scale",100.,20.,0.,500.);
  minuitx_->SetParameter(2,"shape",0.,0.05,-1.0,1.0);
  minuitx_->SetParameter(3,"min",0.,0.05,0,1);
  minuitx_->SetParameter(4,"max",1.,0.05,0,1);
  minuitx_->SetMaxIterations(10000);
  minuitx_->SetPrintLevel(2);
  minuitx_->CreateMinimizer();

  TH1::SetDefaultSumw2(true);
  int nhBins = int((metMax_-metMin_)/50.+0.5);
  TH1* hAllHT = new TH1F("hAllHT","HT (inclusive)",nhBins,metMin_,metMax_);
  TH1* hAll = new TH1F("hAll","MET (inclusive)",nhBins,metMin_,metMax_);
  TH1* hAllS1 = new TH1F("hAllS1","MET (inclusive) s1",nhBins,metMin_,metMax_);
  for ( size_t i=0; i<data_.size(); ++i ) {
    hAllHT->Fill(data_[i].first);
    hAll->Fill(data_[i].second);
    hAllS1->Fill(data_[i].second,data_[i].second);
  }
  hAllS1->Divide(hAll); // mean MET / MET bin

  // char hName[128], hTitle[256];
  char name[128], title[256];
  parGraphs_.clear(); // should delete previous graphs, if any ...
  for ( size_t i=0; i<(size_t)minuitx_->GetNumberTotalParameters(); ++i ) {
    sprintf(name,"gPar%d",(int)i);
    sprintf(title,"par %d (%s)",(int)i,minuitx_->GetParName(i));
    TGraphErrors* graph = new TGraphErrors();
    graph->SetNameTitle(name,title);
    graph->SetMarkerStyle(20);
    parGraphs_.push_back(graph);
  }
  for ( size_t iht=0; iht<htMin.size(); ++iht ) {
    cnv->cd(iht+1);
    fcn_->setHTmin(htMin[iht]);
    fitSingleHT (hAll,hAllS1,htMin[iht]);
    for ( size_t i=0; i<parGraphs_.size(); ++i ) {
      TGraphErrors* graph = parGraphs_[i];
      graph->SetPoint(iht,htMin[iht],minuitx_->GetParameter(i));
      graph->SetPointError(iht,0.,minuitx_->GetParError(i));
    }
  }

  drawGraphs();
  
}

void
LogLErfTest::drawGraphs () const
{
  TCanvas* cpar = new TCanvas("cpar","cpar",800,300);
  cpar->Divide(3,1);
  for ( size_t i=0; i<min(parGraphs_.size(),size_t(3)); ++i ) {
    cpar->cd(i+1);
    parGraphs_[i]->Draw("AP");
    cout << parGraphs_[i]->GetN() << endl;
  }
  cpar->cd();
}

void
LogLErfTest::writeGraphs (string name) const
{
  TFile f(name.c_str(),"recreate");
  for ( size_t i=0; i<parGraphs_.size(); ++i ) {
    TGraphErrors* clone = (TGraphErrors*)parGraphs_[i]->Clone();
    clone->Write();
  }
  f.Write();
  f.Close();
}

