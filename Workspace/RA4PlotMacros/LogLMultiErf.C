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
  if ( res<0.000000001 )  res = 0.000000001;
  if ( res>0.999999999 )  res = 0.999999999;
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
// double-sided response function
//
double doubleSidedResponseFunction (double* xs, double* pars) {
  double effLow = responseFunction(xs,&pars[0]);
  double effHigh = responseFunction(xs,&pars[5]);
  return effLow*(1-effHigh);
}

//
// HT dependence of parameters
//
double htPolynomial (float ht, vector<double>::const_iterator iBegin, vector<double>::const_iterator iEnd) {
  double result(0.);
  double factor(1.);
  for ( vector<double>::const_iterator i=iBegin; i<iEnd; ++i ) {
    result += (*i)*factor;
    factor *= ht;
  }
  return result;
}

double htPolynomial (float ht, const vector<double>& pars) {
  return htPolynomial(ht,pars.begin(),pars.end());
}

double htPolynomialError (float ht, const vector< vector<double> >& covs) {
  double result(0.);
  double factorI(1.);
  for ( size_t i=0; i<covs.size(); ++i ) {
    assert (covs[i].size()==covs.size());
    double sumJ(0.);
    double factorJ(1.);
    for ( size_t j=0; j<covs.size(); ++j ) {
      sumJ += covs[i][j]*factorJ;
      factorJ *= ht;
    }
    result += sumJ*factorI;
    factorI *= ht;
  }
  return sqrt(result);
}

struct FixMask {
  FixMask () : mask(0) {}
  FixMask (unsigned int i) : mask(i) {}
  FixMask& operator|= (unsigned int i) {
    mask |= i;
    return *this;
  }
  // FixMask& operator= (unsigned int i) {
  //   mask = i;
  //   cout << "mask = " << mask << endl;
  //   return *this;
  // }
  bool isFixed (unsigned int i) const {
    return (mask>>i)&1;
  }
  unsigned int mask;
};


struct DataTuple {
  DataTuple () : htBin(-1), ht(0.), met(0.) {}
  DataTuple (int _htBin, float _ht, float _met) : 
    htBin(_htBin), ht(_ht), met(_met) {}
  DataTuple (float _ht, float _met, const vector<float>& htMins) {
    set(_ht,_met,htMins);
  }
  int set (float _ht, float _met, const vector<float>& htMins) {
    ht = _ht;
    met = _met;
    for ( int i=0; i<(int)htMins.size(); ++i ) {
      if ( ht<htMins[i] ) 
	return htBin=(i-1);
    }
    return htBin=(int)htMins.size()-1;
  }
  int htBin;
  float ht;
  float met;
};
typedef vector<DataTuple> DataContainer;
//
// -logL for the unbinned fit to the response function
//
class FcnMultiLogL : public ROOT::Minuit2::FCNBase {
public:
  // constructor from event vector (HT,MET) and HT cut
  FcnMultiLogL (const DataContainer& events, const vector<float>& htMins, FixMask mask) : 
    events_(events), htMins_(htMins), mask_(mask), errorDef_(0.5) {}
  // error definition
  double Up() const {return errorDef_;}
  // (re)set htMin
  void setHTmins (const vector<float>& htMins) {htMins_ = htMins;}
  unsigned int htPolynomialOrder (unsigned int index) const {
    static unsigned int orders[5] = { 1, 1, 0, 0, 0 };
    return orders[index];
  }
  vector<double>::const_iterator retrievePars (unsigned int index, 
					       vector<double>::const_iterator& iBegin,
					       vector<double>::const_iterator iEnd) const {
    vector<double>::const_iterator result(iBegin);
    
    unsigned int length;
    if ( mask_.isFixed(index) ) {
      length = htPolynomialOrder(index) + 1;
    }
    else {
      length = htMins_.size();
    }
    assert ((iEnd-iBegin)>=length);
    iBegin += length;
    return result;
  }

  double getErfPar (unsigned int index, size_t htBin, bool lowerEdge,
		    vector<double>::const_iterator iBegin, 
		    vector<double>::const_iterator iEnd) const {
    if ( mask_.isFixed(index) ) {
      double htRef = lowerEdge ? htMins_[htBin] : htMins_[htBin+1];
      // cout << "index " << index << " " << htRef << " " << *iBegin << " " << *(iEnd-1) << " " << htPolynomial(htRef,iBegin,iEnd) << endl;
      return htPolynomial(htRef,iBegin,iEnd);
    }
    else {
      size_t offset = htBin;
      if ( !lowerEdge )  ++offset;
      return *(iBegin+offset);
    }
  }

  vector< vector<double>::const_iterator > getParIterators (const vector<double>& pars) const {
    vector< vector<double>::const_iterator > parIterators;
    vector<double>::const_iterator ip(pars.begin());
    for ( size_t i=0; i<5; ++i ) {
      parIterators.push_back(retrievePars(i,ip,pars.end()));
    }
    parIterators.push_back(pars.end());
    return parIterators;
  }

  vector< vector<double> >  getResFunPars (const vector< vector<double>::const_iterator >& parIterators,
					   bool lowerEdge) const {
    vector< vector<double> >  result;
    vector<double> resFunParsByBin(5,0);
    for ( size_t i=0; i<htMins_.size(); ++i ) {
      if ( lowerEdge ) {
	for ( size_t j=0; j<5; ++j ) {
	  resFunParsByBin[j] = getErfPar(j,i,true,parIterators[j],parIterators[j+1]);
	}
	result.push_back(resFunParsByBin);
      }
      else if ( (i+1)<htMins_.size() ) {
	for ( size_t j=0; j<5; ++j ) {
	  resFunParsByBin[j] = getErfPar(j,i,false,parIterators[j],parIterators[j+1]);
	}
	result.push_back(resFunParsByBin);
      }
    }
    return result;
  }
  // calculation of the LL
  double operator() (const vector<double>& pars) const
  {
    // initialization
    double result(0.);
    // parameters
    vector< vector<double>::const_iterator > parIterators =
      getParIterators(pars);
    //
    // response function parameters / bin
    //
    vector< vector<double> > resFunParsLow = getResFunPars(parIterators,true);
    vector< vector<double> > resFunParsHigh = getResFunPars(parIterators,false);
    // cout << "iterators " << parIterators[0]-pars.begin() << " ";
    // cout << parIterators[1]-pars.begin() << " ";
    // cout << parIterators[2]-pars.begin() << " ";
    // cout << parIterators[3]-pars.begin() << " ";
    // cout << parIterators[4]-pars.begin() << " ";
    // cout << parIterators[5]-pars.begin() << endl;
    // cout << "Parameters";
    // for ( size_t i=0; i<5; ++i ) {
    //   for ( vector<double>::const_iterator ipp=parIterators[i]; ipp<parIterators[i+1]; ++ipp )
    // 	cout << " " << *ipp;
    //   cout << " ; ";
    // }
    // cout << endl;
    //
    // event loop
    //
    for ( size_t i=0; i<events_.size(); ++i ) {
      // MET value and efficiency at this point
      int ht = events_[i].ht;
      double met = events_[i].met;

      for ( size_t j=0; j<htMins_.size(); ++j ) {
	double eff = responseFunction(met,resFunParsLow[j]);
	// if ( (j+1)<htMins_.size() ) {
	//   double effUp = responseFunction(met,resFunParsHigh[j]);
	//   eff *= (1.-effUp);
	// }
	if ( ht>=htMins_[j] )  result -= log(eff);
	else  result -= log(1.-eff);
      }

    }
    return result;
  }
private:
  const DataContainer& events_;
  vector<float> htMins_;
  FixMask mask_;
  double errorDef_;
};


class LogLMultiErf {
public:
  LogLMultiErf (const string& path, vector<float> htMins = vector<float>(), 
		float metMin = 150, float metMax = 1500) : 
    path_(path), metMin_(metMin), metMax_(metMax), minuitx_(0), fcn_(0) {
    if ( !path.empty() && path[path.size()-1]!='/' )   path_ += "/";
    setHtMins(htMins);
  }
  LogLMultiErf (vector<float> htMins = vector<float>(), float metMin = 150, float metMax = 1500) : 
    metMin_(metMin), metMax_(metMax), minuitx_(0), fcn_(0) {
    setHtMins(htMins);
  }

  // void fitSingleHT (TH1* hAll, TH1* hAllS1, float htMin);

  void fitMultiHT (const vector<string>& filenames, const vector<float>& fileweights, FixMask mask = 28);
  
  void fitMultiHT (const string& filename, FixMask mask = 28) {
    vector<string> filenames(1,filename);
    vector<float> fileweights(1,1.);
    fitMultiHT (filenames,fileweights,mask);
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

  void fitTTJets (FixMask mask = 28) {
    fitMultiHT(fileTTJets(),mask);
  }
  
  void fitWJets (FixMask mask = 28) {
    fitMultiHT(filesWJets(),weightsWJets(),mask);
  }
  
  void drawGraphs () const;

  void writeGraphs (string name="logLErfTest_Graphs.root") const;

private:
  void setHtMins (const vector<float>& htMins) {
    htMins_ = htMins;
    if ( htMins.empty() ) {
      htMins_.push_back(400.);
      htMins_.push_back(450.);
      htMins_.push_back(500.);
      htMins_.push_back(550.);
      htMins_.push_back(600.);
      htMins_.push_back(650);
      htMins_.push_back(700.);
      htMins_.push_back(750.);
      htMins_.push_back(800.);
      htMins_.push_back(1000.);
      htMins_.push_back(1200.);
      htMins_.push_back(1500.);
    }
  }

private:
  string path_;
  vector<float> htMins_;
  float metMin_;
  float metMax_;

  TFitterMinuit* minuitx_;
  FcnMultiLogL* fcn_;

  DataContainer data_;

  vector<TGraphErrors*> parGraphs_;
};

//
// Fit
//
void 
LogLMultiErf::fitMultiHT (const vector<string>& filenames, const vector<float>& fileweights, FixMask mask) {

  mask |= 24;

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
  DataTuple tuple;
  for ( size_t i=0; i<filenames.size(); ++i ) {
    vector< pair<float,float> > dataTmp;
    HtMetTreeReader reader(filenames[i].c_str());
    reader.loop(dataTmp);
    for ( size_t j=0; j<dataTmp.size(); ++j ) {
      if ( dataTmp[j].second<metMin_ )  continue;
      if ( rgen.Uniform()<fileweights[i]/maxFileWgt ) {
	tuple.set(dataTmp[j].first,dataTmp[j].second,htMins_);
	data_.push_back(tuple);
      }
    }
  }
  cout << "Number of events = " << data_.size() << endl;

  TCanvas* cnv = new TCanvas("c","c",800,600);
  int npad = int(sqrt(htMins_.size()/12.)+0.5);
  cnv->Divide(4*npad,3*npad);
  //
  // LL function
  //
  delete minuitx_;
  delete fcn_;
  fcn_ = new FcnMultiLogL(data_,htMins_,mask);
  minuitx_ = new TFitterMinuit();
  minuitx_->SetMinuitFCN(fcn_);
  //
  char parName[128];
  unsigned int index(0);
  vector<unsigned int> parIndices;
  // location
  parIndices.push_back(index);
  if ( mask.isFixed(0) ) {
    minuitx_->SetParameter(index++,"locP0",0.,100.,-metMax_,metMax_);
    minuitx_->SetParameter(index++,"locP1",0.5,0.1,-3.,3.);
  }
  else {
    for ( size_t i=0; i<htMins_.size(); ++i ) {
      sprintf(parName,"locHT%4.4d",int(htMins_[i]+0.5));
      minuitx_->SetParameter(index++,parName,htMins_[i],50.,0.,metMax_);
    }
  }
  // scale
  parIndices.push_back(index);
  if ( mask.isFixed(1) ) {
    minuitx_->SetParameter(index++,"scaleP0",150.,10.,0.,500.);
    minuitx_->SetParameter(index++,"scaleP1",0.,0.1,-1.,1.);
  }
  else {
    for ( size_t i=0; i<htMins_.size(); ++i ) {
      sprintf(parName,"scaleHT%4.4d",int(htMins_[i]+0.5));
      minuitx_->SetParameter(index++,parName,100,20.,0.,500.);
    }
  }
  // shape
  parIndices.push_back(index);
  if ( mask.isFixed(2) ) {
    minuitx_->SetParameter(index++,"shapeP0",0.,0.05,-1.0,1.0);
  }
  else {
    for ( size_t i=0; i<htMins_.size(); ++i ) {
      sprintf(parName,"shapeHT%4.4d",int(htMins_[i]+0.5));
      minuitx_->SetParameter(index++,parName,0.,0.05,-1.0,1.0);
    }
  }
  // emin
  parIndices.push_back(index);
  if ( mask.isFixed(3) ) {
    minuitx_->SetParameter(index++,"eminP0",0.,0.05,0.,1.);
  }
  else {
    for ( size_t i=0; i<htMins_.size(); ++i ) {
      sprintf(parName,"eminHT%4.4d",int(htMins_[i]+0.5));
      minuitx_->SetParameter(index++,parName,0.,0.05,0.,1.);
    }
  }
  // emax
  parIndices.push_back(index);
  if ( mask.isFixed(4) ) {
    minuitx_->SetParameter(index++,"emaxP0",1.,0.05,0.,1.);
  }
  else {
    for ( size_t i=0; i<htMins_.size(); ++i ) {
      sprintf(parName,"emaxHT%4.4d",int(htMins_[i]+0.5));
      minuitx_->SetParameter(index++,parName,1.,0.05,0.,1.);
    }
  }
  //
  parIndices.push_back(index);
  
  minuitx_->SetMaxIterations(10000);
  minuitx_->SetPrintLevel(3);
  minuitx_->CreateMinimizer();

  TH1::SetDefaultSumw2(true);
  int nhBins = int((metMax_-metMin_)/50.+0.5);
  TH1* hAllHT = new TH1F("hAllHT","HT (inclusive)",nhBins,metMin_,metMax_);
  TH1* hAll = new TH1F("hAll","MET (inclusive)",nhBins,metMin_,metMax_);
  TH1* hAllS1 = new TH1F("hAllS1","MET (inclusive) s1",nhBins,metMin_,metMax_);
  vector<TH1*> hCuts;
  char name[128], title[256];
  for ( size_t i=0; i<htMins_.size(); ++i ) {
    sprintf(name,"hCutHT%4.4d",int(htMins_[i]+0.5));
    sprintf(title,"MET (HT > %d)",int(htMins_[i]+0.5));
    hCuts.push_back(new TH1F(name,title,nhBins,metMin_,metMax_));
  }

  for ( size_t i=0; i<data_.size(); ++i ) {
    hAllHT->Fill(data_[i].ht);
    hAll->Fill(data_[i].met);
    hAllS1->Fill(data_[i].met,data_[i].met);
    for ( size_t j=0; j<htMins_.size(); ++j ) {
      if ( data_[i].ht>=htMins_[j] )  hCuts[j]->Fill(data_[i].met);
    }
  }
  hAllS1->Divide(hAll); // mean MET / MET bin

  //
  // first fit: location and scale
  //
  for ( size_t i=parIndices[1]; i<parIndices[5]; ++i )  minuitx_->FixParameter(i);
  int ierr = minuitx_->Minimize();
  for ( size_t i=parIndices[0]; i<parIndices[1]; ++i )  minuitx_->FixParameter(i);
  for ( size_t i=parIndices[1]; i<parIndices[2]; ++i )  minuitx_->ReleaseParameter(i);
  ierr = minuitx_->Minimize();
  for ( size_t i=parIndices[0]; i<parIndices[1]; ++i )  minuitx_->ReleaseParameter(i);
  ierr = minuitx_->Minimize();
  for ( size_t i=parIndices[0]; i<parIndices[2]; ++i )  minuitx_->FixParameter(i);
  for ( size_t i=parIndices[2]; i<parIndices[3]; ++i )  minuitx_->ReleaseParameter(i);
  ierr = minuitx_->Minimize();
  for ( size_t i=parIndices[0]; i<parIndices[2]; ++i )  minuitx_->ReleaseParameter(i);
  ierr = minuitx_->Minimize();

  // char hName[128], hTitle[256];
  parGraphs_.clear(); // should delete previous graphs, if any ...
  for ( size_t i=0; i<3; ++i ) {
    TGraphErrors* graph = new TGraphErrors();
    switch ( i ) {
    case 0: 
      graph->SetNameTitle("Location","Location");
      break;
    case 1: 
      graph->SetNameTitle("Scale","Scale");
      break;
    case 2: 
      graph->SetNameTitle("Shape","Shape");
      break;
    }
    graph->SetMarkerStyle(20);
    parGraphs_.push_back(graph);
  }

  vector<double> resultPars;
  for ( size_t i=0; i<5; ++i ) {
    for ( size_t j=parIndices[i]; j<parIndices[i+1]; ++j ) {
      resultPars.push_back(minuitx_->GetParameter(j));
      // cout << " " << resultPars.back();
    }
  }
  // cout << endl;
  vector< vector<double>::const_iterator > parIterators = fcn_->getParIterators(resultPars);
  //
  // response function parameters / bin
  //
  vector< vector<double> > resFunParsLow = fcn_->getResFunPars(parIterators,true);
  vector< vector<double> > resFunParsHigh = fcn_->getResFunPars(parIterators,false);

  vector<TEfficiency*> efficiencies;
  for ( size_t iht=0; iht<htMins_.size(); ++iht ) {
    cnv->cd(iht+1);
    sprintf(name,"eff%4.4d",int(htMins_[iht]+0.5));
    TEfficiency* eff = new TEfficiency(*hCuts[iht],*hAll);
    eff->SetNameTitle(name,name);
    eff->SetMarkerStyle(24);
    eff->Draw();
    //
    // function corresponding to the result of the unbinned LL fit
    //
    sprintf(name,"fLL%4.4d",int(htMins_[iht]+0.5));
    TF1* f;
    // cout << "Bin " << htMins_[iht] << endl;
    f = new TF1(name,responseFunction,hAll->GetXaxis()->GetXmin(),
		hAll->GetXaxis()->GetXmax(),5);
    for ( size_t k=0; k<5; ++k ) {
      f->SetParameter(k,resFunParsLow[iht][k]);
    }
    f->SetLineColor(2);
    f->SetLineWidth(2);
    f->Draw("same");
  }
  //   fcn_->setHTmin(htMin[iht]);
  //   fitSingleHT (hAll,hAllS1,htMin[iht]);
  for ( size_t i=0; i<parGraphs_.size(); ++i ) {
    TGraphErrors* graph = parGraphs_[i];
    if ( mask.isFixed(i) ) {
      vector<double> polyPars;
      vector< vector<double> > polyCovs;
      for ( size_t ip=parIndices[i]; ip<parIndices[i+1]; ++ip ) {
	polyPars.push_back(minuitx_->GetParameter(ip));
	vector<double> polyCov;
	for ( size_t jp=parIndices[i]; jp<parIndices[i+1]; ++jp ) {
	  polyCov.push_back(minuitx_->GetCovarianceMatrixElement(ip,jp));
	}
	polyCovs.push_back(polyCov);
      }
      for ( size_t iht=0; iht<htMins_.size(); ++iht ) {
	graph->SetPoint(iht,htMins_[iht],htPolynomial(htMins_[iht],polyPars));
	graph->SetPointError(iht,0.,htPolynomialError(htMins_[iht],polyCovs));
      }
    }
    else {
      for ( size_t iht=0; iht<htMins_.size(); ++iht ) {
	graph->SetPoint(iht,htMins_[iht],minuitx_->GetParameter(parIndices[i]+iht));
	graph->SetPointError(iht,0.,minuitx_->GetParError(parIndices[i]+iht));
      }
    }
  }

  drawGraphs();

}

void
LogLMultiErf::drawGraphs () const
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
LogLMultiErf::writeGraphs (string name) const
{
  TFile f(name.c_str(),"recreate");
  for ( size_t i=0; i<parGraphs_.size(); ++i ) {
    TGraphErrors* clone = (TGraphErrors*)parGraphs_[i]->Clone();
    clone->Write();
  }
  f.Write();
  f.Close();
}

