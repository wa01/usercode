#include "TF1.h"
#include "TH1F.h"
#include "TEfficiency.h"
#include "TRandom2.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TKey.h"
#include "TVectorT.h"
#include "TMatrixDSym.h"
#include "TObjString.h"
#include "TObjArray.h"

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

struct FitResults {
  FitResults () {}
  FitResults (FixMask _mask, const vector<unsigned int>& _indices,
	      TFitterMinuit& fitter) : mask(_mask), indices(_indices) {
    cout << "LastIndex: " << indices.back() << endl;
    for ( size_t i=0; i<indices.back(); ++i ) {
      names.push_back(fitter.GetParName(i));
      parameters.push_back(fitter.GetParameter(i));
      vector<double> row;
      for ( size_t j=0; j<indices.back(); ++j ) {
	if ( !fitter.IsFixed(i) && !fitter.IsFixed(j) )
	  row.push_back(fitter.GetCovarianceMatrixElement(i,j));
	else
	  row.push_back(0.);
      }
      covariance.push_back(row);
    }
  }
  TObjArray* namesArray () const {
    TObjArray* result = new TObjArray();
    result->SetOwner(1);
    for ( size_t i=0; i<parameters.size(); ++i ) {
      result->Add(new TObjString(names[i].c_str()));
    }
    return result;
  }
  TVectorD* parTVector () const {
    TVectorD* result = new TVectorD(parameters.size());
    for ( size_t i=0; i<parameters.size(); ++i ) 
      (*result)(i) = parameters[i];
    return result;
  }
  TMatrixDSym* covTMatrix () const {
    TMatrixDSym* result = new TMatrixDSym(parameters.size());
    for ( size_t i=0; i<parameters.size(); ++i )  {
      for ( size_t j=0; j<parameters.size(); ++j )  {
	(*result)(i,j) = covariance[i][j];
      }
    }
    return result;
  }
  FixMask mask;
  vector<unsigned int> indices;
  vector<string> names;
  vector<double> parameters;
  vector< vector<double> > covariance;
};

TCanvas* drawEfficiencies (const vector<TEfficiency*>& efficiencies,
			   const vector<TF1*>& erfFunctions) 
{
  TCanvas* cnv = new TCanvas("cErf","cErf",800,600);
  int npad = int(sqrt(efficiencies.size()/12.)+0.5);
  cnv->Divide(4*npad,3*npad);
  
  for ( size_t iht=0; iht<efficiencies.size(); ++iht ) {
    cnv->cd(iht+1);
    TEfficiency* eff = efficiencies[iht];
    eff->SetMarkerStyle(24);
    eff->Draw();
    
    TF1* f = erfFunctions[iht];
    f->SetLineColor(2);
    f->SetLineWidth(2);
    f->Draw("same");
  } 
  
  return cnv;
}


TCanvas* drawGraphs (const vector<TGraphErrors*>& parGraphsFixed,
		     const vector<TGraphErrors*>& parGraphsFinal) 
{
  TCanvas* cpar = new TCanvas("cpar","cpar",800,300);
  cpar->Divide(3,1);
  for ( size_t i=0; i<min(parGraphsFinal.size(),size_t(3)); ++i ) {
    cpar->cd(i+1);
    parGraphsFinal[i]->Draw("AP");
    cout << parGraphsFinal[i]->GetN() << endl;
    if ( i<parGraphsFixed.size() ) {
      parGraphsFixed[i]->SetMarkerStyle(1);
      parGraphsFixed[i]->SetMarkerColor(2);
      parGraphsFixed[i]->SetLineColor(2);
      parGraphsFixed[i]->Draw("L");
    }
  }
  cpar->cd();
  return cpar;
}

template <class T>
vector<T*> getObjects (TDirectory& dir)
{
  vector<T*> result;

  TIter it(dir.GetListOfKeys());
  TKey* key;
  while ( (key=(TKey*)it.Next()) ) {
    TObject* obj = key->ReadObj();
    if ( obj->InheritsFrom(T::Class()) )  result.push_back((T*)obj);
  }
  return result;
}

template <class T>
T* getObject (TDirectory& dir)
{
  vector<T*> objects = getObjects<T>(dir);
  return objects.empty() ? (T*)0 : objects.front();
}

void drawResults (TFile* file) 
{
  TDirectory* dir = gDirectory;
  file->cd("Efficiencies");
  vector<TEfficiency*> efficiencies = getObjects<TEfficiency>(*gDirectory);
  file->cd("Erfs");
  vector<TF1*> erfs = getObjects<TF1>(*gDirectory);
  file->cd("graphsFixed");
  vector<TGraphErrors*> graphsFixed = getObjects<TGraphErrors>(*gDirectory);
  file->cd("graphsFinal");
  vector<TGraphErrors*> graphsFinal = getObjects<TGraphErrors>(*gDirectory);

  drawEfficiencies(efficiencies,erfs);
  drawGraphs(graphsFixed,graphsFinal);

  file->cd("fixedResults");
  TObjArray* fixedNames = getObject<TObjArray>(*gDirectory);
  TVectorD* fixedPars = getObject<TVectorD>(*gDirectory);
  TMatrixDSym* fixedCov = getObject<TMatrixDSym>(*gDirectory);
  cout << fixedPars->GetNrows() << " " << fixedCov->GetNrows() << " " << fixedCov->GetNcols() << endl;
  for ( int i=0; i<fixedPars->GetNrows(); ++i ) {
    TObjString* str = (TObjString*)(*fixedNames)[i];
    printf("%-20s %12.3f +- %10.3f\n",str->GetString().Data(),(*fixedPars)(i),sqrt((*fixedCov)(i,i)));
  }
  dir->cd();
}


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
  
  void writeResults (string name="logLMultiErf_Results.root") const;

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
  FitResults fitSinglePass (FixMask mask);
  vector<TGraphErrors*> createGraphs (const FitResults& fitResults, string postfix = "");

private:
  string path_;
  vector<float> htMins_;
  float metMin_;
  float metMax_;

  TFitterMinuit* minuitx_;
  FcnMultiLogL* fcn_;
  vector<unsigned int> parIndices_;

  DataContainer data_;
  FitResults fixedResults_;
  FitResults finalResults_;
  vector<TEfficiency*> efficiencies_;
  vector<TF1*> erfFunctions_;
  vector<TGraphErrors*> parGraphsFixed_;
  vector<TGraphErrors*> parGraphsFinal_;
};

vector<TGraphErrors*>
LogLMultiErf::createGraphs (const FitResults& fitResults, string postfix) {
  // char hName[128], hTitle[256];
  string name;
  vector<TGraphErrors*> result;

  for ( size_t i=0; i<3; ++i ) {
    TGraphErrors* graph = new TGraphErrors();
    switch ( i ) {
    case 0: 
      name = "Location" + postfix;
      graph->SetNameTitle(name.c_str(),name.c_str());
      break;
    case 1: 
      name = "Scale" + postfix;
      graph->SetNameTitle(name.c_str(),name.c_str());
      break;
    case 2: 
      name = "Shape" + postfix;
      graph->SetNameTitle(name.c_str(),name.c_str());
      break;
    }
    graph->SetMarkerStyle(20);
    result.push_back(graph);
  }

  // vector<double> resultPars;
  // for ( size_t i=0; i<5; ++i ) {
  //   for ( size_t j=fitResults.indices[i]; j<fitResults.indices[i+1]; ++j ) {
  //     resultPars.push_back(fitResults.parameters[j]);
  //     // cout << " " << resultPars.back();
  //   }
  // }
  // // cout << endl;
  // vector< vector<double>::const_iterator > parIterators = fcn_->getParIterators(resultPars);
  // //
  // // response function parameters / bin
  // //
  // vector< vector<double> > resFunParsLow = fcn_->getResFunPars(parIterators,true);
  // // vector< vector<double> > resFunParsHigh = fcn_->getResFunPars(parIterators,false);

  //   fcn_->setHTmin(htMin[iht]);
  //   fitSingleHT (hAll,hAllS1,htMin[iht]);
  for ( size_t i=0; i<result.size(); ++i ) {
    TGraphErrors* graph = result[i];
    if ( fitResults.mask.isFixed(i) ) {
      vector<double> polyPars;
      vector< vector<double> > polyCovs;
      for ( size_t ip=fitResults.indices[i]; ip<fitResults.indices[i+1]; ++ip ) {
	polyPars.push_back(fitResults.parameters[ip]);
	vector<double> polyCov;
	for ( size_t jp=fitResults.indices[i]; jp<fitResults.indices[i+1]; ++jp ) {
	  polyCov.push_back(fitResults.covariance[ip][jp]);
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
	graph->SetPoint(iht,htMins_[iht],minuitx_->GetParameter(fitResults.indices[i]+iht));
	graph->SetPointError(iht,0.,minuitx_->GetParError(fitResults.indices[i]+iht));
      }
    }
  }
  return result;
}
//
// Fit
//
FitResults
LogLMultiErf::fitSinglePass (FixMask mask) {
  cout << "Starting single pass" << endl;
  //
  // LL function
  //
  delete minuitx_;
  // delete fcn_;
  fcn_ = new FcnMultiLogL(data_,htMins_,mask);
  minuitx_ = new TFitterMinuit();
  minuitx_->SetMinuitFCN(fcn_);
  //
  char parName[128];
  unsigned int index(0);
  parIndices_.clear();
  // location
  parIndices_.push_back(index);
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
  parIndices_.push_back(index);
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
  parIndices_.push_back(index);
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
  parIndices_.push_back(index);
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
  parIndices_.push_back(index);
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
  parIndices_.push_back(index);
  
  minuitx_->SetMaxIterations(10000);
  minuitx_->SetPrintLevel(3);
  minuitx_->CreateMinimizer();

  //
  // first fit: location and scale
  //
  for ( size_t i=parIndices_[1]; i<parIndices_[5]; ++i )  minuitx_->FixParameter(i);
  int ierr = minuitx_->Minimize();
  for ( size_t i=parIndices_[0]; i<parIndices_[1]; ++i )  minuitx_->FixParameter(i);
  for ( size_t i=parIndices_[1]; i<parIndices_[2]; ++i )  minuitx_->ReleaseParameter(i);
  ierr = minuitx_->Minimize();
  for ( size_t i=parIndices_[0]; i<parIndices_[1]; ++i )  minuitx_->ReleaseParameter(i);
  ierr = minuitx_->Minimize();
  for ( size_t i=parIndices_[0]; i<parIndices_[2]; ++i )  minuitx_->FixParameter(i);
  for ( size_t i=parIndices_[2]; i<parIndices_[3]; ++i )  minuitx_->ReleaseParameter(i);
  ierr = minuitx_->Minimize();
  for ( size_t i=parIndices_[0]; i<parIndices_[2]; ++i )  minuitx_->ReleaseParameter(i);
  ierr = minuitx_->Minimize();

  return FitResults(mask,parIndices_,*minuitx_);
}

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

  // TCanvas* cnv = new TCanvas("c","c",800,600);
  // int npad = int(sqrt(htMins_.size()/12.)+0.5);
  // cnv->Divide(4*npad,3*npad);
  // //
  // // LL function
  // //
  // delete minuitx_;
  // delete fcn_;
  // fcn_ = new FcnMultiLogL(data_,htMins_,mask);
  // minuitx_ = new TFitterMinuit();
  // minuitx_->SetMinuitFCN(fcn_);
  // //
  // char parName[128];
  // unsigned int index(0);
  // vector<unsigned int> parIndices_;
  // // location
  // parIndices_.push_back(index);
  // if ( mask.isFixed(0) ) {
  //   minuitx_->SetParameter(index++,"locP0",0.,100.,-metMax_,metMax_);
  //   minuitx_->SetParameter(index++,"locP1",0.5,0.1,-3.,3.);
  // }
  // else {
  //   for ( size_t i=0; i<htMins_.size(); ++i ) {
  //     sprintf(parName,"locHT%4.4d",int(htMins_[i]+0.5));
  //     minuitx_->SetParameter(index++,parName,htMins_[i],50.,0.,metMax_);
  //   }
  // }
  // // scale
  // parIndices_.push_back(index);
  // if ( mask.isFixed(1) ) {
  //   minuitx_->SetParameter(index++,"scaleP0",150.,10.,0.,500.);
  //   minuitx_->SetParameter(index++,"scaleP1",0.,0.1,-1.,1.);
  // }
  // else {
  //   for ( size_t i=0; i<htMins_.size(); ++i ) {
  //     sprintf(parName,"scaleHT%4.4d",int(htMins_[i]+0.5));
  //     minuitx_->SetParameter(index++,parName,100,20.,0.,500.);
  //   }
  // }
  // // shape
  // parIndices_.push_back(index);
  // if ( mask.isFixed(2) ) {
  //   minuitx_->SetParameter(index++,"shapeP0",0.,0.05,-1.0,1.0);
  // }
  // else {
  //   for ( size_t i=0; i<htMins_.size(); ++i ) {
  //     sprintf(parName,"shapeHT%4.4d",int(htMins_[i]+0.5));
  //     minuitx_->SetParameter(index++,parName,0.,0.05,-1.0,1.0);
  //   }
  // }
  // // emin
  // parIndices_.push_back(index);
  // if ( mask.isFixed(3) ) {
  //   minuitx_->SetParameter(index++,"eminP0",0.,0.05,0.,1.);
  // }
  // else {
  //   for ( size_t i=0; i<htMins_.size(); ++i ) {
  //     sprintf(parName,"eminHT%4.4d",int(htMins_[i]+0.5));
  //     minuitx_->SetParameter(index++,parName,0.,0.05,0.,1.);
  //   }
  // }
  // // emax
  // parIndices_.push_back(index);
  // if ( mask.isFixed(4) ) {
  //   minuitx_->SetParameter(index++,"emaxP0",1.,0.05,0.,1.);
  // }
  // else {
  //   for ( size_t i=0; i<htMins_.size(); ++i ) {
  //     sprintf(parName,"emaxHT%4.4d",int(htMins_[i]+0.5));
  //     minuitx_->SetParameter(index++,parName,1.,0.05,0.,1.);
  //   }
  // }
  // //
  // parIndices_.push_back(index);
  
  // minuitx_->SetMaxIterations(10000);
  // minuitx_->SetPrintLevel(3);
  // minuitx_->CreateMinimizer();

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

  fixedResults_ = fitSinglePass(31);
  finalResults_ = fitSinglePass(mask);
  // //
  // // first fit: location and scale
  // //
  // for ( size_t i=parIndices_[1]; i<parIndices_[5]; ++i )  minuitx_->FixParameter(i);
  // int ierr = minuitx_->Minimize();
  // for ( size_t i=parIndices_[0]; i<parIndices_[1]; ++i )  minuitx_->FixParameter(i);
  // for ( size_t i=parIndices_[1]; i<parIndices_[2]; ++i )  minuitx_->ReleaseParameter(i);
  // ierr = minuitx_->Minimize();
  // for ( size_t i=parIndices_[0]; i<parIndices_[1]; ++i )  minuitx_->ReleaseParameter(i);
  // ierr = minuitx_->Minimize();
  // for ( size_t i=parIndices_[0]; i<parIndices_[2]; ++i )  minuitx_->FixParameter(i);
  // for ( size_t i=parIndices_[2]; i<parIndices_[3]; ++i )  minuitx_->ReleaseParameter(i);
  // ierr = minuitx_->Minimize();
  // for ( size_t i=parIndices_[0]; i<parIndices_[2]; ++i )  minuitx_->ReleaseParameter(i);
  // ierr = minuitx_->Minimize();

  // char hName[128], hTitle[256];
  // parGraphsFinal_.clear(); // should delete previous graphs, if any ...
  // for ( size_t i=0; i<3; ++i ) {
  //   TGraphErrors* graph = new TGraphErrors();
  //   switch ( i ) {
  //   case 0: 
  //     graph->SetNameTitle("Location","Location");
  //     break;
  //   case 1: 
  //     graph->SetNameTitle("Scale","Scale");
  //     break;
  //   case 2: 
  //     graph->SetNameTitle("Shape","Shape");
  //     break;
  //   }
  //   graph->SetMarkerStyle(20);
  //   parGraphsFinal_.push_back(graph);
  // }

  vector<double> resultPars;
  for ( size_t i=0; i<5; ++i ) {
    for ( size_t j=parIndices_[i]; j<parIndices_[i+1]; ++j ) {
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

  efficiencies_.clear();
  erfFunctions_.clear();
  for ( size_t iht=0; iht<htMins_.size(); ++iht ) {
    // cnv->cd(iht+1);
    sprintf(name,"eff%4.4d",int(htMins_[iht]+0.5));
    TEfficiency* eff = new TEfficiency(*hCuts[iht],*hAll);
    eff->SetNameTitle(name,name);
    // eff->SetMarkerStyle(24);
    // eff->Draw();
    efficiencies_.push_back(eff);
    //
    // function corresponding to the result of the unbinned LL fit
    //
    sprintf(name,"fLL%4.4d",int(htMins_[iht]+0.5));
    TF1* f = new TF1(name,responseFunction,hAll->GetXaxis()->GetXmin(),
		     hAll->GetXaxis()->GetXmax(),5);
    for ( size_t k=0; k<5; ++k ) {
      f->SetParameter(k,resFunParsLow[iht][k]);
    }
    // f->SetLineColor(2);
    // f->SetLineWidth(2);
    // f->Draw("same");
    erfFunctions_.push_back(f);
  }
  drawEfficiencies(efficiencies_,erfFunctions_);
  //   fcn_->setHTmin(htMin[iht]);
  //   fitSingleHT (hAll,hAllS1,htMin[iht]);
  // for ( size_t i=0; i<parGraphsFinal_.size(); ++i ) {
  //   TGraphErrors* graph = parGraphsFinal_[i];1
  //   if ( mask.isFixed(i) ) {
  //     vector<double> polyPars;
  //     vector< vector<double> > polyCovs;
  //     for ( size_t ip=parIndices_[i]; ip<parIndices_[i+1]; ++ip ) {
  // 	polyPars.push_back(minuitx_->GetParameter(ip));
  // 	vector<double> polyCov;
  // 	for ( size_t jp=parIndices_[i]; jp<parIndices_[i+1]; ++jp ) {
  // 	  polyCov.push_back(minuitx_->GetCovarianceMatrixElement(ip,jp));
  // 	}
  // 	polyCovs.push_back(polyCov);
  //     }
  //     for ( size_t iht=0; iht<htMins_.size(); ++iht ) {
  // 	graph->SetPoint(iht,htMins_[iht],htPolynomial(htMins_[iht],polyPars));
  // 	graph->SetPointError(iht,0.,htPolynomialError(htMins_[iht],polyCovs));
  //     }
  //   }
  //   else {
  //     for ( size_t iht=0; iht<htMins_.size(); ++iht ) {
  // 	graph->SetPoint(iht,htMins_[iht],minuitx_->GetParameter(parIndices_[i]+iht));
  // 	graph->SetPointError(iht,0.,minuitx_->GetParError(parIndices_[i]+iht));
  //     }
  //   }
  // }

  parGraphsFixed_ = createGraphs(fixedResults_,"Fixed");
  parGraphsFinal_ = createGraphs(finalResults_,"Final");
  drawGraphs(parGraphsFixed_,parGraphsFinal_);

}

void
LogLMultiErf::writeResults (string name) const
{
  assert (efficiencies_.size()==htMins_.size());
  assert (erfFunctions_.size()==htMins_.size());

  TFile f(name.c_str(),"recreate");
  TDirectory* dir = f.mkdir("Efficiencies");
  dir->cd();
  for ( size_t i=0; i<efficiencies_.size(); ++i ) {
    TEfficiency* eff = (TEfficiency*)efficiencies_[i]->Clone();
    eff->Write();
  }
  dir = f.mkdir("Erfs");
  dir->cd();
  for ( size_t i=0; i<erfFunctions_.size(); ++i ) {
    TGraphErrors* erf = (TGraphErrors*)erfFunctions_[i]->Clone();
    erf->Write();
  }
  dir = f.mkdir("graphsFixed");
  dir->cd();
  for ( size_t i=0; i<parGraphsFixed_.size(); ++i ) {
    TGraphErrors* clone = (TGraphErrors*)parGraphsFixed_[i]->Clone();
    clone->Write();
  }
  dir = f.mkdir("graphsFinal");
  dir->cd();
  for ( size_t i=0; i<parGraphsFinal_.size(); ++i ) {
    TGraphErrors* clone = (TGraphErrors*)parGraphsFinal_[i]->Clone();
    clone->Write();
  }

  dir = f.mkdir("fixedResults");
  dir->cd();
  TObjArray* fixedNames = fixedResults_.namesArray();
  fixedNames->Write("fixedNames",1);
  TVectorD* fixedPars = fixedResults_.parTVector();
  fixedPars->Write("fixedPars");
  TMatrixDSym* fixedCov = fixedResults_.covTMatrix();
  fixedCov->Write("fixedCov");
  dir = f.mkdir("finalResults");
  dir->cd();
  TObjArray* finalNames = finalResults_.namesArray();
  finalNames->Write("finalNames",1);
  TVectorD* finalPars = finalResults_.parTVector();
  finalPars->Write("finalPars");
  TMatrixDSym* finalCov = finalResults_.covTMatrix();
  finalCov->Write("finalCov");

  f.Write();
  f.Close();
}

