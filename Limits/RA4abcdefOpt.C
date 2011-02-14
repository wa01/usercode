#include "RA4abcd.h"

#include "TStopwatch.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH2F.h"
#include "TTree.h"
#include "RooPlot.h"
#include "RooAbsPdf.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooGlobalFunc.h"
#include "RooFitResult.h"
#include "RooRandom.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/BayesianCalculator.h"
#include "RooStats/MCMCCalculator.h"
#include "RooStats/MCMCInterval.h"
#include "RooStats/MCMCIntervalPlot.h"
#include "RooStats/ProposalHelper.h"
#include "RooStats/SimpleInterval.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include <vector>

#include "RooAbsDataStore.h"

/** Some functions for limit calculation with ABCD.
 *  RA4Single: one test with LM0
 *  RA4Mult:   scan over m0-m1/2 (signal yields in ABCD provided by histograms) */

using namespace RooFit;
using namespace RooStats;

//
// scan over regions
//
double regionContent (TH2* histo, 
		      int iHTbegin, int iMETbegin, 
		      int iHTend=-1, int iMETend=-1)
{
  double sum(0.);
  // sum from begin to end (excluding)
  //  if end<0: sum to edge of histogram, including overflow bin
  int htEnd = iHTend>0 ? iHTend : histo->GetNbinsX()+2;
  int metEnd = iMETend>0 ? iMETend : histo->GetNbinsY()+2;
  for ( int ix=iHTbegin; ix<htEnd; ++ix ) {
    for ( int iy=iMETbegin; iy<metEnd; ++iy ) {
      sum += histo->GetBinContent(ix,iy);
    }
  }
  return sum;
}

double regionContent (TH2* histo, 
		      float HTbegin, float METbegin, 
		      float HTend=-1, float METend=-1)
{
  TAxis* xaxis = histo->GetXaxis();
  TAxis* yaxis = histo->GetYaxis();
  int iHTbegin = xaxis->FindBin(HTbegin);
  int iHTend = HTend>0 ? xaxis->FindBin(HTend) : -1;
  int iMETbegin = yaxis->FindBin(METbegin);
  int iMETend = METend>0 ? yaxis->FindBin(METend) : -1;
  return regionContent(histo,iHTbegin,iMETbegin,iHTend,iMETend);
}


void setupRegions (int iht0, int iht1, int iht2, int imet0, int imet1, int imet2,
		  TH2* hist, double* evts)
{
  evts[0] = regionContent(hist,iht0,imet0,iht1,imet1);
  evts[1] = regionContent(hist,iht1,imet0,iht2,imet1);
  evts[2] = regionContent(hist,iht2,imet0,-1,imet1);
  evts[3] = regionContent(hist,iht0,imet2,iht1,-1);
  evts[4] = regionContent(hist,iht1,imet2,iht2,-1);
  evts[5] = regionContent(hist,iht2,imet2,-1,-1);

//   evts[0] = regionContent(hist,iht0,imet0,iht1,imet1);
//   evts[1] = regionContent(hist,iht0,imet1,iht1,imet2);
//   evts[2] = regionContent(hist,iht0,imet2,iht1,-1);
//   evts[3] = regionContent(hist,iht2,imet0,-1,imet1);
//   evts[4] = regionContent(hist,iht2,imet1,-1,imet2);
//   evts[5] = regionContent(hist,iht2,imet2,-1,-1);
}

int setupRegions (int iht0, int iht1, int iht2, int imet0, int imet1, int imet2,
		  TH2* hBkg, TH2* hTt, TH2* hWjets, TH2* hSig,
		  double* bkgs, double* tt, double* wjets, double* yields)
{
  TAxis* xaxis = hBkg->GetXaxis();
  TAxis* yaxis = hBkg->GetYaxis();
  std::cout << "Limits " 
	    << xaxis->GetBinLowEdge(iht0) << " "
	    << xaxis->GetBinLowEdge(iht1) << " "
	    << xaxis->GetBinLowEdge(iht2) << " "
	    << yaxis->GetBinLowEdge(imet0) << " "
	    << yaxis->GetBinLowEdge(imet1) << " "
	    << yaxis->GetBinLowEdge(imet2) << std::endl;
	      
  setupRegions(iht0,iht1,iht2,imet0,imet1,imet2,hBkg,bkgs);
  setupRegions(iht0,iht1,iht2,imet0,imet1,imet2,hTt,tt);
  setupRegions(iht0,iht1,iht2,imet0,imet1,imet2,hWjets,wjets);
  setupRegions(iht0,iht1,iht2,imet0,imet1,imet2,hSig,yields);

  if ( bkgs[5]<0.001 || yields[5]<0.001 )  return -1;

  double bkgmin(1.e30);
  double ttmin(1.e30);
  double wjetsmin(1.e30);
  for ( unsigned int i=0; i<6; ++i ) {
    bkgmin = min(bkgmin,bkgs[i]);
    ttmin = min(ttmin,tt[i]);
    wjetsmin = min(wjetsmin,wjets[i]);
  }
  if ( bkgmin<0.001 || ttmin<0.001 || wjetsmin<0.001 ) {
    return 1;
  }
  return 0;
}

void RA4Regions (const char* prefix, const char* postfix, const char* sigName,
		 int index, int dHT, int dMET,
		 float HTCut0, float METCut0,
		 float HTCut1, float METCut1,
		 float HTCut2=-1, float METCut2=-1,
		 StatMethod method=ProfileLikelihoodMethod) {

//   if ( dHT>0 && dMET>0 ) {
//     std::cout << "*** Scanning intermediate limits" << std::endl;
//   }
//   else if ( dHT<0 && dMET<0 ) {
//     std::cout << "*** Scanning lower limits" << std::endl;
//   }
//   else {
  if ( index<0 || index>3 ) {
    std::cout << "*** Inconsistency between dHT and dMET" << std::endl;
    return;
  }
  
  std::string spre(prefix);
  std::string spost(postfix);
  TFile* fBkgRegions = new TFile((spre+"mc"+spost).c_str());
  TFile* fTtRegions = new TFile((spre+"ttbar"+spost).c_str());
  TFile* fWjRegions = new TFile((spre+"wjets"+spost).c_str());
  TFile* fSigRegions = new TFile((spre+sigName+spost).c_str());
  if ( fBkgRegions==0 || fBkgRegions->IsZombie() ||
       fTtRegions==0 || fTtRegions->IsZombie() ||
       fWjRegions==0 || fWjRegions->IsZombie() ||
       fSigRegions==0 || fSigRegions->IsZombie() ) {
    std::cout << "Couldn't open one of the files" << std::endl;
    return;
  }
  TH2* hBkg = (TH2*)fBkgRegions->Get("ROOT.c1")->FindObject("ht_vs_kinMetSig");
  TH2* hTt = (TH2*)fTtRegions->Get("ROOT.c1")->FindObject("ht_vs_kinMetSig");
  TH2* hWjets = (TH2*)fWjRegions->Get("ROOT.c1")->FindObject("ht_vs_kinMetSig");
  TH2* hSig = (TH2*)fSigRegions->Get("ROOT.c1")->FindObject("ht_vs_kinMetSig");

  int iHTCut0 = hBkg->GetXaxis()->FindBin(HTCut0);
  int iHTCut1 = hBkg->GetXaxis()->FindBin(HTCut1);
  int iHTCut2 = HTCut2>0 ? hBkg->GetXaxis()->FindBin(HTCut2) : iHTCut1;
  int iMETCut0 = hBkg->GetYaxis()->FindBin(METCut0);
  int iMETCut1 = hBkg->GetYaxis()->FindBin(METCut1);
  int iMETCut2 = METCut2>0 ? hBkg->GetYaxis()->FindBin(METCut2) : iMETCut1;

  gROOT->cd();
  TH2* hExclusion = (TH2*)hBkg->Clone("Exclusion");
  hExclusion->SetTitle("Exclusion");
  hExclusion->Reset();
  hExclusion->SetMinimum(); hExclusion->SetMaximum();
  TH2* hLowerLimit = (TH2*)hExclusion->Clone("LowerLimit");
  hLowerLimit->SetTitle("LowerLimit");
  TH2* hUpperLimit = (TH2*)hExclusion->Clone("UpperLimit");
  hUpperLimit->SetTitle("UpperLimit");
  TH2* hRelUpperLimit = (TH2*)hExclusion->Clone("RelUpperLimit");
  hRelUpperLimit->SetTitle("UpperLimit / Yield");
  TH2* hKappa = (TH2*)hExclusion->Clone("Kappa");
  hKappa->SetTitle("Kappa");
  TH2* hSigKappa = (TH2*)hExclusion->Clone("SigKappa");
  hSigKappa->SetTitle("SigKappa");
  TH2* hBkgD = (TH2*)hExclusion->Clone("BkgD");
  hBkgD->SetTitle("BkgD");
  TH2* hSigD = (TH2*)hExclusion->Clone("SigD");
  hSigD->SetTitle("SigD");
  TH2* hNevD = (TH2*)hExclusion->Clone("NevD");
  hNevD->SetTitle("NevD");
  TH2* hSoB = (TH2*)hExclusion->Clone("SoB");
  hSoB->SetTitle("SoB");
  TH2* hSoRtB = (TH2*)hExclusion->Clone("SoRtB");
  hSoRtB->SetTitle("SoRtB");

  RooWorkspace* wspace = createWorkspace();

  double yields[6];
  double bkgs[6];
  double tt[6];
  double wjets[6];

  int nbx = hBkg->GetNbinsX();
  int nby = hBkg->GetNbinsY();
// //   int iHTbeg = dHT>0 ? iHTCut+dHT : 1;
//   int iHTend = dHT>0 ? nbx+1 : 0;
// //   int iMETbeg = dMET>0 ? iMETCut+dMET : 1;
//   int iMETend = dMET>0 ? nby+1 : 0;

  int iHT0 = iHTCut0;
  int iHT1 = iHTCut1;
  int iHT2 = iHTCut2;
  int dHT0(0), dHT1(0), dHT2(0);
  int dMET0(0), dMET1(0), dMET2(0);
  if ( index==0 ) {
    iHT0 = 1;
    dHT0 = dHT;
    dMET0 = dMET;
  }
  else if ( index==1 ) {
    iHT1 = iHTCut0 + 1;
    dHT1 = dHT;
    dMET1 = dMET;
  }
  else if ( index==2 ) {
    iHT2 = iHTCut1;
    dHT2 = dHT;
    dMET2 = dMET;
  }
  else {
    iHT1 = iHTCut0 + 1;
    iHT2 = iHT1 + (iHTCut2-iHTCut1);
    dHT1 = dHT2 = dHT;
    dMET1 = dMET2 = dMET;
  }
  int ix,iy;
//   for ( int ix=iHTCut+dHT; ix!=iHTend; ix+=dHT ) {
//     for ( int iy=iMETCut+dHT; iy!=iMETend; iy+=dMET ) {

  for ( ; iHT0<iHT1&&iHT1<=iHT2&&iHT2<nbx+1; 
	iHT0+=dHT0,iHT1+=dHT1,iHT2+=dHT2 ) {
    int iMET0 = iMETCut0;
    int iMET1 = iMETCut1;
    int iMET2 = iMETCut2;
    if ( index==0 ) {
      iMET0 = 1;
    }
    else if ( index==1 ) {
      iMET1 = iMETCut0 + 1;
    }
    else if ( index==2 ) {
      iMET2 = iMETCut1;
    }
    else {
      iMET1 = iMETCut0 + 1;
      iMET2 = iMET1 + (iMETCut2-iMETCut1);
    }
    for ( ; iMET0<iMET1&&iMET1<=iMET2&&iMET2<nby+1; 
	  iMET0+=dMET0,iMET1+=dMET1,iMET2+=dMET2 ) {

      if ( index==0 ) {
	ix = iHT0; iy = iMET0;
      }
      else if ( index==1 ) {
	ix = iHT1; iy = iMET1;
      }
      else {
	ix = iHT2; iy = iMET2;
      }

      int err = setupRegions(iHT0,iHT1,iHT2,iMET0,iMET1,iMET2,
			     hBkg,hTt,hWjets,hSig,bkgs,tt,wjets,yields);
      if ( err!=0 ) {
	hExclusion->SetBinContent(ix,iy,-1.);
	hLowerLimit->SetBinContent(ix,iy,-1.);
	hUpperLimit->SetBinContent(ix,iy,-1.);
	hRelUpperLimit->SetBinContent(ix,iy,-1.);
	continue;
      }

      std::cout << "bkgs / yields =";
      for ( unsigned int i=0; i<6; ++i ) 
	std::cout << " ( " << bkgs[i] << " / " << yields[i] << " ) ";
      std::cout << std::endl;
      hBkgD->SetBinContent(ix,iy,bkgs[5]);
      hSigD->SetBinContent(ix,iy,yields[5]);
      hSoB->SetBinContent(ix,iy,yields[5]/bkgs[5]);
      hSoRtB->SetBinContent(ix,iy,yields[5]/sqrt(bkgs[5]));

      double kappa = (bkgs[0]*bkgs[3])/(bkgs[1]*bkgs[2]);
      std::cout << " kappa = " << kappa << std::endl;
      hKappa->SetBinContent(ix,iy,kappa);

      double kappatt = (tt[0]*tt[3])/(tt[1]*tt[2]);
      double kappawjets = (wjets[0]*wjets[3])/(wjets[1]*wjets[2]);
      std::cout << " kappa (tot/tt/w) = " 
		<< kappa << " " << kappatt << " " << kappawjets << std::endl;
      //
      // set uncertainty on kappa:
      //   deviation from 1 (+) diff. tt/W
      //
      double sigma_kappa_abs = kappa - 1.;
      double sigma_kappa_delta = kappatt - kappawjets;
      double sigma_kappa = sqrt(sigma_kappa_abs*sigma_kappa_abs+sigma_kappa_delta*sigma_kappa_delta+0.1*0.1);
      std::cout << "Setting uncertainty on kappa to " 
		<< sigma_kappa_abs << " " << sigma_kappa_delta << " " << sigma_kappa << std::endl;
      hSigKappa->SetBinContent(ix,iy,sigma_kappa);

      if ( fabs(kappa-1.)>0.2 ) {
	hExclusion->SetBinContent(ix,iy,-3.);
	hLowerLimit->SetBinContent(ix,iy,-3.);
	hUpperLimit->SetBinContent(ix,iy,-3.);
	hRelUpperLimit->SetBinContent(ix,iy,-3.);
	continue;	
      }

      setBackgrounds(wspace,bkgs);
      setSignal(wspace,yields);
      setValRange(wspace,"sigmaKappa",sigma_kappa);

      // wspace->Print("v");
      // RooArgSet allVars = wspace->allVars();
      // allVars.printLatex(std::cout,1);

      RooDataSet* data = new RooDataSet("data","data",*wspace->set("obs"));
      data->add(*wspace->set("obs"));
      data->Print("v");
  
      hNevD->SetBinContent(ix,iy,wspace->var("nd")->getVal());

      MyLimit limit = computeLimit(wspace,data,method);
      std::cout << "Checked ( " << hExclusion->GetXaxis()->GetBinCenter(ix) << " , "
		<< hExclusion->GetYaxis()->GetBinCenter(iy) << " ) with signal yield " << yields[5] << std::endl;
      std::cout << "  Limit [ " << limit.lowerLimit << " , "
		<< limit.upperLimit << " ] ; isIn = " << limit.isInInterval << std::endl;
      std::cout << "  yields =" 
		<< " " << yields[0]
		<< " " << yields[1]
		<< " " << yields[2]
		<< " " << yields[3] << std::endl;
      double excl = limit.isInInterval;
      if ( limit.upperLimit<limit.lowerLimit )  excl = -1;
      hExclusion->SetBinContent(ix,iy,excl);
      hLowerLimit->SetBinContent(ix,iy,limit.lowerLimit);
      hUpperLimit->SetBinContent(ix,iy,limit.upperLimit);
      hRelUpperLimit->SetBinContent(ix,iy,limit.upperLimit/yields[5]);

      delete data;

    }
  }

  TFile* out = new TFile("RA4regions.root","RECREATE");
  hExclusion->SetDirectory(out);
  hLowerLimit->SetDirectory(out);
  hUpperLimit->SetDirectory(out);
  hRelUpperLimit->SetDirectory(out);
  hKappa->SetDirectory(out);
  hSigKappa->SetDirectory(out);
  hBkgD->SetDirectory(out);
  hSigD->SetDirectory(out);
  hNevD->SetDirectory(out);
  hSoB->SetDirectory(out);
  hSoRtB->SetDirectory(out);
  out->Write();
  delete out;
}

struct OptResult {
  bool operator< (const OptResult& other) const {
    return upperLimit/yields[5]<other.upperLimit/other.yields[5];
  }
  int htBins[3];
  int metBins[3];
//   float htCuts[3];
//   float metCuts[3];
  int nobs[6];
  float kappa1;
  float kappa1TT;
  float kappa1Wjets;
  float sigKappa1;
  float kappa2;
  float kappa2TT;
  float kappa2Wjets;
  float sigKappa2;
  float bkgs[6];
  float yields[6];
  float upperLimit;
//   float relUpperLimit;
};

void RA4RegionsTot (const char* prefix, const char* postfix, const char* sigName,
		    int dHT=2, int dMET=2,
		    float HTmin=0., float METmin=0.,
		    StatMethod method=ProfileLikelihoodMethod) {
  
  std::string spre(prefix);
  std::string spost(postfix);
  TFile* fBkgRegions = new TFile((spre+"mc"+spost).c_str());
  TFile* fTtRegions = new TFile((spre+"ttbar"+spost).c_str());
  TFile* fWjRegions = new TFile((spre+"wjets"+spost).c_str());
  TFile* fSigRegions = new TFile((spre+sigName+spost).c_str());
  if ( fBkgRegions==0 || fBkgRegions->IsZombie() ||
       fTtRegions==0 || fTtRegions->IsZombie() ||
       fWjRegions==0 || fWjRegions->IsZombie() ||
       fSigRegions==0 || fSigRegions->IsZombie() ) {
    std::cout << "Couldn't open one of the files" << std::endl;
    return;
  }
  TH2* hBkg = (TH2*)fBkgRegions->Get("ROOT.c1")->FindObject("ht_vs_kinMetSig");
  TH2* hTt = (TH2*)fTtRegions->Get("ROOT.c1")->FindObject("ht_vs_kinMetSig");
  TH2* hWjets = (TH2*)fWjRegions->Get("ROOT.c1")->FindObject("ht_vs_kinMetSig");
  TH2* hSig = (TH2*)fSigRegions->Get("ROOT.c1")->FindObject("ht_vs_kinMetSig");

  RooWorkspace* wspace = createWorkspace();
  RooRealVar* vna = wspace->var("na");
  RooRealVar* vnb = wspace->var("nb");
  RooRealVar* vnc = wspace->var("nc");
  RooRealVar* vnd = wspace->var("nd");
  RooRealVar* vne = wspace->var("ne");
  RooRealVar* vnf = wspace->var("nf");

  double yields[6];
  double bkgs[6];
  double tt[6];
  double wjets[6];

  const unsigned int nrmax(100000);
  const unsigned int nrhmax(500);
  std::vector<OptResult> results; results.reserve(nrmax);
  OptResult result;

  int nbx = hBkg->GetNbinsX();
  int nby = hBkg->GetNbinsY();
  TAxis* xaxis = hBkg->GetXaxis();
  TAxis* yaxis = hBkg->GetYaxis();
  int iHTmin = max(xaxis->FindBin(HTmin),1);
  int iMETmin = max(yaxis->FindBin(METmin),1);
  for ( int ix0=iHTmin; ix0<=nbx; ix0+=dHT ) {
    for ( int ix1=ix0+dHT; ix1<=nbx; ix1+=dHT ) {
      for ( int ix2=ix1; ix2<=nbx; ix2+=dHT ) {
	for ( int iy0=iMETmin; iy0<=nby; iy0+=dMET ) {
	  for ( int iy1=iy0+dMET; iy1<=nby; iy1+=dMET ) {
	    for ( int iy2=iy1; iy2<=nby; iy2+=dMET ) {

	      int err = setupRegions(ix0,ix1,ix2,iy0,iy1,iy2,
				     hBkg,hTt,hWjets,hSig,
				     bkgs,tt,wjets,yields);

	      if ( err<0 )  break;
	      if ( err>0 )  continue;

	      double kappa1 = (bkgs[0]*bkgs[5])/(bkgs[2]*bkgs[3]);
	      double kappa1tt = (tt[0]*tt[5])/(tt[2]*tt[3]);
	      double kappa1wjets = (wjets[0]*wjets[5])/(wjets[2]*wjets[3]);
// 	      std::cout << " kappa (tot/tt/w) = " 
// 			<< kappa << " " << kappatt << " " << kappawjets << std::endl;
	      //
	      // set uncertainty on kappa:
	      //   deviation from 1 (+) diff. tt/W
	      //
	      double sigma_kappa1_abs = kappa1 - 1.;
	      double sigma_kappa1_delta = kappa1tt - kappa1wjets;
	      double sigma_kappa1 = sqrt(sigma_kappa1_abs*sigma_kappa1_abs+sigma_kappa1_delta*sigma_kappa1_delta+0.1*0.1);
// 	      std::cout << "Setting uncertainty on kappa to " 
// 			<< sigma_kappa_abs << " " << sigma_kappa_delta << " " << sigma_kappa << std::endl;

	      double kappa2 = (bkgs[1]*bkgs[5])/(bkgs[2]*bkgs[4]);
	      double kappa2tt = (tt[1]*tt[5])/(tt[2]*tt[4]);
	      double kappa2wjets = (wjets[1]*wjets[5])/(wjets[2]*wjets[4]);
// 	      std::cout << " kappa (tot/tt/w) = " 
// 			<< kappa << " " << kappatt << " " << kappawjets << std::endl;
	      //
	      // set uncertainty on kappa:
	      //   deviation from 1 (+) diff. tt/W
	      //
	      double sigma_kappa2_abs = kappa2 - 1.;
	      double sigma_kappa2_delta = kappa2tt - kappa2wjets;
	      double sigma_kappa2 = sqrt(sigma_kappa2_abs*sigma_kappa2_abs+sigma_kappa2_delta*sigma_kappa2_delta+0.1*0.1);
// 	      std::cout << "Setting uncertainty on kappa to " 
// 			<< sigma_kappa_abs << " " << sigma_kappa_delta << " " << sigma_kappa << std::endl;

	      double kappa = (kappa1+kappa2)/2.;
	      double sigma_kappa = max(sigma_kappa1,sigma_kappa2);

	      if ( yields[5]<0.001 ) {
		break;
	      }
	      if ( fabs(kappa1-1.)>0.1 || fabs(kappa2-1.)>0.1 ) {
		continue;	
	      }
	      if ( sigma_kappa>1. ) {
		continue;
	      }

	      setBackgrounds(wspace,bkgs);
	      setSignal(wspace,yields);
	      setValRange(wspace,"sigmaKappa",sigma_kappa);
	      
// 	      wspace->Print("v");
// 	      RooArgSet allVars = wspace->allVars();
// 	      allVars.printLatex(std::cout,1);
	      
	      RooDataSet data("data","data",*wspace->set("obs"));
	      data.add(*wspace->set("obs"));

	      result.nobs[0] = int(vna->getVal()+0.5);
	      result.nobs[1] = int(vnb->getVal()+0.5);
	      result.nobs[2] = int(vnc->getVal()+0.5);
	      result.nobs[3] = int(vnd->getVal()+0.5);
	      result.nobs[4] = int(vne->getVal()+0.5);
	      result.nobs[5] = int(vnf->getVal()+0.5);
// 	      std::cout << result.nobs[0] << " "
// 			<< result.nobs[1] << " "
// 			<< result.nobs[2] << " "
// 			<< result.nobs[3] << std::endl;

	      MyLimit limit = computeLimit(wspace,&data,method);
// 	      std::cout << "  Limit [ " << limit.lowerLimit << " , "
// 			<< limit.upperLimit << " ] ; isIn = " << limit.isInInterval << std::endl;
// 	      std::cout << "  yields =" 
// 			<< " " << yields[0]
// 			<< " " << yields[1]
// 			<< " " << yields[2]
// 			<< " " << yields[3] << std::endl;
// 	      double excl = limit.isInInterval;
	      if ( limit.upperLimit<0.001 )  continue;

	      
	      result.htBins[0] = ix0;
	      result.htBins[1] = ix1;
	      result.htBins[2] = ix2;
	      result.metBins[0] = iy0;
	      result.metBins[1] = iy1;
	      result.metBins[2] = iy2;
	      for ( int i=0; i<6; ++i ) {
		result.bkgs[i] = bkgs[i];
		result.yields[i] = yields[i];
	      }
	      result.kappa1 = kappa1;
	      result.kappa1TT = kappa1tt;
	      result.kappa1Wjets = kappa1wjets;
	      result.sigKappa1 = sigma_kappa1;
	      result.kappa2 = kappa2;
	      result.kappa2TT = kappa2tt;
	      result.kappa2Wjets = kappa2wjets;
	      result.sigKappa2 = sigma_kappa2;
	      result.upperLimit = limit.upperLimit;
// 	      result.relUpperLimit = limit.upperLimit/yields[3];

	      if ( results.size()<(nrmax-1) ) {
		results.push_back(result);
	      }
	      else {
		std::vector<OptResult>::iterator imax = 
		  std::max_element(results.begin(),results.end());
		*imax = result;
	      }

// 	      delete data;
	    }
	  }
	}
      }
    }
  }

  std::sort(results.begin(),results.end());

  float htCuts[3];
  float metCuts[3];

  TFile* out = new TFile("RA4tot.root","RECREATE");
  TTree* tree = new TTree("RA4opt","RA4opt");
  tree->Branch("nev",result.nobs,"nev[6]/I");
  tree->Branch("bkg",result.bkgs,"bkg[6]/F");
  tree->Branch("sig",result.yields,"sig[6]/F");
  tree->Branch("ht",htCuts,"ht[3]/F");
  tree->Branch("met",metCuts,"met[3]/F");
  tree->Branch("kappa1",&result.kappa1,"kappa1/F");
  tree->Branch("kappa1TT",&result.kappa1TT,"kappa1TT/F");
  tree->Branch("kappa1Wjets",&result.kappa1Wjets,"kappa1Wjets/F");
  tree->Branch("sigKappa1",&result.sigKappa1,"sigKappa1/F");
  tree->Branch("kappa2",&result.kappa2,"kappa2/F");
  tree->Branch("kappa2TT",&result.kappa2TT,"kappa2TT/F");
  tree->Branch("kappa2Wjets",&result.kappa2Wjets,"kappa2Wjets/F");
  tree->Branch("sigKappa2",&result.sigKappa2,"sigKappa2/F");
  tree->Branch("upperLimit",&result.upperLimit,"upperLimit/F");

  TH1* hRelUpperLimit = new TH1F("RelUpperLimit","RelUpperLimit",nrhmax,0,nrhmax);
  TH1* hUpperLimit = new TH1F("UpperLimit","UpperLimit",nrhmax,0,nrhmax);
  TH1* hKappa1 = new TH1F("Kappa1","Kappa1",nrhmax,0,nrhmax);
  TH1* hSigKappa1 = new TH1F("sigKappa1","sigKappa1",nrhmax,0,nrhmax);
  TH1* hKappa2 = new TH1F("Kappa2","Kappa2",nrhmax,0,nrhmax);
  TH1* hSigKappa2 = new TH1F("sigKappa2","sigKappa2",nrhmax,0,nrhmax);
  TH1* hHT0 = new TH1F("HT0","HT0",nrhmax,0,nrhmax);
  TH1* hHT1 = new TH1F("HT1","HT1",nrhmax,0,nrhmax);
  TH1* hHT2 = new TH1F("HT2","HT2",nrhmax,0,nrhmax);
  TH1* hMET0 = new TH1F("MET0","MET0",nrhmax,0,nrhmax);
  TH1* hMET1 = new TH1F("MET1","MET1",nrhmax,0,nrhmax);
  TH1* hMET2 = new TH1F("MET2","MET2",nrhmax,0,nrhmax);
  for ( unsigned int ir=0; ir<results.size(); ++ir ) {
    result = results[ir];
    for ( int i=0; i<3; ++i ) {
      htCuts[i] = xaxis->GetBinLowEdge(result.htBins[i]);
      metCuts[i] = yaxis->GetBinLowEdge(result.metBins[i]);
    }
    tree->Fill();
    if ( ir<nrhmax ) {
      std::cout << result.upperLimit/result.yields[5] << " "
		<< result.upperLimit << " "
		<< result.bkgs[5] << " "
		<< result.yields[5] << " "
		<< result.kappa1 << " "
		<< result.sigKappa1 << " "
		<< result.kappa2 << " "
		<< result.sigKappa2 << " ";
      for ( unsigned int i=0; i<3; ++i ) std::cout << result.htBins[i] << " ";
      for ( unsigned int i=0; i<3; ++i ) std::cout << result.metBins[i] << " ";
      std::cout << std::endl;
      hRelUpperLimit->SetBinContent(ir+1,result.upperLimit/result.yields[5]);
      hUpperLimit->SetBinContent(ir+1,result.upperLimit);
      hKappa1->SetBinContent(ir+1,result.kappa1);
      hSigKappa1->SetBinContent(ir+1,result.sigKappa1);
      hKappa2->SetBinContent(ir+1,result.kappa2);
      hSigKappa2->SetBinContent(ir+1,result.sigKappa2);
      hHT0->SetBinContent(ir+1,htCuts[0]);
      hHT1->SetBinContent(ir+1,htCuts[1]);
      hHT2->SetBinContent(ir+1,htCuts[2]);
      hMET0->SetBinContent(ir+1,metCuts[0]);
      hMET1->SetBinContent(ir+1,metCuts[1]);
      hMET2->SetBinContent(ir+1,metCuts[2]);
    }
    
  }
  out->Write();
  delete out;
}

