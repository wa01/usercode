#include "RA4abcd.h"

#include "TStopwatch.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH2F.h"
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
  int metEnd = iMETend>0 ? iMETend : histo->GetNbinsX()+2;
  for ( int ix=iHTbegin; ix<htEnd; ++ix ) {
    for ( int iy=iMETbegin; iy<metEnd; ++iy ) {
      sum += histo->GetBinContent(ix,iy);
    }
  }
  return sum;
}

void RA4Regions (const char* prefix, const char* postfix, const char* sigName,
		 float HTCut, float METCut,
		 int dHT=1, int dMET=1,
		 StatMethod method=ProfileLikelihoodMethod) {

  if ( dHT>0 && dMET>0 ) {
    std::cout << "*** Scanning intermediate limits" << std::endl;
  }
  else if ( dHT<0 && dMET<0 ) {
    std::cout << "*** Scanning lower limits" << std::endl;
  }
  else {
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

  int iHTCut = hBkg->GetXaxis()->FindBin(HTCut);
  int iMETCut = hBkg->GetYaxis()->FindBin(METCut);

  gROOT->cd();
  TH2* hExclusion = (TH2*)hBkg->Clone("Exclusion");
  hExclusion->Reset();
  hExclusion->SetTitle("Exclusion");
  TH2* hLowerLimit = (TH2*)hBkg->Clone("LowerLimit");
  hLowerLimit->Reset();
  hLowerLimit->SetTitle("LowerLimit");
  TH2* hUpperLimit = (TH2*)hBkg->Clone("UpperLimit");
  hUpperLimit->Reset();
  hUpperLimit->SetTitle("UpperLimit");
  TH2* hRelUpperLimit = (TH2*)hBkg->Clone("RelUpperLimit");
  hRelUpperLimit->Reset();
  hRelUpperLimit->SetTitle("UpperLimit / Yield");
  TH2* hKappa = (TH2*)hBkg->Clone("Kappa");
  hKappa->Reset();
  hKappa->SetTitle("Kappa");
  TH2* hSigKappa = (TH2*)hBkg->Clone("SigKappa");
  hSigKappa->Reset();
  hSigKappa->SetTitle("SigKappa");

  RooWorkspace* wspace = createWorkspace();

  double yields[4];
  double bkgs[4];
  double tt[4];
  double wjets[4];

  int nbx = hBkg->GetNbinsX();
  int nby = hBkg->GetNbinsY();
//   int iHTbeg = dHT>0 ? iHTCut+dHT : 1;
  int iHTend = dHT>0 ? nbx+1 : 0;
//   int iMETbeg = dMET>0 ? iMETCut+dMET : 1;
  int iMETend = dMET>0 ? nby+1 : 0;
  for ( int ix=iHTCut+dHT; ix!=iHTend; ix+=dHT ) {
    for ( int iy=iMETCut+dHT; iy!=iMETend; iy+=dMET ) {

      int iHTlow = dHT>0 ? iHTCut : ix;
      int iHTint = dHT>0 ? ix : iHTCut;
      int iMETlow = dMET>0 ? iMETCut : iy;
      int iMETint = dMET>0 ? iy : iMETCut;
      std::cout << "Limits " 
		<< hBkg->GetXaxis()->GetBinLowEdge(iHTlow) << " "
		<< hBkg->GetXaxis()->GetBinLowEdge(iHTint) << " "
		<< hBkg->GetYaxis()->GetBinLowEdge(iMETlow) << " "
		<< hBkg->GetYaxis()->GetBinLowEdge(iMETint) << std::endl;

      bkgs[0] = regionContent(hBkg,iHTlow,iMETlow,iHTint,iMETint);
      bkgs[1] = regionContent(hBkg,iHTint,iMETlow,-1,iMETint);
      bkgs[2] = regionContent(hBkg,iHTlow,iMETint,iHTint,-1);
      bkgs[3] = regionContent(hBkg,iHTint,iMETint,-1,-1);

      tt[0] = regionContent(hTt,iHTlow,iMETlow,iHTint,iMETint);
      tt[1] = regionContent(hTt,iHTint,iMETlow,-1,iMETint);
      tt[2] = regionContent(hTt,iHTlow,iMETint,iHTint,-1);
      tt[3] = regionContent(hTt,iHTint,iMETint,-1,-1);

      wjets[0] = regionContent(hWjets,iHTlow,iMETlow,iHTint,iMETint);
      wjets[1] = regionContent(hWjets,iHTint,iMETlow,-1,iMETint);
      wjets[2] = regionContent(hWjets,iHTlow,iMETint,iHTint,-1);
      wjets[3] = regionContent(hWjets,iHTint,iMETint,-1,-1);

      yields[0] = regionContent(hSig,iHTlow,iMETlow,iHTint,iMETint);
      yields[1] = regionContent(hSig,iHTint,iMETlow,-1,iMETint);
      yields[2] = regionContent(hSig,iHTlow,iMETint,iHTint,-1);
      yields[3] = regionContent(hSig,iHTint,iMETint,-1,-1);

      std::cout << "bkgs / yields =";
      for ( unsigned int i=0; i<4; ++i ) 
	std::cout << " ( " << bkgs[i] << " / " << yields[i] << " ) ";
      std::cout << std::endl;

      double bkgmin(1.e30);
      double ttmin(1.e30);
      double wjetsmin(1.e30);
      for ( unsigned int i=0; i<4; ++i ) {
	bkgmin = min(bkgmin,bkgs[i]);
	ttmin = min(ttmin,tt[i]);
	wjetsmin = min(wjetsmin,wjets[i]);
      }
      if ( bkgmin<0.001 || ttmin<0.001 || wjetsmin<0.001 ) {
	hExclusion->SetBinContent(ix,iy,-1.);
	hLowerLimit->SetBinContent(ix,iy,-1.);
	hUpperLimit->SetBinContent(ix,iy,-1.);
	hRelUpperLimit->SetBinContent(ix,iy,-1.);
	continue;
      }
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

      if ( yields[3]<0.001 ) {
	hExclusion->SetBinContent(ix,iy,-2.);
	hLowerLimit->SetBinContent(ix,iy,-2.);
	hUpperLimit->SetBinContent(ix,iy,-2.);
	hRelUpperLimit->SetBinContent(ix,iy,-2.);
	continue;
      }
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
  
      MyLimit limit = computeLimit(wspace,data,method);
      std::cout << "Checked ( " << hExclusion->GetXaxis()->GetBinCenter(ix) << " , "
		<< hExclusion->GetYaxis()->GetBinCenter(iy) << " ) with signal yield " << yields[3] << std::endl;
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
      hRelUpperLimit->SetBinContent(ix,iy,limit.upperLimit/yields[3]);

      delete data;

    }
  }

  TFile* out = new TFile("RA4regions.root","RECREATE");
  hExclusion->SetDirectory(out);
  hExclusion->SetMinimum(); hExclusion->SetMaximum();
  hLowerLimit->SetDirectory(out);
  hLowerLimit->SetMinimum(); hLowerLimit->SetMaximum();
  hUpperLimit->SetDirectory(out);
  hUpperLimit->SetMinimum(); hUpperLimit->SetMaximum();
  hRelUpperLimit->SetDirectory(out);
  hRelUpperLimit->SetMinimum(); hRelUpperLimit->SetMaximum();
  hKappa->SetDirectory(out);
  hKappa->SetMinimum(); hKappa->SetMaximum();
  hSigKappa->SetDirectory(out);
  hSigKappa->SetMinimum(); hSigKappa->SetMaximum();
  out->Write();
  delete out;
}

struct OptResult {
//   OptResult (TH2* hist, 
// 	     int ht1, int ht2, int ht3,
// 	     int met1, int met2, int met3,
// 	     double k, double sk,
// 	     double b, double y,
// 	     double ul, double rul) :
//     kappa(k), sigKappa(sk), bkg(b), yield(y),
//     upperLimit(ul), relUpperLimit(rul) {
//     htBins[0] = ht1;
//     htBins[1] = ht2;
//     htBins[2] = ht3;
//     metBins[0] = met1;
//     metBins[1] = met2;
//     metBins[2] = met3;
//   }
  bool operator< (const OptResult& other) const {
    return relUpperLimit<other.relUpperLimit;
  }
  int htBins[3];
  int metBins[3];
  double kappa;
  double sigKappa;
  double bkg;
  double yield;
  double upperLimit;
  double relUpperLimit;
};

void RA4RegionsTot (const char* prefix, const char* postfix, const char* sigName,
		    int dHT=2, int dMET=2,
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

  double yields[4];
  double bkgs[4];
  double tt[4];
  double wjets[4];

  const unsigned int nrmax(500);
  std::vector<OptResult> results; results.reserve(2*nrmax);
  OptResult result;

  int nbx = hBkg->GetNbinsX();
  int nby = hBkg->GetNbinsY();
  for ( int ix0=1; ix0<=nbx; ix0+=dHT ) {
    for ( int ix1=ix0+dHT; ix1<=nbx; ix1+=dHT ) {
      for ( int ix2=ix1; ix2<=nbx; ix2+=dHT ) {
	for ( int iy0=1; iy0<=nby; iy0+=dMET ) {
	  for ( int iy1=iy0+dMET; iy1<=nby; iy1+=dMET ) {
	    for ( int iy2=iy1; iy2<=nby; iy2+=dMET ) {

	      std::cout << "Limits " 
			<< hBkg->GetXaxis()->GetBinLowEdge(ix0) << " "
			<< hBkg->GetXaxis()->GetBinLowEdge(ix1) << " "
			<< hBkg->GetXaxis()->GetBinLowEdge(ix2) << " "
			<< hBkg->GetYaxis()->GetBinLowEdge(iy0) << " "
			<< hBkg->GetYaxis()->GetBinLowEdge(iy1) << " "
			<< hBkg->GetYaxis()->GetBinLowEdge(iy2) << std::endl;
	      
	      bkgs[3] = regionContent(hBkg,ix2,iy2,-1,-1);
	      tt[3] = regionContent(hTt,ix2,iy2,-1,-1);
	      wjets[3] = regionContent(hWjets,ix2,iy2,-1,-1);
	      yields[3] = regionContent(hSig,ix2,iy2,-1,-1);
	      if ( bkgs[3]<0.001 || yields[3]<0.001 )  break;

	      bkgs[0] = regionContent(hBkg,ix0,iy0,ix1,iy1);
	      bkgs[1] = regionContent(hBkg,ix2,iy0,-1,iy1);
	      bkgs[2] = regionContent(hBkg,ix0,iy2,ix1,-1);

	      tt[0] = regionContent(hTt,ix0,iy0,ix1,iy1);
	      tt[1] = regionContent(hTt,ix2,iy0,-1,iy1);
	      tt[2] = regionContent(hTt,ix0,iy2,ix1,-1);

	      wjets[0] = regionContent(hWjets,ix0,iy0,ix1,iy1);
	      wjets[1] = regionContent(hWjets,ix2,iy0,-1,iy1);
	      wjets[2] = regionContent(hWjets,ix0,iy2,ix1,-1);

	      yields[0] = regionContent(hSig,ix0,iy0,ix1,iy1);
	      yields[1] = regionContent(hSig,ix2,iy0,-1,iy1);
	      yields[2] = regionContent(hSig,ix0,iy2,ix1,-1);

// 	      std::cout << "bkgs / yields =";
// 	      for ( unsigned int i=0; i<4; ++i ) 
// 		std::cout << " ( " << bkgs[i] << " / " << yields[i] << " ) ";
// 	      std::cout << std::endl;

	      double bkgmin(1.e30);
	      double ttmin(1.e30);
	      double wjetsmin(1.e30);
	      for ( unsigned int i=0; i<4; ++i ) {
		bkgmin = min(bkgmin,bkgs[i]);
		ttmin = min(ttmin,tt[i]);
		wjetsmin = min(wjetsmin,wjets[i]);
	      }
	      if ( bkgmin<0.001 || ttmin<0.001 || wjetsmin<0.001 ) {
		continue;
	      }
	      double kappa = (bkgs[0]*bkgs[3])/(bkgs[1]*bkgs[2]);

	      double kappatt = (tt[0]*tt[3])/(tt[1]*tt[2]);
	      double kappawjets = (wjets[0]*wjets[3])/(wjets[1]*wjets[2]);
// 	      std::cout << " kappa (tot/tt/w) = " 
// 			<< kappa << " " << kappatt << " " << kappawjets << std::endl;
	      //
	      // set uncertainty on kappa:
	      //   deviation from 1 (+) diff. tt/W
	      //
	      double sigma_kappa_abs = kappa - 1.;
	      double sigma_kappa_delta = kappatt - kappawjets;
	      double sigma_kappa = sqrt(sigma_kappa_abs*sigma_kappa_abs+sigma_kappa_delta*sigma_kappa_delta+0.1*0.1);
// 	      std::cout << "Setting uncertainty on kappa to " 
// 			<< sigma_kappa_abs << " " << sigma_kappa_delta << " " << sigma_kappa << std::endl;

	      if ( yields[3]<0.001 ) {
		break;
	      }
	      if ( fabs(kappa-1.)>0.2 ) {
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
  
	      MyLimit limit = computeLimit(wspace,data,method);
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
	      result.kappa = kappa;
	      result.sigKappa = sigma_kappa;
	      result.bkg = bkgs[3];
	      result.yield = yields[3];
	      result.upperLimit = limit.upperLimit;
	      result.relUpperLimit = limit.upperLimit/yields[3];

	      if ( results.size()<nrmax ) {
		results.push_back(result);
	      }
	      else {
		std::vector<OptResult>::iterator imax = 
		  std::max_element(results.begin(),results.end());
		*imax = result;
	      }

	      delete data;
	    }
	  }
	}
      }
    }
  }

  std::sort(results.begin(),results.end());

  TFile* out = new TFile("RA4tot.root","RECREATE");
  TH1* hRelUpperLimit = new TH1F("RelUpperLimit","RelUpperLimit",nrmax,0,nrmax);
  TH1* hUpperLimit = new TH1F("UpperLimit","UpperLimit",nrmax,0,nrmax);
  TH1* hKappa = new TH1F("Kappa","Kappa",nrmax,0,nrmax);
  TH1* hSigKappa = new TH1F("sigKappa","sigKappa",nrmax,0,nrmax);
  TH1* hHT0 = new TH1F("HT0","HT0",nrmax,0,nrmax);
  TH1* hHT1 = new TH1F("HT1","HT1",nrmax,0,nrmax);
  TH1* hHT2 = new TH1F("HT2","HT2",nrmax,0,nrmax);
  TH1* hMET0 = new TH1F("MET0","MET0",nrmax,0,nrmax);
  TH1* hMET1 = new TH1F("MET1","MET1",nrmax,0,nrmax);
  TH1* hMET2 = new TH1F("MET2","MET2",nrmax,0,nrmax);
  for ( unsigned int ir=0; ir<results.size(); ++ir ) {
    std::cout << results[ir].relUpperLimit << " "
	      << results[ir].upperLimit << " "
	      << results[ir].bkg << " "
	      << results[ir].yield << " "
	      << results[ir].kappa << " "
	      << results[ir].sigKappa << " ";
    for ( unsigned int i=0; i<3; ++i ) std::cout << results[ir].htBins[i] << " ";
    for ( unsigned int i=0; i<3; ++i ) std::cout << results[ir].metBins[i] << " ";
    std::cout << std::endl;
    hRelUpperLimit->SetBinContent(ir+1,results[ir].relUpperLimit);
    hUpperLimit->SetBinContent(ir+1,results[ir].upperLimit);
    hKappa->SetBinContent(ir+1,results[ir].kappa);
    hSigKappa->SetBinContent(ir+1,results[ir].sigKappa);;
    hHT0->SetBinContent(ir+1,results[ir].htBins[0]);
    hHT1->SetBinContent(ir+1,results[ir].htBins[1]);
    hHT2->SetBinContent(ir+1,results[ir].htBins[2]);
    hMET0->SetBinContent(ir+1,results[ir].metBins[0]);
    hMET1->SetBinContent(ir+1,results[ir].metBins[1]);
    hMET2->SetBinContent(ir+1,results[ir].metBins[2]);
    
  }
  out->Write();
  delete out;
}
