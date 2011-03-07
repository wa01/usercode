#include "RA4abcd.h"

#include "RA4WorkSpace.h"

#include "TStopwatch.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH2F.h"
#include "TGraph.h"
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
// calculation of the limit: assumes that wspace is set up and observations
//   contained in data
//
MyLimit computeLimit (RooWorkspace* wspace, RooDataSet* data, StatMethod method, bool draw) {

  // let's time this challenging example
  TStopwatch t;

  //
  // get nominal signal
  //
  RooRealVar exp_sig(*wspace->var("s"));
  std::cout << "exp_sig = " << exp_sig.getVal() << std::endl;
  
  /////////////////////////////////////////////////////
  // Now the statistical tests
  // model config
  std::cout << wspace->pdf("model") << " "
	    << wspace->pdf("prior") << " "
	    << wspace->set("poi") << " "
	    << wspace->set("nuis") << std::endl;
  ModelConfig modelConfig("RA4abcd");
  modelConfig.SetWorkspace(*wspace);
  modelConfig.SetPdf(*wspace->pdf("model"));
  modelConfig.SetPriorPdf(*wspace->pdf("prior"));
  modelConfig.SetParametersOfInterest(*wspace->set("poi"));
  modelConfig.SetNuisanceParameters(*wspace->set("nuis"));


  //////////////////////////////////////////////////
  // If you want to see the covariance matrix uncomment
  // wspace->pdf("model")->fitTo(*data);

  // use ProfileLikelihood
  if ( method == ProfileLikelihoodMethod ) {
    ProfileLikelihoodCalculator plc(*data, modelConfig);
    plc.SetConfidenceLevel(0.95);
    RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    LikelihoodInterval* plInt = plc.GetInterval();
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    plInt->LowerLimit( *wspace->var("s") ); // get ugly print out of the way. Fix.
    // RooMsgService::instance().setGlobalKillBelow(RooFit::DEBUG);
    if ( draw ) {
      TCanvas* c = new TCanvas("ProfileLikelihood");
      LikelihoodIntervalPlot* lrplot = new LikelihoodIntervalPlot(plInt);
      lrplot->Draw();
    }
//     RooMsgService::instance().setGlobalKillBelow(msglevel);
    double lowLim = plInt->LowerLimit(*wspace->var("s"));
    double uppLim = plInt->UpperLimit(*wspace->var("s"));
    double exp_sig_val = wspace->var("s")->getVal();
    cout << "Profile Likelihood interval on s = [" << 
      lowLim << ", " <<
      uppLim << "]" << endl; 
//     MyLimit result(plInt->IsInInterval(exp_sig),
    MyLimit result(exp_sig_val>lowLim&&exp_sig_val<uppLim,lowLim,uppLim);
    // std::cout << "isIn " << result << std::endl;
    delete plInt;
//     delete modelConfig;
    return result;
  }

  // use FeldmaCousins (takes ~20 min)  
  if ( method == FeldmanCousinsMethod ) {
    FeldmanCousins fc(*data, modelConfig);
    fc.SetConfidenceLevel(0.95);
    //number counting: dataset always has 1 entry with N events observed
    fc.FluctuateNumDataEntries(false); 
    fc.UseAdaptiveSampling(true);
    fc.SetNBins(100);
    PointSetInterval* fcInt = NULL;
    fcInt = (PointSetInterval*) fc.GetInterval(); // fix cast
    double lowLim = fcInt->LowerLimit(*wspace->var("s"));
    double uppLim = fcInt->UpperLimit(*wspace->var("s"));
    double exp_sig_val = wspace->var("s")->getVal();
    cout << "Feldman Cousins interval on s = [" << lowLim << " " << uppLim << endl;
    // std::cout << "isIn " << result << std::endl;
    MyLimit result(exp_sig_val>lowLim&&exp_sig_val<uppLim,
		   fcInt->LowerLimit(*wspace->var("s")),fcInt->UpperLimit(*wspace->var("s")));
    delete fcInt;
    return result;
  }


  // use BayesianCalculator (only 1-d parameter of interest, slow for this problem)  
  if ( method == BayesianMethod ) {
    BayesianCalculator bc(*data, modelConfig);
    bc.SetConfidenceLevel(0.95);
    bc.SetLeftSideTailFraction(0.5);
    SimpleInterval* bInt = NULL;
    if( wspace->set("poi")->getSize() == 1)   {
      bInt = bc.GetInterval();
      if ( draw ) {
	TCanvas* c = new TCanvas("Bayesian");
	// the plot takes a long time and print lots of error
	// using a scan it is better
	bc.SetScanOfPosterior(50);
	RooPlot* bplot = bc.GetPosteriorPlot();
	bplot->Draw();
      }
      cout << "Bayesian interval on s = [" << 
	bInt->LowerLimit( ) << ", " <<
	bInt->UpperLimit( ) << "]" << endl;
      // std::cout << "isIn " << result << std::endl;
      MyLimit result(bInt->IsInInterval(exp_sig),
		     bInt->LowerLimit(),bInt->UpperLimit());
      delete bInt;
      return result;
    } else {
    cout << "Bayesian Calc. only supports on parameter of interest" << endl;
    return MyLimit();
    }
  }


  // use MCMCCalculator  (takes about 1 min)
  // Want an efficient proposal function, so derive it from covariance
  // matrix of fit
  if ( method == MCMCMethod ) {
    RooFitResult* fit = wspace->pdf("model")->fitTo(*data,Save());
    ProposalHelper ph;
    ph.SetVariables((RooArgSet&)fit->floatParsFinal());
    ph.SetCovMatrix(fit->covarianceMatrix());
    ph.SetUpdateProposalParameters(kTRUE); // auto-create mean vars and add mappings
    ph.SetCacheSize(100);
    ProposalFunction* pf = ph.GetProposalFunction();
    
    MCMCCalculator mc(*data, modelConfig);
    mc.SetConfidenceLevel(0.95);
    mc.SetProposalFunction(*pf);
    mc.SetNumBurnInSteps(100); // first N steps to be ignored as burn-in
    mc.SetNumIters(100000);
    mc.SetLeftSideTailFraction(0.5); // make a central interval
    MCMCInterval* mcInt = NULL;
    mcInt = mc.GetInterval();
    MCMCIntervalPlot mcPlot(*mcInt); 
    mcPlot.Draw();
    cout << "MCMC interval on s = [" << 
      mcInt->LowerLimit(*wspace->var("s") ) << ", " <<
      mcInt->UpperLimit(*wspace->var("s") ) << "]" << endl;
    // std::cout << "isIn " << result << std::endl;
    MyLimit result(mcInt->IsInInterval(exp_sig),
		   mcInt->LowerLimit(*wspace->var("s")),mcInt->UpperLimit(*wspace->var("s")));
    delete mcInt;
    return result;
  }
  

  t.Print();

//   delete modelConfig;
  return MyLimit();

}
//
// set value and range for a variable in the workspace
//
void setValRange (RooWorkspace* workspace, const char* name, double val, double vmin, double vmax)
{
  RooRealVar* var = workspace->var(name);
  if ( var ) {
    if ( vmax>vmin )  var->setRange(vmin,vmax);
    var->setVal(val);
  }
}

//
// single measurement (LM0 or LM1)
//
void RA4Single (StatMethod method, double* sig, double* bkg) {

  // RooWorkspace* wspace = createWorkspace();
  RA4WorkSpace ra4WSpace("wspace",true,true,true);
  ra4WSpace.addChannel(RA4WorkSpace::MuChannel);
  ra4WSpace.finalize();

  double lm0_mc[4] = { 1.09, 7.68, 3.78, 21.13 };
  double lm1_mc[4] = { 0.05 , 0.34 , 1.12 , 3.43 };
  double* lm_mc = sig ? sig : lm0_mc;

  double bkg_mc[4] = {  14.73, 18.20, 8.48, 10.98 };

  ra4WSpace.setBackground(RA4WorkSpace::MuChannel,bkg_mc[0],bkg_mc[1],bkg_mc[2],bkg_mc[3]);
  ra4WSpace.setSignal(RA4WorkSpace::MuChannel,lm_mc[0],lm_mc[1],lm_mc[2],lm_mc[3]);
  
  // setBackgrounds(wspace,bkg);
  // setSignal(wspace,lm_mc);

  RooWorkspace* wspace = ra4WSpace.workspace();
  // wspace->Print("v");
  // RooArgSet allVars = wspace->allVars();
  // // allVars.printLatex(std::cout,1);
  // TIterator* it = allVars.createIterator();
  // RooRealVar* var;
  // while ( var=(RooRealVar*)it->Next() ) {
  //   var->Print("v");
  //   var->printValue(std::cout);
  // }

  ////////////////////////////////////////////////////////////
  // Generate toy data
  // generate toy data assuming current value of the parameters
  // import into workspace. 
  // add Verbose() to see how it's being generated
  // wspace->var("s")->setVal(0.);
  // RooDataSet* data =   wspace->pdf("model")->generate(*wspace->set("obs"),1);
  // data->Print("v");
  // // wspace->import(*data);
  // wspace->var("s")->setVal(lm_mc[3]);
  RooDataSet* data = new RooDataSet("data","data",*wspace->set("obs"));
  data->add(*wspace->set("obs"));
  data->Print("v");

  MyLimit limit = computeLimit(wspace,data,method,true);
  std::cout << "Limit [ " << limit.lowerLimit << " , "
	    << limit.upperLimit << " ] ; isIn = " << limit.isInInterval << std::endl;
}
//
// add one channel
//
bool addChannel (const RA4WorkingPoint& workingPoint, RA4WorkSpace::ChannelType type,
		 RA4WorkSpace& ra4WSpace, 
		 TFile** fYield, TFile** fKFactor, unsigned int& nf, 
		 RA4WorkSpace::ChannelType* channelTypes,
		 const RA4WorkingPoint** workingPoints) {
  if ( workingPoint.valid_ ) {
    std::string fname = workingPoint.prefix_ + workingPoint.yields_;
    fYield[nf] = new TFile(fname.c_str());
    if ( fYield[nf]==0 || fYield[nf]->IsZombie() ) {
      std::cout << "failed to open " << fname << std::endl;
      fYield[nf] = 0;
    }
    fname = workingPoint.prefix_ + workingPoint.kfactors_;
    fKFactor[nf] = new TFile(fname.c_str());
    if ( fKFactor[nf]==0 || fKFactor[nf]->IsZombie() ) {
      std::cout << "failed to open " << fname << std::endl;
      fKFactor[nf] = 0;
    }
    workingPoints[nf] = &workingPoint;
    channelTypes[nf] = type;
    if ( fYield[nf] ) {
      ++nf;
      ra4WSpace.addChannel(type);
      return true;
    }
  }
  return false;
}
//
// scan over parameter space
//
void RA4Mult (const RA4WorkingPoint& muChannel,
	      const RA4WorkingPoint& eleChannel,
	      StatMethod method) {

  //
  // Prepare workspace
  //
  RA4WorkSpace ra4WSpace("wspace",true,false,false);


  TFile* fYield[2];
  TFile* fKFactor[2];
  //
  // Muon channel
  //
  unsigned int nf(0);
  RA4WorkSpace::ChannelType channelTypes[2];
  const RA4WorkingPoint* workingPoints[2];
  addChannel(muChannel,RA4WorkSpace::MuChannel,ra4WSpace,fYield,fKFactor,
	     nf,channelTypes,workingPoints);
  addChannel(eleChannel,RA4WorkSpace::EleChannel,ra4WSpace,fYield,fKFactor,
	     nf,channelTypes,workingPoints);
  if ( nf==0 ) {
    std::cout << "No input file" << std::endl;
    return;
  }
  //
  // finish definition of model
  //
  ra4WSpace.finalize();
  RooWorkspace* wspace = ra4WSpace.workspace();
//   wspace->Print("v");
//   RooArgSet allVars = wspace->allVars();
//   // allVars.printLatex(std::cout,1);
//   TIterator* it = allVars.createIterator();
//   RooRealVar* var;
//   while ( var=(RooRealVar*)it->Next() ) {
//     var->Print("v");
//     var->printValue(std::cout);
//   }
  //
  // preparation of histograms with yields and k-factors
  //
  const char* cRegion = { "ABCD" };
  TH2* hYields[4][2];
  TH2* hYields05[4][2];
  TH2* hYields20[4][2];
  TH2* hYEntries[4][2];
  for ( unsigned int j=0; j<nf; ++j ) {
    for ( unsigned int i=0; i<4; ++i ) {
      hYields[i][j] = 0;
      hYields05[i][j] = 0;
      hYields20[i][j] = 0;
      hYEntries[i][j] = 0;
    }
  }
  TH2* hKF05[2];
  TH2* hKF10[2];
  TH2* hKF20[2];
  for ( unsigned int j=0; j<nf; ++j ) {
    hKF05[j] = 0;
    hKF10[j] = 0;
    hKF20[j] = 0;
  }
  //
  // Retrieval of histograms with k-factors
  //
  for ( unsigned int j=0; j<nf; ++j ) {
    hKF05[j] = (TH2*)fKFactor[j]->Get("hKF05D");
    hKF10[j] = (TH2*)fKFactor[j]->Get("hKF10D");
    hKF20[j] = (TH2*)fKFactor[j]->Get("hKF20D");
    if ( hKF05[j]==0 || hKF10==0 || hKF20==0 ) {
      std::cout << "Missing histogram for kfactor for channel " << j << std::endl;
      return;
    }
  }
  //
  // Retrieval of histograms with yields
  //
  std::string hName;
  for ( unsigned int j=0; j<nf; ++j ) {
    for ( unsigned int i=0; i<4; ++i ) {
      hName = "Events";
      hName += cRegion[i];

      hYields[i][j] = (TH2*)fYield[j]->Get(hName.c_str())->Clone();
      hYields05[i][j] = (TH2*)fYield[j]->Get(hName.c_str())->Clone();
      hYields20[i][j] = (TH2*)fYield[j]->Get(hName.c_str())->Clone();
      if ( hYields[i][j]==0 ) {
	std::cout << "Missing histogram for region " << cRegion[i] << std::endl;
	return;
      }
      hYields[i][j]->Multiply(hYields[i][j],hKF10[j]);
      hYields05[i][j]->Multiply(hYields05[i][j],hKF05[j]);
      hYields20[i][j]->Multiply(hYields20[i][j],hKF20[j]);

      hName = "Entries";
      hName += cRegion[i];
      hYEntries[i][j] = (TH2*)fYield[j]->Get(hName.c_str());
      if ( hYEntries[i][j]==0 ) {
	std::cout << "Missing histogram for region " << cRegion[i] << std::endl;
	return;
      }
      // convert to efficiency (assume 10000 MC events/bin)
      hYEntries[i][j]->Scale(1/10000.);
//       // convert yield to cross section
//       hYields[i][j]->Divide(hYields[i][j],hYEntries[i][j]);
    }
  }
  //
  // histograms with exclusion and limits
  //
  gROOT->cd();
  TH2* hExclusion = (TH2*)hYields[0][0]->Clone("Exclusion");
  hExclusion->Reset();
  hExclusion->SetTitle("Exclusion");
  TH2* hLowerLimit = (TH2*)hYields[0][0]->Clone("LowerLimit");
  hLowerLimit->Reset();
  hLowerLimit->SetTitle("LowerLimit");
  TH2* hUpperLimit = (TH2*)hYields[0][0]->Clone("UpperLimit");
  hUpperLimit->Reset();
  hUpperLimit->SetTitle("UpperLimit");

  double yields[4][2];
  double entries[4][2];

//   double bkgs[4][2];

//   double kappa = (bkgs[0]*bkgs[3])/(bkgs[1]*bkgs[2]);
//   double sigma_kappa_base = 0.10;
//   double delta_kappa_abs = kappa - 1.;
//   double sigma_kappa = sqrt(sigma_kappa_base*sigma_kappa_base+delta_kappa_abs*delta_kappa_abs);
//   sigma_kappa = sqrt(0.129*0.129+0.1*0.1);

//   int nbx = hYields[0][0]->GetNbinsX();
//   int nby = hYields[0][0]->GetNbinsY();
//   for ( int ix=1; ix<=nbx; ++ix ) {
//     for ( int iy=1; iy<=nby; ++iy ) {
  { 
    int ix=40;
    {
      int iy=11;

      bool process(false);
      for ( unsigned int j=0; j<nf; ++j ) {
	ra4WSpace.setBackground(channelTypes[j],
				workingPoints[j]->bkg_[0],workingPoints[j]->bkg_[1],
				workingPoints[j]->bkg_[2],workingPoints[j]->bkg_[3]);
	ra4WSpace.setObserved(channelTypes[j],
			      workingPoints[j]->obs_[0],workingPoints[j]->obs_[1],
			      workingPoints[j]->obs_[2],workingPoints[j]->obs_[3]);
	for ( unsigned int i=0; i<4; ++i ) {
	  yields[i][j] = hYields[i][j]->GetBinContent(ix,iy);
	  entries[i][j] = hYEntries[i][j]->GetBinContent(ix,iy);
	  if ( yields[3][j]>0.01 && entries[i][j]>0.0001 )  process = true;
	}
	double yields05 = hYields05[3][j]->GetBinContent(ix,iy);
	double yields20 = hYields20[3][j]->GetBinContent(ix,iy);
	ra4WSpace.setSignal(channelTypes[j],
			    yields[0][j],yields[1][j],
			    yields[2][j],yields[3][j]);
	std::cout << "yields for channel " << j << " =";
	for ( unsigned int i=0; i<4; ++i )
	  std::cout << " " << yields[i][j];
	std::cout << endl;
	std::cout << "backgrounds for channel " << j << " =";
	for ( unsigned int i=0; i<4; ++i )
	  std::cout << " " << workingPoints[j]->bkg_[i];
	std::cout << endl;
      }

      double sigK(0.);
      for ( unsigned int j=0; j<nf; ++j ) {
	if ( workingPoints[j]->sigKappa_>sigK )
	  sigK = workingPoints[j]->sigKappa_;
      }
      wspace->var("sigmaKappa")->setVal(sigK);
//       wspace->var("sigmaKappa")->setVal(sqrt(0.129*0.129+0.1*0.1)*0.967);
      // for the time being: work with yields
      if ( muChannel.valid_ ) {
	wspace->var("effM")->setVal(1.);
// 	wspace->var("sadM")->setVal(0.);
// 	wspace->var("sbdM")->setVal(0.);
// 	wspace->var("scdM")->setVal(0.);
      }
      if ( eleChannel.valid_ ) {
	wspace->var("effE")->setVal(1.);
// 	wspace->var("sadE")->setVal(0.);
// 	wspace->var("sbdE")->setVal(0.);
// 	wspace->var("scdE")->setVal(0.);
      }

      wspace->Print("v");
      RooArgSet allVars = wspace->allVars();
      // allVars.printLatex(std::cout,1);
      TIterator* it = allVars.createIterator();
      RooRealVar* var;
      while ( var=(RooRealVar*)it->Next() ) {
	var->Print("v");
	var->printValue(std::cout);
	std::cout << std::endl;
      }
      
//       yields[0] =1.52;
//       yields[1] =7.68;
//       yields[2] =5.17;
//       yields[3] =21.12;
//       // *1.3 for NLO
//       for ( unsigned int i=0; i<4; ++i )  yields[i] *= 1.3;

      MyLimit limit(true,0.,999.);
      std::cout << "Checked ( " << hExclusion->GetXaxis()->GetBinCenter(ix) << " , "
		<< hExclusion->GetYaxis()->GetBinCenter(iy) << " ) with signal yield " 
		<< yields[3] << std::endl;
	

      RooDataSet data("data","data",*wspace->set("obs"));
      data.add(*wspace->set("obs"));
      data.Print("v");
      
//       if ( yields[3]>0.01 ) {
      limit = computeLimit(wspace,&data,method);
      std::cout << "  Limit [ " << limit.lowerLimit << " , "
		<< limit.upperLimit << " ] ; isIn = " << limit.isInInterval << std::endl;
      

//       std::cout << "  entries =" 
// 		<< " " << entries[0]
// 		<< " " << entries[1]
// 		<< " " << entries[2]
// 		<< " " << entries[3] << std::endl;
      double excl = limit.isInInterval;
      if ( limit.upperLimit<limit.lowerLimit )  excl = -1;
      hExclusion->SetBinContent(ix,iy,excl);
      hLowerLimit->SetBinContent(ix,iy,limit.lowerLimit);
      hUpperLimit->SetBinContent(ix,iy,limit.upperLimit);
      return;
      
    }
  }
  
  TFile* out = new TFile("RA4abcd.root","RECREATE");
  hExclusion->SetDirectory(out);
  hExclusion->SetMinimum(); hExclusion->SetMaximum();
  hExclusion->SetContour(1); hExclusion->SetContourLevel(0,0.5);
  hLowerLimit->SetDirectory(out);
  hLowerLimit->SetMinimum(); hLowerLimit->SetMaximum();
  hUpperLimit->SetDirectory(out);
  hUpperLimit->SetMinimum(); hUpperLimit->SetMaximum();
  for ( unsigned int j=0; j<nf; ++j ) {
    hYields[3][j]->SetDirectory(out);
    hYields[3][j]->SetMinimum(); hYields[3][j]->SetMaximum();
  }
  out->Write();
  delete out;
}

