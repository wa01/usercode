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

/** Some functions for limit calculation with ABCD.
 *  RA4Single: one test with LM0
 *  RA4Mult:   scan over m0-m1/2 (signal yields in ABCD provided by histograms) */

using namespace RooFit;
using namespace RooStats;

// stat method used
enum StatMethod { ProfileLikelihoodMethod, FeldmanCousinsMethod, BayesianMethod, MCMCMethod };

// simple class to hold the calculated limit
struct MyLimit {
  MyLimit () : isInInterval(false), lowerLimit(999999), upperLimit(-999999) {}
  MyLimit (bool inInterval, double lower, double upper) : 
    isInInterval(inInterval), lowerLimit(lower), upperLimit(upper) {}
  bool isInInterval;
  double lowerLimit;
  double upperLimit;
};

//
// calculation of the limit: assumes that wspace is set up and observations
//   contained in data
//
MyLimit computeLimit (RooWorkspace* wspace, RooDataSet* data, StatMethod method){
  
  // let's time this challenging example
  TStopwatch t;
  // for 2-d plots to inspect correlations:
  //  wspace->defineSet("poi","s,kappa");

  // // test simpler cases where parameters are known.
  // wspace->var("sad")->setConstant(1);
  // wspace->var("sbd")->setConstant(1);
  // wspace->var("scd")->setConstant(1);

  // inspect workspace
  //  wspace->Print();

  //
  // get nominal signal
  //
  RooRealVar exp_sig(*wspace->var("s"));
  // std::cout << "exp_sig = " << exp_sig.getVal() << std::endl;
  
  /////////////////////////////////////////////////////
  // Now the statistical tests
  // model config
  ModelConfig* modelConfig = new ModelConfig("RA4abcd");
  modelConfig->SetWorkspace(*wspace);
  modelConfig->SetPdf(*wspace->pdf("model"));
  modelConfig->SetPriorPdf(*wspace->pdf("prior"));
  modelConfig->SetParametersOfInterest(*wspace->set("poi"));
  modelConfig->SetNuisanceParameters(*wspace->set("nuis"));

  //////////////////////////////////////////////////
  // If you want to see the covariance matrix uncomment
  // wspace->pdf("model")->fitTo(*data);

  // use ProfileLikelihood
  if ( method == ProfileLikelihoodMethod ) {
    ProfileLikelihoodCalculator plc(*data, *modelConfig);
    plc.SetConfidenceLevel(0.95);
    RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    LikelihoodInterval* plInt = plc.GetInterval();
    RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
    plInt->LowerLimit( *wspace->var("s") ); // get ugly print out of the way. Fix.
    // RooMsgService::instance().setGlobalKillBelow(RooFit::DEBUG);
    TCanvas* c = new TCanvas("ProfileLikelihood");
    LikelihoodIntervalPlot* lrplot = new LikelihoodIntervalPlot(plInt);
    lrplot->Draw();
    RooMsgService::instance().setGlobalKillBelow(msglevel);
    cout << "Profile Likelihood interval on s = [" << 
      plInt->LowerLimit( *wspace->var("s") ) << ", " <<
      plInt->UpperLimit( *wspace->var("s") ) << "]" << endl; 
    MyLimit result(plInt->IsInInterval(exp_sig),
		   plInt->LowerLimit(*wspace->var("s")),plInt->UpperLimit(*wspace->var("s")));
    // std::cout << "isIn " << result << std::endl;
    delete plInt;
    return result;
  }

  // use FeldmaCousins (takes ~20 min)  
  if ( method == FeldmanCousinsMethod ) {
    FeldmanCousins fc(*data, *modelConfig);
    fc.SetConfidenceLevel(0.95);
    //number counting: dataset always has 1 entry with N events observed
    fc.FluctuateNumDataEntries(false); 
    fc.UseAdaptiveSampling(true);
    fc.SetNBins(100);
    PointSetInterval* fcInt = NULL;
    fcInt = (PointSetInterval*) fc.GetInterval(); // fix cast
    cout << "Feldman Cousins interval on s = [" << 
      fcInt->LowerLimit( *wspace->var("s") ) << ", " <<
      fcInt->UpperLimit( *wspace->var("s") ) << "]" << endl;
    // std::cout << "isIn " << result << std::endl;
    MyLimit result(fcInt->IsInInterval(exp_sig),
		   fcInt->LowerLimit(*wspace->var("s")),fcInt->UpperLimit(*wspace->var("s")));
    delete fcInt;
    return result;
  }


  // use BayesianCalculator (only 1-d parameter of interest, slow for this problem)  
  if ( method == BayesianMethod ) {
    BayesianCalculator bc(*data, *modelConfig);
    bc.SetConfidenceLevel(0.95);
    bc.SetLeftSideTailFraction(0.5);
    SimpleInterval* bInt = NULL;
    if( wspace->set("poi")->getSize() == 1)   {
      bInt = bc.GetInterval();
      TCanvas* c = new TCanvas("Bayesian");
      // the plot takes a long time and print lots of error
      // using a scan it is better
      bc.SetScanOfPosterior(20);
      RooPlot* bplot = bc.GetPosteriorPlot();
      bplot->Draw();
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
  if ( method = MCMCMethod ) {
    RooFitResult* fit = wspace->pdf("model")->fitTo(*data,Save());
    ProposalHelper ph;
    ph.SetVariables((RooArgSet&)fit->floatParsFinal());
    ph.SetCovMatrix(fit->covarianceMatrix());
    ph.SetUpdateProposalParameters(kTRUE); // auto-create mean vars and add mappings
    ph.SetCacheSize(100);
    ProposalFunction* pf = ph.GetProposalFunction();
    
    MCMCCalculator mc(*data, *modelConfig);
    mc.SetConfidenceLevel(0.95);
    mc.SetProposalFunction(*pf);
    mc.SetNumBurnInSteps(500); // first N steps to be ignored as burn-in
    mc.SetNumIters(50000);
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

  delete modelConfig;

}
//
// set value and range for a variable in the workspace
//
void setValRange (RooWorkspace* workspace, const char* name, double val, double vmin=999999, double vmax=-999999)
{
  RooRealVar* var = workspace->var(name);
  if ( var ) {
    if ( vmax>vmin )  var->setRange(vmin,vmax);
    var->setVal(val);
  }
}
//
// create a workspace with the variables (with dummy values) and the pdfs
//
RooWorkspace* createWorkspace (const char* name = "wspace")
{
  RooWorkspace* wspace = new RooWorkspace(name);
  // observed event counts
  wspace->factory("na[0,0,1000]");
  wspace->factory("nb[0,0,1000]");
  wspace->factory("nc[0,0,1000]");
  wspace->factory("nd[0,0,1000]");
  // signal in D and rel. signal contamination
  wspace->factory("s[1,0,100]");
  wspace->factory("sad[0,0,10]");
  wspace->factory("sbd[0,0,10]");
  wspace->factory("scd[0,0,10]");
  // bkg in A; relative bkg in B&C; kappa
  wspace->factory("bkga[1,0,1000]");
  wspace->factory("bba[0,0,10]");
  wspace->factory("bca[0,0,10]");
  wspace->factory("kappa[1,0,2]");
  // pseudo-measurements for kappa and signal contamination
  wspace->factory("kappanom[0,0,10]");
  wspace->factory("sadnom[0.1,0,10]");
  wspace->factory("sbdnom[0.1,0,10]");
  wspace->factory("scdnom[0.1,0,10]");
  wspace->factory("sigmaKappa[0.1]");
  // uncertainties on pseudo-measurements
  wspace->factory("sigmaSad[0.1]");
  wspace->factory("sigmaSbd[0.1]");
  wspace->factory("sigmaScd[0.1]");
  
  // Poisson distributions in the 4 regions
  wspace->factory("Poisson::a(na, sum::tota(prod::sa(s,sad),bkga))");
  wspace->factory("Poisson::b(nb, sum::totb(prod::sb(s,sbd),prod::bkgb(bkga,bba)))");
  wspace->factory("Poisson::c(nc, sum::totc(prod::sc(s,scd),prod::bkgc(bkga,bca)))");
  wspace->factory("Poisson::d(nd, sum::splusb(s,prod::bkgd(bkga,bba,bca,kappa)))");
  // Pdfs for pseudo-measurements
  wspace->factory("Gaussian::mcKappa(kappanom, kappa, sigmaKappa)");
  wspace->factory("Gaussian::mcSad(sadnom, sad, sigmaSad)");
  wspace->factory("Gaussian::mcSbd(sbdnom, sbd, sigmaSbd)");
  wspace->factory("Gaussian::mcScd(scdnom, scd, sigmaScd)");
  // full model
  wspace->factory("PROD::model(d,c,b,a,mcKappa,mcSad,mcSbd,mcScd)");
  // priors
  wspace->factory("Uniform::prior_poi({s})");
  wspace->factory("Uniform::prior_nuis({bkga,bba,bca,kappa,sad,sbd,scd})");
  wspace->factory("PROD::prior(prior_poi,prior_nuis)"); 
  // sets (observables, POI, nuisance parameters)
  wspace->defineSet("obs","nd,nc,nb,na,kappanom,sadnom,sbdnom,scdnom");
  wspace->defineSet("poi","s");
  wspace->defineSet("nuis","bkga,bba,bca,kappa,sad,sbd,scd");

  return wspace;
}
//
// set background-related variables of the workspace
//
void setBackgrounds (RooWorkspace* wspace) 
{
  // Inputs : expected values
  // double bkg_mc[4] = { 53.26 , 19.48 , 56.13 , 20.22 };
  double bkg_mc[4] = { 43.76 , 38.88 , 23.77 , 22.13 };
  double observed[4];
  for ( unsigned int i=0; i<4; ++i ) 
    observed[i] = int(bkg_mc[i]+0.5);

  // scaling factors
  double sigma_kappa = 0.25;

  // derived quantities
  double bba_mc = bkg_mc[1]/bkg_mc[0];
  double bca_mc = bkg_mc[2]/bkg_mc[0];
  double kappa_mc = (bkg_mc[3]*bkg_mc[0])/(bkg_mc[1]*bkg_mc[2]);

  // Roo variables
  // .. inputs set to expected backgrounds
  setValRange(wspace,"na",observed[0],0,1000);
  setValRange(wspace,"nb",observed[1],0,1000);
  setValRange(wspace,"nc",observed[2],0,1000);
  setValRange(wspace,"nd",observed[3],0,1000);
  // .. pseudo-measurements
  setValRange(wspace,"kappanom",kappa_mc,kappa_mc/10,kappa_mc*10);
  // .. uncertainties on correlation and signal contamination
  setValRange(wspace,"sigmaKappa",sigma_kappa);

  // .. background variables
  setValRange(wspace,"bkga",bkg_mc[0],0,10*bkg_mc[3]);
  setValRange(wspace,"bba",bba_mc,bba_mc/10,bba_mc*10);
  setValRange(wspace,"bca",bca_mc,bca_mc/10,bca_mc*10);

  // .. correlation and signal contamination variables
  setValRange(wspace,"kappa",kappa_mc,kappa_mc/10,kappa_mc*10);
}

void setSignal (RooWorkspace* wspace, double* lm_mc)
{
  // Inputs : expected values

  // relative uncertainties on relative signal contamination
  double sigma_sad_rel = 0.00100;
  double sigma_sbd_rel = 0.00100;
  double sigma_scd_rel = 0.00100;
  // double sigma_sad_rel = 0.0030;
  // double sigma_sbd_rel = 0.0030;
  // double sigma_scd_rel = 0.0030;

  // derived quantities
  double sad_mc = max(lm_mc[0]/lm_mc[3],0.01);
  double sbd_mc = max(lm_mc[1]/lm_mc[3],0.01);
  double scd_mc = max(lm_mc[2]/lm_mc[3],0.01);

  // Roo variables
  // .. pseudo-measurements
  setValRange(wspace,"sadnom",sad_mc,sad_mc/10,sad_mc*10);
  setValRange(wspace,"sbdnom",sbd_mc,sbd_mc/10,sbd_mc*10);
  setValRange(wspace,"scdnom",scd_mc,scd_mc/10,scd_mc*10);
  // .. uncertainties on signal contamination
  setValRange(wspace,"sigmaSad",sad_mc*sigma_sad_rel);
  setValRange(wspace,"sigmaSbd",sbd_mc*sigma_sbd_rel);
  setValRange(wspace,"sigmaScd",scd_mc*sigma_scd_rel);

  // .. background and signal variables
  setValRange(wspace,"s",lm_mc[3],0,10*lm_mc[3]);
  // wspace->var("s")->Print("v");

  // .. correlation and signal contamination variables
  setValRange(wspace,"sad",sad_mc,sad_mc/10,sad_mc*10);
  setValRange(wspace,"sbd",sbd_mc,sbd_mc/10,sbd_mc*10);
  setValRange(wspace,"scd",scd_mc,scd_mc/10,scd_mc*10);

}
//
// single measurement (LM0 or LM1)
//
void RA4Single (StatMethod method) {

  RooWorkspace* wspace = createWorkspace();

  double lm0_mc[4] = { 3.58 , 8.21 , 16.37 , 25.21 };
  double lm1_mc[4] = { 0.05 , 0.34 , 1.12 , 3.43 };

  setBackgrounds(wspace);
  setSignal(wspace,lm0_mc);

  wspace->Print("v");
  // RooArgSet allVars = wspace->allVars();
  // allVars.printLatex(std::cout,1);

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
  
  MyLimit limit = computeLimit(wspace,data,method);
  std::cout << "Limit [ " << limit.lowerLimit << " , "
	    << limit.upperLimit << " ] ; isIn = " << limit.isInInterval << std::endl;
}
//
// scan over parameter space
//
void RA4Mult (const char* file, StatMethod method=ProfileLikelihoodMethod) {

  TFile* fYield = new TFile(file);
  if ( fYield==0 || fYield->IsZombie() ) {
    std::cout << "failed to open " << file << std::endl;
    return;
  }
  const char* cRegion = { "ABCD" };
  TH2* hYields[4];
  TH2* hYEntries[4];

  std::string hName;
  for ( unsigned int i=0; i<4; ++i ) {
    hName = "Events";
    hName += cRegion[i];
    hYields[i] = (TH2*)fYield->Get(hName.c_str());
    hName = "Entries";
    hName += cRegion[i];
    hYEntries[i] = (TH2*)fYield->Get(hName.c_str());
    if ( hYields[i]==0 || hYEntries[i]==0 ) {
      std::cout << "Missing histogram for region " << cRegion[i] << std::endl;
      return;
    }
  }

  gROOT->cd();
  TH2* hExclusion = (TH2*)hYields[0]->Clone("Exclusion");
  hExclusion->Reset();
  hExclusion->SetTitle("Exclusion");
  TH2* hLowerLimit = (TH2*)hYields[0]->Clone("LowerLimit");
  hLowerLimit->Reset();
  hLowerLimit->SetTitle("LowerLimit");
  TH2* hUpperLimit = (TH2*)hYields[0]->Clone("UpperLimit");
  hUpperLimit->Reset();
  hUpperLimit->SetTitle("UpperLimit");

  RooWorkspace* wspace = createWorkspace();

  double yields[4];
  double entries[4];

  int nbx = hYields[0]->GetNbinsX();
  int nby = hYields[0]->GetNbinsY();
  for ( unsigned int ix=1; ix<=nbx; ++ix ) {
    for ( unsigned int iy=1; iy<=nby; ++iy ) {

      for ( unsigned int i=0; i<4; ++i ) {
	yields[i] = hYields[i]->GetBinContent(ix,iy);
	entries[i] = hYEntries[i]->GetBinContent(ix,iy);
      }
      // yields[0] =3.58;
      // yields[1] =8.21;
      // yields[2] =16.37;
      // yields[3] =25.21;

      setBackgrounds(wspace);
      setSignal(wspace,yields);

      // wspace->Print("v");
      // RooArgSet allVars = wspace->allVars();
      // allVars.printLatex(std::cout,1);

      RooDataSet* data = new RooDataSet("data","data",*wspace->set("obs"));
      data->add(*wspace->set("obs"));
  // data->Print("v");
  
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
      std::cout << "  entries =" 
		<< " " << entries[0]
		<< " " << entries[1]
		<< " " << entries[2]
		<< " " << entries[3] << std::endl;
      double excl = limit.isInInterval;
      if ( limit.upperLimit<limit.lowerLimit )  excl = -1;
      hExclusion->SetBinContent(ix,iy,excl);
      hLowerLimit->SetBinContent(ix,iy,limit.lowerLimit);
      hUpperLimit->SetBinContent(ix,iy,limit.upperLimit);

      delete data;

    }
  }

  TFile* out = new TFile("RA4abcd.root","RECREATE");
  hExclusion->SetDirectory(out);
  hExclusion->SetMinimum(); hExclusion->SetMaximum();
  hLowerLimit->SetDirectory(out);
  hLowerLimit->SetMinimum(); hLowerLimit->SetMaximum();
  hUpperLimit->SetDirectory(out);
  hUpperLimit->SetMinimum(); hUpperLimit->SetMaximum();
  hYields[3]->SetDirectory(out);
  hYields[3]->SetMinimum(); hYields[3]->SetMaximum();
  out->Write();
  delete out;
}
