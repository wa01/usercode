#include "RA4abcdef.h"

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
// calculation of the limit: assumes that wspace is set up and observations
//   contained in data
//
MyLimit computeLimit (RooWorkspace* wspace, RooDataSet* data, StatMethod method, bool draw) {

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
  ModelConfig modelConfig("RA4abcdef");
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
    RooMsgService::instance().setGlobalKillBelow(msglevel);
    cout << "Profile Likelihood interval on s = [" << 
      plInt->LowerLimit( *wspace->var("s") ) << ", " <<
      plInt->UpperLimit( *wspace->var("s") ) << "]" << endl; 
    MyLimit result(plInt->IsInInterval(exp_sig),
		   plInt->LowerLimit(*wspace->var("s")),plInt->UpperLimit(*wspace->var("s")));
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
// create a workspace with the variables (with dummy values) and the pdfs
//
RooWorkspace* createWorkspace (const char* name)
{
  RooWorkspace* wspace = new RooWorkspace(name);
  // observed event counts
  wspace->factory("na[0,0,1000]");
  wspace->factory("nb[0,0,1000]");
  wspace->factory("nc[0,0,1000]");
  wspace->factory("nd[0,0,1000]");
  wspace->factory("ne[0,0,1000]");
  wspace->factory("nf[0,0,1000]");
  // signal in D and rel. signal contamination
  wspace->factory("s[1,0,100]");
  wspace->factory("saf[0,0,10]");
  wspace->factory("sbf[0,0,10]");
  wspace->factory("scf[0,0,10]");
  wspace->factory("sdf[0,0,10]");
  wspace->factory("sef[0,0,10]");
//   wspace->factory("eff[1.,0.1,2.]");
  wspace->factory("eff[1.]");
  // bkg in A; relative bkg in B&C; kappa
  wspace->factory("bkgf[1,0,1000]");
  wspace->factory("bef[0,0,10]");
  wspace->factory("bdf[0,0,10]");
  wspace->factory("bcf[0,0,10]");
  wspace->factory("kappa[1,0,2]");
  // pseudo-measurements for kappa and signal contamination
  wspace->factory("kappanom[1]");
  wspace->factory("safnom[0.1]");
  wspace->factory("sbfnom[0.1]");
  wspace->factory("scfnom[0.1]");
  wspace->factory("sdfnom[0.1]");
  wspace->factory("sefnom[0.1]");
  wspace->factory("effnom[1.]");
  wspace->factory("sigmaKappa[0.1]");
  wspace->factory("sigmaSaf[0.1]");
  wspace->factory("sigmaSbf[0.1]");
  wspace->factory("sigmaScf[0.1]");
  wspace->factory("sigmaSdf[0.1]");
  wspace->factory("sigmaSef[0.1]");
  wspace->factory("sigmaEff[0.1]");
  // Poisson distributions in the 4 regions
  wspace->factory("prod::sf(s,eff)");
  wspace->factory("Poisson::a(na, sum::tota(prod::sa(sf,saf),prod::bkga(bkgf,bdf,bcf,kappa)))");
  wspace->factory("Poisson::b(nb, sum::totb(prod::sb(sf,sbf),prod::bkgb(bkgf,bef,bcf,kappa)))");
  wspace->factory("Poisson::c(nc, sum::totc(prod::sc(sf,scf),prod::bkgc(bkgf,bcf)))");
  wspace->factory("Poisson::d(nd, sum::totd(prod::sd(sf,sdf),prod::bkgd(bkgf,bdf)))");
  wspace->factory("Poisson::e(ne, sum::tote(prod::se(sf,sef),prod::bkge(bkgf,bef)))");
  wspace->factory("Poisson::f(nf, sum::splusb(sf,bkgf))");
  // Pdfs for pseudo-measurements
  wspace->factory("Gaussian::mcKappa(kappanom, kappa, sigmaKappa)");
  wspace->factory("Gaussian::mcSaf(safnom, saf, sigmaSaf)");
  wspace->factory("Gaussian::mcSbf(sbfnom, sbf, sigmaSbf)");
  wspace->factory("Gaussian::mcScf(scfnom, scf, sigmaScf)");
  wspace->factory("Gaussian::mcSdf(sdfnom, sdf, sigmaSdf)");
  wspace->factory("Gaussian::mcSef(sefnom, sef, sigmaSef)");
//   wspace->factory("Gaussian::mcEff(effnom, eff, sigmaEff)");
  // full model
//   wspace->factory("PROD::model(d,c,b,a,mcKappa,mcSaf,mcSbf,mcScf,mcEff)");
  wspace->factory("PROD::model(f,e,d,c,b,a,mcKappa,mcSaf,mcSbf,mcScf,mcSdf,mcSef)");
  // priors
  wspace->factory("Uniform::prior_poi({s})");
//   wspace->factory("Uniform::prior_nuis({bkgf,bcf,bcd,kappa,saf,sbf,scf,eff})");
  wspace->factory("Uniform::prior_nuis({bkgf,bcf,bdf,bef,kappa,saf,sbf,scf,sdf,sef})");
  wspace->factory("PROD::prior(prior_poi,prior_nuis)"); 
  // sets (observables, POI, nuisance parameters)
//   wspace->defineSet("obs","nd,nc,nb,na,kappanom,safnom,sbfnom,scfnom,effnom");
  wspace->defineSet("obs","nf,ne,nd,nc,nb,na,kappanom,safnom,sbfnom,scfnom,sdfnom,sefnom");
  wspace->defineSet("poi","s");
//   wspace->defineSet("nuis","bkgf,bcf,bcd,kappa,saf,sbf,scf,eff");
  wspace->defineSet("nuis","bkgf,bcf,bdf,bef,kappa,saf,sbf,scf,sdf,sef");

  return wspace;
}
//
// set background-related variables of the workspace
//
void setBackgrounds (RooWorkspace* wspace, double* bkgs) 
{
  // Inputs : expected values
  // double bkg_mc[4] = { 120. , 12.1 , 18.3 , 1.83 }; // tight settings / HT2
  double bkg_mc[6] = { 14.73, 15.29, 18.20, 8.48, 8.91, 10.98 };
  if ( bkgs ) {
    for ( unsigned int i=0; i<6; ++i )  bkg_mc[i] = bkgs[i];
  }
  //
  // use pessimistic scenario for rounding of expected numbers
  //  
  double observed[6];
//   for ( unsigned int i=0; i<6; ++i ) 
//     observed[i] = int(bkg_mc[i]+0.5);
  observed[0] = int(bkg_mc[0])+1;
  observed[1] = int(bkg_mc[1])+1;
  observed[2] = int(bkg_mc[2]);
  observed[3] = int(bkg_mc[3]);
  observed[4] = int(bkg_mc[4]);
  observed[5] = int(bkg_mc[5])+1;

  // scaling factors
  double sigma_kappa = 0.25;

  // derived quantities
  double bcf_mc = bkg_mc[2]/bkg_mc[5];
  double bdf_mc = bkg_mc[3]/bkg_mc[5];
  double bef_mc = bkg_mc[4]/bkg_mc[5];
  double kappa_mc = (bkg_mc[5]*bkg_mc[0])/(bkg_mc[3]*bkg_mc[2]);

  // Roo variables
  // .. inputs set to expected backgrounds
  setValRange(wspace,"na",observed[0],0,1000);
  setValRange(wspace,"nb",observed[1],0,1000);
  setValRange(wspace,"nc",observed[2],0,1000);
  setValRange(wspace,"nd",observed[3],0,1000);
  setValRange(wspace,"ne",observed[4],0,1000);
  setValRange(wspace,"nf",observed[5],0,1000);
  // .. pseudo-measurements
//   setValRange(wspace,"kappanom",kappa_mc,kappa_mc/10,kappa_mc*10);
  setValRange(wspace,"kappanom",kappa_mc);
  // .. uncertainties on correlation and signal contamination
  setValRange(wspace,"sigmaKappa",sigma_kappa);

  // .. background variables
  setValRange(wspace,"bkgf",bkg_mc[5],0,10*bkg_mc[5]);
  setValRange(wspace,"bcf",bcf_mc,0,bcf_mc*2);
  setValRange(wspace,"bdf",bdf_mc,0,bdf_mc*2);
  setValRange(wspace,"bef",bef_mc,0,bef_mc*2);

  // .. correlation and signal contamination variables
  setValRange(wspace,"kappa",kappa_mc,0,2);
}

void setSignal (RooWorkspace* wspace, double* lm_mc)
{
  // Inputs : expected values

  // relative uncertainties on relative signal contamination
  double sigma_saf_rel = 0.100;
  double sigma_sbf_rel = 0.100;
  double sigma_scf_rel = 0.100;
  double sigma_sdf_rel = 0.100;
  double sigma_sef_rel = 0.100;
  double sigma_eff = 0.001;
  // double sigma_saf_rel = 0.0030;
  // double sigma_sbf_rel = 0.0030;
  // double sigma_scf_rel = 0.0030;

  // derived quantities
  double saf_mc = max(lm_mc[0]/lm_mc[5],0.01);
  double sbf_mc = max(lm_mc[1]/lm_mc[5],0.01);
  double scf_mc = max(lm_mc[2]/lm_mc[5],0.01);
  double sdf_mc = max(lm_mc[3]/lm_mc[5],0.01);
  double sef_mc = max(lm_mc[4]/lm_mc[5],0.01);
  double eff_mc = 1.;

  // Roo variables
  // .. pseudo-measurements
//   setValRange(wspace,"safnom",saf_mc,saf_mc/10,saf_mc*10);
//   setValRange(wspace,"sbfnom",sbf_mc,sbf_mc/10,sbf_mc*10);
//   setValRange(wspace,"scfnom",scf_mc,scf_mc/10,scf_mc*10);
  setValRange(wspace,"safnom",saf_mc);
  setValRange(wspace,"sbfnom",sbf_mc);
  setValRange(wspace,"scfnom",scf_mc);
  setValRange(wspace,"sdfnom",sdf_mc);
  setValRange(wspace,"sefnom",sef_mc);
  // .. uncertainties on signal contamination
  setValRange(wspace,"sigmaSaf",saf_mc*sigma_saf_rel);
  setValRange(wspace,"sigmaSbf",sbf_mc*sigma_sbf_rel);
  setValRange(wspace,"sigmaScf",scf_mc*sigma_scf_rel);
  setValRange(wspace,"sigmaSdf",sdf_mc*sigma_sdf_rel);
  setValRange(wspace,"sigmaSef",sef_mc*sigma_sef_rel);
  setValRange(wspace,"sigmaEff",sigma_eff);

  // .. background and signal variables
  setValRange(wspace,"s",lm_mc[5],0,10*lm_mc[5]);
  // wspace->var("s")->Print("v");

  // .. correlation and signal contamination variables
  setValRange(wspace,"saf",saf_mc,0,saf_mc*3);
  setValRange(wspace,"sbf",sbf_mc,0,sbf_mc*3);
  setValRange(wspace,"scf",scf_mc,0,scf_mc*3);
  setValRange(wspace,"sdf",sdf_mc,0,sdf_mc*3);
  setValRange(wspace,"sef",sef_mc,0,sef_mc*3);
//   setValRange(wspace,"eff",eff_mc,0.1,2);

}
//
// single measurement (LM0 or LM1)
//
void RA4Single (StatMethod method, double* sig, double* bkg) {

  RooWorkspace* wspace = createWorkspace();

  // double lm0_mc[4] = { 16.9 , 13.5 , 16.8 , 7.2 }; // tight settings / HT2
  double lm0_mc[6] = { 1.09, 1.97, 7.68, 3.78, 7.85, 21.13 };
  double* lm_mc = sig ? sig : lm0_mc;

  setBackgrounds(wspace,bkg);
  setSignal(wspace,lm_mc);

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
  
  MyLimit limit = computeLimit(wspace,data,method,true);
  std::cout << "Limit [ " << limit.lowerLimit << " , "
	    << limit.upperLimit << " ] ; isIn = " << limit.isInInterval << std::endl;
}
//
// scan over parameter space
//
void RA4Mult (const char* file, double* bkgs,
	      StatMethod method) {


  TFile* fYield = new TFile(file);
  if ( fYield==0 || fYield->IsZombie() ) {
    std::cout << "failed to open " << file << std::endl;
    return;
  }
  const char* cRegion = { "ABCDEF" };
  TH2* hYields[6];
  TH2* hYEntries[6];

  std::string hName;
  for ( unsigned int i=0; i<6; ++i ) {
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

  double yields[6];
  double entries[6];

//   double bkgs[6];
//   bkgs[0] = bkgA;
//   bkgs[1] = bkgB;
//   bkgs[2] = bkgC;
//   bkgs[3] = bkgD;
//   bkgs[4] = bkgE;
//   bkgs[5] = bkgF;

  // *2 for electrons
  for ( int i=0; i<6; ++i )  bkgs[i] *= 2;

//   // muons : 340 / 400 / 470 ; 2.4 / 4.0 / 5.6
//   bkgs[0] = 18.45;
//   bkgs[1] = 18.20;
//   bkgs[2] = 10.77;
//   bkgs[3] = 10.98;
  double kappa1 = (bkgs[0]*bkgs[5])/(bkgs[2]*bkgs[3]);
  double kappa2 = (bkgs[1]*bkgs[5])/(bkgs[2]*bkgs[4]);
  double sigma_kappa_base = 0.10;
  double delta_kappa1_abs = kappa1 - 1.;
  double delta_kappa2_abs = kappa2 - 1.;
  double sigma_kappa1 = sqrt(sigma_kappa_base*sigma_kappa_base+delta_kappa1_abs*delta_kappa1_abs);
  double sigma_kappa2 = sqrt(sigma_kappa_base*sigma_kappa_base+delta_kappa2_abs*delta_kappa2_abs);
  double sigma_kappa = max(sigma_kappa1,sigma_kappa2);

  int nbx = hYields[0]->GetNbinsX();
  int nby = hYields[0]->GetNbinsY();
  for ( int ix=1; ix<=nbx; ++ix ) {
    for ( int iy=1; iy<=nby; ++iy ) {

      for ( unsigned int i=0; i<6; ++i ) {
	yields[i] = hYields[i]->GetBinContent(ix,iy);
	entries[i] = hYEntries[i]->GetBinContent(ix,iy);
      }
//       yields[0] =1.52;
//       yields[1] =7.68;
//       yields[2] =5.17;
//       yields[3] =21.12;
      // *1.5 for electrons and *1.3 for NLO
      for ( unsigned int i=0; i<6; ++i )  yields[i] *= (1.5*1.3);

      setBackgrounds(wspace,bkgs);
      setSignal(wspace,yields);

      setValRange(wspace,"sigmaKappa",sigma_kappa);
      setValRange(wspace,"s",yields[5],0,100);

      // wspace->Print("v");
      // RooArgSet allVars = wspace->allVars();
      // allVars.printLatex(std::cout,1);

      RooDataSet* data = new RooDataSet("data","data",*wspace->set("obs"));
      data->add(*wspace->set("obs"));
      data->Print("v");
  
      MyLimit limit(false,0.,999.);
      std::cout << "Checked ( " << hExclusion->GetXaxis()->GetBinCenter(ix) << " , "
		<< hExclusion->GetYaxis()->GetBinCenter(iy) << " ) with signal yield " 
		<< yields[5] << std::endl;
      if ( yields[5]>0.01 ) {
	limit = computeLimit(wspace,data,method);
	std::cout << "  Limit [ " << limit.lowerLimit << " , "
		  << limit.upperLimit << " ] ; isIn = " << limit.isInInterval << std::endl;
      }
      std::cout << "  yields =";
      for ( unsigned int i=0; i<4; ++i )   cout << " " << yields[i];
      std::cout << std::endl;
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

      delete data;

    }
  }

  TFile* out = new TFile("RA4abcdef.root","RECREATE");
  hExclusion->SetDirectory(out);
  hExclusion->SetMinimum(); hExclusion->SetMaximum();
  hExclusion->SetContour(1); hExclusion->SetContourLevel(0,0.5);
  hLowerLimit->SetDirectory(out);
  hLowerLimit->SetMinimum(); hLowerLimit->SetMaximum();
  hUpperLimit->SetDirectory(out);
  hUpperLimit->SetMinimum(); hUpperLimit->SetMaximum();
  hYields[5]->SetDirectory(out);
  hYields[5]->SetMinimum(); hYields[5]->SetMaximum();
  out->Write();
  delete out;
}

