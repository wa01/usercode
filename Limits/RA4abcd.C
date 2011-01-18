#include "TStopwatch.h"
#include "TCanvas.h"
#include "TROOT.h"
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


using namespace RooFit;
using namespace RooStats;

void RA4abcd(bool doBayesian=false, bool doFeldmanCousins=false, bool doMCMC=false){
  
  // let's time this challenging example
  TStopwatch t;
  t.Start();

  // Inputs : expected values
  double bkg_mc[4] = { 53.26 , 19.48 , 56.13 , 20.22 };
  double lm0_mc[4] = { 3.58 , 8.21 , 16.37 , 25.21 };
  double lm1_mc[4] = { 0.05 , 0.34 , 1.12 , 3.43 };
  // double kappa_mc = 0.984;

  // choice of signal
  double* lm_mc = lm0_mc;

  // event counts
  double observed[4];
  for ( unsigned int i=0; i<4; ++i ) 
    observed[i] = int(bkg_mc[i]+0.5);

  // relative uncertainty on correlation / scaling factors
  double sigma_kappa = 0.20;
  double sigma_sad_rel = 0.030;
  double sigma_sbd_rel = 0.030;
  double sigma_scd_rel = 0.030;
  // double sigma_kappa = 0.05;       
  // double sigma_sad_rel = 0.0030;
  // double sigma_sbd_rel = 0.0030;
  // double sigma_scd_rel = 0.0030;

  // derived quantities
  double bba_mc = bkg_mc[1]/bkg_mc[0];
  double bca_mc = bkg_mc[2]/bkg_mc[0];
  double kappa_mc = (bkg_mc[3]*bkg_mc[0])/(bkg_mc[1]*bkg_mc[2]);
  double sad_mc = lm_mc[0]/lm_mc[3];
  double sbd_mc = lm_mc[1]/lm_mc[3];
  double scd_mc = lm_mc[2]/lm_mc[3];

  // Roo variables
  // .. inputs
  RooRealVar na_var("na","na",observed[0],0,1000);
  RooRealVar nb_var("nb","nb",observed[1],0,1000);
  RooRealVar nc_var("nc","nc",observed[2],0,1000);
  RooRealVar nd_var("nd","nd",observed[3],0,1000);
  // .. pseudo-measurements
  RooRealVar sad_nom_var("sadnom","sadnom",sad_mc,sad_mc/10,sad_mc*10);
  RooRealVar sbd_nom_var("sbdnom","sbdnom",sbd_mc,sbd_mc/10,sbd_mc*10);
  RooRealVar scd_nom_var("scdnom","scdnom",scd_mc,scd_mc/10,scd_mc*10);
  RooRealVar kappa_nom_var("kappanom","kappanom",kappa_mc,kappa_mc/10,kappa_mc*10);
  // .. uncertainties on correlation and signal contamination
  RooRealVar kappa_sig_var("sigmaKappa","sigmaKappa",sigma_kappa);
  RooRealVar sad_sig_var("sigmaSad","sigmaSad",sad_mc*sigma_sad_rel);
  RooRealVar sbd_sig_var("sigmaSbd","sigmaSbd",sbd_mc*sigma_sbd_rel);
  RooRealVar scd_sig_var("sigmaScd","sigmaScd",scd_mc*sigma_scd_rel);

  // .. background and signal variables
  RooRealVar s_var("s","s",lm_mc[3],0,10*lm_mc[3]);
  RooRealVar bkga_var("bkga","bkga",bkg_mc[0],0,10*bkg_mc[3]);
  RooRealVar bba_var("bba","bba",bba_mc,bba_mc/10,bba_mc*10);
  RooRealVar bca_var("bca","bca",bca_mc,bca_mc/10,bca_mc*10);

  // .. correlation and signal contamination variables
  RooRealVar kappa_var(kappa_nom_var,"kappa");  kappa_var.SetTitle("kappa");
  RooRealVar sad_var(sad_nom_var,"sad");  sad_var.SetTitle("sad");
  RooRealVar sbd_var(sbd_nom_var,"sbd");  sbd_var.SetTitle("sbd");
  RooRealVar scd_var(scd_nom_var,"scd");  scd_var.SetTitle("scd");

  // set RooFit random seed for reproducible results
  RooRandom::randomGenerator()->SetSeed(4357);

  // make model
  RooWorkspace* wspace = new RooWorkspace("wspace");

  wspace->import(na_var);
  wspace->import(nb_var);
  wspace->import(nc_var);
  wspace->import(nd_var);
  wspace->import(s_var);
  wspace->import(sad_var);
  wspace->import(sbd_var);
  wspace->import(scd_var);
  wspace->import(bkga_var);
  wspace->import(bba_var);
  wspace->import(bca_var);
  wspace->import(kappa_var);
  wspace->import(kappa_nom_var);
  wspace->import(sad_nom_var);
  wspace->import(sbd_nom_var);
  wspace->import(scd_nom_var);
  wspace->import(kappa_sig_var);
  wspace->import(sad_sig_var);
  wspace->import(sbd_sig_var);
  wspace->import(scd_sig_var);

  wspace->factory("Poisson::a(na, sum::tota(prod::sa(s,sad),bkga))");
  wspace->factory("Poisson::b(nb, sum::totb(prod::sb(s,sbd),prod::bkgb(bkga,bba)))");
  wspace->factory("Poisson::c(nc, sum::totc(prod::sc(s,scd),prod::bkgc(bkga,bca)))");
  wspace->factory("Poisson::d(nd, sum::splusb(s,prod::bkgd(bkga,bba,bca,kappa)))");
  wspace->factory("Gaussian::mcKappa(kappanom, kappa, sigmaKappa)");
  wspace->factory("Gaussian::mcSad(sadnom, sad, sigmaSad)");
  wspace->factory("Gaussian::mcSbd(sbdnom, sbd, sigmaSbd)");
  wspace->factory("Gaussian::mcScd(scdnom, scd, sigmaScd)");
  wspace->factory("PROD::model(d,c,b,a,mcKappa,mcSad,mcSbd,mcScd)");
  wspace->defineSet("obs","nd,nc,nb,na,kappanom,sadnom,sbdnom,scdnom");

  wspace->factory("Uniform::prior_poi({s})");
  wspace->factory("Uniform::prior_nuis({bkga,bba,bca,kappa,sad,sbd,scd})");
  wspace->factory("PROD::prior(prior_poi,prior_nuis)"); 

  ///////////////////////////////////////////
  // Control some interesting variations
  // define parameers of interest
  // for 1-d plots
  wspace->defineSet("poi","s");
  wspace->defineSet("nuis","bkga,bba,bca,kappa,sad,sbd,scd");
  // for 2-d plots to inspect correlations:
  //  wspace->defineSet("poi","s,kappa");

  // // test simpler cases where parameters are known.
  // wspace->var("sad")->setConstant(1);
  // wspace->var("sbd")->setConstant(1);
  // wspace->var("scd")->setConstant(1);

  // inspect workspace
  //  wspace->Print();

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

  /////////////////////////////////////////////////////
  // Now the statistical tests
  // model config
  ModelConfig* modelConfig = new ModelConfig("RA4abcds");
  modelConfig->SetWorkspace(*wspace);
  modelConfig->SetPdf(*wspace->pdf("model"));
  modelConfig->SetPriorPdf(*wspace->pdf("prior"));
  modelConfig->SetParametersOfInterest(*wspace->set("poi"));
  modelConfig->SetNuisanceParameters(*wspace->set("nuis"));
  // wspace->import(*modelConfig);
  // wspace->writeToFile("RA4abcd.root");

  //////////////////////////////////////////////////
  // If you want to see the covariance matrix uncomment
  wspace->pdf("model")->fitTo(*data);

  // use ProfileLikelihood
  ProfileLikelihoodCalculator plc(*data, *modelConfig);
  plc.SetConfidenceLevel(0.95);
  LikelihoodInterval* plInt = plc.GetInterval();
  RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
  plInt->LowerLimit( *wspace->var("s") ); // get ugly print out of the way. Fix.
  RooMsgService::instance().setGlobalKillBelow(msglevel);
  // RooMsgService::instance().setGlobalKillBelow(RooFit::DEBUG);

  // use FeldmaCousins (takes ~20 min)  
  FeldmanCousins fc(*data, *modelConfig);
  fc.SetConfidenceLevel(0.95);
  //number counting: dataset always has 1 entry with N events observed
  fc.FluctuateNumDataEntries(false); 
  fc.UseAdaptiveSampling(true);
  fc.SetNBins(40);
  PointSetInterval* fcInt = NULL;
  if(doFeldmanCousins){ // takes 7 minutes
    fcInt = (PointSetInterval*) fc.GetInterval(); // fix cast
  }


  // use BayesianCalculator (only 1-d parameter of interest, slow for this problem)  
  BayesianCalculator bc(*data, *modelConfig);
  bc.SetConfidenceLevel(0.95);
  SimpleInterval* bInt = NULL;
  if(doBayesian && wspace->set("poi")->getSize() == 1)   {
    bInt = bc.GetInterval();
  } else{
    cout << "Bayesian Calc. only supports on parameter of interest" << endl;
  }


  // use MCMCCalculator  (takes about 1 min)
  // Want an efficient proposal function, so derive it from covariance
  // matrix of fit
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
  if(doMCMC)
    mcInt = mc.GetInterval();

  //////////////////////////////////////
  // Make some  plots
  TCanvas* c1 = (TCanvas*) gROOT->Get("c1");  
  if(!c1)
    c1 = new TCanvas("c1");

  if(doBayesian && doMCMC){
    c1->Divide(3);
    c1->cd(1);
  }
  else if (doBayesian || doMCMC){
    c1->Divide(2);
    c1->cd(1);
  }

  LikelihoodIntervalPlot* lrplot = new LikelihoodIntervalPlot(plInt);
  lrplot->Draw();

  if(doBayesian && wspace->set("poi")->getSize() == 1)   {
    c1->cd(2);
    // the plot takes a long time and print lots of error
    // using a scan it is better
    bc.SetScanOfPosterior(20);
    RooPlot* bplot = bc.GetPosteriorPlot();
    bplot->Draw();
  } 

  if(doMCMC){
    if(doBayesian && wspace->set("poi")->getSize() == 1) 
      c1->cd(3);
    else 
      c1->cd(2);
    MCMCIntervalPlot mcPlot(*mcInt); 
    mcPlot.Draw();
  }

  ////////////////////////////////////
  // querry intervals
  cout << "Profile Likelihood interval on s = [" << 
    plInt->LowerLimit( *wspace->var("s") ) << ", " <<
    plInt->UpperLimit( *wspace->var("s") ) << "]" << endl; 
  //Profile Likelihood interval on s = [12.1902, 88.6871]

   
  if(doBayesian && wspace->set("poi")->getSize() == 1)   {
    cout << "Bayesian interval on s = [" << 
      bInt->LowerLimit( ) << ", " <<
      bInt->UpperLimit( ) << "]" << endl;
  }  
  
  if(doFeldmanCousins){    
    cout << "Feldman Cousins interval on s = [" << 
      fcInt->LowerLimit( *wspace->var("s") ) << ", " <<
      fcInt->UpperLimit( *wspace->var("s") ) << "]" << endl;
    //Feldman Cousins interval on s = [18.75 +/- 2.45, 83.75 +/- 2.45]
  }

  if(doMCMC){
    cout << "MCMC interval on s = [" << 
      mcInt->LowerLimit(*wspace->var("s") ) << ", " <<
      mcInt->UpperLimit(*wspace->var("s") ) << "]" << endl;
    //MCMC interval on s = [15.7628, 84.7266]

  }


  t.Print();
   

}
