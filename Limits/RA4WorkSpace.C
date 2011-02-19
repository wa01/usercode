#include "RA4WorkSpace.h"

// #include "TStopwatch.h"
// #include "TCanvas.h"
#include "TROOT.h"
// #include "TFile.h"
// #include "TH2F.h"
// #include "RooPlot.h"
#include "RooAbsPdf.h"
#include "RooProdPdf.h"
#include "RooUniform.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
// #include "RooGlobalFunc.h"
// #include "RooFitResult.h"
// #include "RooRandom.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooPoisson.h"
// #include "RooStats/ProfileLikelihoodCalculator.h"
// #include "RooStats/LikelihoodInterval.h"
// #include "RooStats/LikelihoodIntervalPlot.h"
// #include "RooStats/BayesianCalculator.h"
// #include "RooStats/MCMCCalculator.h"
// #include "RooStats/MCMCInterval.h"
// #include "RooStats/MCMCIntervalPlot.h"
// #include "RooStats/ProposalHelper.h"
// #include "RooStats/SimpleInterval.h"
// #include "RooStats/FeldmanCousins.h"
// #include "RooStats/PointSetInterval.h"
// #include "RooStats/SimpleLikelihoodRatioTestStat.h"
// #include <vector>

using namespace RooFit;
// using namespace RooStats;

RA4WorkSpace::RA4WorkSpace (const char* name) :
  wspace_(new RooWorkspace(name)), finalized_(false), hasEle_(false), hasMu_(false) 
{
  // total signal in D
  wspace_->factory("s[1,0,100]");
  vS_ = wspace_->var("s");
  // systematic uncertainty on efficiency
  wspace_->factory("effsys[1.,0,2]");
  vEffsys_ = wspace_->var("effsys");
  // systematic uncertainty on kappa
  wspace_->factory("kappasys[1,0,2]");
  vKappasys_ = wspace_->var("kappasys");
  // systematic uncertainty on signal contamination
  wspace_->factory("scontsys[1,0,2]");
  vSContsys_ = wspace_->var("scontsys");
  // pseudo-measurements for kappa and signal contamination syst
  wspace_->factory("effnom[1]");
  wspace_->factory("kappanom[1]");
  wspace_->factory("scontnom[1]");
  wspace_->factory("sigmaEff[0.15]");
  wspace_->factory("sigmaKappa[0.15]");
  wspace_->factory("sigmaScont[0.15]");
  // Pdfs for pseudo-measurements
  wspace_->factory("Gaussian::mcEff(effnom, effsys, sigmaEff)");
  pdfMcEff_ = wspace_->pdf("mcEff");
  wspace_->factory("Gaussian::mcKappa(kappanom, kappasys, sigmaKappa)");
  pdfMcKappa_ = wspace_->pdf("mcKappa");
  wspace_->factory("Gaussian::mcScont(scontnom, scontsys, sigmaScont)");
  pdfMcSCont_ = wspace_->pdf("mcScont");
  // // model
  // wspace_->factory("PROD::model(mcKappa,mcScont)");
  // pdfModel_ = wspace_->pdf("model");
  // // priors
  // wspace_->factory("Uniform::prior_poi({s})");
  wspace_->defineSet("obs","effnom,kappanom,scontnom");
  setObs_ = wspace_->set("obs");
  wspace_->defineSet("poi","s");
  wspace_->defineSet("nuis","effsys,kappasys,scontsys");
  setNuis_ = wspace_->set("nuis");

  wspace_->Print("v");
}

void 
RA4WorkSpace::setBackground (float bkgA, float bkgB, float bkgC, float bkgD)
{

}

void
RA4WorkSpace::addChannel (ChannelType channel)
{
  if ( finalized_ ) {
    std::cout << "Workspace has already been finalized" << std::endl;
    return;
  }
  if ( channel!=EleChannel && channel!=MuChannel ) {
    std::cout << "No such channel type" << std::endl;
    return;
  }
  if ( hasEle_ && channel==EleChannel ) {
    std::cout << "Electron channel is already defined" << std::endl;
    return;
  }
  if ( hasMu_ && channel==MuChannel ) {
    std::cout << "Muon channel is already defined" << std::endl;
    return;
  }
  
  std::string suffix(channel==EleChannel?"E":"M");
  const char regions[] = { "abcd" };
  //
  // observation
  //
  for ( unsigned int i=0; i<4; ++i ) {
    std::string name("n");
    name = name + regions[i] + suffix;
    RooRealVar tmp(name.c_str(),name.c_str(),0,0,1000);
    wspace_->import(tmp);
    vObs_[i][channel] = wspace_->var(name.c_str());
    wspace_->extendSet("obs",name.c_str());
  }
  // wspace_->Print("v");
  //
  // efficiency
  //
  {
    std::string name("eff");
    name += suffix;
    RooRealVar tmp(name.c_str(),name.c_str(),1.);
    wspace_->import(tmp);
    vEff_[channel] = wspace_->var(name.c_str());
  }
  //
  // kappa
  //
  {
    std::string name("kappa");
    name += suffix;
    RooRealVar tmp(name.c_str(),name.c_str(),1.);
    wspace_->import(tmp);
    vKappa_[channel] = wspace_->var(name.c_str());
  }
  // std::cout << "Eff" << std::endl;
  // wspace_->Print("v");
  // std::cout << "Eff2" << std::endl;
  //
  // signal contamination
  //
  for ( unsigned int i=0; i<3; ++i ) {
    // std::cout << i << std::endl;
    std::string name("s");
    name = name + regions[i] + "d" + suffix;
    RooRealVar tmp(name.c_str(),name.c_str(),0,0,10);
    wspace_->import(tmp);
    vSCont_[i][channel] = wspace_->var(name.c_str());
    // wspace_->extendSet("nuis",name.c_str());
  }
  // wspace_->Print("v");

  // bkg in A; relative bkg in B&C; kappa
  for ( unsigned int i=0; i<2; ++i ) {
    std::string name("b");
    name = name + regions[i+1] + "d" + suffix;
    RooRealVar tmp(name.c_str(),name.c_str(),0,0,10);
    wspace_->import(tmp);
    vBxd_[i][channel] = wspace_->var(name.c_str());
    wspace_->extendSet("nuis",name.c_str());
  }
  // wspace_->Print("v");
  {
    std::string name("bkgd");
    name = name + suffix;
    RooRealVar tmp(name.c_str(),name.c_str(),1,0,1000);
    wspace_->import(tmp);
    vBkgd_[channel] = wspace_->var(name.c_str());
    wspace_->extendSet("nuis",name.c_str());
  }
  //
  // pdfs for the 4 regions
  //
  // signal in D
  //
  {
    std::string name("sd");
    name += suffix;
    RooArgSet tmpSet(*wspace_->var("s"),*vEff_[channel],*wspace_->var("effsys"));
    RooProduct tmp(name.c_str(),name.c_str(),tmpSet);
    wspace_->import(tmp);
    pSd_[channel] = (RooProduct*)wspace_->var(name.c_str());
  }
  //
  // region A
  //
  {
    std::string name1("sa");
    name1 += suffix;
    RooArgSet tmpSetSa(*wspace_->var("s"),*vEff_[channel],*vSCont_[0][channel],*wspace_->var("scontsys"));
    RooProduct tmpSa(name1.c_str(),name1.c_str(),tmpSetSa);
    wspace_->import(tmpSa);
    std::string name2("bkga");
    name2 += suffix;
    RooArgSet tmpSetBkga(*vBkgd_[channel],*vBxd_[0][channel],*vBxd_[1][channel],*vKappa_[channel],*wspace_->var("kappasys"));
    RooProduct tmpBkga(name2.c_str(),name2.c_str(),tmpSetBkga);
    wspace_->import(tmpBkga);
    std::string name3("tota");
    name3 += suffix;
    RooArgSet tmpSetTota(*wspace_->function(name1.c_str()),*wspace_->function(name2.c_str()));
    RooAddition tmpTota(name3.c_str(),name3.c_str(),tmpSetTota);
    wspace_->import(tmpTota);
    std::string name("a");
    name += suffix;
    RooPoisson tmpPoisA(name.c_str(),name.c_str(),*vObs_[0][channel],*wspace_->function(name3.c_str()));
    wspace_->import(tmpPoisA);
    pdfReg_[0][channel] = wspace_->pdf(name.c_str());
  }
  //
  // region B
  //
  {
    std::string name1("sb");
    name1 += suffix;
    RooArgSet tmpSetSb(*wspace_->var("s"),*vEff_[channel],*vSCont_[1][channel],*wspace_->var("scontsys"));
    RooProduct tmpSb(name1.c_str(),name1.c_str(),tmpSetSb);
    wspace_->import(tmpSb);
    std::string name2("bkgb");
    name2 += suffix;
    RooArgSet tmpSetBkgb(*vBkgd_[channel],*vBxd_[0][channel]);
    RooProduct tmpBkgb(name2.c_str(),name2.c_str(),tmpSetBkgb);
    wspace_->import(tmpBkgb);
    std::string name3("totb");
    name3 += suffix;
    RooArgSet tmpSetTotb(*wspace_->function(name1.c_str()),*wspace_->function(name2.c_str()));
    RooAddition tmpTotb(name3.c_str(),name3.c_str(),tmpSetTotb);
    wspace_->import(tmpTotb);
    std::string name("b");
    name += suffix;
    RooPoisson tmpPoisB(name.c_str(),name.c_str(),*vObs_[1][channel],*wspace_->function(name3.c_str()));
    wspace_->import(tmpPoisB);
    pdfReg_[1][channel] = wspace_->pdf(name.c_str());
  }
  //
  // region C
  //
  {
    std::string name1("sc");
    name1 += suffix;
    RooArgSet tmpSetSc(*wspace_->var("s"),*vEff_[channel],*vSCont_[2][channel],*wspace_->var("scontsys"));
    RooProduct tmpSc(name1.c_str(),name1.c_str(),tmpSetSc);
    wspace_->import(tmpSc);
    std::string name2("bkgc");
    name2 += suffix;
    RooArgSet tmpSetBkgc(*vBkgd_[channel],*vBxd_[1][channel]);
    RooProduct tmpBkgc(name2.c_str(),name2.c_str(),tmpSetBkgc);
    wspace_->import(tmpBkgc);
    std::string name3("totc");
    name3 += suffix;
    RooArgSet tmpSetTotc(*wspace_->function(name1.c_str()),*wspace_->function(name2.c_str()));
    RooAddition tmpTotc(name3.c_str(),name3.c_str(),tmpSetTotc);
    wspace_->import(tmpTotc);
    std::string name("c");
    name += suffix;
    RooPoisson tmpPoisC(name.c_str(),name.c_str(),*vObs_[2][channel],*wspace_->function(name3.c_str()));
    wspace_->import(tmpPoisC);
    pdfReg_[2][channel] = wspace_->pdf(name.c_str());
  }
  //
  // region D
  //
  {
    std::string name1("sd");
    name1 += suffix;
    std::string name2("bkgd");
    name2 += suffix;
    std::string name3("totd");
    name3 += suffix;
    RooArgSet tmpSetTotd(*wspace_->function(name1.c_str()),*wspace_->function(name2.c_str()));
    RooAddition tmpTotd(name3.c_str(),name3.c_str(),tmpSetTotd);
    wspace_->import(tmpTotd);
    std::string name("d");
    name += suffix;
    RooPoisson tmpPoisD(name.c_str(),name.c_str(),*vObs_[3][channel],*wspace_->function(name3.c_str()));
    wspace_->import(tmpPoisD);
    pdfReg_[3][channel] = wspace_->pdf(name.c_str());
  }

  switch ( channel ) {
  case EleChannel:
    hasEle_ = true;
    break;
  case MuChannel:
    hasMu_ = true;
    break;
  }
}

void 
RA4WorkSpace::finalize ()
{
  //
  // check
  //
  if ( !hasEle_ && !hasMu_ ) {
    std::cout << "no channel has been defined" << std::endl;
    return;
  }
  //
  // combined model
  RooArgSet modelSet(*pdfMcKappa_,*pdfMcSCont_);
  if ( hasEle_ ) {
    for ( unsigned int i=0; i<4; ++i )  modelSet.add(*pdfReg_[i][0]);
  }
  if ( hasMu_ ) {
    for ( unsigned int i=0; i<4; ++i )  modelSet.add(*pdfReg_[i][1]);
  }
  RooProdPdf model("model","model",modelSet);
  wspace_->import(model);
  //
  // priors
  //
  wspace_->factory("Uniform::prior_poi({s})");
  RooUniform priorNuis("prior_nuis","prior_nuis",*setNuis_);
  wspace_->import(priorNuis);
  wspace_->Print("v");

  finalized_ = true;
}

  // {
  //   std::string name("sd");
  //   name += suffix;
  //   RooArgSet tmpSet(*wspace_->var("s"),*vEff_[channel],*wspace_->var("effsys"));
  //   RooProduct tmp(name.c_str(),name.c_str(),tmpSet);
  //   wspace_->import(tmp);
  //   pSd_[channel] = (RooProduct*)wspace_->var(name.c_str());
  //   RooArgList pdfList(
  // }

//   wspace->factory("bbd[0,0,10]");
//   wspace->factory("bcd[0,0,10]");
//   wspace->factory("bkgd[1,0,1000]");
//   wspace->factory("kappa[1,0,2]");
//   // pseudo-measurements for kappa and signal contamination
//   wspace->factory("kappanom[1]");
//   wspace->factory("sadnom[0.1]");
//   wspace->factory("sbdnom[0.1]");
//   wspace->factory("scdnom[0.1]");
//   wspace->factory("effnom[1.]");
//   wspace->factory("sigmaKappa[0.1]");
//   wspace->factory("sigmaSad[0.1]");
//   wspace->factory("sigmaSbd[0.1]");
//   wspace->factory("sigmaScd[0.1]");
//   wspace->factory("sigmaEff[0.1]");
//   // Poisson distributions in the 4 regions
//   wspace->factory("prod::sd(s,eff)");
//   wspace->factory("Poisson::a(na, sum::tota(prod::sa(sd,sad),prod::bkga(bkgd,bbd,bcd,kappa)))");
//   wspace->factory("Poisson::b(nb, sum::totb(prod::sb(sd,sbd),prod::bkgb(bkgd,bbd)))");
//   wspace->factory("Poisson::c(nc, sum::totc(prod::sc(sd,scd),prod::bkgc(bkgd,bcd)))");
//   wspace->factory("Poisson::d(nd, sum::splusb(sd,bkgd))");
//   // Pdfs for pseudo-measurements
//   wspace->factory("Gaussian::mcKappa(kappanom, kappa, sigmaKappa)");
//   wspace->factory("Gaussian::mcSad(sadnom, sad, sigmaSad)");
//   wspace->factory("Gaussian::mcSbd(sbdnom, sbd, sigmaSbd)");
//   wspace->factory("Gaussian::mcScd(scdnom, scd, sigmaScd)");
// //   wspace->factory("Gaussian::mcEff(effnom, eff, sigmaEff)");
//   // full model
// //   wspace->factory("PROD::model(d,c,b,a,mcKappa,mcSad,mcSbd,mcScd,mcEff)");
//   wspace->factory("PROD::model(d,c,b,a,mcKappa,mcSad,mcSbd,mcScd)");
//   // priors
//   wspace->factory("Uniform::prior_poi({s})");
// //   wspace->factory("Uniform::prior_nuis({bkgd,bbd,bcd,kappa,sad,sbd,scd,eff})");
//   wspace->factory("Uniform::prior_nuis({bkgd,bbd,bcd,kappa,sad,sbd,scd})");
//   wspace->factory("PROD::prior(prior_poi,prior_nuis)"); 
//   // sets (observables, POI, nuisance parameters)
// //   wspace->defineSet("obs","nd,nc,nb,na,kappanom,sadnom,sbdnom,scdnom,effnom");
//   wspace->defineSet("obs","nd,nc,nb,na,kappanom,sadnom,sbdnom,scdnom");
//   wspace->defineSet("poi","s");
// //   wspace->defineSet("nuis","bkgd,bbd,bcd,kappa,sad,sbd,scd,eff");
//   wspace->defineSet("nuis","bkgd,bbd,bcd,kappa,sad,sbd,scd");

//   return wspace;

// }

void 
RA4WorkSpace::setSignal (float sigA, float sigB, float sigC, float sigD)
{
}

// //
// // calculation of the limit: assumes that wspace is set up and observations
// //   contained in data
// //
// MyLimit computeLimit (RooWorkspace* wspace, RooDataSet* data, StatMethod method, bool draw) {

//   // let's time this challenging example
//   TStopwatch t;
//   // for 2-d plots to inspect correlations:
//   //  wspace->defineSet("poi","s,kappa");

//   // // test simpler cases where parameters are known.
//   // wspace->var("sad")->setConstant(1);
//   // wspace->var("sbd")->setConstant(1);
//   // wspace->var("scd")->setConstant(1);

//   // inspect workspace
//   //  wspace->Print();

//   //
//   // get nominal signal
//   //
//   RooRealVar exp_sig(*wspace->var("s"));
//   // std::cout << "exp_sig = " << exp_sig.getVal() << std::endl;
  
//   /////////////////////////////////////////////////////
//   // Now the statistical tests
//   // model config
//   ModelConfig modelConfig("RA4abcd");
//   modelConfig.SetWorkspace(*wspace);
//   modelConfig.SetPdf(*wspace->pdf("model"));
//   modelConfig.SetPriorPdf(*wspace->pdf("prior"));
//   modelConfig.SetParametersOfInterest(*wspace->set("poi"));
//   modelConfig.SetNuisanceParameters(*wspace->set("nuis"));

//   //////////////////////////////////////////////////
//   // If you want to see the covariance matrix uncomment
//   // wspace->pdf("model")->fitTo(*data);

//   // use ProfileLikelihood
//   if ( method == ProfileLikelihoodMethod ) {
//     ProfileLikelihoodCalculator plc(*data, modelConfig);
//     plc.SetConfidenceLevel(0.95);
//     RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();
//     RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
//     LikelihoodInterval* plInt = plc.GetInterval();
//     RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
//     plInt->LowerLimit( *wspace->var("s") ); // get ugly print out of the way. Fix.
//     // RooMsgService::instance().setGlobalKillBelow(RooFit::DEBUG);
//     if ( draw ) {
//       TCanvas* c = new TCanvas("ProfileLikelihood");
//       LikelihoodIntervalPlot* lrplot = new LikelihoodIntervalPlot(plInt);
//       lrplot->Draw();
//     }
//     RooMsgService::instance().setGlobalKillBelow(msglevel);
//     cout << "Profile Likelihood interval on s = [" << 
//       plInt->LowerLimit( *wspace->var("s") ) << ", " <<
//       plInt->UpperLimit( *wspace->var("s") ) << "]" << endl; 
//     MyLimit result(plInt->IsInInterval(exp_sig),
// 		   plInt->LowerLimit(*wspace->var("s")),plInt->UpperLimit(*wspace->var("s")));
//     // std::cout << "isIn " << result << std::endl;
//     delete plInt;
// //     delete modelConfig;
//     return result;
//   }

//   // use FeldmaCousins (takes ~20 min)  
//   if ( method == FeldmanCousinsMethod ) {
//     FeldmanCousins fc(*data, modelConfig);
//     fc.SetConfidenceLevel(0.95);
//     //number counting: dataset always has 1 entry with N events observed
//     fc.FluctuateNumDataEntries(false); 
//     fc.UseAdaptiveSampling(true);
//     fc.SetNBins(100);
//     PointSetInterval* fcInt = NULL;
//     fcInt = (PointSetInterval*) fc.GetInterval(); // fix cast
//     cout << "Feldman Cousins interval on s = [" << 
//       fcInt->LowerLimit( *wspace->var("s") ) << ", " <<
//       fcInt->UpperLimit( *wspace->var("s") ) << "]" << endl;
//     // std::cout << "isIn " << result << std::endl;
//     MyLimit result(fcInt->IsInInterval(exp_sig),
// 		   fcInt->LowerLimit(*wspace->var("s")),fcInt->UpperLimit(*wspace->var("s")));
//     delete fcInt;
//     return result;
//   }


//   // use BayesianCalculator (only 1-d parameter of interest, slow for this problem)  
//   if ( method == BayesianMethod ) {
//     BayesianCalculator bc(*data, modelConfig);
//     bc.SetConfidenceLevel(0.95);
//     bc.SetLeftSideTailFraction(0.5);
//     SimpleInterval* bInt = NULL;
//     if( wspace->set("poi")->getSize() == 1)   {
//       bInt = bc.GetInterval();
//       if ( draw ) {
// 	TCanvas* c = new TCanvas("Bayesian");
// 	// the plot takes a long time and print lots of error
// 	// using a scan it is better
// 	bc.SetScanOfPosterior(50);
// 	RooPlot* bplot = bc.GetPosteriorPlot();
// 	bplot->Draw();
//       }
//       cout << "Bayesian interval on s = [" << 
// 	bInt->LowerLimit( ) << ", " <<
// 	bInt->UpperLimit( ) << "]" << endl;
//       // std::cout << "isIn " << result << std::endl;
//       MyLimit result(bInt->IsInInterval(exp_sig),
// 		     bInt->LowerLimit(),bInt->UpperLimit());
//       delete bInt;
//       return result;
//     } else {
//     cout << "Bayesian Calc. only supports on parameter of interest" << endl;
//     return MyLimit();
//     }
//   }


//   // use MCMCCalculator  (takes about 1 min)
//   // Want an efficient proposal function, so derive it from covariance
//   // matrix of fit
//   if ( method == MCMCMethod ) {
//     RooFitResult* fit = wspace->pdf("model")->fitTo(*data,Save());
//     ProposalHelper ph;
//     ph.SetVariables((RooArgSet&)fit->floatParsFinal());
//     ph.SetCovMatrix(fit->covarianceMatrix());
//     ph.SetUpdateProposalParameters(kTRUE); // auto-create mean vars and add mappings
//     ph.SetCacheSize(100);
//     ProposalFunction* pf = ph.GetProposalFunction();
    
//     MCMCCalculator mc(*data, modelConfig);
//     mc.SetConfidenceLevel(0.95);
//     mc.SetProposalFunction(*pf);
//     mc.SetNumBurnInSteps(100); // first N steps to be ignored as burn-in
//     mc.SetNumIters(100000);
//     mc.SetLeftSideTailFraction(0.5); // make a central interval
//     MCMCInterval* mcInt = NULL;
//     mcInt = mc.GetInterval();
//     MCMCIntervalPlot mcPlot(*mcInt); 
//     mcPlot.Draw();
//     cout << "MCMC interval on s = [" << 
//       mcInt->LowerLimit(*wspace->var("s") ) << ", " <<
//       mcInt->UpperLimit(*wspace->var("s") ) << "]" << endl;
//     // std::cout << "isIn " << result << std::endl;
//     MyLimit result(mcInt->IsInInterval(exp_sig),
// 		   mcInt->LowerLimit(*wspace->var("s")),mcInt->UpperLimit(*wspace->var("s")));
//     delete mcInt;
//     return result;
//   }
  

//   t.Print();

// //   delete modelConfig;
//   return MyLimit();

// }
// //
// // set value and range for a variable in the workspace
// //
// void setValRange (RooWorkspace* workspace, const char* name, double val, double vmin, double vmax)
// {
//   RooRealVar* var = workspace->var(name);
//   if ( var ) {
//     if ( vmax>vmin )  var->setRange(vmin,vmax);
//     var->setVal(val);
//   }
// }
// //
// // create a workspace with the variables (with dummy values) and the pdfs
// //
// RooWorkspace* createWorkspace (const char* name)
// {
//   RooWorkspace* wspace = new RooWorkspace(name);
//   // observed event counts
//   wspace->factory("na[0,0,1000]");
//   wspace->factory("nb[0,0,1000]");
//   wspace->factory("nc[0,0,1000]");
//   wspace->factory("nd[0,0,1000]");
//   // signal in D and rel. signal contamination
//   wspace->factory("s[1,0,100]");
//   wspace->factory("sad[0,0,10]");
//   wspace->factory("sbd[0,0,10]");
//   wspace->factory("scd[0,0,10]");
// //   wspace->factory("eff[1.,0.1,2.]");
//   wspace->factory("eff[1.]");
//   // bkg in A; relative bkg in B&C; kappa
//   wspace->factory("bbd[0,0,10]");
//   wspace->factory("bcd[0,0,10]");
//   wspace->factory("bkgd[1,0,1000]");
//   wspace->factory("kappa[1,0,2]");
//   // pseudo-measurements for kappa and signal contamination
//   wspace->factory("kappanom[1]");
//   wspace->factory("sadnom[0.1]");
//   wspace->factory("sbdnom[0.1]");
//   wspace->factory("scdnom[0.1]");
//   wspace->factory("effnom[1.]");
//   wspace->factory("sigmaKappa[0.1]");
//   wspace->factory("sigmaSad[0.1]");
//   wspace->factory("sigmaSbd[0.1]");
//   wspace->factory("sigmaScd[0.1]");
//   wspace->factory("sigmaEff[0.1]");
//   // Poisson distributions in the 4 regions
//   wspace->factory("prod::sd(s,eff)");
//   wspace->factory("Poisson::a(na, sum::tota(prod::sa(sd,sad),prod::bkga(bkgd,bbd,bcd,kappa)))");
//   wspace->factory("Poisson::b(nb, sum::totb(prod::sb(sd,sbd),prod::bkgb(bkgd,bbd)))");
//   wspace->factory("Poisson::c(nc, sum::totc(prod::sc(sd,scd),prod::bkgc(bkgd,bcd)))");
//   wspace->factory("Poisson::d(nd, sum::splusb(sd,bkgd))");
//   // Pdfs for pseudo-measurements
//   wspace->factory("Gaussian::mcKappa(kappanom, kappa, sigmaKappa)");
//   wspace->factory("Gaussian::mcSad(sadnom, sad, sigmaSad)");
//   wspace->factory("Gaussian::mcSbd(sbdnom, sbd, sigmaSbd)");
//   wspace->factory("Gaussian::mcScd(scdnom, scd, sigmaScd)");
// //   wspace->factory("Gaussian::mcEff(effnom, eff, sigmaEff)");
//   // full model
// //   wspace->factory("PROD::model(d,c,b,a,mcKappa,mcSad,mcSbd,mcScd,mcEff)");
//   wspace->factory("PROD::model(d,c,b,a,mcKappa,mcSad,mcSbd,mcScd)");
//   // priors
//   wspace->factory("Uniform::prior_poi({s})");
// //   wspace->factory("Uniform::prior_nuis({bkgd,bbd,bcd,kappa,sad,sbd,scd,eff})");
//   wspace->factory("Uniform::prior_nuis({bkgd,bbd,bcd,kappa,sad,sbd,scd})");
//   wspace->factory("PROD::prior(prior_poi,prior_nuis)"); 
//   // sets (observables, POI, nuisance parameters)
// //   wspace->defineSet("obs","nd,nc,nb,na,kappanom,sadnom,sbdnom,scdnom,effnom");
//   wspace->defineSet("obs","nd,nc,nb,na,kappanom,sadnom,sbdnom,scdnom");
//   wspace->defineSet("poi","s");
// //   wspace->defineSet("nuis","bkgd,bbd,bcd,kappa,sad,sbd,scd,eff");
//   wspace->defineSet("nuis","bkgd,bbd,bcd,kappa,sad,sbd,scd");

//   return wspace;
// }
// //
// // set background-related variables of the workspace
// //
// void setBackgrounds (RooWorkspace* wspace, double* bkgs) 
// {
//   // Inputs : expected values
//   // double bkg_mc[4] = { 120. , 12.1 , 18.3 , 1.83 }; // tight settings / HT2
//   double bkg_mc[4] = {  14.73, 18.20, 8.48, 10.98 };
//   if ( bkgs ) {
//     for ( unsigned int i=0; i<4; ++i )  bkg_mc[i] = bkgs[i];
//   }
//   //
//   // use pessimistic scenario for rounding of expected numbers
//   //  
//   double observed[4];
// //   for ( unsigned int i=0; i<4; ++i ) 
// //     observed[i] = int(bkg_mc[i]+0.5);
//   observed[0] = int(bkg_mc[0])+1;
//   observed[1] = int(bkg_mc[1]);
//   observed[2] = int(bkg_mc[2]);
//   observed[3] = int(bkg_mc[3])+1;

//   // scaling factors
//   double sigma_kappa = 0.25;

//   // derived quantities
//   double bbd_mc = bkg_mc[1]/bkg_mc[3];
//   double bcd_mc = bkg_mc[2]/bkg_mc[3];
//   double kappa_mc = (bkg_mc[3]*bkg_mc[0])/(bkg_mc[1]*bkg_mc[2]);

//   // Roo variables
//   // .. inputs set to expected backgrounds
//   setValRange(wspace,"na",observed[0],0,1000);
//   setValRange(wspace,"nb",observed[1],0,1000);
//   setValRange(wspace,"nc",observed[2],0,1000);
//   setValRange(wspace,"nd",observed[3],0,1000);
//   // .. pseudo-measurements
// //   setValRange(wspace,"kappanom",kappa_mc,kappa_mc/10,kappa_mc*10);
//   setValRange(wspace,"kappanom",kappa_mc);
//   // .. uncertainties on correlation and signal contamination
//   setValRange(wspace,"sigmaKappa",sigma_kappa);

//   // .. background variables
//   setValRange(wspace,"bkgd",bkg_mc[3],0,10*bkg_mc[3]);
//   setValRange(wspace,"bbd",bbd_mc,0,bbd_mc*2);
//   setValRange(wspace,"bcd",bcd_mc,0,bcd_mc*2);

//   // .. correlation and signal contamination variables
//   setValRange(wspace,"kappa",kappa_mc,0,2);
// }

// void setSignal (RooWorkspace* wspace, double* lm_mc)
// {
//   // Inputs : expected values

//   // relative uncertainties on relative signal contamination
//   double sigma_sad_rel = 0.100;
//   double sigma_sbd_rel = 0.100;
//   double sigma_scd_rel = 0.100;
//   double sigma_eff = 0.001;
//   // double sigma_sad_rel = 0.0030;
//   // double sigma_sbd_rel = 0.0030;
//   // double sigma_scd_rel = 0.0030;

//   // derived quantities
//   double sad_mc = max(lm_mc[0]/lm_mc[3],0.01);
//   double sbd_mc = max(lm_mc[1]/lm_mc[3],0.01);
//   double scd_mc = max(lm_mc[2]/lm_mc[3],0.01);
//   double eff_mc = 1.;

//   // Roo variables
//   // .. pseudo-measurements
// //   setValRange(wspace,"sadnom",sad_mc,sad_mc/10,sad_mc*10);
// //   setValRange(wspace,"sbdnom",sbd_mc,sbd_mc/10,sbd_mc*10);
// //   setValRange(wspace,"scdnom",scd_mc,scd_mc/10,scd_mc*10);
//   setValRange(wspace,"sadnom",sad_mc);
//   setValRange(wspace,"sbdnom",sbd_mc);
//   setValRange(wspace,"scdnom",scd_mc);
//   // .. uncertainties on signal contamination
//   setValRange(wspace,"sigmaSad",sad_mc*sigma_sad_rel);
//   setValRange(wspace,"sigmaSbd",sbd_mc*sigma_sbd_rel);
//   setValRange(wspace,"sigmaScd",scd_mc*sigma_scd_rel);
//   setValRange(wspace,"sigmaEff",sigma_eff);

//   // .. background and signal variables
//   setValRange(wspace,"s",lm_mc[3],0,10*lm_mc[3]);
//   // wspace->var("s")->Print("v");

//   // .. correlation and signal contamination variables
//   setValRange(wspace,"sad",sad_mc,0,sad_mc*3);
//   setValRange(wspace,"sbd",sbd_mc,0,sbd_mc*3);
//   setValRange(wspace,"scd",scd_mc,0,scd_mc*3);
// //   setValRange(wspace,"eff",eff_mc,0.1,2);

// }
// //
// // single measurement (LM0 or LM1)
// //
// void RA4Single (StatMethod method, double* sig, double* bkg) {

//   RooWorkspace* wspace = createWorkspace();

//   // double lm0_mc[4] = { 16.9 , 13.5 , 16.8 , 7.2 }; // tight settings / HT2
//   double lm0_mc[4] = { 1.09, 7.68, 3.78, 21.13 };
//   double lm1_mc[4] = { 0.05 , 0.34 , 1.12 , 3.43 };
//   double* lm_mc = sig ? sig : lm0_mc;

//   setBackgrounds(wspace,bkg);
//   setSignal(wspace,lm_mc);

//   wspace->Print("v");
//   // RooArgSet allVars = wspace->allVars();
//   // allVars.printLatex(std::cout,1);

//   ////////////////////////////////////////////////////////////
//   // Generate toy data
//   // generate toy data assuming current value of the parameters
//   // import into workspace. 
//   // add Verbose() to see how it's being generated
//   // wspace->var("s")->setVal(0.);
//   // RooDataSet* data =   wspace->pdf("model")->generate(*wspace->set("obs"),1);
//   // data->Print("v");
//   // // wspace->import(*data);
//   // wspace->var("s")->setVal(lm_mc[3]);
//   RooDataSet* data = new RooDataSet("data","data",*wspace->set("obs"));
//   data->add(*wspace->set("obs"));
//   data->Print("v");
  
//   MyLimit limit = computeLimit(wspace,data,method,true);
//   std::cout << "Limit [ " << limit.lowerLimit << " , "
// 	    << limit.upperLimit << " ] ; isIn = " << limit.isInInterval << std::endl;
// }
// //
// // scan over parameter space
// //
// void RA4Mult (const char* file, float bkgA, float bkgB, float bkgC, float bkgD, 
// 	      StatMethod method) {


//   TFile* fYield = new TFile(file);
//   if ( fYield==0 || fYield->IsZombie() ) {
//     std::cout << "failed to open " << file << std::endl;
//     return;
//   }
//   const char* cRegion = { "ABCD" };
//   TH2* hYields[4];
//   TH2* hYEntries[4];

//   std::string hName;
//   for ( unsigned int i=0; i<4; ++i ) {
//     hName = "Events";
//     hName += cRegion[i];
//     hYields[i] = (TH2*)fYield->Get(hName.c_str());
//     hName = "Entries";
//     hName += cRegion[i];
//     hYEntries[i] = (TH2*)fYield->Get(hName.c_str());
//     if ( hYields[i]==0 || hYEntries[i]==0 ) {
//       std::cout << "Missing histogram for region " << cRegion[i] << std::endl;
//       return;
//     }
//   }

//   gROOT->cd();
//   TH2* hExclusion = (TH2*)hYields[0]->Clone("Exclusion");
//   hExclusion->Reset();
//   hExclusion->SetTitle("Exclusion");
//   TH2* hLowerLimit = (TH2*)hYields[0]->Clone("LowerLimit");
//   hLowerLimit->Reset();
//   hLowerLimit->SetTitle("LowerLimit");
//   TH2* hUpperLimit = (TH2*)hYields[0]->Clone("UpperLimit");
//   hUpperLimit->Reset();
//   hUpperLimit->SetTitle("UpperLimit");

//   RooWorkspace* wspace = createWorkspace();

//   double yields[4];
//   double entries[4];

//   double bkgs[4];
//   bkgs[0] = bkgA;
//   bkgs[1] = bkgB;
//   bkgs[2] = bkgC;
//   bkgs[3] = bkgD;

//   // *2 for electrons
//   for ( int i=0; i<4; ++i )  bkgs[i] *= 2;

// //   // muons : 340 / 400 / 470 ; 2.4 / 4.0 / 5.6
// //   bkgs[0] = 18.45;
// //   bkgs[1] = 18.20;
// //   bkgs[2] = 10.77;
// //   bkgs[3] = 10.98;
//   double kappa = (bkgs[0]*bkgs[3])/(bkgs[1]*bkgs[2]);
//   double sigma_kappa_base = 0.10;
//   double delta_kappa_abs = kappa - 1.;
//   double sigma_kappa = sqrt(sigma_kappa_base*sigma_kappa_base+delta_kappa_abs*delta_kappa_abs);

//   int nbx = hYields[0]->GetNbinsX();
//   int nby = hYields[0]->GetNbinsY();
//   for ( int ix=1; ix<=nbx; ++ix ) {
//     for ( int iy=1; iy<=nby; ++iy ) {

//       for ( unsigned int i=0; i<4; ++i ) {
// 	yields[i] = hYields[i]->GetBinContent(ix,iy);
// 	entries[i] = hYEntries[i]->GetBinContent(ix,iy);
//       }
// //       yields[0] =1.52;
// //       yields[1] =7.68;
// //       yields[2] =5.17;
// //       yields[3] =21.12;
//       // *1.5 for electrons and *1.3 for NLO
//       for ( unsigned int i=0; i<4; ++i )  yields[i] *= (1.5*1.3);

//       setBackgrounds(wspace,bkgs);
//       setSignal(wspace,yields);

//       setValRange(wspace,"sigmaKappa",sigma_kappa);
//       setValRange(wspace,"s",yields[3],0,100);

//       // wspace->Print("v");
//       // RooArgSet allVars = wspace->allVars();
//       // allVars.printLatex(std::cout,1);

//       RooDataSet* data = new RooDataSet("data","data",*wspace->set("obs"));
//       data->add(*wspace->set("obs"));
//       data->Print("v");
  
//       MyLimit limit(false,0.,999.);
//       std::cout << "Checked ( " << hExclusion->GetXaxis()->GetBinCenter(ix) << " , "
// 		<< hExclusion->GetYaxis()->GetBinCenter(iy) << " ) with signal yield " 
// 		<< yields[3] << std::endl;
//       if ( yields[3]>0.01 ) {
// 	limit = computeLimit(wspace,data,method);
// 	std::cout << "  Limit [ " << limit.lowerLimit << " , "
// 		  << limit.upperLimit << " ] ; isIn = " << limit.isInInterval << std::endl;
//       }
//       std::cout << "  yields =" 
// 		<< " " << yields[0]
// 		<< " " << yields[1]
// 		<< " " << yields[2]
// 		<< " " << yields[3] << std::endl;
// //       std::cout << "  entries =" 
// // 		<< " " << entries[0]
// // 		<< " " << entries[1]
// // 		<< " " << entries[2]
// // 		<< " " << entries[3] << std::endl;
//       double excl = limit.isInInterval;
//       if ( limit.upperLimit<limit.lowerLimit )  excl = -1;
//       hExclusion->SetBinContent(ix,iy,excl);
//       hLowerLimit->SetBinContent(ix,iy,limit.lowerLimit);
//       hUpperLimit->SetBinContent(ix,iy,limit.upperLimit);

//       delete data;

//     }
//   }

//   TFile* out = new TFile("RA4abcd.root","RECREATE");
//   hExclusion->SetDirectory(out);
//   hExclusion->SetMinimum(); hExclusion->SetMaximum();
//   hExclusion->SetContour(1); hExclusion->SetContourLevel(0,0.5);
//   hLowerLimit->SetDirectory(out);
//   hLowerLimit->SetMinimum(); hLowerLimit->SetMaximum();
//   hUpperLimit->SetDirectory(out);
//   hUpperLimit->SetMinimum(); hUpperLimit->SetMaximum();
//   hYields[3]->SetDirectory(out);
//   hYields[3]->SetMinimum(); hYields[3]->SetMaximum();
//   out->Write();
//   delete out;
// }

