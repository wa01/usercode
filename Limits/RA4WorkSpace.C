#include "RA4WorkSpace.h"

#include "TROOT.h"
#include "RooAbsPdf.h"
#include "RooProdPdf.h"
#include "RooUniform.h"
#include "RooWorkspace.h"
#include "RooDataSet.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "RooPoisson.h"

using namespace RooFit;
// using namespace RooStats;

RA4WorkSpace::RA4WorkSpace (const char* name, bool constEff, bool constSCont, bool constKappa) :
  wspace_(new RooWorkspace(name)), 
  constEff_(constEff), constSCont_(constSCont), constKappa_(constKappa),
  finalized_(false), hasEle_(false), hasMu_(false) 
{
  //
  // definition of common variables, sets and pdfs
  //
  // total signal in D
  wspace_->factory("s[1,0,100]");
  vS_ = wspace_->var("s");
  // systematic uncertainty on efficiency
  wspace_->factory("effSys[1.,0,2]");
  vEffSys_ = wspace_->var("effSys");
  if ( constEff_ )  vEffSys_->setConstant(true);
  // systematic uncertainty on kappa
  wspace_->factory("kappaSys[1,0,2]");
  vKappaSys_ = wspace_->var("kappaSys");
  if ( constKappa_ )  vKappaSys_->setConstant(true);
  // systematic uncertainty on signal contamination
  wspace_->factory("scontSys[1,0,2]");
  vSContSys_ = wspace_->var("scontSys");
  if ( constSCont_ )  vSContSys_->setConstant(true);
  // pseudo-measurements for systematics on efficiency,
  // kappa and signal contamination syst
  if ( !constEff_ ) {
    wspace_->factory("effScale[1]");
    wspace_->factory("sigmaEff[0.15]");
  }
  if ( !constKappa_ ) {
    wspace_->factory("kappaScale[1]");
    wspace_->factory("sigmaKappa[0.15]");
  }
  if ( !constSCont_ ) {
    wspace_->factory("scontScale[1]");
    wspace_->factory("sigmaScont[0.15]");
  }
  // Pdfs for pseudo-measurements
  pdfMcEff_ = 0;
  if ( !constEff_ ) {
    wspace_->factory("Gaussian::mcEff(effScale, effSys, sigmaEff)");
    pdfMcEff_ = wspace_->pdf("mcEff");
  }
  pdfMcKappa_ = 0;
  if ( !constKappa_ ) {
    wspace_->factory("Gaussian::mcKappa(kappaScale, kappaSys, sigmaKappa)");
    pdfMcKappa_ = wspace_->pdf("mcKappa");
  }
  pdfMcSCont_ = 0;
  if ( !constSCont_ ) {
    wspace_->factory("Gaussian::mcScont(scontScale, scontSys, sigmaScont)");
    pdfMcSCont_ = wspace_->pdf("mcScont");
  }
  // sets
  wspace_->defineSet("poi","s");
  wspace_->defineSet("obs","");
  if ( !constEff_ )  wspace_->extendSet("obs","effScale");
  if ( !constKappa_ )  wspace_->extendSet("obs","kappaScale");
  if ( !constSCont_ )  wspace_->extendSet("obs","scontScale");
  setObs_ = wspace_->set("obs");
  wspace_->defineSet("nuis","");
  if ( !constEff_ )  wspace_->extendSet("nuis","effSys");
  if ( !constKappa_ )  wspace_->extendSet("nuis","kappaSys");
  if ( !constSCont_ )  wspace_->extendSet("nuis","scontSys");
  setNuis_ = wspace_->set("nuis");

  // wspace_->Print("v");
}

//
// set value and range for a variable in the workspace
//
void 
RA4WorkSpace::setValRange (const char* name, double val, double vmin, double vmax)
{
  RooRealVar* var = wspace_->var(name);
  setValRange(var,val,vmin,vmax);
}

//
// set value and range for a variable in the workspace
//
void 
RA4WorkSpace::setValRange (RooRealVar* var, double val, double vmin, double vmax)
{
  if ( var ) {
    if ( vmax>vmin )  var->setRange(vmin,vmax);
    var->setVal(val);
  }
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
    vEff_[channel]->setConstant(true);
  }
  //
  // kappa
  //
  {
    std::string name("kappa");
    name += suffix;
//     RooRealVar tmp(name.c_str(),name.c_str(),1.,0,2);
    RooRealVar tmp(name.c_str(),name.c_str(),1.);
    wspace_->import(tmp);
    vKappa_[channel] = wspace_->var(name.c_str());
    vKappa_[channel]->setConstant(true);
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
//     RooRealVar tmp(name.c_str(),name.c_str(),0,0,10);
    RooRealVar tmp(name.c_str(),name.c_str(),0);
    wspace_->import(tmp);
    vSCont_[i][channel] = wspace_->var(name.c_str());
    vSCont_[i][channel]->setConstant(true);
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
    RooArgSet tmpSet(*wspace_->var("s"),*vEff_[channel],*wspace_->var("effSys"));
    RooProduct tmp(name.c_str(),name.c_str(),tmpSet);
    wspace_->import(tmp);
    // pSd_[channel] = (RooProduct*)wspace_->var(name.c_str());
  }
  //
  // region A
  //
  {
    std::string name1("sa");
    name1 += suffix;
    RooArgSet tmpSetSa(*wspace_->var("s"),*vEff_[channel],*vSCont_[0][channel],*wspace_->var("scontSys"));
    RooProduct tmpSa(name1.c_str(),name1.c_str(),tmpSetSa);
    wspace_->import(tmpSa);
    std::string name2("bkga");
    name2 += suffix;
    RooArgSet tmpSetBkga(*vBkgd_[channel],*vBxd_[0][channel],*vBxd_[1][channel],*vKappa_[channel],*wspace_->var("kappaSys"));
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
    RooArgSet tmpSetSb(*wspace_->var("s"),*vEff_[channel],*vSCont_[1][channel],*wspace_->var("scontSys"));
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
    RooArgSet tmpSetSc(*wspace_->var("s"),*vEff_[channel],*vSCont_[2][channel],*wspace_->var("scontSys"));
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
  RooArgSet modelSet;
  if ( !constEff_ )  modelSet.add(*pdfMcEff_);
  if ( !constKappa_ )  modelSet.add(*pdfMcKappa_);
  if ( !constSCont_ )  modelSet.add(*pdfMcSCont_);
  // RooArgSet modelSet(*pdfMcEff_,*pdfMcKappa_,*pdfMcSCont_);
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
  wspace_->factory("PROD::prior(prior_poi,prior_nuis)"); 
  wspace_->Print("v");

  finalized_ = true;
}


void 
RA4WorkSpace::setBackground (ChannelType channel,
			     float bkgA, float bkgB, float bkgC, float bkgD)
{
  //
  // set background-related variables of the workspace
  //
  if ( (channel==EleChannel && !hasEle_) ||
       (channel==MuChannel && !hasMu_) ) {
    std::cout << "Channel " << channel << " has not been defined" << std::endl;
    return;
  }
  //
  std::string suffix;
  switch ( channel ) {
  case EleChannel:
    suffix = "E";
    break;
  case MuChannel:
    suffix = "M";
    break;
  }
  //
  // background numbers
  //
  // vBkgd_[channel]->setVal(bkgD);
  // vBxd_[0][channel]->setVal(bkgB/bkgD);
  // vBxd_[1][channel]->setVal(bkgC/bkgD);
  setValRange(vBkgd_[channel],bkgD,0,10*bkgD);
  setValRange(vBxd_[0][channel],bkgB/bkgD,0,5*bkgB/bkgD);
  setValRange(vBxd_[1][channel],bkgC/bkgD,0,5*bkgC/bkgD);
  //
  // kappa
  //
  // vKappa_[channel]->setVal((bkgA*bkgD)/(bkgB*bkgC));
//   setValRange(vKappa_[channel],(bkgA*bkgD)/(bkgB*bkgC),0,2);
  setValRange(vKappa_[channel],(bkgA*bkgD)/(bkgB*bkgC));
  vKappa_[channel]->Print("v");
  vKappa_[channel]->printValue(std::cout);
  //
  // set observed values to expectations
  // use pessimistic scenario for rounding of expected numbers
  //  
  setObserved(channel,int(bkgA)+1,int(bkgB),int(bkgC),int(bkgD)+1);
}

void 
RA4WorkSpace::setSignal (ChannelType channel,
			 float sigA, float sigB, float sigC, float sigD)
{
 //
  // set signal-related variables of the workspace
  //
  if ( (channel==EleChannel && !hasEle_) ||
       (channel==MuChannel && !hasMu_) ) {
    std::cout << "Channel " << channel << " has not been defined" << std::endl;
    return;
  }
  //
  std::string suffix;
  switch ( channel ) {
  case EleChannel:
    suffix = "E";
    break;
  case MuChannel:
    suffix = "M";
    break;
  }

  // Inputs : expected values

  // relative uncertainties on relative signal contamination
  // add explicit k-factor uncertainty
  // double sigma_sad_rel = 0.10;
  // double sigma_sbd_rel = 0.10;
  // double sigma_scd_rel = 0.10;
  // double sigma_eff = 0.05;

  // derived quantities
  double sad_mc = max(sigA/sigD,float(0.01));
  double sbd_mc = max(sigB/sigD,float(0.01));
  double scd_mc = max(sigC/sigD,float(0.01));
  vSCont_[0][channel]->setVal(sad_mc);
  vSCont_[1][channel]->setVal(sbd_mc);
  vSCont_[2][channel]->setVal(scd_mc);
  vSCont_[0][channel]->setConstant(true);
  vSCont_[1][channel]->setConstant(true);
  vSCont_[2][channel]->setConstant(true);

  // .. background and signal variables
  // vS_->setRange(0,10*sigD);
  // vS_->setVal(sigD);
//   setValRange(vS_,sigD,0,10*sigD);
  setValRange(vS_,sigD,0,100);
 }

void 
RA4WorkSpace::setObserved (ChannelType channel,
			   int obsA, int obsB, int obsC, int obsD)
{
  //
  // set background-related variables of the workspace
  //
  if ( (channel==EleChannel && !hasEle_) ||
       (channel==MuChannel && !hasMu_) ) {
    std::cout << "Channel " << channel << " has not been defined" << std::endl;
    return;
  }
  //
  for ( unsigned int i=0; i<4; ++i ) 
    vObs_[i][channel]->setRange(0,1000);
  vObs_[0][channel]->setVal(obsA);
  vObs_[1][channel]->setVal(obsB);
  vObs_[2][channel]->setVal(obsC);
  vObs_[3][channel]->setVal(obsD);
}
