#ifndef RA4WorkSpace_H_
#define RA4WorkSpace_H_

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooAbsPdf.h"


class RA4WorkSpace {
public:
  enum ChannelType { EleChannel, MuChannel };
  RA4WorkSpace (const char* name = "wspace");
  ~RA4WorkSpace () {delete wspace_;}

  RooWorkspace* workspace () {return wspace_;}

  void addChannel (ChannelType channel);

  void setBackground (float bkgA, float bkgB, float bkgC, float bkgD);
  void setSignal (float sigA, float sigB, float sigC, float sigD);


private:

private:
  RooWorkspace* wspace_;

  bool hasEle_;
  bool hasMu_;

  RooRealVar* vS_;
  RooRealVar* vEffsys_;
  RooRealVar* vKappasys_;
  RooRealVar* vSContsys_;
  RooAbsPdf* pdfMcKappa_;
  RooAbsPdf* pdfMcSCont_;
  RooAbsPdf* pdfModel_;
  const RooArgSet* setObs_;
  const RooArgSet* setNuis_;
  RooRealVar* vObs_[4][2];
  RooRealVar* vEff_[2];
  RooRealVar* vKappa_[2];
  RooRealVar* vSCont_[3][2];
  RooRealVar* vBxd_[2][2];
  RooRealVar* vBkgd_[2];
  RooAbsReal* pSd_[2];
  RooAbsPdf* pdfReg_[4][2];
};


#endif

