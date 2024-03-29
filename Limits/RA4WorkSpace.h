#ifndef RA4WorkSpace_H_
#define RA4WorkSpace_H_

/** Class to build and hold the RooStats work space for the two channels 
 *  of the RA4 ABCD method */

#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooAbsPdf.h"

class RA4WorkSpace {
public:
  /// channel definitions
  enum ChannelType { EleChannel, MuChannel };
  /// constructor
  RA4WorkSpace (const char* name = "wspace", bool constEff = false, bool constSCont = false,
		bool constKappa = false);
  /// destructor
  ~RA4WorkSpace () {delete wspace_;}
  //
  // creation of the work space
  //
  /// definition of a new channel (only 1/type and before finalization)
  void addChannel (ChannelType channel);
  /// definition of the final model (no addChannel after using this method)
  void finalize ();

  /// setting background-related numbers
  void setBackground (ChannelType channel, float bkgA, float bkgB, float bkgC, float bkgD);
  /// setting signal-related numbers
  void setSignal (ChannelType channel, 
		  float sigA, float sigB, float sigC, float sigD,
		  float effA, float effB, float effC, float effD);
  /// setting observations
  void setObserved (ChannelType channel, int obsA, int obsB, int obsC, int obsD);

  /// access to the work space
  RooWorkspace* workspace () {return wspace_;}
  
 private:
  void setValRange (const char* name, double val, double vmin=999999, double vmax=-999999);
  void setValRange (RooRealVar* var, double val, double vmin=999999, double vmax=-999999);

 private:
  /// the work space
  RooWorkspace* wspace_;

  bool constEff_;   //< no systematics on efficiency
  bool constSCont_; //< no systematics on signal contamination
  bool constKappa_; //< no systematics on kappa

  bool finalized_;  //< finalization flag
  bool hasEle_;     //< electron channel defined?
  bool hasMu_;      //< muon channel defined?

  //
  // variables in work space
  //
  // common variables
  RooRealVar* vS_;                      //< signal
  RooRealVar* vEffSys_;                 //< rel. uncertainty on the efficiency
  RooRealVar* vKappaSys_;               //< rel. uncertainty on kappa
  RooRealVar* vSContSys_;               //< rel. uncertainty on the signal contamination
  // variables / channel
  RooRealVar* vObs_[4][2];              //< observed values / region / channel
  RooRealVar* vEff_[2];                 //< efficiency*acceptance in D / channel
  RooRealVar* vKappa_[2];               //< kappa / channel
  RooRealVar* vSCont_[3][2];            //< rel. signal contamination for A,B,C / channel
  RooRealVar* vBxd_[2][2];              //< bkg. rel. to D for B,C / channel
  RooRealVar* vBkgd_[2];                //< background in D / channel
  //
  // sets of variables in work space
  //
  const RooArgSet* setObs_;             //< observables
  const RooArgSet* setNuis_;            //< nuisance parameters
  //
  // pdfs in work space
  //
  // common
  RooAbsPdf* pdfMcEff_;                 //< pseudo-measurement uncertainty efficiency
  RooAbsPdf* pdfMcKappa_;               //< pseudo-measurement uncertainty kappa
  RooAbsPdf* pdfMcSCont_;               //< pseudo-measurement uncertainty sig. cont.
  RooAbsPdf* pdfModel_;                 //< the full model
  // pdfs / channel
  RooAbsPdf* pdfReg_[4][2];             //< pdf / region / channel
};


#endif

