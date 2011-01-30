#ifndef RA4abcd_h_
#define RA4abcd_h_

#include "RooWorkspace.h"
#include "RooDataSet.h"

/** Some functions for limit calculation with ABCD.
 *  RA4Single: one test with LM0
 *  RA4Mult:   scan over m0-m1/2 (signal yields in ABCD provided by histograms) */

// using namespace RooFit;
// using namespace RooStats;

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
MyLimit computeLimit (RooWorkspace* wspace, RooDataSet* data, StatMethod method, bool draw=false);
//
// set value and range for a variable in the workspace
//
void setValRange (RooWorkspace* workspace, const char* name, double val, double vmin=999999, double vmax=-999999);
//
// create a workspace with the variables (with dummy values) and the pdfs
//
RooWorkspace* createWorkspace (const char* name = "wspace");
//
// set background-related variables of the workspace
//
void setBackgrounds (RooWorkspace* wspace, double* bkgs=0);

void setSignal (RooWorkspace* wspace, double* lm_mc);
//
// single measurement (LM0 or LM1)
//
void RA4Single (StatMethod method, double* sig=0, double* bkg=0);
//
// scan over parameter space
//
void RA4Mult (const char* file, StatMethod method=ProfileLikelihoodMethod);

#endif
