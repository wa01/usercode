#ifndef RA4WorkingPoint_h_
#define RA4WorkingPoint_h_

#include <string>

struct RA4WorkingPoint {
  RA4WorkingPoint () : valid_(false) {}
  RA4WorkingPoint (float bkgA, float bkgB, float bkgC, float bkgD, float sigKappa) :
    valid_(true), sigKappa_(sigKappa) {
    bkg_[0] = bkgA; bkg_[1] = bkgB; bkg_[2] = bkgC; bkg_[3] = bkgD;
    obs_[0] = obs_[1] = obs_[2] = obs_[3] = -1;
  }
  RA4WorkingPoint (float bkgA, float bkgB, float bkgC, float bkgD, float sigKappa,
		   float effA, float effB, float effC, float effD, float sigEff) :
    valid_(true), sigKappa_(sigKappa), sigEff_(sigEff) {
    bkg_[0] = bkgA; bkg_[1] = bkgB; bkg_[2] = bkgC; bkg_[3] = bkgD;
    eff_[0] = effA; eff_[1] = effB; eff_[2] = effC; eff_[3] = effD;
    obs_[0] = obs_[1] = obs_[2] = obs_[3] = -1;
  }
  RA4WorkingPoint (float bkgA, float bkgB, float bkgC, float bkgD, float sigKappa,
		   const char* prefix, const char* yields, const char* kfactors,
		   float sigEff) :
    valid_(true), sigKappa_(sigKappa), sigEff_(sigEff),
    prefix_(prefix), yields_(yields), kfactors_(kfactors) {
    bkg_[0] = bkgA; bkg_[1] = bkgB; bkg_[2] = bkgC; bkg_[3] = bkgD;
    obs_[0] = obs_[1] = obs_[2] = obs_[3] = -1;
  }
  void setExpectedPessimistic () {
    obs_[0] = int(bkg_[0]) + 1;
    obs_[1] = int(bkg_[1]);
    obs_[2] = int(bkg_[2]);
    obs_[3] = int(bkg_[3]) + 1;
  }
  void setExpectedOptimistic () {
    obs_[0] = int(bkg_[0]);
    obs_[1] = int(bkg_[1]) + 1;
    obs_[2] = int(bkg_[2]) + 1;
    obs_[3] = int(bkg_[3]);
  }
  void setObserved (int obsA, int obsB, int obsC, int obsD) {
    obs_[0] = obsA; obs_[1] = obsB; obs_[2] = obsC; obs_[3] = obsD;
  }
  bool valid_;
  int obs_[4];
  double bkg_[4];
  double sigKappa_;
  double eff_[4];
  double sigEff_;
  std::string prefix_;
  std::string yields_;
  std::string kfactors_;
};
#endif
