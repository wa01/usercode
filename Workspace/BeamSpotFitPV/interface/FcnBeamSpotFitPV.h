#ifndef FcnBeamSpotFitPV_h_
#define FcnBeamSpotFitPV_h_
#include "Workspace/BeamSpotFitPV/interface/BeamSpotFitPVData.h"
#include "Minuit2/FCNBase.h"

#include <vector> 


class FcnBeamSpotFitPV : public ROOT::Minuit2::FCNBase { 
public: 
  FcnBeamSpotFitPV(const std::vector<BeamSpotFitPVData>& data);
  ~FcnBeamSpotFitPV() {} 
  void setLimits (float xmin, float xmax,
		  float ymin, float ymax,
		  float zmin, float zmax) ;
  double Up() const {return errorDef_;} 
  double operator() (const std::vector<double>&) const; 
  unsigned int nrOfVerticesUsed () const;
private: 
  const std::vector<BeamSpotFitPVData>& data_;
  double errorDef_; 

  float lowerLimits_[3];
  float upperLimits_[3];
}; 
#endif
