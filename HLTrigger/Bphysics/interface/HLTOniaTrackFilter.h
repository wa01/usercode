#ifndef HLTOniaTrackFilter_h
#define HLTOniaTrackFilter_h

#include "HLTrigger/HLTcore/interface/HLTFilter.h"

class HLTOniaTrackFilter : public HLTFilter {
  
public:
  explicit HLTOniaTrackFilter(const edm::ParameterSet&);
  ~HLTOniaTrackFilter() {}
  virtual bool filter(edm::Event&, const edm::EventSetup&);

private:
  edm::InputTag candTag_; 
  bool saveTag_;
  double minPt_;
  double minP_;
  double maxEta_;
};

#endif 
