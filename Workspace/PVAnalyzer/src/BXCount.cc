#include <memory>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <map>

class BXCount : public edm::EDAnalyzer {

public:
  
  explicit BXCount(const edm::ParameterSet &) {}
  ~BXCount() {}
  
private:
  
  virtual void beginJob() {}
  virtual void analyze(const edm::Event &, const edm::EventSetup &);
  virtual void endJob();
  
  std::map<int, unsigned int> bxCounts_;
};



void 
BXCount::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {
  int bx = iEvent.bunchCrossing();
  ++bxCounts_[bx];
}

void
BXCount::endJob() 
{
  std::cout << "BX count summary:" << std::endl;
  for ( std::map<int,unsigned int>::const_iterator i=bxCounts_.begin();
	i!=bxCounts_.end(); ++i )  std::cout << (*i).first << " " << (*i).second << std::endl;
}

#include "FWCore/Framework/interface/MakerMacros.h"

//define this as a plug-in
DEFINE_FWK_MODULE(BXCount);
