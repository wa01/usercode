#ifndef HtMetTreeReader_h_
#define HtMetTreeReader_h_

#include "TFile.h"
#include "TTree.h"

#include <vector>
#include <algorithm>

class HtMetTreeReader {
public:
  HtMetTreeReader (const char* filename);
  HtMetTreeReader (TFile* file);
  ~HtMetTreeReader ();
  void loop (std::vector< std::pair<float,float> >& data);
private:
  void init (TFile* file);
  
private:
  TFile* file_;
  TTree* tree_;

  TBranch* b_htmet_;
  float ht_;
  float met_;
};
#endif
