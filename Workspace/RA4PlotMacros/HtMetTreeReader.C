#include "HtMetTreeReader.h"

#include <iostream>
using namespace std;

HtMetTreeReader::HtMetTreeReader (const char* filename)
{
  file_ = new TFile(filename);
  init(file_);
}

HtMetTreeReader::HtMetTreeReader (TFile* file) : file_(0) {
  init(file);
}

HtMetTreeReader::~HtMetTreeReader ()
{
  delete tree_;
  delete file_;
}

void
HtMetTreeReader::init (TFile* file)
{
  tree_ = (TTree*)file->FindObjectAny("HTMETtree");
  if ( tree_ ) {
    tree_->SetBranchAddress("htmet",&ht_,&b_htmet_);
  }
}

void
HtMetTreeReader::loop (vector< pair<float,float> >& data) 
{
  if ( tree_ == 0 )  return;

  data.clear();
  pair<float,float> p;

  Long64_t nevt = tree_->GetEntriesFast();
  for ( Long64_t i=0; i<nevt; ++i ) {
    tree_->GetEntry(i);
    p.first = ht_; p.second = met_;
    data.push_back(p);
  }
}
