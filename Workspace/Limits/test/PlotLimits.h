//////////////////////////////////////////////////////////
// Extract observed and expected limits from combined 
//   higgs limits output file
// Standard sequence:
//   PlotLimits::Loop to fill the histograms
//   PlotLimits::drawHistograms() to draw histograms and
//     contour lines (automatically saves the Canvas
//     as pdf)
//  PlotLimits::saveContours() to save the contour lines
//     to a root file
//////////////////////////////////////////////////////////


#ifndef PlotLimits_h
#define PlotLimits_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TH2F.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include <string>
#include <iostream>
#include <vector>

class PlotLimits {
public:
  //
  // constructor from TFile and level (0.05 for "single points"
  //   runs where CLs is stored, 1.0 for all other cases
  //   where the limit on the signal strength is stored)
  //
  PlotLimits(TFile* file, float level=0.05, float relmax=5.);
  virtual ~PlotLimits();
  // 
  // loop over the tree and filling of histograms
  //
  virtual void     Loop();
  virtual void     Show(Long64_t entry = -1);
  //
  // drawing of all available histograms and contours
  //
  void drawHistograms ();
  //
  // saving contours to file (for use with limit plot program)
  //
  void saveContours();

private:
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual Bool_t   Notify();
  void drawHistogram (TH2* histogram, TGraph** graph = 0);
  std::vector<TGraph*> getContours (const char* name);
  std::vector<TGraph*> getContours (TH2* histogram);
  void saveCanvas () {
    if ( canvas_ ) {
      std::string cname = name_;
      cname += ".pdf";
      canvas_->SaveAs(cname.c_str());
    }
  }
  TH2* getHistogram (const char* name) {return (TH2*)gROOT->Get(name);}
  bool acceptGraph(TGraph* graph, TH2* histo);

private:
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        limit;
   Double_t        limitErr;
   Double_t        mh;
   Int_t           syst;
   Int_t           iToy;
   Int_t           iSeed;
   Int_t           iChannel;
   Float_t         t_cpu;
   Float_t         t_real;
   Float_t         quantileExpected;

   // List of branches
   TBranch        *b_limit;   //!
   TBranch        *b_limitErr;   //!
   TBranch        *b_mh;   //!
   TBranch        *b_syst;   //!
   TBranch        *b_iToy;   //!
   TBranch        *b_iSeed;   //!
   TBranch        *b_iChannel;   //!
   TBranch        *b_t_cpu;   //!
   TBranch        *b_t_real;   //!
   TBranch        *b_quantileExpected;   //!

  TCanvas* canvas_;

  double level_;
  double relmax_;

  TH2* hExist;
  TH2* hObs;
  TH2* hExpMinus2;
  TH2* hExpMinus1;
  TH2* hExpMedian;
  TH2* hExpPlus1;
  TH2* hExpPlus2;

  TGraph* gObs;
  TGraph* gExpMinus2;
  TGraph* gExpMinus1;
  TGraph* gExpMedian;
  TGraph* gExpPlus1;
  TGraph* gExpPlus2;

  std::string name_;

};

#endif

#ifdef PlotLimits_cxx
PlotLimits::PlotLimits(TFile* file, float level, float relmax) : 
  level_(level), relmax_(relmax) {
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.

  gStyle->SetOptStat(0);
  TTree* tree = (TTree*)file->Get("limit");
  Init(tree);
  
  TDirectory* curDir = gDirectory;
  gROOT->cd();

  //
  // histogram pointers (histograms will be created in Loop)
  //
  hExist = 0;
  hObs = 0;
  hExpMinus2 = 0;
  hExpMinus1 = 0;
  hExpMedian = 0;
  hExpPlus1 = 0;
  hExpPlus2 = 0;
  //
  // graphs holding the contours
  //
  gObs = new TGraph();
  gExpMinus2 = new TGraph();
  gExpMinus1 = new TGraph();
  gExpMedian = new TGraph();
  gExpPlus1 = new TGraph();
  gExpPlus2 = new TGraph();
  curDir->cd();

  canvas_ = 0;

  name_ = file->GetName();
  size_t ipos = name_.rfind('/');
  if ( ipos!=std::string::npos )  name_.erase(0,ipos+1);
  ipos = name_.find('.');
  if ( ipos!=std::string::npos )  name_.erase(ipos);
  std::cout << ":" << name_ << ":" << std::endl;

}

PlotLimits::~PlotLimits()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PlotLimits::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PlotLimits::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (!fChain->InheritsFrom(TChain::Class()))  return centry;
   TChain *chain = (TChain*)fChain;
   if (chain->GetTreeNumber() != fCurrent) {
      fCurrent = chain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void PlotLimits::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("limit", &limit, &b_limit);
   fChain->SetBranchAddress("limitErr", &limitErr, &b_limitErr);
   fChain->SetBranchAddress("mh", &mh, &b_mh);
   fChain->SetBranchAddress("syst", &syst, &b_syst);
   fChain->SetBranchAddress("iToy", &iToy, &b_iToy);
   fChain->SetBranchAddress("iSeed", &iSeed, &b_iSeed);
   fChain->SetBranchAddress("iChannel", &iChannel, &b_iChannel);
   fChain->SetBranchAddress("t_cpu", &t_cpu, &b_t_cpu);
   fChain->SetBranchAddress("t_real", &t_real, &b_t_real);
   fChain->SetBranchAddress("quantileExpected", &quantileExpected, &b_quantileExpected);
   Notify();
}

Bool_t PlotLimits::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PlotLimits::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PlotLimits::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PlotLimits_cxx
