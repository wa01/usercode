#include "TFile.h"
#include "TDirectory.h"
#include "TKey.h"
#include "TGraphErrors.h"
#include "TH1.h"

#include <iostream>
#include <vector>
#include <map>
#include <string>

struct FitResult {
  FitResult () : ls(0.), els(0.) {
    for ( unsigned int i=0; i<9; ++i ) {
      values[i] = errors[i] = 0.;
    }
  }
  bool operator< (const FitResult& other) const {
    return ls<other.ls;
  }
  float ls;
  float els;
  float values[9];
  float errors[9];
};
typedef std::vector<FitResult> FitResults;

struct RunResult {
  FitResults fitResults;
  std::map<unsigned int, unsigned int> pvCountMap;
};

const unsigned int VariableSize = 9;
char* GraphNames[VariableSize] = { "x", "y", "z", "ex", "ey", "ez", "corrxy", "dxdz", "dydz" };

std::vector<std::string> histoFiles()
{
  std::vector<std::string> result;
  result.push_back("crab_0_100407_122314/res/histo_10_0.root");
  result.push_back("crab_0_100407_122314/res/histo_11_0.root");
  result.push_back("crab_0_100407_122314/res/histo_12_0.root");
  result.push_back("crab_0_100407_122314/res/histo_13_0.root");
  result.push_back("crab_0_100407_122314/res/histo_14_0.root");
  result.push_back("crab_0_100407_122314/res/histo_15_0.root");
  result.push_back("crab_0_100407_122314/res/histo_16_0.root");
  result.push_back("crab_0_100407_122314/res/histo_17_0.root");
  result.push_back("crab_0_100407_122314/res/histo_18_0.root");
  result.push_back("crab_0_100407_122314/res/histo_19_0.root");
  result.push_back("crab_0_100407_122314/res/histo_1_1.root");
  result.push_back("crab_0_100407_122314/res/histo_20_0.root");
  result.push_back("crab_0_100407_122314/res/histo_21_0.root");
  result.push_back("crab_0_100407_122314/res/histo_22_0.root");
  result.push_back("crab_0_100407_122314/res/histo_23_0.root");
  result.push_back("crab_0_100407_122314/res/histo_24_0.root");
  result.push_back("crab_0_100407_122314/res/histo_25_0.root");
  result.push_back("crab_0_100407_122314/res/histo_27_0.root");
  result.push_back("crab_0_100407_122314/res/histo_28_0.root");
  result.push_back("crab_0_100407_122314/res/histo_29_0.root");
  result.push_back("crab_0_100407_122314/res/histo_2_1.root");
  result.push_back("crab_0_100407_122314/res/histo_30_0.root");
  result.push_back("crab_0_100407_122314/res/histo_31_0.root");
  result.push_back("crab_0_100407_122314/res/histo_32_0.root");
  result.push_back("crab_0_100407_122314/res/histo_33_0.root");
  result.push_back("crab_0_100407_122314/res/histo_34_0.root");
  result.push_back("crab_0_100407_122314/res/histo_35_0.root");
  result.push_back("crab_0_100407_122314/res/histo_36_0.root");
  result.push_back("crab_0_100407_122314/res/histo_37_0.root");
  result.push_back("crab_0_100407_122314/res/histo_38_0.root");
  result.push_back("crab_0_100407_122314/res/histo_39_0.root");
  result.push_back("crab_0_100407_122314/res/histo_3_1.root");
  result.push_back("crab_0_100407_122314/res/histo_40_0.root");
  result.push_back("crab_0_100407_122314/res/histo_42_0.root");
  result.push_back("crab_0_100407_122314/res/histo_43_0.root");
  result.push_back("crab_0_100407_122314/res/histo_44_0.root");
  result.push_back("crab_0_100407_122314/res/histo_45_0.root");
  result.push_back("crab_0_100407_122314/res/histo_46_0.root");
  result.push_back("crab_0_100407_122314/res/histo_47_0.root");
  result.push_back("crab_0_100407_122314/res/histo_48_0.root");
  result.push_back("crab_0_100407_122314/res/histo_4_1.root");
  result.push_back("crab_0_100407_122314/res/histo_50_0.root");
  result.push_back("crab_0_100407_122314/res/histo_51_0.root");
  result.push_back("crab_0_100407_122314/res/histo_52_0.root");
  result.push_back("crab_0_100407_122314/res/histo_53_0.root");
  result.push_back("crab_0_100407_122314/res/histo_54_0.root");
  result.push_back("crab_0_100407_122314/res/histo_55_0.root");
  result.push_back("crab_0_100407_122314/res/histo_56_0.root");
  result.push_back("crab_0_100407_122314/res/histo_57_0.root");
  result.push_back("crab_0_100407_122314/res/histo_58_0.root");
  result.push_back("crab_0_100407_122314/res/histo_59_0.root");
  result.push_back("crab_0_100407_122314/res/histo_5_1.root");
  result.push_back("crab_0_100407_122314/res/histo_60_0.root");
  result.push_back("crab_0_100407_122314/res/histo_61_0.root");
  result.push_back("crab_0_100407_122314/res/histo_62_0.root");
  result.push_back("crab_0_100407_122314/res/histo_63_0.root");
  result.push_back("crab_0_100407_122314/res/histo_64_0.root");
  result.push_back("crab_0_100407_122314/res/histo_65_0.root");
  result.push_back("crab_0_100407_122314/res/histo_66_0.root");
  result.push_back("crab_0_100407_122314/res/histo_67_0.root");
  result.push_back("crab_0_100407_122314/res/histo_68_0.root");
  result.push_back("crab_0_100407_122314/res/histo_69_0.root");
  result.push_back("crab_0_100407_122314/res/histo_6_1.root");
  result.push_back("crab_0_100407_122314/res/histo_70_0.root");
  result.push_back("crab_0_100407_122314/res/histo_71_0.root");
  result.push_back("crab_0_100407_122314/res/histo_72_0.root");
  result.push_back("crab_0_100407_122314/res/histo_8_1.root");
  return result;
}

class MergeBeamSpotGraphs {
public:
  MergeBeamSpotGraphs (std::vector<std::string>& fileNames) :
    fileNames_(fileNames) {}
  ~MergeBeamSpotGraphs () {
    for ( std::map<unsigned int, RunResult>::iterator ir=runResultMap_.begin();
	  ir!=runResultMap_.end(); ++ir ) {
      std::cout << "Results for run " << ir->first 
		<< " : found " << ir->second.fitResults.size() 
		<< " fit points / " << ir->second.pvCountMap.size() 
		<< " pv counts " << std::endl;
    }
  }

  void merge ();
  void readFile (TFile& file);
  void readResults(TDirectory* dir,unsigned int run,
		   RunResult& result);
  void writeResults (const char* file);

private:
  std::vector<std::string> fileNames_;
  std::map<unsigned int, RunResult> runResultMap_;
};

void
MergeBeamSpotGraphs::merge ()
{
  for ( unsigned int i=0; i<fileNames_.size(); ++i ) {
    std::cout << "Opening " << fileNames_[i] << std::endl;
    TFile file(fileNames_[i].c_str());
    if ( file.IsZombie() || !file.IsOpen() )  continue;
    readFile(file);
    std::cout << "Back in merge" << std::endl;
  }
}

void
MergeBeamSpotGraphs::readFile (TFile& file)
{
  std::cout << "Read file " << file.GetName() << std::endl;
  TDirectory* module = (TDirectory*)file.FindObjectAny("beamSpotFitPV");
  if ( module==0 )  return;

  std::cout << "module = " << module << endl;
  module->ls();
  TIter iter(module->GetListOfKeys());
  std::cout << "Nr. of keys = " << module->GetListOfKeys()->GetSize() << std::endl;
  TKey* key;
  unsigned int irun;
  while ( (key=(TKey*)iter()) ) {
    std::cout << "  Using dir " << key->GetName() << std::endl;
    sscanf(key->GetName(),"Run%d",&irun);
    RunResult& result = runResultMap_[irun];
    readResults((TDirectory*)module->Get(key->GetName()),irun,result);
    std::cout << "Back in readFile" << std::endl;
  }
} 

void 
MergeBeamSpotGraphs::readResults(TDirectory* dir, unsigned int run,
				 RunResult& result)
{
  std::cout << "Getting pvcounts" << std::endl;
  TH1* hcount = (TH1*)dir->Get("pvcounts");
  if ( hcount==0 ) {
    std::cout << "No pvcounts histogram!!!" << endl;
    return;
  }
  unsigned int ls;
  TAxis* axis = hcount->GetXaxis();
  for ( unsigned int i=1; i<=hcount->GetNbinsX(); ++i ) {
    sscanf(axis->GetBinLabel(i),"%d",&ls);
//     if ( result.pvCountMap.find(ls)!=result.pvCountMap.end() )  
//       std::cout << "Found duplicate LS " << ls << " for run " << run << std::endl;
    result.pvCountMap[ls] += hcount->GetBinContent(i);
  }

  int np(-1);
  TGraphErrors* graphs[VariableSize];
  for ( unsigned int iv=0; iv<VariableSize; ++iv ) {
    graphs[iv] = (TGraphErrors*)dir->Get(GraphNames[iv]);
    if ( graphs[iv]==0 ) {
      std::cout << "No graph for  " << GraphNames[iv] << std::endl;
      return;
    }
    if ( iv==0 ) {
      np = graphs[iv]->GetN();
    }
    else {
      if ( np!=graphs[iv]->GetN() ) {
	std::cout << "Inconsistent number of points between graphs" << std::endl;
	return;
      }
    }
  }

  double x,ex,y,ey;
  for ( unsigned int ip=0; ip<np; ++ip ) {
    std::cout << "Getting results for point " << ip << std::endl;
    FitResult fitResult;
    for ( unsigned int iv=0; iv<VariableSize; ++iv ) {
      graphs[iv]->GetPoint(ip,x,y);
      if ( iv==0 ) {
	fitResult.ls = x;
	fitResult.els = graphs[iv]->GetErrorX(ip);
      }
      fitResult.values[iv] = y;
      fitResult.errors[iv] = graphs[iv]->GetErrorY(ip);
    }
    result.fitResults.push_back(fitResult);
  }
  
  std::cout << "Done reading results" << std::endl;
}

void
MergeBeamSpotGraphs::writeResults (const char* filename)
{
  TFile* file = new TFile(filename,"NEW");

  char title[64];
  for ( std::map<unsigned int, RunResult>::iterator ir=runResultMap_.begin();
	ir!=runResultMap_.end(); ++ir ) {
    unsigned int irun = ir->first;
    RunResult& runResult = ir->second;
    std::cout << "Results for run " << irun 
	      << " : found " << runResult.fitResults.size() 
	      << " fit points / " << runResult.pvCountMap.size() 
	      << " pv counts " << std::endl;
    sprintf(title,"Run%d",irun);
    TDirectory* dir = file->mkdir(title);
    dir->cd();

    unsigned int npv = runResult.pvCountMap.size();
    unsigned int ls1 = runResult.pvCountMap.begin()->first;
    unsigned int ls2 = runResult.pvCountMap.rbegin()->first;
    TH1* h_count = new TH1F("pvcounts","Nr. of selected primary vertices",ls2-ls1+1,ls1-0.5,ls2+0.5);
    for ( std::map<unsigned int, unsigned int>::iterator i=runResult.pvCountMap.begin();
	  i!=runResult.pvCountMap.end(); ++i )  h_count->Fill(i->first,i->second);
    h_count->Write();

    TGraphErrors* graphs[VariableSize];
    for ( unsigned int iv=0; iv<VariableSize; ++iv ) {
      graphs[iv] = new TGraphErrors();
      graphs[iv]->SetName(GraphNames[iv]);
    }

    std::vector<FitResult>& fitResults = runResult.fitResults;
    for ( unsigned int i=0; i<fitResults.size(); ++i ) {
      for ( unsigned int iv=0; iv<VariableSize; ++iv ) {
	graphs[iv]->SetPoint(i,fitResults[i].ls,fitResults[i].values[iv]);
	graphs[iv]->SetPointError(i,fitResults[i].els,fitResults[i].errors[iv]);
      }
    }
    for ( unsigned int iv=0; iv<VariableSize; ++iv )  graphs[iv]->Write();

    dir->Write();

  }
  file->Write();
  
}
