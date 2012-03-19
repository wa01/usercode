#define PlotLimits_cxx
#include "PlotLimits.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TGraph.h"
#include "TMarker.h"

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <set>
#include <cmath>
#include <string>

using namespace std;


void PlotLimits::Loop()
{
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();
   cout << nentries << endl;
   //
   // first loop to determine range and spacing of the
   //   bins in the two coordinates
   //
   set<int> m0s;
   set<int> m12s;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //
      // mass pair is encoded in the seed!
      //
      m0s.insert(iSeed/10000);
      m12s.insert(iSeed%10000);
   }
   //
   // extraction of ranges and bin widths
   //
   if ( m0s.size()<2 || m12s.size()<2 ) {
     cout << "too few points" << endl;
     return;
   }
   // for M0
   int dm0(0);
   int m0 = *m0s.begin();
   int m0min(m0); int m0max(m0);
   set<int>::const_iterator im = m0s.begin(); ++im;
   for ( ; im!=m0s.end(); ++im ) {
     int dm = *im - m0;
     if ( dm0==0 || abs(dm)<dm0 )  dm0 = abs(dm);
     m0 = *im;
     if ( m0<m0min )  m0min = m0;
     if ( m0>m0max )  m0max = m0;
   }
   // and for M12
   int dm12(0);
   int m12 = *m12s.begin();
   int m12min(m12); int m12max(m12);
   im = m12s.begin(); ++im;
   for ( ; im!=m12s.end(); ++im ) {
     int dm = *im - m12;
     if ( dm12==0 || abs(dm)<dm12 )  dm12 = abs(dm);
     m12 = *im;
     if ( m12<m12min )  m12min = m12;
     if ( m12>m12max )  m12max = m12;
   }
   if ( dm0==0 || dm12==0 ) {
     cout << "dm0 or dm12 == 0 ??" << endl;
     return;
   }
   //
   // creating histograms (grid points at bin centre)
   //
   int nb0 = (m0max-m0min)/dm0 + 1;
   if ( (m0max-m0min)%dm0 )  ++nb0;
   float fm0min = m0min - dm0/2.;
   float fm0max = fm0min + nb0*dm0;
   int nb12 = (m12max-m12min)/dm12 + 1;
   if ( (m12max-m12min)%dm12 )  ++nb12;
   float fm12min = m12min - dm12/2.;
   float fm12max = fm12min + nb12*dm12;
   cout << "m0 " << " " << dm0 << " " << m0min << " " << m0max << " ; "
	<< nb0 << " " << fm0min << " " << fm0max << endl;
   cout << "m12 " << " " << dm12 << " " << m12min << " " << m12max << " ; "
	<< nb0 << " " << fm12min << " " << fm12max << endl;
   hObs.first = new TH2F("hObs","hObs",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hExpMinus2.first = new TH2F("hExpMinus2","hExpMinus2",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hExpMinus1.first = new TH2F("hExpMinus1","hExpMinus1",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hExpMedian.first = new TH2F("hExpMedian","hExpMedian",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hExpPlus1.first = new TH2F("hExpPlus1","hExpPlus1",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hExpPlus2.first = new TH2F("hExpPlus2","hExpPlus2",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hObs.second = new TH2F("hObsExist","hObsExist",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hExpMinus2.second = new TH2F("hExpMinus2Exist","hExpMinus2Exist",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hExpMinus1.second = new TH2F("hExpMinus1Exist","hExpMinus1Exist",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hExpMedian.second = new TH2F("hExpMedianExist","hExpMedianExist",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hExpPlus1.second = new TH2F("hExpPlus1Exist","hExpPlus1Exist",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hExpPlus2.second = new TH2F("hExpPlus2Exist","hExpPlus2Exist",nb0,fm0min,fm0max,nb12,fm12min,fm12max);

   //
   // second loop: filling of histograms
   //
   nbytes = 0;
   nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      float m0 = iSeed/10000;
      float m12 = iSeed%10000;

      TH2* hLim(0);
      TH2* hExist(0);
      // observed
      if ( quantileExpected<0.) {
	hLim = hObs.first;
	hExist = hObs.second;
      }
      // -2 sigma
      else if ( fabs(quantileExpected-0.025)<0.001 ) {
	hLim = hExpMinus2.first;
	hExist = hExpMinus2.second;
      }
      // -1 sigma
      else if ( fabs(quantileExpected-0.16)<0.01 ) {
	hLim = hExpMinus1.first;
	hExist = hExpMinus1.second;
      }
      // median
      else if ( fabs(quantileExpected-0.50)<0.01 ) {
	hLim = hExpMedian.first;
	hExist = hExpMedian.second;
      }
      // +1 sigma
      else if ( fabs(quantileExpected-0.84)<0.01 ) {
	hLim = hExpPlus1.first;
	hExist = hExpPlus1.second;
      }
      // +2 sigma
      else if ( fabs(quantileExpected-0.975)<0.001 ) {
	hLim = hExpPlus2.first;
	hExist = hExpPlus2.second;
      }
      if ( hExist ) {
	int ibin = hExist->FindBin(m0,m12);
	if ( hExist->GetBinContent(ibin) > 0.5 ) {
	  std::cout << "***** Duplicate entry for " << m0 << " "
		    << m12 << " " << hExist->GetTitle() << std::endl;
	}
	hLim->Fill(m0,m12,limit);
	hExist->Fill(m0,m12);
      }
   }
}

void 
PlotLimits::drawHistograms()
{
  vector<TGraph*> contours;

  canvas_ = new TCanvas("c","c",1000,1000);
  canvas_->Divide(2,3);
  canvas_->cd(1);
  drawHistogram(hExpMinus2,&gExpMinus2);
  canvas_->cd(2);
  drawHistogram(hExpMinus1,&gExpMinus1);
  canvas_->cd(3);
  drawHistogram(hObs,&gObs);
  canvas_->cd(4);
  drawHistogram(hExpMedian,&gExpMedian);
  canvas_->cd(5);
  drawHistogram(hExpPlus1,&gExpPlus1);
  canvas_->cd(6);
  drawHistogram(hExpPlus2,&gExpPlus2);

  string cname = name_;
  cname += "-histos.pdf";
  canvas_->SaveAs(cname.c_str());

//   TCanvas* cExists = new TCanvas("cexist","cexist");
//   hExist->Draw("zcol");

}


void 
PlotLimits::drawHistogram (std::pair<TH2*,TH2*> histos, TGraph** graph)
{
  if ( histos.first->GetEntries()==0 )  return;

  double levels[1] = { level_ };

  histos.first->SetMaximum(relmax_*level_);
  histos.first->Draw("zcol");
  //
  // need to draw and update the pad in order to have access
  //   to the contour lines ...
  //
  TH2* hc = (TH2*)histos.first->Clone();
  hc->SetContour(1,levels);
  hc->Draw("cont3 list same");
  gPad->Update();
  //
  // extract contour lines
  //
  TVirtualPad* curPad = gPad;
  vector<TGraph*> contours = getContours(histos.first,histos.second);
  curPad->cd();
//   cout << "nr. contours = " << contours.size();
  for ( unsigned int i=0; i<contours.size(); ++i )  contours[i]->Draw();
  if ( graph && !contours.empty() )  *graph = contours.back();
  cout << "Nr. of points = " << (*graph)->GetN() << endl;

//   TH1* hc1 = new TH1F("hc1","hc1",hc->GetNbinsX(),hc->GetXaxis()->GetXmin(),hc->GetXaxis()->GetXmax());
//   double* xc1 = (**graph).GetX();
//   double* yc1 = (**graph).GetY();
//   for ( unsigned int i=0; i<(**graph).GetN(); ++i ) {
//     int ibin = hc1->FindBin(xc1[i]);
//     if ( hc1->GetBinContent(ibin) > 1.e-10 ) {
//       cout << "Two points for bin " << ibin << "????" << endl;
//       continue;
//     }
//     hc1->SetBinContent(ibin,yc1[i]);
//   }
//   hc1->Draw("same");
}


vector<TGraph*>
PlotLimits::getContours (TH2* histo, TH2* hExist)
{
  vector<TGraph*> result;
  if ( histo == 0 )  return result;

  double levels[] = { level_ };
  //
  // work on clone (avoid modifying
  //   a histogram which has already been drawn)
  //
  string newName(histo->GetName());
  newName += "Cont";
  TH2* hc = (TH2*)histo->Clone(newName.c_str());
  hc->SetContour(1,levels);
  hc->SetLineColor(1);
  TCanvas* c = new TCanvas();
  hc->Draw("cont list");
  c->Update();

  //
  // extract all contours
  //
  TObjArray* conts = (TObjArray*)gROOT->GetListOfSpecials()->FindObject("contours");
  if ( conts ) {
    // loop over levels (here: 1 level)
    for ( int i=0; i<conts->GetSize(); i++ ) {
      TList* contLevel = (TList*)conts->At(i);
      cout << "Found cont level of size " << contLevel->GetSize() << endl;
      // loop over contours at this level
      for(int j = 0; j < contLevel->GetSize(); j++){
	// keep a clone of the graph
	TGraph* curv = (TGraph*)contLevel->At(j)->Clone();
	cout << "Found graph with " << curv->GetN() << " points" << endl;
// 	{
// 	  double x,y;
// 	  cout << "First / last point = ";
// 	  curv->GetPoint(0,x,y);
// 	  cout << x << "," << y << " / ";
// 	  curv->GetPoint(curv->GetN()-1,x,y);
// 	  cout << x << "," << y << endl;
// 	}
// // 	cout << "Found graph" << endl;
//  	curv->SetLineColor(4);
// 	curv->Draw();
        // only use graphs with at least 1 point and passing "quality" cuts
	if ( curv->GetN()>0 && acceptGraph(curv,hc) ) {
// 	  TGraph* curv1 = (TGraph*)curv->Clone();
// 	  curv1->SetLineColor(2);
// 	  curv1->Draw();
	  curv->SetLineColor(3);
	  curv->SetLineWidth(3);
	  curv->Draw();
	  result.push_back(curv);
	}
      }
    }
  }
//   c->Update();
  delete c;
  
  //
  // now try to combine in one graph
  //
  // index: m0 bin number
  typedef map<int,pair<float,float> > BestPointMap;
  BestPointMap bestPointMap;
  TAxis* xaxis = hc->GetXaxis();
  TAxis* yaxis = hc->GetYaxis();
  // loop over graphs
  for ( unsigned int i=0; i<result.size(); ++i ) {
    double* xvals = result[i]->GetX();
    double* yvals = result[i]->GetY();
    double x,y;
    // loop over points
    for ( unsigned j=0; j<result[i]->GetN(); ++j ) {
      int ibx = xaxis->FindBin(xvals[j]);
      int iby = yaxis->FindBin(yvals[j]);
      if ( hExist->GetBinContent(ibx,iby)<1.e-6 || 
	   (iby>1&&hExist->GetBinContent(ibx,iby-1)<1.e-6) ||
	   (iby<hc->GetYaxis()->GetNbins()&&hExist->GetBinContent(ibx,iby+1)<1.e-6) ) continue;
      int ibin = xaxis->FindBin(xvals[j]);
      BestPointMap::iterator it = bestPointMap.find(ibin);
      // take lowest point in each column
      if ( it==bestPointMap.end() ) {
	bestPointMap[ibin] = pair<float,float>(xvals[j],yvals[j]);
      }
      else {
	if ( yvals[j]<it->second.second )  
	  it->second = pair<float,float>(xvals[j],yvals[j]);
      }
    }
  }

  //
  // create new graph from the "best" points and
  //   add it as last in the result list
  //
  TGraph* combGraph = new TGraph();
  int nCombGraph(0);
  for ( BestPointMap::const_iterator it=bestPointMap.begin(); it!=bestPointMap.end(); ++it ) {
    int ibx = it->first;
    double y = it->second.second;
    int iby = yaxis->FindBin(y);
    // skip point if the bin or the bin directly below are empty
    //  (most probably a missing job)
//     if ( hc->GetBinContent(ibx,iby)<1.e-6 || 
// 	 (iby>1&&hc->GetBinContent(ibx,iby-1)<1.e-6) ) continue;
//     if ( hExist->GetBinContent(ibx,iby)<1.e-6 || 
// 	 (iby>1&&hExist->GetBinContent(ibx,iby-1)<1.e-6) ||
// 	 (iby<hc->GetYaxis()->GetNbins()&&hExist->GetBinContent(ibx,iby+1)<1.e-6) ) continue;
    combGraph->SetPoint(nCombGraph++,it->second.first,it->second.second);
  }
  combGraph->SetLineColor(2);
  combGraph->SetLineWidth(3);
  result.push_back(combGraph);

  return result;
}

bool
PlotLimits::acceptGraph (TGraph* graph, TH2* histo)
{
  // maps holding maximum width in x or y
  typedef map< int, vector<int> > MinMaxMap;
  MinMaxMap xMap;
  MinMaxMap yMap;

  double* xvalues = graph->GetX();
  double* yvalues = graph->GetY();
  TAxis* xaxis = histo->GetXaxis();
  TAxis* yaxis = histo->GetYaxis();
  // store bin numbers of points:
  //   row numbers by column
  //   column numbers by row
  for ( int ip=0; ip<graph->GetN(); ++ip ) {
    int ix = xaxis->FindBin(xvalues[ip]);
    int iy = yaxis->FindBin(yvalues[ip]);
    xMap[ix].push_back(iy);
    yMap[iy].push_back(ix);
  }
//   if ( graph->GetN()>100 ) {
//     for ( MinMaxMap::iterator ix=xMap.begin(); ix!=xMap.end(); ++ix ) {
//       cout << ix->first << " :";
//       for ( size_t i=0; i<ix->second.size(); ++i )  cout << " " << ix->second[i];
//       cout << endl;
//     }
//   }
  //
  // get maximum spread in each column
  // 
  int dyMax(0);
  for ( MinMaxMap::iterator it=xMap.begin(); it!=xMap.end(); ++it ) {
    vector<int> ybins = it->second;
    if ( ybins.size()<2 )  continue;
    sort(ybins.begin(),ybins.end());
    if ( (ybins.back()-ybins.front())>dyMax )  dyMax = ybins.back() - ybins.front();
  }
  //
  // get maximum spread in each row
  // 
  int dxMax(0);
  for ( MinMaxMap::iterator it=yMap.begin(); it!=yMap.end(); ++it ) {
    vector<int> xbins = it->second;
    if ( xbins.size()<2 )  continue;
    sort(xbins.begin(),xbins.end());
    if ( (xbins.back()-xbins.front())>dxMax )  dxMax = xbins.back() - xbins.front();
  }
//   if ( graph->GetN()>100 )  cout << "dxMax / dyMax = " << dxMax << " " << dyMax << endl;
  //
  // require a maximum spread of >1 bin both in x and y
  //   (excludes 1 bin "islands" or parts of single rows / columns
  //    selected due to fluctuations or missing jobs)
  return dxMax>1 && dyMax>1;
}

void
PlotLimits::saveContours ()
{
  string cname = name_;
  cname += "-contours.root";
  TDirectory* curDir = gDirectory;
  
  TFile* f = new TFile(cname.c_str(),"recreate");
  if ( gObs->GetN()>0 ) gObs->Write("gObs");
  if ( gExpMinus2->GetN()>0 ) gExpMinus2->Write("gExpMinus2");
  if ( gExpMinus1->GetN()>0 ) gExpMinus1->Write("gExpMinus1");
  if ( gExpMedian->GetN()>0 ) gExpMedian->Write("gExpMedian");
  if ( gExpPlus1->GetN()>0 ) gExpPlus1->Write("gExpPlus1");
  if ( gExpPlus2->GetN()>0 ) gExpPlus2->Write("gExpPlus2");
  f->Write();
  delete f;

  curDir->cd();
}

void PlotLimits::drawSlices (float m0) {
  if ( m0 < 0. )  return;

  int ibx = hObs.first->GetXaxis()->FindBin(m0);
  if ( ibx<1 || ibx>hObs.first->GetNbinsX() )  return;

  TCanvas* c = new TCanvas("cSlice","cSlice",1000,1000);
  c->Divide(2,3);
  c->cd(1);
  drawSlice(hExpMinus2.first,gExpMinus2,ibx,m0);
  c->cd(2);
  drawSlice(hExpMinus1.first,gExpMinus1,ibx,m0);
  c->cd(3);
  drawSlice(hObs.first,gObs,ibx,m0);
  c->cd(4);
  drawSlice(hExpMedian.first,gExpMedian,ibx,m0);
  c->cd(5);
  drawSlice(hExpPlus1.first,gExpPlus1,ibx,m0);
  c->cd(6);
  drawSlice(hExpPlus2.first,gExpPlus2,ibx,m0);
}

void PlotLimits::drawSlice (TH2* histo, TGraph* graph, int ibx, float m0) {
  string hName(histo->GetName());
  TH1* histo1D = histo->ProjectionY((hName+"py").c_str(),ibx,ibx);
  gPad->SetLogy(1); gPad->SetGridx(1); gPad->SetGridy(1);
  histo1D->SetMinimum(0.0001); 
  histo1D->SetMarkerStyle(20);
  histo1D->SetLineWidth(2);
  histo1D->Draw("PL");

  int np = graph->GetN();
  if ( np > 1 ) {
    double* xg = graph->GetX();
    double* yg = graph->GetY();
    double x0 = xg[0];
    for ( int i=1; i<np; ++i ) {
      double x1 = xg[i];
      if ( m0>=x0 && m0<=x1 ) {
	double y = (yg[i]*(m0-x0)-yg[i-1]*(m0-x1))/(x1-x0);
	TMarker* marker = new TMarker(y,0.05,29);
	marker->SetMarkerColor(2);
	marker->Draw();
	break;
      }
      x0 = x1;
    }
  }
}
