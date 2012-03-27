#define PlotLimits2_cxx
#include "PlotLimits2.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TGraph.h"
#include "TMarker.h"

#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <algorithm>
#include <set>
#include <cmath>
#include <string>

using namespace std;


void PlotLimits2::Loop()
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

   hObs.hExist = new TH2F("hObsExist","hObsExist",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hExpMinus2.hExist = new TH2F("hExpMinus2Exist","hExpMinus2Exist",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hExpMinus1.hExist = new TH2F("hExpMinus1Exist","hExpMinus1Exist",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hExpMedian.hExist = new TH2F("hExpMedianExist","hExpMedianExist",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hExpPlus1.hExist = new TH2F("hExpPlus1Exist","hExpPlus1Exist",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hExpPlus2.hExist = new TH2F("hExpPlus2Exist","hExpPlus2Exist",nb0,fm0min,fm0max,nb12,fm12min,fm12max);

   hObs.hLimit = new TH2F("hObs","hObs",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hExpMinus2.hLimit = new TH2F("hExpMinus2","hExpMinus2",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hExpMinus1.hLimit = new TH2F("hExpMinus1","hExpMinus1",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hExpMedian.hLimit = new TH2F("hExpMedian","hExpMedian",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hExpPlus1.hLimit = new TH2F("hExpPlus1","hExpPlus1",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hExpPlus2.hLimit = new TH2F("hExpPlus2","hExpPlus2",nb0,fm0min,fm0max,nb12,fm12min,fm12max);

   hObs.hLimitErr = new TH2F("hObsErr","hObsErr",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hExpMinus2.hLimitErr = new TH2F("hExpMinus2Err","hExpMinus2Err",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hExpMinus1.hLimitErr = new TH2F("hExpMinus1Err","hExpMinus1Err",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hExpMedian.hLimitErr = new TH2F("hExpMedianErr","hExpMedianErr",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hExpPlus1.hLimitErr = new TH2F("hExpPlus1Err","hExpPlus1Err",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
   hExpPlus2.hLimitErr = new TH2F("hExpPlus2Err","hExpPlus2Err",nb0,fm0min,fm0max,nb12,fm12min,fm12max);
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

      LimitHistograms* histos;
//       TH2* hLim(0);
//       TH2* hExist(0);
      // observed
      if ( quantileExpected<0.) {
	histos = &hObs;
// 	hLim = hObs.hLimit;
// 	hExist = hObs.hExist;
      }
      // -2 sigma
      else if ( fabs(quantileExpected-0.025)<0.001 ) {
	histos = &hExpMinus2;
// 	hLim = hExpMinus2.hLimit;
// 	hExist = hExpMinus2.hExist;
      }
      // -1 sigma
      else if ( fabs(quantileExpected-0.16)<0.01 ) {
	histos = &hExpMinus1;
// 	hLim = hExpMinus1.hLimit;
// 	hExist = hExpMinus1.hExist;
      }
      // median
      else if ( fabs(quantileExpected-0.50)<0.01 ) {
	histos = &hExpMedian;
// 	hLim = hExpMedian.hLimit;
// 	hExist = hExpMedian.hExist;
      }
      // +1 sigma
      else if ( fabs(quantileExpected-0.84)<0.01 ) {
	histos = &hExpPlus1;
// 	hLim = hExpPlus1.hLimit;
// 	hExist = hExpPlus1.hExist;
      }
      // +2 sigma
      else if ( fabs(quantileExpected-0.975)<0.001 ) {
	histos = &hExpPlus2;
// 	hLim = hExpPlus2.hLimit;
// 	hExist = hExpPlus2.hExist;
      }
      if ( histos->hExist ) {
	int ibin = histos->hExist->FindBin(m0,m12);
	if ( histos->hExist->GetBinContent(ibin) > 0.5 ) {
	  std::cout << "***** Duplicate entry for " << m0 << " "
		    << m12 << " " << histos->hExist->GetTitle() << std::endl;
	}
	histos->hLimit->Fill(m0,m12,limit);
	histos->hLimitErr->Fill(m0,m12,limitErr);
	histos->hExist->Fill(m0,m12);
      }
   }
}

void 
PlotLimits2::drawHistograms()
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

bool
PlotLimits2::interpolateExp (double ys[], double zs[], 
			    double level, double& limit) const 
{
  if ( zs[1]>zs[0] )  return false;

  double logLev = log(level);
  double y0 = ys[0]; double y1 = ys[1];
  double z0 = log(zs[0]); double z1 = log(zs[1]);
  double y = ((logLev-z0)*y1-(logLev-z1)*y0) / (z1-z0);
  if ( y>=min(y0,y1) && y<=max(y0,y1) ) {
    limit = y;
    return true;
  }
  return false;
}

// TGraph*
TH1*
PlotLimits2::scanLimit(LimitHistograms& histos)
{
  int nbx = histos.hExist->GetNbinsX();
  int nby = histos.hExist->GetNbinsY();

  TAxis* xAxis = histos.hExist->GetXaxis();
  TAxis* yAxis = histos.hExist->GetYaxis();

//   TGraph* result = new TGraph();
  TH1* result = new TH1F("lim","lim",nbx,xAxis->GetXmin(),xAxis->GetXmax());

  int nGraph(0);
  for ( int ix=1; ix<=nbx; ++ix ) {
    int nEmpty(0);
    int nUp(0);
    bool found(false);
    double yLim;
    int np(0);
    double ys[2];
    double zs[2];
    double ezs[2];
    for ( int iy=nby; iy>0; --iy ) {
      if ( histos.hExist->GetBinContent(ix,iy)<1.e-6 || 
	   histos.hLimit->GetBinContent(ix,iy)<=0. ) {
	if ( ++nEmpty>1 )  np = 0;
	continue;
      }
      nEmpty = 0;
      double y = yAxis->GetBinCenter(iy);
      double z = histos.hLimit->GetBinContent(ix,iy);
      double ez = histos.hLimitErr->GetBinContent(ix,iy);
      if ( np>0 && z>zs[np] ) {
	double dz = z - zs[np];
	double edz = sqrt(ez*ez+ezs[np]*ezs[np]);
	if ( dz>edz ) {
	  ++nUp;
	  if ( nUp>1 )  break;
	}
      }
      ys[np] = y;
      zs[np] = z;
      ezs[np] = ez;
      ++np;
      if ( np>1 ) {
	found = interpolateExp(ys,zs,level_,yLim);
	if ( found )  break;
	np = 1;
	ys[0] = ys[1];
	zs[0] = zs[1];
	ezs[0] = ezs[1];
      }
    }
    if ( found ) {
      nGraph++;
      result->SetBinContent(ix,yLim);
      std::cout << ix << " " << yLim << " " << result->GetBinContent(ix) << std::endl;
//       result->SetPoint(nGraph++,xAxis->GetBinCenter(ix),yLim);
    }
  }
  //
  // interpolate empty bins
  //
  int ifirst(0);
  int ilast(0);
  int np = 0;
  int ixs[2];
  double ys[2];
  for ( int ix=1; ix<=nbx; ++ix ) {
    double y = result->GetBinContent(ix);
    if ( y<1.e-6 )  continue;
    if ( ifirst == 0 )  ifirst = ix;
    ilast = ix;
    ixs[np] = ix;
    ys[np] = y;
    ++np;
    if ( np>1 ) {
      int dix = ixs[1] - ixs[0];
      if ( dix>1 ) {
	double a = (ys[1]-ys[0]) / dix;
	double b = (ys[0]*ixs[1]-ys[1]*ixs[0]) / dix;
	for ( int i=ixs[0]+1; i<ixs[1]; ++i )  result->SetBinContent(i,a*i+b);
      }
      ixs[0] = ixs[1];
      ys[0] = ys[1];
      np = 1;
    }
  }
  if ( ifirst!=0 )  result->GetXaxis()->SetRange(ifirst,ilast);
  result->Smooth(6,"R");
  return result;
}


void 
PlotLimits2::drawHistogram (LimitHistograms histos, TGraph** graph)
{
  if ( histos.hLimit->GetEntries()==0 )  return;

  double levels[1] = { level_ };

  histos.hLimit->SetMaximum(relmax_*level_);
  histos.hLimit->Draw("zcol");
  //
  // need to draw and update the pad in order to have access
  //   to the contour lines ...
  //
  TH2* hc = (TH2*)histos.hLimit->Clone();
  hc->SetContour(1,levels);
  hc->Draw("cont3 list same");
  *graph = 0;
//   *graph = scanLimit(histos);
//   (**graph).SetLineColor(2);
//   (**graph).SetLineWidth(2);  
//   (**graph).Draw("same");
  TH1* hlim = scanLimit(histos);
  hlim->SetLineColor(2);
  hlim->SetLineWidth(2);
  hlim->Draw("same C");
  gPad->Update();
  return;
}



void
PlotLimits2::saveContours ()
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

void PlotLimits2::drawSlices (float m0) {
  if ( m0 < 0. )  return;

  int ibx = hObs.hLimit->GetXaxis()->FindBin(m0);
  if ( ibx<1 || ibx>hObs.hLimit->GetNbinsX() )  return;

  TCanvas* c = new TCanvas("cSlice","cSlice",1000,1000);
  c->Divide(2,3);
  c->cd(1);
  drawSlice(hExpMinus2,gExpMinus2,ibx,m0);
  c->cd(2);
  drawSlice(hExpMinus1,gExpMinus1,ibx,m0);
  c->cd(3);
  drawSlice(hObs,gObs,ibx,m0);
  c->cd(4);
  drawSlice(hExpMedian,gExpMedian,ibx,m0);
  c->cd(5);
  drawSlice(hExpPlus1,gExpPlus1,ibx,m0);
  c->cd(6);
  drawSlice(hExpPlus2,gExpPlus2,ibx,m0);
}

void PlotLimits2::drawSlice (LimitHistograms& histos, TGraph* graph, int ibx, float m0) {
  string hName(histos.hLimit->GetName());
  TH1* histo1D = histos.hLimit->ProjectionY((hName+"py").c_str(),ibx,ibx);
  hName = histos.hLimitErr->GetName();
  TH1* histoErr1D = histos.hLimitErr->ProjectionY((hName+"py").c_str(),ibx,ibx);
  for ( int i=1; i<=histo1D->GetNbinsX(); ++i ) {
    histo1D->SetBinError(i,histoErr1D->GetBinContent(i));
  }
  gPad->SetLogy(1); gPad->SetGridx(1); gPad->SetGridy(1);
  histo1D->SetMinimum(0.0001); 
  histo1D->SetMarkerStyle(20);
  histo1D->SetLineWidth(2);
  histo1D->Draw("HIST L");
  histo1D->Draw("EX0 same");

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
