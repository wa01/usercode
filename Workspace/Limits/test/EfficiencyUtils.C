#include "TH2.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TCanvas.h"

#include <iostream>

using namespace std;

//
// Smoothing of efficiency in msugra plane
//

//
// fill missing bins in histogram h fitting a plane to the 
//   surrounding bins (minimum number nmin of non-empty bins) 
//
TH2* fillMissing (TH2* h, int nmin=5) {
  TH2* hNew = h;
  int nbx = h->GetNbinsX();
  int nby = h->GetNbinsX();

  TMatrixD mat(3,3);
  TVectorD cvec(3);
  // loop over histogram (excluding a 1-bin wide margin)
  for ( int ix=2; ix<nbx; ++ix ) {
    for ( int iy=2; iy<nby; ++iy ) {
      // check only empty bins
      if ( h->GetBinContent(ix,iy)>1.e-10 )  continue;
      // clear matrix and vector used for fit
      int nn(0);
      for ( int i=0; i<3; ++i ) {
	cvec(i) = 0.;
	for ( int j=0; j<3; ++j ) mat(i,j) = 0.;
      }
      // loop over neighbours
      for ( int jx=-1; jx<2; ++jx ) {
	for ( int jy=-1; jy<2; ++jy ) {
	  // skip the central bin (the one to be filled)
	  if ( jx==0 && jy==0 )  continue;
	  // skip empty neighbours
	  double z =  h->GetBinContent(ix+jx,iy+jy);
	  if ( z<1.e-10 )  continue;
	  // update matrix and vector
	  ++nn;
	  mat(0,0) += jx*jx; mat(0,1) += jx*jy; mat(0,2) += jx;
	  mat(1,0) += jx*jy; mat(1,1) += jy*jy; mat(1,2) += jy;
	  mat(2,0) += jx;    mat(2,1) += jy;    mat(2,2) += 1;
	  cvec(0)  += jx*z;  cvec(1)  += jy*z;  cvec(2)  += z;
	}
      }
      // don't change empty bin if <nmin non-empty neighbours
      if ( nn<nmin ) continue;
//       cout << "x / y = " << h->GetXaxis()->GetBinCenter(ix) << " " << h->GetYaxis()->GetBinCenter(iy) << endl;
      // 
      // linear 2D fit to neighbours (in units of bin number): 
      //   par(0)*(x-ix)+par(1)*(y-iy)+par(2)
      //
      double det;
      TMatrixD mat1(mat);
      mat.Invert(&det);
      TVectorD par = mat*cvec;
      TVectorD tmp = mat1*par;
      hNew->SetBinContent(ix,iy,par(2));
    }
  }
  return hNew;
}
//
// clear bins with more than 1 empty neighbour
//   (remove artefacts of smoothing at the edges of
//    the physical region)
//
void clearBins (TH2* histo, TH2* refHisto) {
  int nbx = refHisto->GetNbinsX();
  int nby = refHisto->GetNbinsY();
  // loop over bins
  for ( int ix=1; ix<=nbx; ++ix ) {
    for ( int iy=1; iy<=nby; ++iy ) {
      int nempty(0);
      // loop over neighbours
      for ( int jx=ix-1; jx<ix+2; ++jx ) {
	if ( jx<1 || jx>nbx )  continue;
	for ( int jy=iy-1; jy<iy+2; ++jy ) {
	  if ( jy<1 || jy>nby )  continue;
	  if ( refHisto->GetBinContent(jx,jy)<1.e-10 )  ++nempty;
	}
      }
      // clear bin if >1 empty neighbour
      if ( nempty>1 )  histo->SetBinContent(ix,iy,0.);
    }
  }
}

//
// perform filling of (isolated) empty bins and smoothing
//   on an efficiency histogram in the msugra plane
//
TH1* doAll (TH2* hRaw) {
  // fill isolated empty bins
  TH2* hFilled = fillMissing(hRaw);
  TCanvas* c = new TCanvas("cFilled","cFilled");
  hFilled->Draw("ZCOL");
  // smooth histogram (2 steps, 3x3)
  TH2* hSmooth = (TH2*)hFilled->Clone("hSmooth");
  hSmooth->SetTitle("hSmooth");
  hSmooth->Smooth(1,"k3a");
  hSmooth->Smooth(1,"k3a");
  clearBins(hSmooth,hFilled);
  c = new TCanvas("cSmooth","cSmooth");
  hSmooth->Draw("ZCOL");
  // relative difference smoothed/raw (cross check)
  TH2* hRelDiff = (TH2*)hFilled->Clone("hRelDiff");
  hRelDiff->SetTitle("hRelDiff");
  hRelDiff->Add(hSmooth,-1);
  hRelDiff->Divide(hSmooth);
  hRelDiff->SetMinimum(-0.15);
  hRelDiff->SetMaximum(0.15);
  hRelDiff->Smooth();
  hRelDiff->Smooth();
  c = new TCanvas("cRelDiff","cRelDiff");
  hRelDiff->Draw("ZCOL");
  // ratio smoothed / raw
  TH2* hRatio = (TH2*)hSmooth->Clone("hRatio");
  hRatio->SetTitle("hRatio");
  hRatio->Divide(hFilled);
  hRatio->SetMinimum(0.5);
  hRatio->SetMaximum(1.5);
//   hRatio->Smooth();
//   hRatio->Smooth();
  c = new TCanvas("cRatio","cRatio");
  hRatio->Draw("ZCOL");

  return hRatio;
}

