#include <iostream>
#include <fstream>
#include <TROOT.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TGraphErrors.h>
#include <TH2F.h>
#include <TF1.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TText.h>
#include <TString.h>

#include "TMatrixDSym.h"
#include "TVectorD.h"

enum { NBINS = 20 };

float sf[NBINS], sf_stat[NBINS], sf_syst[NBINS], sf_eror[NBINS];

float eff_PT[NBINS], eff_stat_PT[NBINS], effMC_PT[NBINS], effMC_stat_PT[NBINS];
float sf_PT[NBINS], sf_stat_PT[NBINS], sf_syst_PT[NBINS], sf_eror_PT[NBINS];
float mupt_syst_PT[NBINS], gluon_syst_PT[NBINS], pu_syst_PT[NBINS], total_syst_PT[NBINS];
float away_syst_PT[NBINS];
float lc_syst_PT[NBINS];

float eff_SY[NBINS], eff_stat_SY[NBINS], effMC_SY[NBINS], effMC_stat_SY[NBINS];
float sf_SY[NBINS], sf_stat_SY[NBINS], sf_syst_SY[NBINS], sf_eror_SY[NBINS];
float mupt_syst_SY[NBINS], gluon_syst_SY[NBINS], pu_syst_SY[NBINS], total_syst_SY[NBINS];
float away_syst_SY[NBINS];
float ptrel_syst_SY[NBINS], clos_syst_SY[NBINS];

float eff_IP[NBINS], eff_stat_IP[NBINS], effMC_IP[NBINS], effMC_stat_IP[NBINS];
float sf_IP[NBINS], sf_stat_IP[NBINS], sf_syst_IP[NBINS], sf_eror_IP[NBINS];
float mupt_syst_IP[NBINS], gluon_syst_IP[NBINS], pu_syst_IP[NBINS], total_syst_IP[NBINS];
float away_syst_IP[NBINS];

float eff_JP[NBINS], eff_stat_JP[NBINS], effMC_JP[NBINS], effMC_stat_JP[NBINS];
float sf_JP[NBINS], sf_stat_JP[NBINS], sf_syst_JP[NBINS], sf_eror_JP[NBINS];
float mupt_syst_JP[NBINS], gluon_syst_JP[NBINS], pu_syst_JP[NBINS], total_syst_JP[NBINS];
float cor_syst_JP[NBINS], inc_syst_JP[NBINS], bias_syst_JP[NBINS];

float frac_PTJP[NBINS], frac_SYJP[NBINS], frac_PTSY[NBINS];

float away[4];
float mu[4];
float ptrel[4];
float gsplit[4];
float closure[4];
float pu[4];

//
// correlations (order: PT/SY, PT/IP, PT/JP, SY/IP, SY/JP, IP/JP)
//
float corrStat[] =      { 1,     1,     1,     1,     1,    1 };
float corrSystPU[] =    { 1,     1,     1,     1,     1,    1 };
float corrSystMuPt[] =  { 1,    -1,    -1,    -1,    -1,    1 };
float corrSystGluon[] = { 1,    -1,    -1,    -1,    -1,    1 };
float corrSystAway[] =  { 0,     1,     0,     0,     0,    0 };
 


// debug: print sym. matrix
void printMatrix (TMatrixDSym& mat, const char* tit = "matrix")
{
  cout << tit << endl;
  for ( size_t i=0; i<mat.GetNrows(); ++i ) {
    for ( size_t j=0; j<mat.GetNrows(); ++j ) 
      printf("%12.9f ",mat(i,j));
    printf("\n");
  }
}
// debug: print matrix
void printMatrix (TMatrixD& mat, const char* tit = "matrix")
{
  cout << tit << endl;
  for ( size_t i=0; i<mat.GetNrows(); ++i ) {
    for ( size_t j=0; j<mat.GetNrows(); ++j ) 
      printf("%12.9f ",mat(i,j));
    printf("\n");
  }
}

// clear matrix contents
void clearMatrix (TMatrixDSym& mat)
{
  for ( size_t i=0; i<mat.GetNrows(); ++i ) {
    for ( size_t j=0; j<mat.GetNrows(); ++j )  mat(i,j) = 0;
  }  
}

// square
inline double sqr (double arg) {return arg*arg;}

// combined variance for a list of uncertainties
double combVar (double e0, double e1=0, double e2=0, double e3=0, double e4=0, double e5=0, double e6=0, double e7=0, double e8=0, double e9=0)
{
  double sum(0.);
  sum += sqr(e0);
  sum += sqr(e1);
  sum += sqr(e2);
  sum += sqr(e3);
  sum += sqr(e4);
  sum += sqr(e5);
  sum += sqr(e6);
  sum += sqr(e7);
  sum += sqr(e8);
  sum += sqr(e9);
  return sum;
}

// combined s.d. for a list of uncertainties
double combError (double e0, double e1=0, double e2=0, double e3=0, double e4=0, double e5=0, double e6=0, double e7=0, double e8=0, double e9=0)
{
  return sqrt(combVar(e0, e1, e2, e3, e4, e5, e6, e7, e8, e9));
}

// fill off-diagonal terms of a covariance matrix assuming 100% (positive) correlation
void correlate (TMatrixDSym& matrix, float* correlations)
{
  size_t ind = 0;
  for ( size_t i=0; i<matrix.GetNrows(); ++i ) {
    for ( size_t j=i+1; j<matrix.GetNrows(); ++j ) {
      matrix(i,j) = matrix(j,i) = correlations[ind++]*sqrt(matrix(i,i)*matrix(j,j));
    }
  }
}


// remove unused rows & columns from matrix
TMatrixDSym compactify (TMatrixDSym& cov, bool usePT, bool useSY, bool useIP, bool useJP)
{
  size_t n(0);
  vector<int> flags(4,-1);
  if ( usePT )  flags[0] = n++;
  if ( useSY )  flags[1] = n++;
  if ( useIP )  flags[2] = n++;
  if ( useJP )  flags[3] = n++;

  TMatrixDSym newCov(n);
  for ( size_t i=0; i<4; ++i ) {
    if ( flags[i]<0 )  continue;
    for ( size_t j=0; j<4; ++j ) {
      if ( flags[j]<0 )  continue;
      newCov(flags[i],flags[j]) = cov(i,j);
    }
  }
  return newCov;
}

// remove unused rows & columns from vector
TVectorD compactify (TVectorD& vec, bool usePT, bool useSY, bool useIP, bool useJP)
{
  size_t n(0);
  vector<int> flags(4,-1);
  if ( usePT )  flags[0] = n++;
  if ( useSY )  flags[1] = n++;
  if ( useIP )  flags[2] = n++;
  if ( useJP )  flags[3] = n++;

  TVectorD newVec(n);
  for ( size_t i=0; i<4; ++i ) {
    if ( flags[i]<0 )  continue;
    newVec(flags[i]) = vec[i];
  }
  return newVec;
}

//
// LSQ combination of several correlated variables
//
pair<double, double> combine (const TVectorD& vec, const TMatrixDSym& cov, double& chi2)
{
  pair<double,double> result(0.,0.);
  chi2 = -1.;

  // weight matrix
  double det;
  TMatrixDSym wgt(cov);
  wgt.Invert(&det);
  // cout << "Det = " << det << std::endl;
  if ( det<0 )  return result;

 // for ( size_t i=0; i<cov.GetNrows(); ++i ) {
 //   printf("%8.3f ; ",vec(i));
 //   for ( size_t j=0; j<cov.GetNrows(); ++j )  printf("%9.6f ",cov(i,j));
 //   printf("\n");
 // }

  // auxiliary matrix (Nx1)
  TVectorD atg(wgt.GetNrows());
  for ( size_t i=0; i<wgt.GetNrows(); ++i ) {
    for ( size_t j=0; j<wgt.GetNrows(); ++j )  atg(i) += wgt(i,j);
  }

  double atgm(0);
  double atga(0);
  for ( size_t i=0; i<wgt.GetNrows(); ++i ) {
    atgm += atg(i)*vec(i);
    atga += atg(i);
  }
  result.first = atgm/atga;
  result.second = sqrt(1./atga);
  
  TVectorD dm(vec);
  for ( size_t i=0; i<wgt.GetNrows(); ++i ) dm(i) -= result.first;
  chi2 = wgt.Similarity(dm);

  double wgtMin(999999);
  for ( size_t i=0; i<wgt.GetNrows(); ++i ) {
    double wi = atg(i)/atga;
    if ( wi<wgtMin )  wgtMin = wi;
  }
  if ( wgtMin<0 )  cout << "**** at least one negative weight ****" << endl;
  return result;
}


// (result,error) for combination of bin "k" & tagger "tagger" for the methods indicated in the list
pair<double,double> combine (int k, const string& tagger, bool usePT, bool useSY, bool useIP, bool useJP)
{
 //
 // syst. error components for PT, SY and IP are assumed to be relative !!!
 //
 TMatrixDSym totalCov(4);
 clearMatrix(totalCov);
 TMatrixDSym currCov(4);

 // start with stat errors
 // start with stat error (fully correlated)
 // totalCov *= 0.;
 // totalCov(0,0) = sqr(sf_stat_PT[k]);
 // totalCov(1,1) = sqr(sf_stat_SY[k]);
 // totalCov(2,2) = sqr(sf_stat_IP[k]);
 // totalCov(3,3) = sqr(sf_stat_JP[k]);
 // correlate(totalCov,corrStat);
 // assume PT/IP contained in SY, SY contained in JP
 // total stat variance is divided according to the fraction of the subsample
 // for each tagger
 // events used for PT/IP, SY and JP
 clearMatrix(currCov);
 currCov(0,0) = sqr(sf_stat_PT[k]);
 currCov(1,1) = sqr(frac_PTSY[k]*sf_stat_SY[k]);
 currCov(2,2) = sqr(sf_stat_IP[k]);
 currCov(3,3) = sqr(frac_PTSY[k]*frac_SYJP[k]*sf_stat_JP[k]);
 correlate(totalCov,corrStat);
 totalCov += currCov;
 // events used for SY and JP
 clearMatrix(currCov);
 currCov(1,1) = sqr((1-frac_PTSY[k])*sf_stat_SY[k]);
 currCov(3,3) = sqr(frac_SYJP[k]*(1-frac_PTSY[k])*sf_stat_JP[k]);
 correlate(totalCov,corrStat);
 totalCov += currCov;
 // remaining (JP only)
 clearMatrix(currCov);
 currCov(3,3) = sqr((1-frac_SYJP[k])*sf_stat_JP[k]);
 correlate(totalCov,corrStat);
 totalCov += currCov;

 // pileup (fully correlated)
 currCov *= 0.;
 currCov(0,0) = sqr(pu_syst_PT[k]);
 currCov(1,1) = sqr(pu_syst_SY[k]*sf_SY[k]);
 currCov(2,2) = sqr(pu_syst_IP[k]);
 currCov(3,3) = sqr(pu_syst_JP[k]*sf_JP[k]);
 correlate(currCov,corrSystPU);
 totalCov += currCov;
 // mupt (fully correlated)
 currCov *= 0.;
 currCov(0,0) = sqr(mupt_syst_PT[k]);
 currCov(1,1) = sqr(mupt_syst_SY[k]*sf_SY[k]);
 currCov(2,2) = sqr(mupt_syst_IP[k]);
 currCov(3,3) = sqr(mupt_syst_JP[k]*sf_JP[k]);
 correlate(currCov,corrSystMuPt);
 totalCov += currCov;
 // gluon (fully correlated)
 currCov *= 0.;
 currCov(0,0) = sqr(gluon_syst_PT[k]);
 currCov(1,1) = sqr(gluon_syst_SY[k]*sf_SY[k]);
 currCov(2,2) = sqr(gluon_syst_IP[k]);
 currCov(3,3) = sqr(gluon_syst_JP[k]*sf_JP[k]);
 correlate(currCov,corrSystGluon);
 totalCov += currCov;
 // away (fully correlated between PT and IP; uncorrelated for SY; none for JP)
 currCov *= 0.;
 currCov(0,0) = sqr(away_syst_PT[k]);
 currCov(1,1) = sqr(away_syst_SY[k]*sf_SY[k]);
 currCov(2,2) = sqr(away_syst_IP[k]);
 correlate(currCov,corrSystAway);
 totalCov += currCov;
 // uncorr syst: l/c (PT), ptrel and closure (SY), Cb, inclusive and bias (JP)
 currCov *= 0.;
 currCov(0,0) = sqr(lc_syst_PT[k]);
 currCov(1,1) = combVar(ptrel_syst_SY[k]*sf_SY[k],clos_syst_SY[k]*sf_SY[k]);
 currCov(3,3) = combVar(cor_syst_JP[k]*sf_JP[k],inc_syst_JP[k]*sf_JP[k]
 					       ,bias_syst_JP[k]*sf_JP[k]);
 totalCov += currCov;

 // vector of measurements
 TVectorD sfs(4);
 sfs[0] = sf_PT[k];
 sfs[1] = sf_SY[k];
 sfs[2] = sf_IP[k];
 sfs[3] = sf_JP[k];
 
 // removal of unused taggers
 TMatrixDSym cov = compactify(totalCov,usePT,useSY,useIP,useJP);
 TVectorD vec = compactify(sfs,usePT,useSY,useIP,useJP);

 double chi2;
 pair<double,double> result = combine(vec,cov,chi2);

 cout << k << " " << result.first << " " << result.second << endl;
 double chi2norm(chi2);
 if ( cov.GetNrows()>1 )  chi2norm /= (cov.GetNrows()-1);
 cout << "NDF = " << cov.GetNrows()-1 << " chi2 = " << chi2 << " " << chi2norm << endl;
 if ( chi2norm > 1. ) result.second *= sqrt(chi2norm);
//$$

 return result;
}

void plot(const char* tit)
{

// TString title = "CSVM";
TString title = tit;

int stati=0;
bool fit= 0;
bool logx=1;

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

TCanvas *c1 = new TCanvas("c1", "plots",200,10,700,750);
c1->SetFillColor(10);
c1->SetFillStyle(4000);
c1->SetBorderSize(2);

TPad* pad1 = new TPad("pad1","This is pad1",0.02,0.52,0.98,0.98,21);
TPad* pad2 = new TPad("pad2","This is pad2",0.02,0.02,0.98,0.48,21);

pad1->SetFillColor(0);
pad1->SetBorderMode(0);
pad1->SetFrameFillColor(10);
pad1->Draw();
pad1->SetLogx(logx);
// pad1->SetGrid();
   pad1->SetTopMargin(0.05);
   pad1->SetBottomMargin(0.15);
   pad1->SetRightMargin(0.04);
   pad1->SetLeftMargin(0.16);

pad2->SetFillColor(0);
pad2->SetBorderMode(0);
pad2->SetFrameFillColor(10);
pad2->Draw();
pad2->SetLogx(logx);
// pad2->SetGrid();
   pad2->SetTopMargin(0.05);
   pad2->SetBottomMargin(0.15);
   pad2->SetRightMargin(0.04);
   pad2->SetLeftMargin(0.16);

gStyle->SetOptDate(0);
gStyle->SetStatColor(0);
gStyle->SetTitleFont(42);
gStyle->SetTitleColor(1);
gStyle->SetTitleTextColor(1);
gStyle->SetTitleFillColor(10);
gStyle->SetTitleFontSize(0.05);
gStyle->SetTitleW(0.4);
gStyle->SetTitleH(0.09);
// gStyle->SetTitleX(0); // Set the position of the title box
// gStyle->SetTitleY(0.985); // Set the position of the title box
// gStyle->SetTitleStyle(Style_t style = 1001);
// gStyle->SetTitleBorderSize(2);
gStyle->SetOptStat(stati);
// gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
// gStyle->SetPadGridX(true); gStyle->SetPadGridY(true);
gStyle->SetPadGridX(false); gStyle->SetPadGridY(false);

if (fit) {
gStyle->SetStatW(0.2);
gStyle->SetStatH(0.1);
gStyle->SetOptFit(111);
} else {
gStyle->SetStatW(0.3);
gStyle->SetStatH(0.2);
gStyle->SetOptFit(0);
}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

 int nbin = 1000;
 float x[20], ex[20], x1[20], ex1[20], x2[20], ex2[20], x3[20], exstat[20];

// for JP and the full combination
  x[0]  =  25.; ex[0]  =  5.;
  x[1]  =  35.; ex[1]  =  5.;
  x[2]  =  45.; ex[2]  =  5.;
  x[3]  =  55.; ex[3]  =  5.;
  x[4]  =  65.; ex[4]  =  5.;
  x[5]  =  75.; ex[5]  =  5.;
  x[6]  =  90.; ex[6]  = 10.;
  x[7]  = 110.; ex[7]  = 10.;
  x[8]  = 140.; ex[8]  = 20.;
  x[9]  = 185.; ex[9]  = 25.;
  x[10] = 235.; ex[10] = 25.;
  x[11] = 290.; ex[11] = 30.;
  x[12] = 360.; ex[12] = 40.;
  x[13] = 450.; ex[13] = 50.;
  x[14] = 585.; ex[14] = 85.;

// for PtRel and IP3D
  x1[0]  =  24.0; x3[0]  =  25.5; ex1[0]  =  5.;
  x1[1]  =  34.0; x3[1]  =  35.5; ex1[1]  =  5.;
  x1[2]  =  44.0; x3[2]  =  45.5; ex1[2]  =  5.;
  x1[3]  =  54.0; x3[3]  =  55.5; ex1[3]  =  5.;
  x1[4]  =  64.0; x3[4]  =  65.5; ex1[4]  =  5.;
  x1[5]  =  74.0; x3[5]  =  75.5; ex1[5]  =  5.;
  x1[6]  =  88.0; x3[6]  =  91.0; ex1[6]  = 10.;
  x1[7]  = 108.0; x3[7]  = 111.0; ex1[7]  = 10.;
  x1[8]  = 138.0; x3[8]  = 142.0; ex1[8]  = 20.;
  x1[9]  = 240.0; x3[9]  = 240.0; ex1[9]  = 80.;
  x1[10] = 495.0; x3[10] = 495.0; ex1[10] =175.;

// for System8
  x2[0]  =  26.0; ex2[0]  =  5.;
  x2[1]  =  36.0; ex2[1]  =  5.;
  x2[2]  =  46.0; ex2[2]  =  5.;
  x2[3]  =  56.0; ex2[3]  =  5.;
  x2[4]  =  70.0; ex2[4]  = 10.;
  x2[5]  =  92.0; ex2[5]  = 10.;
  x2[6]  = 112.0; ex2[6]  = 10.;
  x2[7]  = 142.0; ex2[7]  = 20.;

  double par[10];

  float ttbar_pt[15] = {0., 205400., 228400., 228400., 212500., 186300.,
  286500., 183900., 174000., 65840., 19960., 8084., 3790., 2630., 0. };
  float ttbar_tot = 1806000.;

  for (int k = 1; k <= 13; k++) {
    ttbar_pt[k] = ttbar_pt[k] / ttbar_tot;
  }


  frac_PTJP[0] = 0.117946; frac_SYJP[0] = 0.193411; frac_PTSY[0] = 0.609818;
  frac_PTJP[1] = 0.0964905; frac_SYJP[1] = 0.206224; frac_PTSY[1] = 0.467892;
  frac_PTJP[2] = 0.112438; frac_SYJP[2] = 0.230399; frac_PTSY[2] = 0.488014;
  frac_PTJP[3] = 0.128465; frac_SYJP[3] = 0.310995; frac_PTSY[3] = 0.413078;
  frac_PTJP[4] = 0.144013; frac_SYJP[4] = 0.331386; frac_PTSY[4] = 0.434579;
  frac_PTJP[5] = 0.159494; frac_SYJP[5] = 0.35607; frac_PTSY[5] = 0.44793;
  frac_PTJP[6] = 0.162212; frac_SYJP[6] = 0.403861; frac_PTSY[6] = 0.401653;
  frac_PTJP[7] = 0.178839; frac_SYJP[7] = 0.432426; frac_PTSY[7] = 0.413572;
  frac_PTJP[8] = 0.182612; frac_SYJP[8] = 0.476164; frac_PTSY[8] = 0.383506;
  frac_PTJP[9] = 0.203323;
  frac_PTJP[10] = 0.223828;
  frac_PTJP[11] = 0.235041;
  frac_PTJP[12] = 0.227775;
  frac_PTJP[13] = 0.242758;
  frac_PTJP[14] = 0.239796;

// #############################################################################
if ( title == "TCHEL" ) { 
// TCHEL : PTREL
// float eff_PT[11];
eff_PT[0] = 0.445746; eff_PT[1] = 0.611461; eff_PT[2] = 0.705049; 
eff_PT[3] = 0.742701; eff_PT[4] = 0.773479; eff_PT[5] = 0.80028; 
eff_PT[6] = 0.830316; eff_PT[7] = 0.849354; eff_PT[8] = 0.851167; 
eff_PT[9] = 0.858751; eff_PT[10] = 0.803269; 
// float eff_stat_PT[11];
eff_stat_PT[0] = 0.0105682; eff_stat_PT[1] = 0.00775012; eff_stat_PT[2] = 0.00832163; 
eff_stat_PT[3] = 0.00828347; eff_stat_PT[4] = 0.00996871; eff_stat_PT[5] = 0.0123276; 
eff_stat_PT[6] = 0.00971215; eff_stat_PT[7] = 0.0109529; eff_stat_PT[8] = 0.0116127; 
eff_stat_PT[9] = 0.011897; eff_stat_PT[10] = 0.0286867; 
// float effMC_PT[11];
effMC_PT[0] = 0.518431; effMC_PT[1] = 0.659265; effMC_PT[2] = 0.732966; 
effMC_PT[3] = 0.778953; effMC_PT[4] = 0.804013; effMC_PT[5] = 0.828024; 
effMC_PT[6] = 0.843949; effMC_PT[7] = 0.850492; effMC_PT[8] = 0.864351; 
effMC_PT[9] = 0.872484; effMC_PT[10] = 0.891032; 
// float effMC_stat_PT[11];
effMC_stat_PT[0] = 0.00337147; effMC_stat_PT[1] = 0.0019906; effMC_stat_PT[2] = 0.00167447; 
effMC_stat_PT[3] = 0.00217041; effMC_stat_PT[4] = 0.00203545; effMC_stat_PT[5] = 0.00206263; 
effMC_stat_PT[6] = 0.00213585; effMC_stat_PT[7] = 0.00197953; effMC_stat_PT[8] = 0.00247783; 
effMC_stat_PT[9] = 0.00250313; effMC_stat_PT[10] = 0.00350537; 
// float sf_PT[11];
sf_PT[0] = 0.859798; sf_PT[1] = 0.927489; sf_PT[2] = 0.961913; 
sf_PT[3] = 0.953461; sf_PT[4] = 0.962024; sf_PT[5] = 0.966494; 
sf_PT[6] = 0.983846; sf_PT[7] = 0.998662; sf_PT[8] = 0.984747; 
sf_PT[9] = 0.98426; sf_PT[10] = 0.901505; 
// float sf_stat_PT[11];
sf_stat_PT[0] = 0.0211379; sf_stat_PT[1] = 0.0120847; sf_stat_PT[2] = 0.0115641; 
sf_stat_PT[3] = 0.0109609; sf_stat_PT[4] = 0.0126356; sf_stat_PT[5] = 0.0150813; 
sf_stat_PT[6] = 0.0117743; sf_stat_PT[7] = 0.0130864; sf_stat_PT[8] = 0.0137286; 
sf_stat_PT[9] = 0.0139251; sf_stat_PT[10] = 0.0323897; 
// float away_syst_PT[11];
away_syst_PT[0] = 0.00922501; away_syst_PT[1] = 0.0250744; away_syst_PT[2] = 0.0260051; 
away_syst_PT[3] = 0.015838; away_syst_PT[4] = 0.0159802; away_syst_PT[5] = 0.0160545; 
away_syst_PT[6] = 0.0143073; away_syst_PT[7] = 0.0145228; away_syst_PT[8] = 0.0373829; 
away_syst_PT[9] = 0.0373645; away_syst_PT[10] = 0.070444; 
// float mupt_syst_PT[11];
mupt_syst_PT[0] = 0.023069; mupt_syst_PT[1] = 0.00750977; mupt_syst_PT[2] = 0.0077885; 
mupt_syst_PT[3] = 0.0166044; mupt_syst_PT[4] = 0.0167536; mupt_syst_PT[5] = 0.0168314; 
mupt_syst_PT[6] = 0.0152467; mupt_syst_PT[7] = 0.0154763; mupt_syst_PT[8] = 0.004975; 
mupt_syst_PT[9] = 0.00497254; mupt_syst_PT[10] = 0.00835502; 
// float lc_syst_PT[11];
lc_syst_PT[0] = 0.00501704; lc_syst_PT[1] = 0.00440174; lc_syst_PT[2] = 0.00456511; 
lc_syst_PT[3] = 0.00532418; lc_syst_PT[4] = 0.00537199; lc_syst_PT[5] = 0.00539695; 
lc_syst_PT[6] = 0.0043277; lc_syst_PT[7] = 0.00439287; lc_syst_PT[8] = 0.0104152; 
lc_syst_PT[9] = 0.01041; lc_syst_PT[10] = 0.008021; 
// float gluon_syst_PT[11];
gluon_syst_PT[0] = 0.000407994; gluon_syst_PT[1] = 0.00185256; gluon_syst_PT[2] = 0.00192132; 
gluon_syst_PT[3] = 0.00318729; gluon_syst_PT[4] = 0.00321592; gluon_syst_PT[5] = 0.00323086; 
gluon_syst_PT[6] = 0.00359866; gluon_syst_PT[7] = 0.00365286; gluon_syst_PT[8] = 0.0218643; 
gluon_syst_PT[9] = 0.0218535; gluon_syst_PT[10] = 0.056732; 
// float pu_syst_PT[11];
pu_syst_PT[0] = 0.020003; pu_syst_PT[1] = 0.0125264; pu_syst_PT[2] = 0.0129914; 
pu_syst_PT[3] = 0.00836997; pu_syst_PT[4] = 0.00844514; pu_syst_PT[5] = 0.00848438; 
pu_syst_PT[6] = 0.00670927; pu_syst_PT[7] = 0.00681031; pu_syst_PT[8] = 0.00243768; 
pu_syst_PT[9] = 0.00243648; pu_syst_PT[10] = 0.01741; 
// float total_syst_PT[11];
total_syst_PT[0] = 0.0322914; total_syst_PT[1] = 0.0294082; total_syst_PT[2] = 0.0304997; 
total_syst_PT[3] = 0.0252014; total_syst_PT[4] = 0.0254277; total_syst_PT[5] = 0.0255459; 
total_syst_PT[6] = 0.0226683; total_syst_PT[7] = 0.0230097; total_syst_PT[8] = 0.0448854; 
total_syst_PT[9] = 0.0448632; total_syst_PT[10] = 0.0928339; 

// TCHEL : S8
eff_SY[0]=0.4612; eff_SY[1]=0.6004; eff_SY[2]=0.6838; eff_SY[3]=0.7431; eff_SY[4]=0.7678; eff_SY[5]=0.8344; eff_SY[6]=0.8389; eff_SY[7]=0.8280; 
eff_stat_SY[0]=0.0102; eff_stat_SY[1]=0.0074; eff_stat_SY[2]=0.0090; eff_stat_SY[3]=0.0096; eff_stat_SY[4]=0.0105;eff_stat_SY[5]=0.0181; eff_stat_SY[6]=0.0203; eff_stat_SY[7]=0.0172; 
effMC_SY[0]=0.4945; effMC_SY[1]=0.6361; effMC_SY[2]=0.7240; effMC_SY[3]=0.7801; effMC_SY[4]=0.8195; effMC_SY[5]=0.8529; effMC_SY[6]=0.8666; effMC_SY[7]=0.8802; 
effMC_stat_SY[0]=0.0007; effMC_stat_SY[1]=0.0007; effMC_stat_SY[2]=0.0007; effMC_stat_SY[3]=0.0006; effMC_stat_SY[4]=0.0008;effMC_stat_SY[5]=0.0007; effMC_stat_SY[6]=0.0008; effMC_stat_SY[7]=0.0007; 
sf_SY[0]=0.9327; sf_SY[1]=0.9439; sf_SY[2]=0.9445; sf_SY[3]=0.9526; sf_SY[4]=0.9369; sf_SY[5]=0.9783; sf_SY[6]=0.968; sf_SY[7]=0.9407; 
sf_stat_SY[0]=0.0207; sf_stat_SY[1]=0.0117; sf_stat_SY[2]=0.0125; sf_stat_SY[3]=0.0128; sf_stat_SY[4]=0.0139; sf_stat_SY[5]=0.0212; sf_stat_SY[6]=0.0234; sf_stat_SY[7]=0.0196; 
away[0]=2.7678; away[1]=2.3401; away[2]=1.5661; away[3]=5.044; 
mu[0]=0.1538; mu[1]=2.3534; mu[2]=4.8789; mu[3]=1.8463; 
ptrel[0]=1.4352; ptrel[1]=1.4227; ptrel[2]=1.7347; ptrel[3]=1.3032; 
gsplit[0]=0.2573; gsplit[1]=0.5751; gsplit[2]=1.027; gsplit[3]=0.7552; 
closure[0]=0.3056; closure[1]=0.1125; closure[2]=0.1517; closure[3]=0.3945; 
pu[0]=0.2221; pu[1]=0.1728; pu[2]=0.7589; pu[3]=0.1327; 

// TCHEL : IP3D
// float eff_IP[11];
eff_IP[0] = 0.431179; eff_IP[1] = 0.587223; eff_IP[2] = 0.679693; 
eff_IP[3] = 0.70567; eff_IP[4] = 0.783285; eff_IP[5] = 0.801354; 
eff_IP[6] = 0.805512; eff_IP[7] = 0.811242; eff_IP[8] = 0.847788; 
eff_IP[9] = 0.867091; eff_IP[10] = 0.917731; 
// float eff_stat_IP[11];
eff_stat_IP[0] = 0.0135925; eff_stat_IP[1] = 0.012543; eff_stat_IP[2] = 0.012834; 
eff_stat_IP[3] = 0.0106791; eff_stat_IP[4] = 0.0123747; eff_stat_IP[5] = 0.0146699; 
eff_stat_IP[6] = 0.00773698; eff_stat_IP[7] = 0.00900054; eff_stat_IP[8] = 0.00689342; 
eff_stat_IP[9] = 0.00744557; eff_stat_IP[10] = 0.0118243; 
// float effMC_IP[11];
effMC_IP[0] = 0.520958; effMC_IP[1] = 0.66147; effMC_IP[2] = 0.732438; 
effMC_IP[3] = 0.781923; effMC_IP[4] = 0.805375; effMC_IP[5] = 0.825966; 
effMC_IP[6] = 0.841829; effMC_IP[7] = 0.849809; effMC_IP[8] = 0.863054; 
effMC_IP[9] = 0.871349; effMC_IP[10] = 0.892757; 
// float effMC_stat_IP[11];
effMC_stat_IP[0] = 0.00335517; effMC_stat_IP[1] = 0.00213307; effMC_stat_IP[2] = 0.00176827; 
effMC_stat_IP[3] = 0.00215911; effMC_stat_IP[4] = 0.00200808; effMC_stat_IP[5] = 0.00203269; 
effMC_stat_IP[6] = 0.00208478; effMC_stat_IP[7] = 0.0019227; effMC_stat_IP[8] = 0.00239944; 
effMC_stat_IP[9] = 0.00240765; effMC_stat_IP[10] = 0.00350666; 
// float sf_IP[11];
sf_IP[0] = 0.827664; sf_IP[1] = 0.887755; sf_IP[2] = 0.927987; 
sf_IP[3] = 0.90248; sf_IP[4] = 0.972571; sf_IP[5] = 0.970202; 
sf_IP[6] = 0.956859; sf_IP[7] = 0.954617; sf_IP[8] = 0.982312; 
sf_IP[9] = 0.995113; sf_IP[10] = 1.02797; 
// float sf_stat_IP[11];
sf_stat_IP[0] = 0.0266303; sf_stat_IP[1] = 0.0191773; sf_stat_IP[2] = 0.017665; 
sf_stat_IP[3] = 0.013883; sf_stat_IP[4] = 0.0155553; sf_stat_IP[5] = 0.0179207; 
sf_stat_IP[6] = 0.00949125; sf_stat_IP[7] = 0.0108092; sf_stat_IP[8] = 0.00844123; 
sf_stat_IP[9] = 0.00897638; sf_stat_IP[10] = 0.0138465; 
// float away_syst_IP[11];
away_syst_IP[0] = 0.0365002; away_syst_IP[1] = 0.00375857; away_syst_IP[2] = 0.0039289; 
away_syst_IP[3] = 0.0118989; away_syst_IP[4] = 0.012823; away_syst_IP[5] = 0.0127918; 
away_syst_IP[6] = 0.00545709; away_syst_IP[7] = 0.0054443; away_syst_IP[8] = 0.0136526; 
away_syst_IP[9] = 0.0138306; away_syst_IP[10] = 0.00482795; 
// float mupt_syst_IP[11];
mupt_syst_IP[0] = 0.0556241; mupt_syst_IP[1] = 0.0260505; mupt_syst_IP[2] = 0.0272311; 
mupt_syst_IP[3] = 0.00537861; mupt_syst_IP[4] = 0.00579634; mupt_syst_IP[5] = 0.00578222; 
mupt_syst_IP[6] = 0.00590676; mupt_syst_IP[7] = 0.00589292; mupt_syst_IP[8] = 0.00101022; 
mupt_syst_IP[9] = 0.00102338; mupt_syst_IP[10] = 0.0123791; 
// float gluon_syst_IP[11];
gluon_syst_IP[0] = 0.00116053; gluon_syst_IP[1] = 0.000782418; gluon_syst_IP[2] = 0.000817876; 
gluon_syst_IP[3] = 0.000974364; gluon_syst_IP[4] = 0.00105004; gluon_syst_IP[5] = 0.00104748; 
gluon_syst_IP[6] = 0.00344615; gluon_syst_IP[7] = 0.00343807; gluon_syst_IP[8] = 0.0089984; 
gluon_syst_IP[9] = 0.00911566; gluon_syst_IP[10] = 0.0292369; 
// float pu_syst_IP[11];
pu_syst_IP[0] = 0.0130108; pu_syst_IP[1] = 0.00744548; pu_syst_IP[2] = 0.0077829; 
pu_syst_IP[3] = 0.00341745; pu_syst_IP[4] = 0.00368287; pu_syst_IP[5] = 0.0036739; 
pu_syst_IP[6] = 0.00593577; pu_syst_IP[7] = 0.00592186; pu_syst_IP[8] = 0.00135831; 
pu_syst_IP[9] = 0.00137601; pu_syst_IP[10] = 0.0142345; 
// float total_syst_IP[11];
total_syst_IP[0] = 0.0678007; total_syst_IP[1] = 0.0273643; total_syst_IP[2] = 0.0286044; 
total_syst_IP[3] = 0.013533; total_syst_IP[4] = 0.014584; total_syst_IP[5] = 0.0145485; 
total_syst_IP[6] = 0.0105726; total_syst_IP[7] = 0.0105478; total_syst_IP[8] = 0.0164387; 
total_syst_IP[9] = 0.0166529; total_syst_IP[10] = 0.0351279; 


// TCHEL : JP
eff_JP[0] = 0.573951; eff_stat_JP[0] = 0.0109688; effMC_JP[0] = 0.584754; effMC_stat_JP[0] = 0.00114153; 
sf_JP[0] = 0.981526; sf_stat_JP[0] = 0.0187579; sf_syst_JP[0] = 0.0981526; sf_eror_JP[0] = 0.0999289; 
mupt_syst_JP[0] = 0.005; gluon_syst_JP[0] = 0.002; pu_syst_JP[0] = 0.001; total_syst_JP[0] = 0.1; 
cor_syst_JP[0] = 0.044; inc_syst_JP[0] = 0.09; bias_syst_JP[0] = 0; 
eff_JP[1] = 0.693381; eff_stat_JP[1] = 0.00813674; effMC_JP[1] = 0.702562; effMC_stat_JP[1] = 0.000800859; 
sf_JP[1] = 0.986933; sf_stat_JP[1] = 0.0115815; sf_syst_JP[1] = 0.043425; sf_eror_JP[1] = 0.0449429; 
mupt_syst_JP[1] = 0.005; gluon_syst_JP[1] = 0; pu_syst_JP[1] = 0; total_syst_JP[1] = 0.044; 
cor_syst_JP[1] = 0.022; inc_syst_JP[1] = 0.038; bias_syst_JP[1] = 0; 
eff_JP[2] = 0.757924; eff_stat_JP[2] = 0.00668911; effMC_JP[2] = 0.764837; effMC_stat_JP[2] = 0.000799126; 
sf_JP[2] = 0.990961; sf_stat_JP[2] = 0.00874579; sf_syst_JP[2] = 0.0396384; sf_eror_JP[2] = 0.0405918; 
mupt_syst_JP[2] = 0.005; gluon_syst_JP[2] = 0.001; pu_syst_JP[2] = 0; total_syst_JP[2] = 0.04; 
cor_syst_JP[2] = 0.013; inc_syst_JP[2] = 0.038; bias_syst_JP[2] = 0; 
eff_JP[3] = 0.809754; eff_stat_JP[3] = 0.00619074; effMC_JP[3] = 0.81293; effMC_stat_JP[3] = 0.00101703; 
sf_JP[3] = 0.996093; sf_stat_JP[3] = 0.00761535; sf_syst_JP[3] = 0.0268945; sf_eror_JP[3] = 0.0279519; 
mupt_syst_JP[3] = 0.005; gluon_syst_JP[3] = 0.001; pu_syst_JP[3] = 0.001; total_syst_JP[3] = 0.027; 
cor_syst_JP[3] = 0.008; inc_syst_JP[3] = 0.025; bias_syst_JP[3] = 0; 
eff_JP[4] = 0.829075; eff_stat_JP[4] = 0.00607968; effMC_JP[4] = 0.833165; effMC_stat_JP[4] = 0.000973417; 
sf_JP[4] = 0.99509; sf_stat_JP[4] = 0.00729708; sf_syst_JP[4] = 0.0258723; sf_eror_JP[4] = 0.0268817; 
mupt_syst_JP[4] = 0.005; gluon_syst_JP[4] = 0.001; pu_syst_JP[4] = 0; total_syst_JP[4] = 0.026; 
cor_syst_JP[4] = 0.006; inc_syst_JP[4] = 0.025; bias_syst_JP[4] = 0; 
eff_JP[5] = 0.842056; eff_stat_JP[5] = 0.0051331; effMC_JP[5] = 0.84539; effMC_stat_JP[5] = 0.000983613; 
sf_JP[5] = 0.996057; sf_stat_JP[5] = 0.00607187; sf_syst_JP[5] = 0.0258975; sf_eror_JP[5] = 0.0265998; 
mupt_syst_JP[5] = 0.005; gluon_syst_JP[5] = 0.001; pu_syst_JP[5] = 0; total_syst_JP[5] = 0.026; 
cor_syst_JP[5] = 0.005; inc_syst_JP[5] = 0.025; bias_syst_JP[5] = 0; 
eff_JP[6] = 0.864526; eff_stat_JP[6] = 0.00291262; effMC_JP[6] = 0.864687; effMC_stat_JP[6] = 0.000891088; 
sf_JP[6] = 0.999814; sf_stat_JP[6] = 0.00336841; sf_syst_JP[6] = 0.0119978; sf_eror_JP[6] = 0.0124616; 
mupt_syst_JP[6] = 0.001; gluon_syst_JP[6] = 0.001; pu_syst_JP[6] = 0; total_syst_JP[6] = 0.012; 
cor_syst_JP[6] = 0.004; inc_syst_JP[6] = 0.011; bias_syst_JP[6] = 0; 
eff_JP[7] = 0.873498; eff_stat_JP[7] = 0.00555282; effMC_JP[7] = 0.874211; effMC_stat_JP[7] = 0.000866809; 
sf_JP[7] = 0.999185; sf_stat_JP[7] = 0.0063518; sf_syst_JP[7] = 0.010991; sf_eror_JP[7] = 0.0126944; 
mupt_syst_JP[7] = 0.001; gluon_syst_JP[7] = 0; pu_syst_JP[7] = 0; total_syst_JP[7] = 0.011; 
cor_syst_JP[7] = 0.003; inc_syst_JP[7] = 0.011; bias_syst_JP[7] = 0; 
eff_JP[8] = 0.884927; eff_stat_JP[8] = 0.00499561; effMC_JP[8] = 0.884691; effMC_stat_JP[8] = 0.000803766; 
sf_JP[8] = 1.00027; sf_stat_JP[8] = 0.00564673; sf_syst_JP[8] = 0.011003; sf_eror_JP[8] = 0.0123673; 
mupt_syst_JP[8] = 0.001; gluon_syst_JP[8] = 0.001; pu_syst_JP[8] = 0; total_syst_JP[8] = 0.011; 
cor_syst_JP[8] = 0.003; inc_syst_JP[8] = 0.011; bias_syst_JP[8] = 0; 
eff_JP[9] = 0.893959; eff_stat_JP[9] = 0.0155169; effMC_JP[9] = 0.892257; effMC_stat_JP[9] = 0.00094935; 
sf_JP[9] = 1.00191; sf_stat_JP[9] = 0.0173906; sf_syst_JP[9] = 0.0200382; sf_eror_JP[9] = 0.0265323; 
mupt_syst_JP[9] = 0.001; gluon_syst_JP[9] = 0.001; pu_syst_JP[9] = 0; total_syst_JP[9] = 0.02; 
cor_syst_JP[9] = 0.002; inc_syst_JP[9] = 0.02; bias_syst_JP[9] = 0; 
eff_JP[10] = 0.894889; eff_stat_JP[10] = 0.00785028; effMC_JP[10] = 0.897759; effMC_stat_JP[10] = 0.00165286; 
sf_JP[10] = 0.996804; sf_stat_JP[10] = 0.00874431; sf_syst_JP[10] = 0.0199361; sf_eror_JP[10] = 0.0217695; 
mupt_syst_JP[10] = 0.002; gluon_syst_JP[10] = 0.001; pu_syst_JP[10] = 0.001; total_syst_JP[10] = 0.02; 
cor_syst_JP[10] = 0.002; inc_syst_JP[10] = 0.02; bias_syst_JP[10] = 0; 
eff_JP[11] = 0.901045; eff_stat_JP[11] = 0.00265045; effMC_JP[11] = 0.901718; effMC_stat_JP[11] = 0.00273248; 
sf_JP[11] = 0.999253; sf_stat_JP[11] = 0.00293934; sf_syst_JP[11] = 0.0199851; sf_eror_JP[11] = 0.0202001; 
mupt_syst_JP[11] = 0.002; gluon_syst_JP[11] = 0.003; pu_syst_JP[11] = 0.001; total_syst_JP[11] = 0.02; 
cor_syst_JP[11] = 0.001; inc_syst_JP[11] = 0.02; bias_syst_JP[11] = 0; 
eff_JP[12] = 0.916023; eff_stat_JP[12] = 0.00798108; effMC_JP[12] = 0.907837; effMC_stat_JP[12] = 0.00451277; 
sf_JP[12] = 1.00902; sf_stat_JP[12] = 0.00879131; sf_syst_JP[12] = 0.0171533; sf_eror_JP[12] = 0.019275; 
mupt_syst_JP[12] = 0.002; gluon_syst_JP[12] = 0.005; pu_syst_JP[12] = 0.001; total_syst_JP[12] = 0.017; 
cor_syst_JP[12] = 0.001; inc_syst_JP[12] = 0.016; bias_syst_JP[12] = 0; 
eff_JP[13] = 0.939532; eff_stat_JP[13] = 0.05357; effMC_JP[13] = 0.917675; effMC_stat_JP[13] = 0.00376296; 
sf_JP[13] = 1.02382; sf_stat_JP[13] = 0.0583758; sf_syst_JP[13] = 0.0215002; sf_eror_JP[13] = 0.0622093; 
mupt_syst_JP[13] = 0.002; gluon_syst_JP[13] = 0.013; pu_syst_JP[13] = 0.003; total_syst_JP[13] = 0.021; 
cor_syst_JP[13] = 0.001; inc_syst_JP[13] = 0.016; bias_syst_JP[13] = 0; 
eff_JP[14] = 0.92444; eff_stat_JP[14] = 0.0114244; effMC_JP[14] = 0.903687; effMC_stat_JP[14] = 0.00341195; 
sf_JP[14] = 1.02296; sf_stat_JP[14] = 0.012642; sf_syst_JP[14] = 0.0276199; sf_eror_JP[14] = 0.0303757; 
mupt_syst_JP[14] = 0.002; gluon_syst_JP[14] = 0.019; pu_syst_JP[14] = 0.005; total_syst_JP[14] = 0.027; 
cor_syst_JP[14] = 0.008; inc_syst_JP[14] = 0.016; bias_syst_JP[14] = 0; 
}

// #############################################################################
if ( title == "TCHEM" ) { 
// TCHEM : PTREL
// float eff_PT[11];
eff_PT[0] = 0.34077; eff_PT[1] = 0.489312; eff_PT[2] = 0.571932; 
eff_PT[3] = 0.602859; eff_PT[4] = 0.639993; eff_PT[5] = 0.661936; 
eff_PT[6] = 0.681971; eff_PT[7] = 0.69443; eff_PT[8] = 0.70993; 
eff_PT[9] = 0.699618; eff_PT[10] = 0.650965; 
// float eff_stat_PT[11];
eff_stat_PT[0] = 0.00953315; eff_stat_PT[1] = 0.00755767; eff_stat_PT[2] = 0.00857883; 
eff_stat_PT[3] = 0.00877433; eff_stat_PT[4] = 0.0106868; eff_stat_PT[5] = 0.013988; 
eff_stat_PT[6] = 0.0110117; eff_stat_PT[7] = 0.012113; eff_stat_PT[8] = 0.0146398; 
eff_stat_PT[9] = 0.0129272; eff_stat_PT[10] = 0.0254389; 
// float effMC_PT[11];
effMC_PT[0] = 0.396447; effMC_PT[1] = 0.522805; effMC_PT[2] = 0.599218; 
effMC_PT[3] = 0.65822; effMC_PT[4] = 0.681362; effMC_PT[5] = 0.707244; 
effMC_PT[6] = 0.734159; effMC_PT[7] = 0.733808; effMC_PT[8] = 0.754462; 
effMC_PT[9] = 0.754612; effMC_PT[10] = 0.766522; 
// float effMC_stat_PT[11];
effMC_stat_PT[0] = 0.00332124; effMC_stat_PT[1] = 0.00213572; effMC_stat_PT[2] = 0.00189972; 
effMC_stat_PT[3] = 0.00256284; effMC_stat_PT[4] = 0.00246368; effMC_stat_PT[5] = 0.00253506; 
effMC_stat_PT[6] = 0.00267056; effMC_stat_PT[7] = 0.00253602; effMC_stat_PT[8] = 0.00325604; 
effMC_stat_PT[9] = 0.00341624; effMC_stat_PT[10] = 0.00508199; 
// float sf_PT[11];
sf_PT[0] = 0.859561; sf_PT[1] = 0.935937; sf_PT[2] = 0.954464; 
sf_PT[3] = 0.915893; sf_PT[4] = 0.939284; sf_PT[5] = 0.935937; 
sf_PT[6] = 0.928915; sf_PT[7] = 0.946338; sf_PT[8] = 0.940975; 
sf_PT[9] = 0.927122; sf_PT[10] = 0.849245; 
// float sf_stat_PT[11];
sf_stat_PT[0] = 0.0251015; sf_stat_PT[1] = 0.0149531; sf_stat_PT[2] = 0.014633; 
sf_stat_PT[3] = 0.0137992; sf_stat_PT[4] = 0.0160479; sf_stat_PT[5] = 0.0200606; 
sf_stat_PT[6] = 0.015375; sf_stat_PT[7] = 0.016828; sf_stat_PT[8] = 0.0198247; 
sf_stat_PT[9] = 0.0176376; sf_stat_PT[10] = 0.0336617; 
// float away_syst_PT[11];
away_syst_PT[0] = 0.019694; away_syst_PT[1] = 0.0373128; away_syst_PT[2] = 0.0380514; 
away_syst_PT[3] = 0.0201917; away_syst_PT[4] = 0.0207073; away_syst_PT[5] = 0.0206336; 
away_syst_PT[6] = 0.0147936; away_syst_PT[7] = 0.0150711; away_syst_PT[8] = 0.0547567; 
away_syst_PT[9] = 0.0539506; away_syst_PT[10] = 0.079519; 
// float mupt_syst_PT[11];
mupt_syst_PT[0] = 0.017437; mupt_syst_PT[1] = 0.00761788; mupt_syst_PT[2] = 0.00776868; 
mupt_syst_PT[3] = 0.0105092; mupt_syst_PT[4] = 0.0107776; mupt_syst_PT[5] = 0.0107392; 
mupt_syst_PT[6] = 0.0126441; mupt_syst_PT[7] = 0.0128813; mupt_syst_PT[8] = 0.0689595; 
mupt_syst_PT[9] = 0.0679442; mupt_syst_PT[10] = 0.025804; 
// float lc_syst_PT[11];
lc_syst_PT[0] = 0.005247; lc_syst_PT[1] = 0.00499526; lc_syst_PT[2] = 0.00509415; 
lc_syst_PT[3] = 0.00633918; lc_syst_PT[4] = 0.00650108; lc_syst_PT[5] = 0.00647791; 
lc_syst_PT[6] = 0.00411728; lc_syst_PT[7] = 0.00419451; lc_syst_PT[8] = 0.0178712; 
lc_syst_PT[9] = 0.0176081; lc_syst_PT[10] = 0.037425; 
// float gluon_syst_PT[11];
gluon_syst_PT[0] = 0.00197804; gluon_syst_PT[1] = 0.000385834; gluon_syst_PT[2] = 0.000393471; 
gluon_syst_PT[3] = 0.00193052; gluon_syst_PT[4] = 0.00197983; gluon_syst_PT[5] = 0.00197277; 
gluon_syst_PT[6] = 0.004763; gluon_syst_PT[7] = 0.00485233; gluon_syst_PT[8] = 0.0132615; 
gluon_syst_PT[9] = 0.0130663; gluon_syst_PT[10] = 0.055487; 
// float pu_syst_PT[11];
pu_syst_PT[0] = 0.01818; pu_syst_PT[1] = 0.0126419; pu_syst_PT[2] = 0.0128922; 
pu_syst_PT[3] = 0.00602235; pu_syst_PT[4] = 0.00617616; pu_syst_PT[5] = 0.00615415; 
pu_syst_PT[6] = 0.0140817; pu_syst_PT[7] = 0.0143458; pu_syst_PT[8] = 0.008517; 
pu_syst_PT[9] = 0.00839161; pu_syst_PT[10] = 0.016941; 
// float total_syst_PT[11];
total_syst_PT[0] = 0.0324632; total_syst_PT[1] = 0.0404375; total_syst_PT[2] = 0.041238; 
total_syst_PT[3] = 0.0244607; total_syst_PT[4] = 0.0250854; total_syst_PT[5] = 0.024996; 
total_syst_PT[6] = 0.0248325; total_syst_PT[7] = 0.0252983; total_syst_PT[8] = 0.0912222; 
total_syst_PT[9] = 0.0898793; total_syst_PT[10] = 0.108423; 

// TCHEM : S8
eff_SY[0]=0.3436; eff_SY[1]=0.4739; eff_SY[2]=0.5506; eff_SY[3]=0.6163; eff_SY[4]=0.6368; eff_SY[5]=0.7352; eff_SY[6]=0.7300; eff_SY[7]=0.6709; 
eff_stat_SY[0]=0.0089; eff_stat_SY[1]=0.0072; eff_stat_SY[2]=0.0092; eff_stat_SY[3]=0.0106; eff_stat_SY[4]=0.0122; eff_stat_SY[5]=0.0257; eff_stat_SY[6]=0.0373; eff_stat_SY[7]=0.0315; 
effMC_SY[0]=0.3692; effMC_SY[1]=0.4985; effMC_SY[2]=0.5876; effMC_SY[3]=0.6496; effMC_SY[4]=0.6965; effMC_SY[5]=0.7382; effMC_SY[6]=0.7558; effMC_SY[7]=0.7716; 
effMC_stat_SY[0]=0.0006; effMC_stat_SY[1]=0.0007; effMC_stat_SY[2]=0.0008; effMC_stat_SY[3]=0.0009; effMC_stat_SY[4]=0.0007; effMC_stat_SY[5]=0.0009; effMC_stat_SY[6]=0.0010; effMC_stat_SY[7]=0.0009; 
sf_SY[0]=0.9307; sf_SY[1]=0.9507; sf_SY[2]=0.937; sf_SY[3]=0.9487; sf_SY[4]=0.9143; sf_SY[5]=0.9959; sf_SY[6]=0.9659; sf_SY[7]=0.8695; 
sf_stat_SY[0]=0.0242; sf_stat_SY[1]=0.0145; sf_stat_SY[2]=0.0157; sf_stat_SY[3]=0.0175; sf_stat_SY[4]=0.0182; sf_stat_SY[5]=0.0348; sf_stat_SY[6]=0.0494; sf_stat_SY[7]=0.0408; 
away[0]=2.8995; away[1]=1.2697; away[2]=2.1035; away[3]=13.0984; 
mu[0]=0.3924; mu[1]=0.7554; mu[2]=1.4254; mu[3]=1.5462; 
ptrel[0]=0.7412; ptrel[1]=3.1019; ptrel[2]=2.3803; ptrel[3]=2.8992; 
gsplit[0]=0.1637; gsplit[1]=0.5642; gsplit[2]=1.3025; gsplit[3]=1.0283; 
closure[0]=0.2456; closure[1]=0.0891; closure[2]=0.0269; closure[3]=0.4499; 
pu[0]=0.0436; pu[1]=0.0804; pu[2]=0.8995; pu[3]=0.1784; 

// TCHEM : IP3D
// float eff_IP[11];
eff_IP[0] = 0.318364; eff_IP[1] = 0.45281; eff_IP[2] = 0.527138; 
eff_IP[3] = 0.550776; eff_IP[4] = 0.630089; eff_IP[5] = 0.651658; 
eff_IP[6] = 0.662487; eff_IP[7] = 0.666643; eff_IP[8] = 0.714593; 
eff_IP[9] = 0.719321; eff_IP[10] = 0.78074; 
// float eff_stat_IP[11];
eff_stat_IP[0] = 0.0133702; eff_stat_IP[1] = 0.0134002; eff_stat_IP[2] = 0.0141022; 
eff_stat_IP[3] = 0.0127196; eff_stat_IP[4] = 0.0160889; eff_stat_IP[5] = 0.0190044; 
eff_stat_IP[6] = 0.0103553; eff_stat_IP[7] = 0.0116975; eff_stat_IP[8] = 0.0098051; 
eff_stat_IP[9] = 0.0110933; eff_stat_IP[10] = 0.0204628; 
// float effMC_IP[11];
effMC_IP[0] = 0.399137; effMC_IP[1] = 0.525989; effMC_IP[2] = 0.599272; 
effMC_IP[3] = 0.659541; effMC_IP[4] = 0.679442; effMC_IP[5] = 0.707365; 
effMC_IP[6] = 0.731063; effMC_IP[7] = 0.731752; effMC_IP[8] = 0.751682; 
effMC_IP[9] = 0.751315; effMC_IP[10] = 0.76939; 
// float effMC_stat_IP[11];
effMC_stat_IP[0] = 0.00330275; effMC_stat_IP[1] = 0.00228512; effMC_stat_IP[2] = 0.00200025; 
effMC_stat_IP[3] = 0.00254801; effMC_stat_IP[4] = 0.00242878; effMC_stat_IP[5] = 0.00249005; 
effMC_stat_IP[6] = 0.00260406; effMC_stat_IP[7] = 0.00246057; effMC_stat_IP[8] = 0.00315256; 
effMC_stat_IP[9] = 0.00329064; effMC_stat_IP[10] = 0.00508227; 
// float sf_IP[11];
sf_IP[0] = 0.797631; sf_IP[1] = 0.860874; sf_IP[2] = 0.87963; 
sf_IP[3] = 0.83509; sf_IP[4] = 0.927363; sf_IP[5] = 0.921247; 
sf_IP[6] = 0.906196; sf_IP[7] = 0.911022; sf_IP[8] = 0.950659; 
sf_IP[9] = 0.957416; sf_IP[10] = 1.01475; 
// float sf_stat_IP[11];
sf_stat_IP[0] = 0.0341418; sf_stat_IP[1] = 0.0257493; sf_stat_IP[2] = 0.0237147; 
sf_stat_IP[3] = 0.0195535; sf_stat_IP[4] = 0.0239105; sf_stat_IP[5] = 0.0270615; 
sf_stat_IP[6] = 0.0145278; sf_stat_IP[7] = 0.0162764; sf_stat_IP[8] = 0.0136399; 
sf_stat_IP[9] = 0.0153491; sf_stat_IP[10] = 0.0274278; 
// float away_syst_IP[11];
away_syst_IP[0] = 0.0151572; away_syst_IP[1] = 0.0050607; away_syst_IP[2] = 0.00517096; 
away_syst_IP[3] = 0.032665; away_syst_IP[4] = 0.0362744; away_syst_IP[5] = 0.0360351; 
away_syst_IP[6] = 0.0177862; away_syst_IP[7] = 0.0178809; away_syst_IP[8] = 0.0231102; 
away_syst_IP[9] = 0.0232744; away_syst_IP[10] = 0.029016; 
// float mupt_syst_IP[11];
mupt_syst_IP[0] = 0.0244011; mupt_syst_IP[1] = 0.0230457; mupt_syst_IP[2] = 0.0235478; 
mupt_syst_IP[3] = 0.00568876; mupt_syst_IP[4] = 0.00631734; mupt_syst_IP[5] = 0.00627568; 
mupt_syst_IP[6] = 0.0103142; mupt_syst_IP[7] = 0.0103691; mupt_syst_IP[8] = 0.0042173; 
mupt_syst_IP[9] = 0.00424727; mupt_syst_IP[10] = 0.00829601; 
// float gluon_syst_IP[11];
gluon_syst_IP[0] = 0.0273083; gluon_syst_IP[1] = 0.000445379; gluon_syst_IP[2] = 0.000455083; 
gluon_syst_IP[3] = 0.000199258; gluon_syst_IP[4] = 0.000221275; gluon_syst_IP[5] = 0.000219816; 
gluon_syst_IP[6] = 0.00369742; gluon_syst_IP[7] = 0.00371711; gluon_syst_IP[8] = 0.0110566; 
gluon_syst_IP[9] = 0.0111352; gluon_syst_IP[10] = 0.064061; 
// float pu_syst_IP[11];
pu_syst_IP[0] = 0.0413443; pu_syst_IP[1] = 0.0280925; pu_syst_IP[2] = 0.0287046; 
pu_syst_IP[3] = 0.0126666; pu_syst_IP[4] = 0.0140662; pu_syst_IP[5] = 0.0139734; 
pu_syst_IP[6] = 0.00495234; pu_syst_IP[7] = 0.00497871; pu_syst_IP[8] = 0.000860651; 
pu_syst_IP[9] = 0.000866768; pu_syst_IP[10] = 0.0023264; 
// float total_syst_IP[11];
total_syst_IP[0] = 0.0572734; total_syst_IP[1] = 0.0366893; total_syst_IP[2] = 0.0374886; 
total_syst_IP[3] = 0.0354944; total_syst_IP[4] = 0.0394163; total_syst_IP[5] = 0.0391563; 
total_syst_IP[6] = 0.0214692; total_syst_IP[7] = 0.0215836; total_syst_IP[8] = 0.025978; 
total_syst_IP[9] = 0.0261626; total_syst_IP[10] = 0.0708518; 

// TCHEM : JP
eff_JP[0] = 0.41382; eff_stat_JP[0] = 0.0130397; effMC_JP[0] = 0.453084; effMC_stat_JP[0] = 0.00115318; 
sf_JP[0] = 0.91334; sf_stat_JP[0] = 0.0287798; sf_syst_JP[0] = 0.122388; sf_eror_JP[0] = 0.125726; 
mupt_syst_JP[0] = 0.006; gluon_syst_JP[0] = 0.003; pu_syst_JP[0] = 0.001; total_syst_JP[0] = 0.134; 
cor_syst_JP[0] = 0.044; inc_syst_JP[0] = 0.126; bias_syst_JP[0] = 0; 
eff_JP[1] = 0.534968; eff_stat_JP[1] = 0.00988017; effMC_JP[1] = 0.563854; effMC_stat_JP[1] = 0.000868789; 
sf_JP[1] = 0.94877; sf_stat_JP[1] = 0.0175225; sf_syst_JP[1] = 0.0521824; sf_eror_JP[1] = 0.0550458; 
mupt_syst_JP[1] = 0.006; gluon_syst_JP[1] = 0.003; pu_syst_JP[1] = 0.001; total_syst_JP[1] = 0.055; 
cor_syst_JP[1] = 0.022; inc_syst_JP[1] = 0.05; bias_syst_JP[1] = 0; 
eff_JP[2] = 0.602582; eff_stat_JP[2] = 0.00876745; effMC_JP[2] = 0.629788; effMC_stat_JP[2] = 0.000909848; 
sf_JP[2] = 0.956801; sf_stat_JP[2] = 0.0139213; sf_syst_JP[2] = 0.0497537; sf_eror_JP[2] = 0.0516646; 
mupt_syst_JP[2] = 0.006; gluon_syst_JP[2] = 0.003; pu_syst_JP[2] = 0.001; total_syst_JP[2] = 0.052; 
cor_syst_JP[2] = 0.013; inc_syst_JP[2] = 0.05; bias_syst_JP[2] = 0; 
eff_JP[3] = 0.662443; eff_stat_JP[3] = 0.00912239; effMC_JP[3] = 0.684347; effMC_stat_JP[3] = 0.00121212; 
sf_JP[3] = 0.967993; sf_stat_JP[3] = 0.0133301; sf_syst_JP[3] = 0.0338798; sf_eror_JP[3] = 0.0364078; 
mupt_syst_JP[3] = 0.006; gluon_syst_JP[3] = 0.005; pu_syst_JP[3] = 0.001; total_syst_JP[3] = 0.035; 
cor_syst_JP[3] = 0.008; inc_syst_JP[3] = 0.033; bias_syst_JP[3] = 0; 
eff_JP[4] = 0.685272; eff_stat_JP[4] = 0.00956307; effMC_JP[4] = 0.710022; effMC_stat_JP[4] = 0.0011847; 
sf_JP[4] = 0.965141; sf_stat_JP[4] = 0.0134687; sf_syst_JP[4] = 0.0337799; sf_eror_JP[4] = 0.0363661; 
mupt_syst_JP[4] = 0.006; gluon_syst_JP[4] = 0.006; pu_syst_JP[4] = 0; total_syst_JP[4] = 0.035; 
cor_syst_JP[4] = 0.006; inc_syst_JP[4] = 0.033; bias_syst_JP[4] = 0; 
eff_JP[5] = 0.702908; eff_stat_JP[5] = 0.00949777; effMC_JP[5] = 0.724984; effMC_stat_JP[5] = 0.00121484; 
sf_JP[5] = 0.969549; sf_stat_JP[5] = 0.0131007; sf_syst_JP[5] = 0.0339342; sf_eror_JP[5] = 0.0363753; 
mupt_syst_JP[5] = 0.006; gluon_syst_JP[5] = 0.007; pu_syst_JP[5] = 0.001; total_syst_JP[5] = 0.035; 
cor_syst_JP[5] = 0.005; inc_syst_JP[5] = 0.033; bias_syst_JP[5] = 0; 
eff_JP[6] = 0.735104; eff_stat_JP[6] = 0.00788351; effMC_JP[6] = 0.752813; effMC_stat_JP[6] = 0.00112377; 
sf_JP[6] = 0.976476; sf_stat_JP[6] = 0.0104721; sf_syst_JP[6] = 0.0283178; sf_eror_JP[6] = 0.0301921; 
mupt_syst_JP[6] = 0.003; gluon_syst_JP[6] = 0.007; pu_syst_JP[6] = 0; total_syst_JP[6] = 0.029; 
cor_syst_JP[6] = 0.004; inc_syst_JP[6] = 0.028; bias_syst_JP[6] = 0; 
eff_JP[7] = 0.743925; eff_stat_JP[7] = 0.00967973; effMC_JP[7] = 0.762839; effMC_stat_JP[7] = 0.00111181; 
sf_JP[7] = 0.975206; sf_stat_JP[7] = 0.0126891; sf_syst_JP[7] = 0.028281; sf_eror_JP[7] = 0.0309972; 
mupt_syst_JP[7] = 0.003; gluon_syst_JP[7] = 0.008; pu_syst_JP[7] = 0; total_syst_JP[7] = 0.029; 
cor_syst_JP[7] = 0.003; inc_syst_JP[7] = 0.028; bias_syst_JP[7] = 0; 
eff_JP[8] = 0.762221; eff_stat_JP[8] = 0.0103064; effMC_JP[8] = 0.777947; effMC_stat_JP[8] = 0.00104594; 
sf_JP[8] = 0.979784; sf_stat_JP[8] = 0.0132481; sf_syst_JP[8] = 0.0284137; sf_eror_JP[8] = 0.0313505; 
mupt_syst_JP[8] = 0.003; gluon_syst_JP[8] = 0.008; pu_syst_JP[8] = 0.001; total_syst_JP[8] = 0.029; 
cor_syst_JP[8] = 0.003; inc_syst_JP[8] = 0.028; bias_syst_JP[8] = 0; 
eff_JP[9] = 0.763985; eff_stat_JP[9] = 0.00984478; effMC_JP[9] = 0.781977; effMC_stat_JP[9] = 0.00126425; 
sf_JP[9] = 0.976992; sf_stat_JP[9] = 0.0125896; sf_syst_JP[9] = 0.0302868; sf_eror_JP[9] = 0.0327992; 
mupt_syst_JP[9] = 0.003; gluon_syst_JP[9] = 0.013; pu_syst_JP[9] = 0; total_syst_JP[9] = 0.031; 
cor_syst_JP[9] = 0.002; inc_syst_JP[9] = 0.028; bias_syst_JP[9] = 0; 
eff_JP[10] = 0.761731; eff_stat_JP[10] = 0.0136492; effMC_JP[10] = 0.787451; effMC_stat_JP[10] = 0.00223194; 
sf_JP[10] = 0.967337; sf_stat_JP[10] = 0.0173334; sf_syst_JP[10] = 0.0309548; sf_eror_JP[10] = 0.0354774; 
mupt_syst_JP[10] = 0.003; gluon_syst_JP[10] = 0.016; pu_syst_JP[10] = 0.001; total_syst_JP[10] = 0.032; 
cor_syst_JP[10] = 0.002; inc_syst_JP[10] = 0.028; bias_syst_JP[10] = 0; 
eff_JP[11] = 0.751687; eff_stat_JP[11] = 0.0172777; effMC_JP[11] = 0.785563; effMC_stat_JP[11] = 0.00376725; 
sf_JP[11] = 0.956876; sf_stat_JP[11] = 0.0219941; sf_syst_JP[11] = 0.0325338; sf_eror_JP[11] = 0.0392707; 
mupt_syst_JP[11] = 0.003; gluon_syst_JP[11] = 0.019; pu_syst_JP[11] = 0.003; total_syst_JP[11] = 0.034; 
cor_syst_JP[11] = 0.001; inc_syst_JP[11] = 0.028; bias_syst_JP[11] = 0; 
eff_JP[12] = 0.793312; eff_stat_JP[12] = 0.0168; effMC_JP[12] = 0.799752; effMC_stat_JP[12] = 0.00624341; 
sf_JP[12] = 0.991947; sf_stat_JP[12] = 0.0210066; sf_syst_JP[12] = 0.0238067; sf_eror_JP[12] = 0.0317496; 
mupt_syst_JP[12] = 0.003; gluon_syst_JP[12] = 0.006; pu_syst_JP[12] = 0.003; total_syst_JP[12] = 0.024; 
cor_syst_JP[12] = 0.001; inc_syst_JP[12] = 0.023; bias_syst_JP[12] = 0; 
eff_JP[13] = 0.783641; eff_stat_JP[13] = 0.0241804; effMC_JP[13] = 0.784733; effMC_stat_JP[13] = 0.0056269; 
sf_JP[13] = 0.998608; sf_stat_JP[13] = 0.0308135; sf_syst_JP[13] = 0.0259638; sf_eror_JP[13] = 0.0402938; 
mupt_syst_JP[13] = 0.003; gluon_syst_JP[13] = 0.004; pu_syst_JP[13] = 0.01; total_syst_JP[13] = 0.026; 
cor_syst_JP[13] = 0.001; inc_syst_JP[13] = 0.023; bias_syst_JP[13] = 0; 
eff_JP[14] = 0.734029; eff_stat_JP[14] = 0.0451813; effMC_JP[14] = 0.777839; effMC_stat_JP[14] = 0.00480761; 
sf_JP[14] = 0.943677; sf_stat_JP[14] = 0.0580857; sf_syst_JP[14] = 0.0235919; sf_eror_JP[14] = 0.0626939; 
mupt_syst_JP[14] = 0.003; gluon_syst_JP[14] = 0.003; pu_syst_JP[14] = 0.006; total_syst_JP[14] = 0.025; 
cor_syst_JP[14] = 0.008; inc_syst_JP[14] = 0.023; bias_syst_JP[14] = 0; 
}

// #############################################################################
if ( title == "TCHPM" ) { 
// TCHPM : PTREL
// float eff_PT[11];
eff_PT[0] = 0.163318; eff_PT[1] = 0.294264; eff_PT[2] = 0.402482; 
eff_PT[3] = 0.456872; eff_PT[4] = 0.485777; eff_PT[5] = 0.51893; 
eff_PT[6] = 0.572201; eff_PT[7] = 0.577089; eff_PT[8] = 0.606924; 
eff_PT[9] = 0.646563; eff_PT[10] = 0.718029; 
// float eff_stat_PT[11];
eff_stat_PT[0] = 0.00674078; eff_stat_PT[1] = 0.00639266; eff_stat_PT[2] = 0.00795556; 
eff_stat_PT[3] = 0.0084692; eff_stat_PT[4] = 0.0102179; eff_stat_PT[5] = 0.01178; 
eff_stat_PT[6] = 0.00872078; eff_stat_PT[7] = 0.0104439; eff_stat_PT[8] = 0.0122929; 
eff_stat_PT[9] = 0.0119566; eff_stat_PT[10] = 0.0255568; 
// float effMC_PT[11];
effMC_PT[0] = 0.197494; effMC_PT[1] = 0.332986; effMC_PT[2] = 0.423897; 
effMC_PT[3] = 0.495089; effMC_PT[4] = 0.52745; effMC_PT[5] = 0.567049; 
effMC_PT[6] = 0.600704; effMC_PT[7] = 0.605322; effMC_PT[8] = 0.643303; 
effMC_PT[9] = 0.660855; effMC_PT[10] = 0.712199; 
// float effMC_stat_PT[11];
effMC_stat_PT[0] = 0.00263915; effMC_stat_PT[1] = 0.00200112; effMC_stat_PT[2] = 0.00192065; 
effMC_stat_PT[3] = 0.00272354; effMC_stat_PT[4] = 0.00267906; effMC_stat_PT[5] = 0.00278821; 
effMC_stat_PT[6] = 0.00300079; effMC_stat_PT[7] = 0.00285985; effMC_stat_PT[8] = 0.00370808; 
effMC_stat_PT[9] = 0.00384331; effMC_stat_PT[10] = 0.00542885; 
// float sf_PT[11];
sf_PT[0] = 0.826952; sf_PT[1] = 0.883714; sf_PT[2] = 0.94948; 
sf_PT[3] = 0.922809; sf_PT[4] = 0.92099; sf_PT[5] = 0.915142; 
sf_PT[6] = 0.952551; sf_PT[7] = 0.953359; sf_PT[8] = 0.943449; 
sf_PT[9] = 0.978374; sf_PT[10] = 1.00819; 
// float sf_stat_PT[11];
sf_stat_PT[0] = 0.035876; sf_stat_PT[1] = 0.019919; sf_stat_PT[2] = 0.0192544; 
sf_stat_PT[3] = 0.0178438; sf_stat_PT[4] = 0.0199291; sf_stat_PT[5] = 0.021256; 
sf_stat_PT[6] = 0.0152775; sf_stat_PT[7] = 0.0178317; sf_stat_PT[8] = 0.0198677; 
sf_stat_PT[9] = 0.0189662; sf_stat_PT[10] = 0.036698; 
// float away_syst_PT[11];
away_syst_PT[0] = 0.00730699; away_syst_PT[1] = 0.0330805; away_syst_PT[2] = 0.0355423; 
away_syst_PT[3] = 0.0179172; away_syst_PT[4] = 0.0178819; away_syst_PT[5] = 0.0177684; 
away_syst_PT[6] = 0.0104786; away_syst_PT[7] = 0.0104875; away_syst_PT[8] = 0.0677157; 
away_syst_PT[9] = 0.0702224; away_syst_PT[10] = 0.14203; 
// float mupt_syst_PT[11];
mupt_syst_PT[0] = 0.056975; mupt_syst_PT[1] = 0.00898395; mupt_syst_PT[2] = 0.00965254; 
mupt_syst_PT[3] = 0.0178842; mupt_syst_PT[4] = 0.0178489; mupt_syst_PT[5] = 0.0177356; 
mupt_syst_PT[6] = 0.02162; mupt_syst_PT[7] = 0.0216384; mupt_syst_PT[8] = 0.0636975; 
mupt_syst_PT[9] = 0.0660555; mupt_syst_PT[10] = 0.060555; 
// float lc_syst_PT[11];
lc_syst_PT[0] = 0.00468296; lc_syst_PT[1] = 0.00441343; lc_syst_PT[2] = 0.00474187; 
lc_syst_PT[3] = 0.00647476; lc_syst_PT[4] = 0.00646199; lc_syst_PT[5] = 0.00642096; 
lc_syst_PT[6] = 0.00535911; lc_syst_PT[7] = 0.00536366; lc_syst_PT[8] = 0.0129067; 
lc_syst_PT[9] = 0.0133845; lc_syst_PT[10] = 0.03894; 
// float gluon_syst_PT[11];
gluon_syst_PT[0] = 0.00205201; gluon_syst_PT[1] = 0.00116751; gluon_syst_PT[2] = 0.0012544; 
gluon_syst_PT[3] = 0.00745109; gluon_syst_PT[4] = 0.0074364; gluon_syst_PT[5] = 0.00738918; 
gluon_syst_PT[6] = 0.000887365; gluon_syst_PT[7] = 0.000888117; gluon_syst_PT[8] = 0.00183609; 
gluon_syst_PT[9] = 0.00190406; gluon_syst_PT[10] = 0.0964099; 
// float pu_syst_PT[11];
pu_syst_PT[0] = 0.00785202; pu_syst_PT[1] = 0.016921; pu_syst_PT[2] = 0.0181803; 
pu_syst_PT[3] = 0.0126174; pu_syst_PT[4] = 0.0125925; pu_syst_PT[5] = 0.0125126; 
pu_syst_PT[6] = 0.0102141; pu_syst_PT[7] = 0.0102228; pu_syst_PT[8] = 0.0105742; 
pu_syst_PT[9] = 0.0109656; pu_syst_PT[10] = 0.037313; 
// float total_syst_PT[11];
total_syst_PT[0] = 0.0582008; total_syst_PT[1] = 0.0384992; total_syst_PT[2] = 0.0413643; 
total_syst_PT[3] = 0.0299585; total_syst_PT[4] = 0.0298994; total_syst_PT[5] = 0.0297096; 
total_syst_PT[6] = 0.0266658; total_syst_PT[7] = 0.0266884; total_syst_PT[8] = 0.0944698; 
total_syst_PT[9] = 0.097967; total_syst_PT[10] = 0.18985; 

// TCHPM : S8
eff_SY[0]=0.1667; eff_SY[1]=0.2835; eff_SY[2]=0.3869; eff_SY[3]=0.4609; eff_SY[4]=0.5040;eff_SY[5]=0.5733; eff_SY[6]=0.6215; eff_SY[7]=0.5821; 
eff_stat_SY[0]=0.0065; eff_stat_SY[1]=0.0061; eff_stat_SY[2]=0.0084; eff_stat_SY[3]=0.0095; eff_stat_SY[4]=0.0106; eff_stat_SY[5]=0.0171; eff_stat_SY[6]=0.0267; eff_stat_SY[7]=0.0209; 
effMC_SY[0]=0.1984; effMC_SY[1]=0.3241; effMC_SY[2]=0.4229; effMC_SY[3]=0.4976; effMC_SY[4]=0.5599; effMC_SY[5]=0.6177; effMC_SY[6]=0.6462; effMC_SY[7]=0.6750; 
effMC_stat_SY[0]=0.0005; effMC_stat_SY[1]=0.0006; effMC_stat_SY[2]=0.0007; effMC_stat_SY[3]=0.0009; effMC_stat_SY[4]=0.0008; effMC_stat_SY[5]=0.0010; effMC_stat_SY[6]=0.0011; effMC_stat_SY[7]=0.0010; 
sf_SY[0]=0.8402; sf_SY[1]=0.8747; sf_SY[2]=0.9149; sf_SY[3]=0.9262; sf_SY[4]=0.9002; sf_SY[5]=0.9281; sf_SY[6]=0.9618; sf_SY[7]=0.8624; 
sf_stat_SY[0]=0.0328; sf_stat_SY[1]=0.0189; sf_stat_SY[2]=0.0199; sf_stat_SY[3]=0.0192; sf_stat_SY[4]=0.0190; sf_stat_SY[5]=0.0277; sf_stat_SY[6]=0.0414; sf_stat_SY[7]=0.031; 
away[0]=2.2743; away[1]=1.8538; away[2]=0.5296; away[3]=16.1023; 
mu[0]=0.5686; mu[1]=1.3539; mu[2]=1.008; mu[3]=0.332; 
ptrel[0]=0.7818; ptrel[1]=1.0831; ptrel[2]=3.9638; ptrel[3]=0.415; 
gsplit[0]=0.7217; gsplit[1]=1.5478; gsplit[2]=2.655; gsplit[3]=2.4135; 
closure[0]=0.1569; closure[1]=0.0566; closure[2]=0.0636; closure[3]=0.8287; 
pu[0]=0.1421; pu[1]=0.1875; pu[2]=1.1447; pu[3]=0.5976; 

// TCHPM : IP3D
// float eff_IP[11];
eff_IP[0] = 0.16432; eff_IP[1] = 0.300007; eff_IP[2] = 0.392573; 
eff_IP[3] = 0.432895; eff_IP[4] = 0.485408; eff_IP[5] = 0.518872; 
eff_IP[6] = 0.537421; eff_IP[7] = 0.561785; eff_IP[8] = 0.598967; 
eff_IP[9] = 0.63815; eff_IP[10] = 0.738518; 
// float eff_stat_IP[11];
eff_stat_IP[0] = 0.00677816; eff_stat_IP[1] = 0.0088027; eff_stat_IP[2] = 0.00943088; 
eff_stat_IP[3] = 0.00993591; eff_stat_IP[4] = 0.0120955; eff_stat_IP[5] = 0.0151593; 
eff_stat_IP[6] = 0.00879021; eff_stat_IP[7] = 0.0103056; eff_stat_IP[8] = 0.00935007; 
eff_stat_IP[9] = 0.0104884; eff_stat_IP[10] = 0.0205236; 
// float effMC_IP[11];
effMC_IP[0] = 0.200753; effMC_IP[1] = 0.335528; effMC_IP[2] = 0.424378; 
effMC_IP[3] = 0.500054; effMC_IP[4] = 0.52891; effMC_IP[5] = 0.565159; 
effMC_IP[6] = 0.599601; effMC_IP[7] = 0.606688; effMC_IP[8] = 0.640417; 
effMC_IP[9] = 0.658947; effMC_IP[10] = 0.715725; 
// float effMC_stat_IP[11];
effMC_stat_IP[0] = 0.00261952; effMC_stat_IP[1] = 0.00213704; effMC_stat_IP[2] = 0.00201594; 
effMC_stat_IP[3] = 0.0027038; effMC_stat_IP[4] = 0.00263637; effMC_stat_IP[5] = 0.00273577; 
effMC_stat_IP[6] = 0.00292415; effMC_stat_IP[7] = 0.00277751; effMC_stat_IP[8] = 0.00359106; 
effMC_stat_IP[9] = 0.00370817; effMC_stat_IP[10] = 0.00542852; 
// float sf_IP[11];
sf_IP[0] = 0.818519; sf_IP[1] = 0.894132; sf_IP[2] = 0.925056; 
sf_IP[3] = 0.865698; sf_IP[4] = 0.917753; sf_IP[5] = 0.918099; 
sf_IP[6] = 0.896297; sf_IP[7] = 0.925988; sf_IP[8] = 0.935277; 
sf_IP[9] = 0.968439; sf_IP[10] = 1.03185; 
// float sf_stat_IP[11];
sf_stat_IP[0] = 0.0354127; sf_stat_IP[1] = 0.0268463; sf_stat_IP[2] = 0.0226531; 
sf_stat_IP[3] = 0.0204136; sf_stat_IP[4] = 0.0233218; sf_stat_IP[5] = 0.0271887; 
sf_stat_IP[6] = 0.0152979; sf_stat_IP[7] = 0.0175077; sf_stat_IP[8] = 0.0155133; 
sf_stat_IP[9] = 0.016824; sf_stat_IP[10] = 0.0297241; 
// float away_syst_IP[11];
away_syst_IP[0] = 0.0229973; away_syst_IP[1] = 0.0242812; away_syst_IP[2] = 0.025121; 
away_syst_IP[3] = 0.0143937; away_syst_IP[4] = 0.0152592; away_syst_IP[5] = 0.0152649; 
away_syst_IP[6] = 0.00519076; away_syst_IP[7] = 0.00536271; away_syst_IP[8] = 0.0207214; 
away_syst_IP[9] = 0.0214561; away_syst_IP[10] = 0.0274859; 
// float mupt_syst_IP[11];
mupt_syst_IP[0] = 0.119146; mupt_syst_IP[1] = 0.0441533; mupt_syst_IP[2] = 0.0456804; 
mupt_syst_IP[3] = 0.00524259; mupt_syst_IP[4] = 0.00555783; mupt_syst_IP[5] = 0.00555993; 
mupt_syst_IP[6] = 0.0102529; mupt_syst_IP[7] = 0.0105926; mupt_syst_IP[8] = 0.0063483; 
mupt_syst_IP[9] = 0.00657339; mupt_syst_IP[10] = 0.0284092; 
// float gluon_syst_IP[11];
gluon_syst_IP[0] = 0.000959338; gluon_syst_IP[1] = 0.00758048; gluon_syst_IP[2] = 0.00784265; 
gluon_syst_IP[3] = 0.00346811; gluon_syst_IP[4] = 0.00367665; gluon_syst_IP[5] = 0.00367804; 
gluon_syst_IP[6] = 0.00573653; gluon_syst_IP[7] = 0.00592656; gluon_syst_IP[8] = 0.021614; 
gluon_syst_IP[9] = 0.0223803; gluon_syst_IP[10] = 0.0823072; 
// float pu_syst_IP[11];
pu_syst_IP[0] = 0.0191957; pu_syst_IP[1] = 0.00898904; pu_syst_IP[2] = 0.00929993; 
pu_syst_IP[3] = 0.00639865; pu_syst_IP[4] = 0.0067834; pu_syst_IP[5] = 0.00678596; 
pu_syst_IP[6] = 0.00695866; pu_syst_IP[7] = 0.00718917; pu_syst_IP[8] = 0.00215987; 
pu_syst_IP[9] = 0.00223645; pu_syst_IP[10] = 0.00338173; 
// float total_syst_IP[11];
total_syst_IP[0] = 0.122858; total_syst_IP[1] = 0.0517432; total_syst_IP[2] = 0.0535327; 
total_syst_IP[3] = 0.0169597; total_syst_IP[4] = 0.0179795; total_syst_IP[5] = 0.0179863; 
total_syst_IP[6] = 0.0146081; total_syst_IP[7] = 0.015092; total_syst_IP[8] = 0.030684; 
total_syst_IP[9] = 0.0317719; total_syst_IP[10] = 0.0913699; 

// TCHPM : JP
eff_JP[0] = 0.227018; eff_stat_JP[0] = 0.00947889; effMC_JP[0] = 0.265318; effMC_stat_JP[0] = 0.00102278; 
sf_JP[0] = 0.855643; sf_stat_JP[0] = 0.0357265; sf_syst_JP[0] = 0.169417; sf_eror_JP[0] = 0.173143; 
mupt_syst_JP[0] = 0.014; gluon_syst_JP[0] = 0.007; pu_syst_JP[0] = 0.004; total_syst_JP[0] = 0.198; 
cor_syst_JP[0] = 0.044; inc_syst_JP[0] = 0.192; bias_syst_JP[0] = 0; 
eff_JP[1] = 0.365419; eff_stat_JP[1] = 0.00789366; effMC_JP[1] = 0.393618; effMC_stat_JP[1] = 0.000855906; 
sf_JP[1] = 0.928358; sf_stat_JP[1] = 0.0200541; sf_syst_JP[1] = 0.0742686; sf_eror_JP[1] = 0.0769285; 
mupt_syst_JP[1] = 0.014; gluon_syst_JP[1] = 0.004; pu_syst_JP[1] = 0.001; total_syst_JP[1] = 0.08; 
cor_syst_JP[1] = 0.022; inc_syst_JP[1] = 0.076; bias_syst_JP[1] = 0; 
eff_JP[2] = 0.445764; eff_stat_JP[2] = 0.00755917; effMC_JP[2] = 0.478127; effMC_stat_JP[2] = 0.000941241; 
sf_JP[2] = 0.932313; sf_stat_JP[2] = 0.0158099; sf_syst_JP[2] = 0.0736527; sf_eror_JP[2] = 0.0753305; 
mupt_syst_JP[2] = 0.014; gluon_syst_JP[2] = 0.007; pu_syst_JP[2] = 0.001; total_syst_JP[2] = 0.079; 
cor_syst_JP[2] = 0.013; inc_syst_JP[2] = 0.076; bias_syst_JP[2] = 0; 
eff_JP[3] = 0.515454; eff_stat_JP[3] = 0.00830062; effMC_JP[3] = 0.544021; effMC_stat_JP[3] = 0.00129892; 
sf_JP[3] = 0.94749; sf_stat_JP[3] = 0.0152579; sf_syst_JP[3] = 0.0454795; sf_eror_JP[3] = 0.0479707; 
mupt_syst_JP[3] = 0.014; gluon_syst_JP[3] = 0.005; pu_syst_JP[3] = 0; total_syst_JP[3] = 0.048; 
cor_syst_JP[3] = 0.008; inc_syst_JP[3] = 0.045; bias_syst_JP[3] = 0; 
eff_JP[4] = 0.55456; eff_stat_JP[4] = 0.00892759; effMC_JP[4] = 0.581867; effMC_stat_JP[4] = 0.00128783; 
sf_JP[4] = 0.95307; sf_stat_JP[4] = 0.015343; sf_syst_JP[4] = 0.0457474; sf_eror_JP[4] = 0.0482517; 
mupt_syst_JP[4] = 0.014; gluon_syst_JP[4] = 0.006; pu_syst_JP[4] = 0.001; total_syst_JP[4] = 0.048; 
cor_syst_JP[4] = 0.006; inc_syst_JP[4] = 0.045; bias_syst_JP[4] = 0; 
eff_JP[5] = 0.576402; eff_stat_JP[5] = 0.00912072; effMC_JP[5] = 0.604206; effMC_stat_JP[5] = 0.00133047; 
sf_JP[5] = 0.953983; sf_stat_JP[5] = 0.0150954; sf_syst_JP[5] = 0.0457912; sf_eror_JP[5] = 0.0482152; 
mupt_syst_JP[5] = 0.014; gluon_syst_JP[5] = 0.004; pu_syst_JP[5] = 0; total_syst_JP[5] = 0.048; 
cor_syst_JP[5] = 0.005; inc_syst_JP[5] = 0.045; bias_syst_JP[5] = 0; 
eff_JP[6] = 0.619944; eff_stat_JP[6] = 0.00783631; effMC_JP[6] = 0.640752; effMC_stat_JP[6] = 0.00124987; 
sf_JP[6] = 0.967526; sf_stat_JP[6] = 0.0122299; sf_syst_JP[6] = 0.020318; sf_eror_JP[6] = 0.0237148; 
mupt_syst_JP[6] = 0; gluon_syst_JP[6] = 0.003; pu_syst_JP[6] = 0.001; total_syst_JP[6] = 0.021; 
cor_syst_JP[6] = 0.004; inc_syst_JP[6] = 0.02; bias_syst_JP[6] = 0; 
eff_JP[7] = 0.638988; eff_stat_JP[7] = 0.00928929; effMC_JP[7] = 0.661531; effMC_stat_JP[7] = 0.00123688; 
sf_JP[7] = 0.965922; sf_stat_JP[7] = 0.0140421; sf_syst_JP[7] = 0.0193184; sf_eror_JP[7] = 0.0238827; 
mupt_syst_JP[7] = 0; gluon_syst_JP[7] = 0.001; pu_syst_JP[7] = 0.001; total_syst_JP[7] = 0.02; 
cor_syst_JP[7] = 0.003; inc_syst_JP[7] = 0.02; bias_syst_JP[7] = 0; 
eff_JP[8] = 0.667675; eff_stat_JP[8] = 0.00998369; effMC_JP[8] = 0.68612; effMC_stat_JP[8] = 0.00116784; 
sf_JP[8] = 0.973118; sf_stat_JP[8] = 0.0145509; sf_syst_JP[8] = 0.0194624; sf_eror_JP[8] = 0.0243005; 
mupt_syst_JP[8] = 0; gluon_syst_JP[8] = 0.002; pu_syst_JP[8] = 0; total_syst_JP[8] = 0.02; 
cor_syst_JP[8] = 0.003; inc_syst_JP[8] = 0.02; bias_syst_JP[8] = 0; 
eff_JP[9] = 0.683131; eff_stat_JP[9] = 0.00971955; effMC_JP[9] = 0.70153; effMC_stat_JP[9] = 0.00140107; 
sf_JP[9] = 0.973773; sf_stat_JP[9] = 0.0138548; sf_syst_JP[9] = 0.0165541; sf_eror_JP[9] = 0.0215869; 
mupt_syst_JP[9] = 0; gluon_syst_JP[9] = 0.002; pu_syst_JP[9] = 0; total_syst_JP[9] = 0.017; 
cor_syst_JP[9] = 0.002; inc_syst_JP[9] = 0.017; bias_syst_JP[9] = 0; 
eff_JP[10] = 0.687939; eff_stat_JP[10] = 0.0134546; effMC_JP[10] = 0.710355; effMC_stat_JP[10] = 0.00247465; 
sf_JP[10] = 0.968444; sf_stat_JP[10] = 0.0189407; sf_syst_JP[10] = 0.0193689; sf_eror_JP[10] = 0.0270907; 
mupt_syst_JP[10] = 0.008; gluon_syst_JP[10] = 0.007; pu_syst_JP[10] = 0.001; total_syst_JP[10] = 0.02; 
cor_syst_JP[10] = 0.002; inc_syst_JP[10] = 0.017; bias_syst_JP[10] = 0; 
eff_JP[11] = 0.693183; eff_stat_JP[11] = 0.0176899; effMC_JP[11] = 0.720702; effMC_stat_JP[11] = 0.00411809; 
sf_JP[11] = 0.961817; sf_stat_JP[11] = 0.0245453; sf_syst_JP[11] = 0.0192363; sf_eror_JP[11] = 0.0311851; 
mupt_syst_JP[11] = 0.008; gluon_syst_JP[11] = 0.006; pu_syst_JP[11] = 0.004; total_syst_JP[11] = 0.02; 
cor_syst_JP[11] = 0.001; inc_syst_JP[11] = 0.017; bias_syst_JP[11] = 0; 
eff_JP[12] = 0.752064; eff_stat_JP[12] = 0.0229799; effMC_JP[12] = 0.754246; effMC_stat_JP[12] = 0.00671688; 
sf_JP[12] = 0.997108; sf_stat_JP[12] = 0.0304674; sf_syst_JP[12] = 0.027919; sf_eror_JP[12] = 0.0413247; 
mupt_syst_JP[12] = 0.008; gluon_syst_JP[12] = 0.007; pu_syst_JP[12] = 0.001; total_syst_JP[12] = 0.028; 
cor_syst_JP[12] = 0.001; inc_syst_JP[12] = 0.026; bias_syst_JP[12] = 0; 
eff_JP[13] = 0.754782; eff_stat_JP[13] = 0.0250309; effMC_JP[13] = 0.737526; effMC_stat_JP[13] = 0.00602354; 
sf_JP[13] = 1.0234; sf_stat_JP[13] = 0.033939; sf_syst_JP[13] = 0.0296786; sf_eror_JP[13] = 0.0450852; 
mupt_syst_JP[13] = 0.008; gluon_syst_JP[13] = 0.007; pu_syst_JP[13] = 0.006; total_syst_JP[13] = 0.029; 
cor_syst_JP[13] = 0.001; inc_syst_JP[13] = 0.026; bias_syst_JP[13] = 0; 
eff_JP[14] = 0.714846; eff_stat_JP[14] = 0.047877; effMC_JP[14] = 0.740927; effMC_stat_JP[14] = 0.00506698; 
sf_JP[14] = 0.9648; sf_stat_JP[14] = 0.0646178; sf_syst_JP[14] = 0.0299088; sf_eror_JP[14] = 0.0712039; 
mupt_syst_JP[14] = 0.008; gluon_syst_JP[14] = 0.002; pu_syst_JP[14] = 0.013; total_syst_JP[14] = 0.031; 
cor_syst_JP[14] = 0.008; inc_syst_JP[14] = 0.026; bias_syst_JP[14] = 0; 
}

// #############################################################################
if ( title == "TCHPT" ) { 
// TCHPT : PTREL
// float eff_PT[11];
eff_PT[0] = 0.117491; eff_PT[1] = 0.217876; eff_PT[2] = 0.309282; 
eff_PT[3] = 0.348859; eff_PT[4] = 0.370343; eff_PT[5] = 0.390878; 
eff_PT[6] = 0.42702; eff_PT[7] = 0.419973; eff_PT[8] = 0.429038; 
eff_PT[9] = 0.419457; eff_PT[10] = 0.478172; 
// float eff_stat_PT[11];
eff_stat_PT[0] = 0.00567006; eff_stat_PT[1] = 0.00560229; eff_stat_PT[2] = 0.00727244; 
eff_stat_PT[3] = 0.00786886; eff_stat_PT[4] = 0.0090658; eff_stat_PT[5] = 0.0107545; 
eff_stat_PT[6] = 0.00806656; eff_stat_PT[7] = 0.00954718; eff_stat_PT[8] = 0.0112236; 
eff_stat_PT[9] = 0.0102366; eff_stat_PT[10] = 0.0231888; 
// float effMC_PT[11];
effMC_PT[0] = 0.143114; effMC_PT[1] = 0.247508; effMC_PT[2] = 0.323431; 
effMC_PT[3] = 0.39104; effMC_PT[4] = 0.419548; effMC_PT[5] = 0.451149; 
effMC_PT[6] = 0.488094; effMC_PT[7] = 0.484964; effMC_PT[8] = 0.521938; 
effMC_PT[9] = 0.52245; effMC_PT[10] = 0.547825; 
// float effMC_stat_PT[11];
effMC_stat_PT[0] = 0.00226117; effMC_stat_PT[1] = 0.0018035; effMC_stat_PT[2] = 0.00179336; 
effMC_stat_PT[3] = 0.00262516; effMC_stat_PT[4] = 0.00263001; effMC_stat_PT[5] = 0.00276787; 
effMC_stat_PT[6] = 0.00303277; effMC_stat_PT[7] = 0.00292358; effMC_stat_PT[8] = 0.00388447; 
effMC_stat_PT[9] = 0.00410593; effMC_stat_PT[10] = 0.00605899; 
// float sf_PT[11];
sf_PT[0] = 0.820962; sf_PT[1] = 0.880281; sf_PT[2] = 0.956254; 
sf_PT[3] = 0.892131; sf_PT[4] = 0.88272; sf_PT[5] = 0.866405; 
sf_PT[6] = 0.874872; sf_PT[7] = 0.865988; sf_PT[8] = 0.822011; 
sf_PT[9] = 0.802866; sf_PT[10] = 0.872856; 
// float sf_stat_PT[11];
sf_stat_PT[0] = 0.0416884; sf_stat_PT[1] = 0.0235261; sf_stat_PT[2] = 0.023102; 
sf_stat_PT[3] = 0.0209953; sf_stat_PT[4] = 0.0223057; sf_stat_PT[5] = 0.0244236; 
sf_stat_PT[6] = 0.0173977; sf_stat_PT[7] = 0.0203668; sf_stat_PT[8] = 0.0223571; 
sf_stat_PT[9] = 0.0205843; sf_stat_PT[10] = 0.0434159; 
// float away_syst_PT[11];
away_syst_PT[0] = 0.007447; away_syst_PT[1] = 0.0464228; away_syst_PT[2] = 0.0504293; 
away_syst_PT[3] = 0.0237832; away_syst_PT[4] = 0.0235323; away_syst_PT[5] = 0.0230974; 
away_syst_PT[6] = 0.0333633; away_syst_PT[7] = 0.0330245; away_syst_PT[8] = 0.0806159; 
away_syst_PT[9] = 0.0787383; away_syst_PT[10] = 0.068767; 
// float mupt_syst_PT[11];
mupt_syst_PT[0] = 0.057506; mupt_syst_PT[1] = 0.0178181; mupt_syst_PT[2] = 0.0193559; 
mupt_syst_PT[3] = 0.0199758; mupt_syst_PT[4] = 0.0197651; mupt_syst_PT[5] = 0.0193997; 
mupt_syst_PT[6] = 0.0318362; mupt_syst_PT[7] = 0.031513; mupt_syst_PT[8] = 0.0852835; 
mupt_syst_PT[9] = 0.0832972; mupt_syst_PT[10] = 0.099666; 
// float lc_syst_PT[11];
lc_syst_PT[0] = 0.00485402; lc_syst_PT[1] = 0.00443317; lc_syst_PT[2] = 0.00481577; 
lc_syst_PT[3] = 0.00639947; lc_syst_PT[4] = 0.00633196; lc_syst_PT[5] = 0.00621493; 
lc_syst_PT[6] = 0.00485453; lc_syst_PT[7] = 0.00480523; lc_syst_PT[8] = 0.017225; 
lc_syst_PT[9] = 0.0168238; lc_syst_PT[10] = 0.060546; 
// float gluon_syst_PT[11];
gluon_syst_PT[0] = 0.00203103; gluon_syst_PT[1] = 0.00322459; gluon_syst_PT[2] = 0.00350289; 
gluon_syst_PT[3] = 0.00383744; gluon_syst_PT[4] = 0.00379696; gluon_syst_PT[5] = 0.00372678; 
gluon_syst_PT[6] = 0.00781219; gluon_syst_PT[7] = 0.00773286; gluon_syst_PT[8] = 0.0209076; 
gluon_syst_PT[9] = 0.0204207; gluon_syst_PT[10] = 0.047471; 
// float pu_syst_PT[11];
pu_syst_PT[0] = 0.019413; pu_syst_PT[1] = 0.0160923; pu_syst_PT[2] = 0.0174811; 
pu_syst_PT[3] = 0.0172089; pu_syst_PT[4] = 0.0170274; pu_syst_PT[5] = 0.0167127; 
pu_syst_PT[6] = 0.0133225; pu_syst_PT[7] = 0.0131872; pu_syst_PT[8] = 0.0118266; 
pu_syst_PT[9] = 0.0115512; pu_syst_PT[10] = 0.032467; 
// float total_syst_PT[11];
total_syst_PT[0] = 0.0613754; total_syst_PT[1] = 0.0525506; total_syst_PT[2] = 0.057086; 
total_syst_PT[3] = 0.0362836; total_syst_PT[4] = 0.0359008; total_syst_PT[5] = 0.0352373; 
total_syst_PT[6] = 0.0488747; total_syst_PT[7] = 0.0483784; total_syst_PT[8] = 0.12102; 
total_syst_PT[9] = 0.118202; total_syst_PT[10] = 0.147091; 

// TCHPT : S8
eff_SY[0]=0.1160; eff_SY[1]=0.2047; eff_SY[2]=0.2920; eff_SY[3]=0.3531; eff_SY[4]=0.3821;eff_SY[5]=0.4462; eff_SY[6]=0.4899; eff_SY[7]=0.4403; 
eff_stat_SY[0]=0.0059; eff_stat_SY[1]=0.0054; eff_stat_SY[2]=0.0074; eff_stat_SY[3]=0.0087; eff_stat_SY[4]=0.0095;eff_stat_SY[5]=0.0171; eff_stat_SY[6]=0.0283; eff_stat_SY[7]=0.0252; 
effMC_SY[0]=0.1403; effMC_SY[1]=0.2372; effMC_SY[2]=0.3196; effMC_SY[3]=0.3855; effMC_SY[4]=0.4420; effMC_SY[5]=0.4970; effMC_SY[6]=0.5224; effMC_SY[7]=0.5474; 
effMC_stat_SY[0]=0.0004; effMC_stat_SY[1]=0.0005; effMC_stat_SY[2]=0.0007; effMC_stat_SY[3]=0.0009; effMC_stat_SY[4]=0.0008; effMC_stat_SY[5]=0.0010; effMC_stat_SY[6]=0.0012; effMC_stat_SY[7]=0.0011; 
sf_SY[0]=0.8268; sf_SY[1]=0.863; sf_SY[2]=0.9136; sf_SY[3]=0.916; sf_SY[4]=0.8645; sf_SY[5]=0.8978; sf_SY[6]=0.9378; sf_SY[7]=0.8043; 
sf_stat_SY[0]=0.0421; sf_stat_SY[1]=0.0228; sf_stat_SY[2]=0.0232; sf_stat_SY[3]=0.0227; sf_stat_SY[4]=0.0216; sf_stat_SY[5]=0.0345; sf_stat_SY[6]=0.0542; sf_stat_SY[7]=0.0461; 
away[0]=2.5304; away[1]=0.876; away[2]=1.4471; away[3]=17.2527; 
mu[0]=1.5085; mu[1]=1.2045; mu[2]=2.1706; mu[3]=0.0655; 
ptrel[0]=0.2433; ptrel[1]=1.533; ptrel[2]=2.1706; ptrel[3]=2.6425; 
gsplit[0]=0.7237; gsplit[1]=1.4454; gsplit[2]=2.7384; gsplit[3]=0.0536; 
closure[0]=0.1703; closure[1]=0.0241; closure[2]=0.1576; closure[3]=1.1981; 
pu[0]=0.1946; pu[1]=0.3559; pu[2]=1.0963; pu[3]=1.1356; 

// TCHPT : IP3D
// float eff_IP[11];
eff_IP[0] = 0.1123; eff_IP[1] = 0.214309; eff_IP[2] = 0.28853; 
eff_IP[3] = 0.32587; eff_IP[4] = 0.367936; eff_IP[5] = 0.403269; 
eff_IP[6] = 0.41848; eff_IP[7] = 0.428879; eff_IP[8] = 0.462201; 
eff_IP[9] = 0.463824; eff_IP[10] = 0.535648; 
// float eff_stat_IP[11];
eff_stat_IP[0] = 0.0072884; eff_stat_IP[1] = 0.00695647; eff_stat_IP[2] = 0.0637083; 
eff_stat_IP[3] = 0.0686587; eff_stat_IP[4] = 0.0730022; eff_stat_IP[5] = 0.0146543; 
eff_stat_IP[6] = 0.00880579; eff_stat_IP[7] = 0.0108191; eff_stat_IP[8] = 0.00957529; 
eff_stat_IP[9] = 0.0110482; eff_stat_IP[10] = 0.0212214; 
// float effMC_IP[11];
effMC_IP[0] = 0.145044; effMC_IP[1] = 0.250331; effMC_IP[2] = 0.323934; 
effMC_IP[3] = 0.392816; effMC_IP[4] = 0.419709; effMC_IP[5] = 0.450517; 
effMC_IP[6] = 0.485088; effMC_IP[7] = 0.487569; effMC_IP[8] = 0.516393; 
effMC_IP[9] = 0.519442; effMC_IP[10] = 0.541952; 
// float effMC_stat_IP[11];
effMC_stat_IP[0] = 0.00224388; effMC_stat_IP[1] = 0.00192261; effMC_stat_IP[2] = 0.00187784; 
effMC_stat_IP[3] = 0.00260274; effMC_stat_IP[4] = 0.00258402; effMC_stat_IP[5] = 0.00271225; 
effMC_stat_IP[6] = 0.00295242; effMC_stat_IP[7] = 0.00283801; effMC_stat_IP[8] = 0.00376142; 
effMC_stat_IP[9] = 0.00396156; effMC_stat_IP[10] = 0.00605773; 
// float sf_IP[11];
sf_IP[0] = 0.774249; sf_IP[1] = 0.856103; sf_IP[2] = 0.890708; 
sf_IP[3] = 0.829574; sf_IP[4] = 0.876645; sf_IP[5] = 0.895124; 
sf_IP[6] = 0.862688; sf_IP[7] = 0.879627; sf_IP[8] = 0.895058; 
sf_IP[9] = 0.892928; sf_IP[10] = 0.988369; 
// float sf_stat_IP[11];
sf_stat_IP[0] = 0.0516575; sf_stat_IP[1] = 0.0285564; sf_stat_IP[2] = 0.196739; 
sf_stat_IP[3] = 0.174872; sf_stat_IP[4] = 0.174019; sf_stat_IP[5] = 0.0329712; 
sf_stat_IP[6] = 0.0188971; sf_stat_IP[7] = 0.022773; sf_stat_IP[8] = 0.0196554; 
sf_stat_IP[9] = 0.022333; sf_stat_IP[10] = 0.0406859; 
// float away_syst_IP[11];
away_syst_IP[0] = 0.0380992; away_syst_IP[1] = 0.000291936; away_syst_IP[2] = 0.000303736; 
away_syst_IP[3] = 0.0166121; away_syst_IP[4] = 0.0175547; away_syst_IP[5] = 0.0179248; 
away_syst_IP[6] = 0.0197747; away_syst_IP[7] = 0.020163; away_syst_IP[8] = 0.013535; 
away_syst_IP[9] = 0.0135028; away_syst_IP[10] = 0.0493564; 
// float mupt_syst_IP[11];
mupt_syst_IP[0] = 0.00920665; mupt_syst_IP[1] = 0.0273128; mupt_syst_IP[2] = 0.0284168; 
mupt_syst_IP[3] = 0.00263853; mupt_syst_IP[4] = 0.00278824; mupt_syst_IP[5] = 0.00284702; 
mupt_syst_IP[6] = 0.0200954; mupt_syst_IP[7] = 0.02049; mupt_syst_IP[8] = 0.00778511; 
mupt_syst_IP[9] = 0.00776659; mupt_syst_IP[10] = 0.0169748; 
// float gluon_syst_IP[11];
gluon_syst_IP[0] = 0; gluon_syst_IP[1] = 0.00609836; gluon_syst_IP[2] = 0.00634487; 
gluon_syst_IP[3] = 0.00754543; gluon_syst_IP[4] = 0.00797357; gluon_syst_IP[5] = 0.00814164; 
gluon_syst_IP[6] = 0.0119507; gluon_syst_IP[7] = 0.0121853; gluon_syst_IP[8] = 0.0195174; 
gluon_syst_IP[9] = 0.019471; gluon_syst_IP[10] = 0.13394; 
// float pu_syst_IP[11];
pu_syst_IP[0] = 0.012987; pu_syst_IP[1] = 0.00805748; pu_syst_IP[2] = 0.00838318; 
pu_syst_IP[3] = 0.00242644; pu_syst_IP[4] = 0.00256412; pu_syst_IP[5] = 0.00261817; 
pu_syst_IP[6] = 0.00305587; pu_syst_IP[7] = 0.00311587; pu_syst_IP[8] = 0.00781822; 
pu_syst_IP[9] = 0.00779961; pu_syst_IP[10] = 0.02398; 
// float total_syst_IP[11];
total_syst_IP[0] = 0.0412913; total_syst_IP[1] = 0.0291237; total_syst_IP[2] = 0.0303009; 
total_syst_IP[3] = 0.0185942; total_syst_IP[4] = 0.0196493; total_syst_IP[5] = 0.0200635; 
total_syst_IP[6] = 0.0307737; total_syst_IP[7] = 0.031378; total_syst_IP[8] = 0.0261889; 
total_syst_IP[9] = 0.0261266; total_syst_IP[10] = 0.145737; 

// TCHPT : JP
eff_JP[0] = 0.157606; eff_stat_JP[0] = 0.00773439; effMC_JP[0] = 0.191742; effMC_stat_JP[0] = 0.000911971; 
sf_JP[0] = 0.821969; sf_stat_JP[0] = 0.0403374; sf_syst_JP[0] = 0.22111; sf_eror_JP[0] = 0.224759; 
mupt_syst_JP[0] = 0.033; gluon_syst_JP[0] = 0.024; pu_syst_JP[0] = 0.009; total_syst_JP[0] = 0.269; 
cor_syst_JP[0] = 0.044; inc_syst_JP[0] = 0.262; bias_syst_JP[0] = 0; 
eff_JP[1] = 0.247634; eff_stat_JP[1] = 0.012424; effMC_JP[1] = 0.294012; effMC_stat_JP[1] = 0.000798172; 
sf_JP[1] = 0.842256; sf_stat_JP[1] = 0.0422567; sf_syst_JP[1] = 0.207195; sf_eror_JP[1] = 0.21146; 
mupt_syst_JP[1] = 0.033; gluon_syst_JP[1] = 0.011; pu_syst_JP[1] = 0.011; total_syst_JP[1] = 0.246; 
cor_syst_JP[1] = 0.022; inc_syst_JP[1] = 0.242; bias_syst_JP[1] = 0; 
eff_JP[2] = 0.323346; eff_stat_JP[2] = 0.00743845; effMC_JP[2] = 0.364902; effMC_stat_JP[2] = 0.000907099; 
sf_JP[2] = 0.88612; sf_stat_JP[2] = 0.0203848; sf_syst_JP[2] = 0.217099; sf_eror_JP[2] = 0.218054; 
mupt_syst_JP[2] = 0.033; gluon_syst_JP[2] = 0.01; pu_syst_JP[2] = 0.003; total_syst_JP[2] = 0.245; 
cor_syst_JP[2] = 0.013; inc_syst_JP[2] = 0.242; bias_syst_JP[2] = 0; 
eff_JP[3] = 0.393391; eff_stat_JP[3] = 0.00798533; effMC_JP[3] = 0.424827; effMC_stat_JP[3] = 0.00128916; 
sf_JP[3] = 0.926002; sf_stat_JP[3] = 0.0187967; sf_syst_JP[3] = 0.0527821; sf_eror_JP[3] = 0.0560292; 
mupt_syst_JP[3] = 0.033; gluon_syst_JP[3] = 0.005; pu_syst_JP[3] = 0.003; total_syst_JP[3] = 0.057; 
cor_syst_JP[3] = 0.008; inc_syst_JP[3] = 0.045; bias_syst_JP[3] = 0; 
eff_JP[4] = 0.425323; eff_stat_JP[4] = 0.00784317; effMC_JP[4] = 0.460156; effMC_stat_JP[4] = 0.0013013; 
sf_JP[4] = 0.924302; sf_stat_JP[4] = 0.0170446; sf_syst_JP[4] = 0.0526852; sf_eror_JP[4] = 0.0553737; 
mupt_syst_JP[4] = 0.033; gluon_syst_JP[4] = 0.008; pu_syst_JP[4] = 0.001; total_syst_JP[4] = 0.057; 
cor_syst_JP[4] = 0.006; inc_syst_JP[4] = 0.045; bias_syst_JP[4] = 0; 
eff_JP[5] = 0.447169; eff_stat_JP[5] = 0.00816136; effMC_JP[5] = 0.480883; effMC_stat_JP[5] = 0.00135934; 
sf_JP[5] = 0.929891; sf_stat_JP[5] = 0.0169716; sf_syst_JP[5] = 0.0530038; sf_eror_JP[5] = 0.0556546; 
mupt_syst_JP[5] = 0.033; gluon_syst_JP[5] = 0.011; pu_syst_JP[5] = 0.002; total_syst_JP[5] = 0.057; 
cor_syst_JP[5] = 0.005; inc_syst_JP[5] = 0.045; bias_syst_JP[5] = 0; 
eff_JP[6] = 0.48967; eff_stat_JP[6] = 0.00697142; effMC_JP[6] = 0.516891; effMC_stat_JP[6] = 0.0013018; 
sf_JP[6] = 0.947338; sf_stat_JP[6] = 0.0134872; sf_syst_JP[6] = 0.0407355; sf_eror_JP[6] = 0.0429102; 
mupt_syst_JP[6] = 0.005; gluon_syst_JP[6] = 0.01; pu_syst_JP[6] = 0.001; total_syst_JP[6] = 0.043; 
cor_syst_JP[6] = 0.004; inc_syst_JP[6] = 0.041; bias_syst_JP[6] = 0; 
eff_JP[7] = 0.500771; eff_stat_JP[7] = 0.00809994; effMC_JP[7] = 0.533689; effMC_stat_JP[7] = 0.001304; 
sf_JP[7] = 0.938319; sf_stat_JP[7] = 0.0151773; sf_syst_JP[7] = 0.041286; sf_eror_JP[7] = 0.0439874; 
mupt_syst_JP[7] = 0.005; gluon_syst_JP[7] = 0.014; pu_syst_JP[7] = 0.001; total_syst_JP[7] = 0.044; 
cor_syst_JP[7] = 0.003; inc_syst_JP[7] = 0.041; bias_syst_JP[7] = 0; 
eff_JP[8] = 0.528852; eff_stat_JP[8] = 0.00872623; effMC_JP[8] = 0.558016; effMC_stat_JP[8] = 0.00124977; 
sf_JP[8] = 0.947735; sf_stat_JP[8] = 0.0156379; sf_syst_JP[8] = 0.0407526; sf_eror_JP[8] = 0.04365; 
mupt_syst_JP[8] = 0.005; gluon_syst_JP[8] = 0.013; pu_syst_JP[8] = 0.001; total_syst_JP[8] = 0.043; 
cor_syst_JP[8] = 0.003; inc_syst_JP[8] = 0.041; bias_syst_JP[8] = 0; 
eff_JP[9] = 0.53659; eff_stat_JP[9] = 0.00834745; effMC_JP[9] = 0.56622; effMC_stat_JP[9] = 0.00151745; 
sf_JP[9] = 0.94767; sf_stat_JP[9] = 0.0147424; sf_syst_JP[9] = 0.0350638; sf_eror_JP[9] = 0.0380369; 
mupt_syst_JP[9] = 0.005; gluon_syst_JP[9] = 0.021; pu_syst_JP[9] = 0; total_syst_JP[9] = 0.037; 
cor_syst_JP[9] = 0.002; inc_syst_JP[9] = 0.03; bias_syst_JP[9] = 0; 
eff_JP[10] = 0.531424; eff_stat_JP[10] = 0.0118709; effMC_JP[10] = 0.567096; effMC_stat_JP[10] = 0.00270313; 
sf_JP[10] = 0.937097; sf_stat_JP[10] = 0.0209328; sf_syst_JP[10] = 0.0356097; sf_eror_JP[10] = 0.0413066; 
mupt_syst_JP[10] = 0.001; gluon_syst_JP[10] = 0.024; pu_syst_JP[10] = 0.001; total_syst_JP[10] = 0.038; 
cor_syst_JP[10] = 0.002; inc_syst_JP[10] = 0.03; bias_syst_JP[10] = 0; 
eff_JP[11] = 0.528413; eff_stat_JP[11] = 0.0171281; effMC_JP[11] = 0.566349; effMC_stat_JP[11] = 0.0045488; 
sf_JP[11] = 0.933017; sf_stat_JP[11] = 0.030243; sf_syst_JP[11] = 0.0363877; sf_eror_JP[11] = 0.0473149; 
mupt_syst_JP[11] = 0.001; gluon_syst_JP[11] = 0.025; pu_syst_JP[11] = 0.005; total_syst_JP[11] = 0.039; 
cor_syst_JP[11] = 0.001; inc_syst_JP[11] = 0.03; bias_syst_JP[11] = 0; 
eff_JP[12] = 0.581113; eff_stat_JP[12] = 0.0255928; effMC_JP[12] = 0.595302; effMC_stat_JP[12] = 0.00765764; 
sf_JP[12] = 0.976165; sf_stat_JP[12] = 0.0429912; sf_syst_JP[12] = 0.0409989; sf_eror_JP[12] = 0.0594067; 
mupt_syst_JP[12] = 0.001; gluon_syst_JP[12] = 0.009; pu_syst_JP[12] = 0.007; total_syst_JP[12] = 0.042; 
cor_syst_JP[12] = 0.001; inc_syst_JP[12] = 0.04; bias_syst_JP[12] = 0; 
eff_JP[13] = 0.559159; eff_stat_JP[13] = 0.0286416; effMC_JP[13] = 0.568927; effMC_stat_JP[13] = 0.0067799; 
sf_JP[13] = 0.98283; sf_stat_JP[13] = 0.0503431; sf_syst_JP[13] = 0.040296; sf_eror_JP[13] = 0.0644841; 
mupt_syst_JP[13] = 0.001; gluon_syst_JP[13] = 0.006; pu_syst_JP[13] = 0.007; total_syst_JP[13] = 0.041; 
cor_syst_JP[13] = 0.001; inc_syst_JP[13] = 0.04; bias_syst_JP[13] = 0; 
eff_JP[14] = 0.52804; eff_stat_JP[14] = 0.0495154; effMC_JP[14] = 0.573268; effMC_stat_JP[14] = 0.00572014; 
sf_JP[14] = 0.921105; sf_stat_JP[14] = 0.086374; sf_syst_JP[14] = 0.0386864; sf_eror_JP[14] = 0.094642; 
mupt_syst_JP[14] = 0.001; gluon_syst_JP[14] = 0.004; pu_syst_JP[14] = 0.007; total_syst_JP[14] = 0.042; 
cor_syst_JP[14] = 0.008; inc_syst_JP[14] = 0.04; bias_syst_JP[14] = 0; 
}

// #############################################################################
if ( title == "SSVHEM" ) { 
// SSVHEM : PTREL
// float eff_PT[11];
eff_PT[0] = 0.299331; eff_PT[1] = 0.447391; eff_PT[2] = 0.547259; 
eff_PT[3] = 0.586314; eff_PT[4] = 0.617608; eff_PT[5] = 0.631873; 
eff_PT[6] = 0.648677; eff_PT[7] = 0.654867; eff_PT[8] = 0.610588; 
eff_PT[9] = 0.569635; eff_PT[10] = 0.475878; 
// float eff_stat_PT[11];
eff_stat_PT[0] = 0.00914445; eff_stat_PT[1] = 0.00742553; eff_stat_PT[2] = 0.00853115; 
eff_stat_PT[3] = 0.00876623; eff_stat_PT[4] = 0.0106228; eff_stat_PT[5] = 0.0139637; 
eff_stat_PT[6] = 0.0110392; eff_stat_PT[7] = 0.0115739; eff_stat_PT[8] = 0.014732; 
eff_stat_PT[9] = 0.0125974; eff_stat_PT[10] = 0.0224132; 
// float effMC_PT[11];
effMC_PT[0] = 0.34822; effMC_PT[1] = 0.491137; effMC_PT[2] = 0.580258; 
effMC_PT[3] = 0.649577; effMC_PT[4] = 0.670212; effMC_PT[5] = 0.695759; 
effMC_PT[6] = 0.712357; effMC_PT[7] = 0.698603; effMC_PT[8] = 0.687535; 
effMC_PT[9] = 0.650637; effMC_PT[10] = 0.586097; 
// float effMC_stat_PT[11];
effMC_stat_PT[0] = 0.00322019; effMC_stat_PT[1] = 0.00214142; effMC_stat_PT[2] = 0.00191719; 
effMC_stat_PT[3] = 0.00258722; effMC_stat_PT[4] = 0.00249013; effMC_stat_PT[5] = 0.00256671; 
effMC_stat_PT[6] = 0.00274277; effMC_stat_PT[7] = 0.00265601; effMC_stat_PT[8] = 0.00353238; 
effMC_stat_PT[9] = 0.00386501; effMC_stat_PT[10] = 0.00603279; 
// float sf_PT[11];
sf_PT[0] = 0.859605; sf_PT[1] = 0.910929; sf_PT[2] = 0.94313; 
sf_PT[3] = 0.902609; sf_PT[4] = 0.921512; sf_PT[5] = 0.908178; 
sf_PT[6] = 0.910606; sf_PT[7] = 0.937394; sf_PT[8] = 0.888083; 
sf_PT[9] = 0.875504; sf_PT[10] = 0.811945; 
// float sf_stat_PT[11];
sf_stat_PT[0] = 0.0274374; sf_stat_PT[1] = 0.015632; sf_stat_PT[2] = 0.0150289; 
sf_stat_PT[3] = 0.0139659; sf_stat_PT[4] = 0.0162156; sf_stat_PT[5] = 0.0203474; 
sf_stat_PT[6] = 0.0158884; sf_stat_PT[7] = 0.0169461; sf_stat_PT[8] = 0.0219077; 
sf_stat_PT[9] = 0.020048; sf_stat_PT[10] = 0.039144; 
// float away_syst_PT[11];
away_syst_PT[0] = 0.0195757; away_syst_PT[1] = 0.0361345; away_syst_PT[2] = 0.0374118; 
away_syst_PT[3] = 0.0187467; away_syst_PT[4] = 0.0191393; away_syst_PT[5] = 0.0188624; 
away_syst_PT[6] = 0.00335917; away_syst_PT[7] = 0.00345799; away_syst_PT[8] = 0.0806116; 
away_syst_PT[9] = 0.0794698; away_syst_PT[10] = 0.0763538; 
// float mupt_syst_PT[11];
mupt_syst_PT[0] = 0.0485843; mupt_syst_PT[1] = 0.0050476; mupt_syst_PT[2] = 0.00522604; 
mupt_syst_PT[3] = 0.0106317; mupt_syst_PT[4] = 0.0108544; mupt_syst_PT[5] = 0.0106973; 
mupt_syst_PT[6] = 0.0177885; mupt_syst_PT[7] = 0.0183118; mupt_syst_PT[8] = 0.0682057; 
mupt_syst_PT[9] = 0.0672396; mupt_syst_PT[10] = 0.0564646; 
// float lc_syst_PT[11];
lc_syst_PT[0] = 0.00540993; lc_syst_PT[1] = 0.00482728; lc_syst_PT[2] = 0.00499792; 
lc_syst_PT[3] = 0.00601429; lc_syst_PT[4] = 0.00614025; lc_syst_PT[5] = 0.0060514; 
lc_syst_PT[6] = 0.00382939; lc_syst_PT[7] = 0.00394204; lc_syst_PT[8] = 0.0198273; 
lc_syst_PT[9] = 0.0195465; lc_syst_PT[10] = 0.0579066; 
// float gluon_syst_PT[11];
gluon_syst_PT[0] = 0.00860291; gluon_syst_PT[1] = 0.00190872; gluon_syst_PT[2] = 0.00197619; 
gluon_syst_PT[3] = 0.000868678; gluon_syst_PT[4] = 0.00088687; gluon_syst_PT[5] = 0.000874038; 
gluon_syst_PT[6] = 0.012159; gluon_syst_PT[7] = 0.0125167; gluon_syst_PT[8] = 0.0539326; 
gluon_syst_PT[9] = 0.0531687; gluon_syst_PT[10] = 0.0704088; 
// float pu_syst_PT[11];
pu_syst_PT[0] = 0.0261086; pu_syst_PT[1] = 0.0160553; pu_syst_PT[2] = 0.0166228; 
pu_syst_PT[3] = 0.00741708; pu_syst_PT[4] = 0.00757241; pu_syst_PT[5] = 0.00746284; 
pu_syst_PT[6] = 0.00647737; pu_syst_PT[7] = 0.00666792; pu_syst_PT[8] = 0.00845338; 
pu_syst_PT[9] = 0.00833364; pu_syst_PT[10] = 0.0335934; 
// float total_syst_PT[11];
total_syst_PT[0] = 0.0594018; total_syst_PT[1] = 0.0401982; total_syst_PT[2] = 0.0416192; 
total_syst_PT[3] = 0.0235884; total_syst_PT[4] = 0.0240824; total_syst_PT[5] = 0.0237339; 
total_syst_PT[6] = 0.0230689; total_syst_PT[7] = 0.0237476; total_syst_PT[8] = 0.120514; 
total_syst_PT[9] = 0.118807; total_syst_PT[10] = 0.135857; 

// SSVHEM : S8
eff_SY[0]=0.2986; eff_SY[1]=0.4286; eff_SY[2]=0.5156; eff_SY[3]=0.5798; eff_SY[4]=0.5991;eff_SY[5]=0.6379; eff_SY[6]=0.6707; eff_SY[7]=0.5568; eff_SY[9]=0.7120; 
eff_stat_SY[0]=0.0089; eff_stat_SY[1]=0.0070; eff_stat_SY[2]=0.0089; eff_stat_SY[3]=0.0101; eff_stat_SY[4]=0.0116; eff_stat_SY[5]=0.0201; eff_stat_SY[6]=0.0368; eff_stat_SY[7]=0.0303; eff_stat_SY[9]=0.1032; 
effMC_SY[0]=0.2979; effMC_SY[1]=0.4449; effMC_SY[2]=0.5517; effMC_SY[3]=0.6232; effMC_SY[4]=0.6714; effMC_SY[5]=0.7005; effMC_SY[6]=0.7049; effMC_SY[7]=0.6963; effMC_SY[9]=0.6733; 
effMC_stat_SY[0]=0.0006; effMC_stat_SY[1]=0.0006; effMC_stat_SY[2]=0.0008; effMC_stat_SY[3]=0.0009; effMC_stat_SY[4]=0.0007; effMC_stat_SY[5]=0.0009; effMC_stat_SY[6]=0.0011; effMC_stat_SY[7]=0.0010; effMC_stat_SY[9]=0.0013; 
sf_SY[0]=1.0023; sf_SY[1]=0.9634; sf_SY[2]=0.9346; sf_SY[3]=0.9304; sf_SY[4]=0.8923; sf_SY[5]=0.9106; sf_SY[6]=0.9515; sf_SY[7]=0.7997; sf_SY[9]=1.0575; 
sf_stat_SY[0]=0.0299; sf_stat_SY[1]=0.0158; sf_stat_SY[2]=0.0162; sf_stat_SY[3]=0.0163; sf_stat_SY[4]=0.0173; sf_stat_SY[5]=0.0287; sf_stat_SY[6]=0.0522; sf_stat_SY[7]=0.0435; sf_stat_SY[9]=0.1533; 
away[0]=2.1472; away[1]=0.7353; away[2]=2.6403; away[3]=14.7069; 
mu[0]=0.2654; mu[1]=0.2394; mu[2]=0.22; mu[3]=1.1327; 
ptrel[0]=0.0483; ptrel[1]=2.3598; ptrel[2]=2.4202; ptrel[3]=1.7979; 
gsplit[0]=0.7805; gsplit[1]=0.6172; gsplit[2]=0.0854; gsplit[3]=4.405; 
closure[0]=0.2296; closure[1]=0.0463; closure[2]=0.0427; closure[3]=0.6788; 
pu[0]=0.0483; pu[1]=0.1197; pu[2]=1.3201; pu[3]=1.4923; 

// SSVHEM : IP3D
// float eff_IP[11];
eff_IP[0] = 0.289183; eff_IP[1] = 0.407842; eff_IP[2] = 0.50178; 
eff_IP[3] = 0.548837; eff_IP[4] = 0.606369; eff_IP[5] = 0.638588; 
eff_IP[6] = 0.649208; eff_IP[7] = 0.639621; eff_IP[8] = 0.672128; 
eff_IP[9] = 0.639138; eff_IP[10] = 0.606494; 
// float eff_stat_IP[11];
eff_stat_IP[0] = 0.0110715; eff_stat_IP[1] = 0.0784516; eff_stat_IP[2] = 0.012801; 
eff_stat_IP[3] = 0.012486; eff_stat_IP[4] = 0.0169324; eff_stat_IP[5] = 0.0190931; 
eff_stat_IP[6] = 0.0105445; eff_stat_IP[7] = 0.0115332; eff_stat_IP[8] = 0.00966974; 
eff_stat_IP[9] = 0.0107587; eff_stat_IP[10] = 0.0208013; 
// float effMC_IP[11];
effMC_IP[0] = 0.351116; effMC_IP[1] = 0.491105; effMC_IP[2] = 0.579997; 
effMC_IP[3] = 0.647574; effMC_IP[4] = 0.670725; effMC_IP[5] = 0.691532; 
effMC_IP[6] = 0.709763; effMC_IP[7] = 0.699042; effMC_IP[8] = 0.686029; 
effMC_IP[9] = 0.649153; effMC_IP[10] = 0.584128; 
// float effMC_stat_IP[11];
effMC_stat_IP[0] = 0.0032007; effMC_stat_IP[1] = 0.0022906; effMC_stat_IP[2] = 0.00201815; 
effMC_stat_IP[3] = 0.00257232; effMC_stat_IP[4] = 0.00245307; effMC_stat_IP[5] = 0.00252223; 
effMC_stat_IP[6] = 0.00267198; effMC_stat_IP[7] = 0.00257783; effMC_stat_IP[8] = 0.0034205; 
effMC_stat_IP[9] = 0.00372542; effMC_stat_IP[10] = 0.00603134; 
// float sf_IP[11];
sf_IP[0] = 0.823612; sf_IP[1] = 0.830457; sf_IP[2] = 0.865142; 
sf_IP[3] = 0.847528; sf_IP[4] = 0.904049; sf_IP[5] = 0.92344; 
sf_IP[6] = 0.914683; sf_IP[7] = 0.914996; sf_IP[8] = 0.979738; 
sf_IP[9] = 0.984573; sf_IP[10] = 1.03829; 
// float sf_stat_IP[11];
sf_stat_IP[0] = 0.0324138; sf_stat_IP[1] = 0.159792; sf_stat_IP[2] = 0.0222751; 
sf_stat_IP[3] = 0.0195729; sf_stat_IP[4] = 0.0254605; sf_stat_IP[5] = 0.0278146; 
sf_stat_IP[6] = 0.0152502; sf_stat_IP[7] = 0.01684; sf_stat_IP[8] = 0.0149177; 
sf_stat_IP[9] = 0.0175102; sf_stat_IP[10] = 0.0371896; 
// float away_syst_IP[11];
away_syst_IP[0] = 0.0296969; away_syst_IP[1] = 0.0229649; away_syst_IP[2] = 0.023924; 
away_syst_IP[3] = 0.0279377; away_syst_IP[4] = 0.0298008; away_syst_IP[5] = 0.03044; 
away_syst_IP[6] = 0.0172755; away_syst_IP[7] = 0.0172814; away_syst_IP[8] = 0.0279393; 
away_syst_IP[9] = 0.0280772; away_syst_IP[10] = 0.0208318; 
// float mupt_syst_IP[11];
mupt_syst_IP[0] = 0.196623; mupt_syst_IP[1] = 0.0461427; mupt_syst_IP[2] = 0.0480699; 
mupt_syst_IP[3] = 0.00329016; mupt_syst_IP[4] = 0.00350957; mupt_syst_IP[5] = 0.00358485; 
mupt_syst_IP[6] = 0.00532199; mupt_syst_IP[7] = 0.00532381; mupt_syst_IP[8] = 0.00318931; 
mupt_syst_IP[9] = 0.00320505; mupt_syst_IP[10] = 0.0249443; 
// float gluon_syst_IP[11];
gluon_syst_IP[0] = 0.052487; gluon_syst_IP[1] = 0.00379905; gluon_syst_IP[2] = 0.00395772; 
gluon_syst_IP[3] = 0.00146173; gluon_syst_IP[4] = 0.00155921; gluon_syst_IP[5] = 0.00159266; 
gluon_syst_IP[6] = 0.000330382; gluon_syst_IP[7] = 0.000330495; gluon_syst_IP[8] = 0.0213213; 
gluon_syst_IP[9] = 0.0214265; gluon_syst_IP[10] = 0.108211; 
// float pu_syst_IP[11];
pu_syst_IP[0] = 0.160844; pu_syst_IP[1] = 0.0217011; pu_syst_IP[2] = 0.0226075; 
pu_syst_IP[3] = 0.00198832; pu_syst_IP[4] = 0.00212092; pu_syst_IP[5] = 0.00216641; 
pu_syst_IP[6] = 0.00381831; pu_syst_IP[7] = 0.00381962; pu_syst_IP[8] = 0.00570899; 
pu_syst_IP[9] = 0.00573717; pu_syst_IP[10] = 0.0267666; 
// float total_syst_IP[11];
total_syst_IP[0] = 0.26109; total_syst_IP[1] = 0.0560527; total_syst_IP[2] = 0.0583938; 
total_syst_IP[3] = 0.0282388; total_syst_IP[4] = 0.030122; total_syst_IP[5] = 0.0307681; 
total_syst_IP[6] = 0.0184785; total_syst_IP[7] = 0.0184848; total_syst_IP[8] = 0.0357487; 
total_syst_IP[9] = 0.0359251; total_syst_IP[10] = 0.116114; 

// SSVHEM : JP
eff_JP[0] = 0.379733; eff_stat_JP[0] = 0.0079811; effMC_JP[0] = 0.397712; effMC_stat_JP[0] = 0.00113379; 
sf_JP[0] = 0.954796; sf_stat_JP[0] = 0.0200676; sf_syst_JP[0] = 0.11744; sf_eror_JP[0] = 0.119142; 
mupt_syst_JP[0] = 0; gluon_syst_JP[0] = 0.003; pu_syst_JP[0] = 0.001; total_syst_JP[0] = 0.123; 
cor_syst_JP[0] = 0.044; inc_syst_JP[0] = 0.115; bias_syst_JP[0] = 0; 
eff_JP[1] = 0.504939; eff_stat_JP[1] = 0.00683955; effMC_JP[1] = 0.525279; effMC_stat_JP[1] = 0.000874842; 
sf_JP[1] = 0.961276; sf_stat_JP[1] = 0.0130208; sf_syst_JP[1] = 0.0547927; sf_eror_JP[1] = 0.0563186; 
mupt_syst_JP[1] = 0; gluon_syst_JP[1] = 0.003; pu_syst_JP[1] = 0.001; total_syst_JP[1] = 0.057; 
cor_syst_JP[1] = 0.022; inc_syst_JP[1] = 0.053; bias_syst_JP[1] = 0; 
eff_JP[2] = 0.581634; eff_stat_JP[2] = 0.00624044; effMC_JP[2] = 0.601729; effMC_stat_JP[2] = 0.000922436; 
sf_JP[2] = 0.966606; sf_stat_JP[2] = 0.0103709; sf_syst_JP[2] = 0.0531633; sf_eror_JP[2] = 0.0541654; 
mupt_syst_JP[2] = 0; gluon_syst_JP[2] = 0.001; pu_syst_JP[2] = 0; total_syst_JP[2] = 0.055; 
cor_syst_JP[2] = 0.013; inc_syst_JP[2] = 0.053; bias_syst_JP[2] = 0; 
eff_JP[3] = 0.638933; eff_stat_JP[3] = 0.00655197; effMC_JP[3] = 0.657423; effMC_stat_JP[3] = 0.00123767; 
sf_JP[3] = 0.971874; sf_stat_JP[3] = 0.00996613; sf_syst_JP[3] = 0.0349875; sf_eror_JP[3] = 0.0363792; 
mupt_syst_JP[3] = 0; gluon_syst_JP[3] = 0.001; pu_syst_JP[3] = 0; total_syst_JP[3] = 0.036; 
cor_syst_JP[3] = 0.008; inc_syst_JP[3] = 0.035; bias_syst_JP[3] = 0; 
eff_JP[4] = 0.660592; eff_stat_JP[4] = 0.00686323; effMC_JP[4] = 0.68076; effMC_stat_JP[4] = 0.00121715; 
sf_JP[4] = 0.970374; sf_stat_JP[4] = 0.0100817; sf_syst_JP[4] = 0.0349335; sf_eror_JP[4] = 0.0363591; 
mupt_syst_JP[4] = 0; gluon_syst_JP[4] = 0.001; pu_syst_JP[4] = 0.001; total_syst_JP[4] = 0.036; 
cor_syst_JP[4] = 0.006; inc_syst_JP[4] = 0.035; bias_syst_JP[4] = 0; 
eff_JP[5] = 0.66869; eff_stat_JP[5] = 0.00688899; effMC_JP[5] = 0.694303; effMC_stat_JP[5] = 0.00125342; 
sf_JP[5] = 0.963109; sf_stat_JP[5] = 0.00992217; sf_syst_JP[5] = 0.0337088; sf_eror_JP[5] = 0.0351388; 
mupt_syst_JP[5] = 0; gluon_syst_JP[5] = 0.001; pu_syst_JP[5] = 0; total_syst_JP[5] = 0.035; 
cor_syst_JP[5] = 0.005; inc_syst_JP[5] = 0.035; bias_syst_JP[5] = 0; 
eff_JP[6] = 0.691503; eff_stat_JP[6] = 0.00545777; effMC_JP[6] = 0.709932; effMC_stat_JP[6] = 0.00118217; 
sf_JP[6] = 0.974041; sf_stat_JP[6] = 0.00768774; sf_syst_JP[6] = 0.0204549; sf_eror_JP[6] = 0.0218518; 
mupt_syst_JP[6] = 0.003; gluon_syst_JP[6] = 0.002; pu_syst_JP[6] = 0; total_syst_JP[6] = 0.021; 
cor_syst_JP[6] = 0.004; inc_syst_JP[6] = 0.02; bias_syst_JP[6] = 0; 
eff_JP[7] = 0.687905; eff_stat_JP[7] = 0.00673792; effMC_JP[7] = 0.706042; effMC_stat_JP[7] = 0.00119084; 
sf_JP[7] = 0.974312; sf_stat_JP[7] = 0.00954323; sf_syst_JP[7] = 0.0204606; sf_eror_JP[7] = 0.0225767; 
mupt_syst_JP[7] = 0.003; gluon_syst_JP[7] = 0.005; pu_syst_JP[7] = 0; total_syst_JP[7] = 0.021; 
cor_syst_JP[7] = 0.003; inc_syst_JP[7] = 0.02; bias_syst_JP[7] = 0; 
eff_JP[8] = 0.682763; eff_stat_JP[8] = 0.00740096; effMC_JP[8] = 0.695967; effMC_stat_JP[8] = 0.0011576; 
sf_JP[8] = 0.981028; sf_stat_JP[8] = 0.0106341; sf_syst_JP[8] = 0.0215826; sf_eror_JP[8] = 0.0240602; 
mupt_syst_JP[8] = 0.003; gluon_syst_JP[8] = 0.009; pu_syst_JP[8] = 0; total_syst_JP[8] = 0.022; 
cor_syst_JP[8] = 0.003; inc_syst_JP[8] = 0.02; bias_syst_JP[8] = 0; 
eff_JP[9] = 0.664132; eff_stat_JP[9] = 0.00650113; effMC_JP[9] = 0.675225; effMC_stat_JP[9] = 0.00143384; 
sf_JP[9] = 0.983572; sf_stat_JP[9] = 0.0096281; sf_syst_JP[9] = 0.0668829; sf_eror_JP[9] = 0.0675723; 
mupt_syst_JP[9] = 0.003; gluon_syst_JP[9] = 0.015; pu_syst_JP[9] = 0; total_syst_JP[9] = 0.068; 
cor_syst_JP[9] = 0.002; inc_syst_JP[9] = 0.066; bias_syst_JP[9] = 0; 
eff_JP[10] = 0.630499; eff_stat_JP[10] = 0.00859246; effMC_JP[10] = 0.646984; effMC_stat_JP[10] = 0.00260727; 
sf_JP[10] = 0.974521; sf_stat_JP[10] = 0.0132808; sf_syst_JP[10] = 0.069191; sf_eror_JP[10] = 0.0704541; 
mupt_syst_JP[10] = 0.01; gluon_syst_JP[10] = 0.023; pu_syst_JP[10] = 0; total_syst_JP[10] = 0.071; 
cor_syst_JP[10] = 0.002; inc_syst_JP[10] = 0.066; bias_syst_JP[10] = 0; 
eff_JP[11] = 0.615207; eff_stat_JP[11] = 0.0124828; effMC_JP[11] = 0.633456; effMC_stat_JP[11] = 0.00442289; 
sf_JP[11] = 0.971192; sf_stat_JP[11] = 0.0197058; sf_syst_JP[11] = 0.070897; sf_eror_JP[11] = 0.0735847; 
mupt_syst_JP[11] = 0.01; gluon_syst_JP[11] = 0.03; pu_syst_JP[11] = 0.003; total_syst_JP[11] = 0.073; 
cor_syst_JP[11] = 0.001; inc_syst_JP[11] = 0.066; bias_syst_JP[11] = 0; 
eff_JP[12] = 0.601794; eff_stat_JP[12] = 0.0174858; effMC_JP[12] = 0.611595; effMC_stat_JP[12] = 0.00760388; 
sf_JP[12] = 0.983975; sf_stat_JP[12] = 0.0285905; sf_syst_JP[12] = 0.0669103; sf_eror_JP[12] = 0.0727627; 
mupt_syst_JP[12] = 0.01; gluon_syst_JP[12] = 0.006; pu_syst_JP[12] = 0.004; total_syst_JP[12] = 0.068; 
cor_syst_JP[12] = 0.001; inc_syst_JP[12] = 0.067; bias_syst_JP[12] = 0; 
eff_JP[13] = 0.555429; eff_stat_JP[13] = 0.0194048; effMC_JP[13] = 0.587998; effMC_stat_JP[13] = 0.00673841; 
sf_JP[13] = 0.94461; sf_stat_JP[13] = 0.0330015; sf_syst_JP[13] = 0.0651781; sf_eror_JP[13] = 0.0730567; 
mupt_syst_JP[13] = 0.01; gluon_syst_JP[13] = 0.005; pu_syst_JP[13] = 0.01; total_syst_JP[13] = 0.069; 
cor_syst_JP[13] = 0.001; inc_syst_JP[13] = 0.067; bias_syst_JP[13] = 0; 
eff_JP[14] = 0.545121; eff_stat_JP[14] = 0.034643; effMC_JP[14] = 0.552303; effMC_stat_JP[14] = 0.00575084; 
sf_JP[14] = 0.986996; sf_stat_JP[14] = 0.0627245; sf_syst_JP[14] = 0.0690897; sf_eror_JP[14] = 0.0933153; 
mupt_syst_JP[14] = 0.01; gluon_syst_JP[14] = 0.005; pu_syst_JP[14] = 0.013; total_syst_JP[14] = 0.07; 
cor_syst_JP[14] = 0.008; inc_syst_JP[14] = 0.067; bias_syst_JP[14] = 0; 
}

// #############################################################################
if ( title == "SSVHPT" ) { 
// SSVHPT : PTREL
// float eff_PT[11];
eff_PT[0] = 0.102624; eff_PT[1] = 0.213138; eff_PT[2] = 0.308854; 
eff_PT[3] = 0.36449; eff_PT[4] = 0.39078; eff_PT[5] = 0.418914; 
eff_PT[6] = 0.437801; eff_PT[7] = 0.435632; eff_PT[8] = 0.408376; 
eff_PT[9] = 0.356415; eff_PT[10] = 0.301485; 
// float eff_stat_PT[11];
eff_stat_PT[0] = 0.00532656; eff_stat_PT[1] = 0.0056841; eff_stat_PT[2] = 0.00733428; 
eff_stat_PT[3] = 0.00790654; eff_stat_PT[4] = 0.00951848; eff_stat_PT[5] = 0.010992; 
eff_stat_PT[6] = 0.00797643; eff_stat_PT[7] = 0.0094494; eff_stat_PT[8] = 0.0109527; 
eff_stat_PT[9] = 0.00908193; eff_stat_PT[10] = 0.0177026; 
// float effMC_PT[11];
effMC_PT[0] = 0.129133; effMC_PT[1] = 0.245807; effMC_PT[2] = 0.340827; 
effMC_PT[3] = 0.415292; effMC_PT[4] = 0.448122; effMC_PT[5] = 0.483838; 
effMC_PT[6] = 0.505311; effMC_PT[7] = 0.495813; effMC_PT[8] = 0.495881; 
effMC_PT[9] = 0.457647; effMC_PT[10] = 0.397521; 
// float effMC_stat_PT[11];
effMC_stat_PT[0] = 0.00212199; effMC_stat_PT[1] = 0.00178595; effMC_stat_PT[2] = 0.00181921; 
effMC_stat_PT[3] = 0.00265782; effMC_stat_PT[4] = 0.0026562; effMC_stat_PT[5] = 0.00278807; 
effMC_stat_PT[6] = 0.00304071; effMC_stat_PT[7] = 0.00292558; effMC_stat_PT[8] = 0.00388578; 
effMC_stat_PT[9] = 0.00409412; effMC_stat_PT[10] = 0.00573646; 
// float sf_PT[11];
sf_PT[0] = 0.794712; sf_PT[1] = 0.867097; sf_PT[2] = 0.906191; 
sf_PT[3] = 0.877672; sf_PT[4] = 0.872038; sf_PT[5] = 0.865815; 
sf_PT[6] = 0.866398; sf_PT[7] = 0.878623; sf_PT[8] = 0.823537; 
sf_PT[9] = 0.778799; sf_PT[10] = 0.758413; 
// float sf_stat_PT[11];
sf_stat_PT[0] = 0.0432664; sf_stat_PT[1] = 0.0239671; sf_stat_PT[2] = 0.022056; 
sf_stat_PT[3] = 0.0198498; sf_stat_PT[4] = 0.0218607; sf_stat_PT[5] = 0.0232597; 
sf_stat_PT[6] = 0.0166239; sf_stat_PT[7] = 0.019751; sf_stat_PT[8] = 0.0230109; 
sf_stat_PT[9] = 0.0210323; sf_stat_PT[10] = 0.0458575; 
// float away_syst_PT[11];
away_syst_PT[0] = 0.023229; away_syst_PT[1] = 0.033902; away_syst_PT[2] = 0.0354305; 
away_syst_PT[3] = 0.0330297; away_syst_PT[4] = 0.0328177; away_syst_PT[5] = 0.0325835; 
away_syst_PT[6] = 0.0017637; away_syst_PT[7] = 0.00178859; away_syst_PT[8] = 0.0868279; 
away_syst_PT[9] = 0.0821111; away_syst_PT[10] = 0.129279; 
// float mupt_syst_PT[11];
mupt_syst_PT[0] = 0.037992; mupt_syst_PT[1] = 0.00657429; mupt_syst_PT[2] = 0.0068707; 
mupt_syst_PT[3] = 0.0190099; mupt_syst_PT[4] = 0.0188879; mupt_syst_PT[5] = 0.0187531; 
mupt_syst_PT[6] = 0.0202688; mupt_syst_PT[7] = 0.0205548; mupt_syst_PT[8] = 0.0817739; 
mupt_syst_PT[9] = 0.0773316; mupt_syst_PT[10] = 0.094705; 
// float lc_syst_PT[11];
lc_syst_PT[0] = 0.004807; lc_syst_PT[1] = 0.00425123; lc_syst_PT[2] = 0.0044429; 
lc_syst_PT[3] = 0.00599044; lc_syst_PT[4] = 0.00595198; lc_syst_PT[5] = 0.00590951; 
lc_syst_PT[6] = 0.00422187; lc_syst_PT[7] = 0.00428144; lc_syst_PT[8] = 0.0185321; 
lc_syst_PT[9] = 0.0175254; lc_syst_PT[10] = 0.061316; 
// float gluon_syst_PT[11];
gluon_syst_PT[0] = 0.010237; gluon_syst_PT[1] = 0.00184402; gluon_syst_PT[2] = 0.00192716; 
gluon_syst_PT[3] = 0.000185646; gluon_syst_PT[4] = 0.000184454; gluon_syst_PT[5] = 0.000183138; 
gluon_syst_PT[6] = 0.0177311; gluon_syst_PT[7] = 0.0179813; gluon_syst_PT[8] = 0.0563654; 
gluon_syst_PT[9] = 0.0533034; gluon_syst_PT[10] = 0.025405; 
// float pu_syst_PT[11];
pu_syst_PT[0] = 0.038406; pu_syst_PT[1] = 0.0197884; pu_syst_PT[2] = 0.0206806; 
pu_syst_PT[3] = 0.00806448; pu_syst_PT[4] = 0.00801271; pu_syst_PT[5] = 0.00795553; 
pu_syst_PT[6] = 0.0124233; pu_syst_PT[7] = 0.0125986; pu_syst_PT[8] = 0.0111525; 
pu_syst_PT[9] = 0.0105467; pu_syst_PT[10] = 0.090738; 
// float total_syst_PT[11];
total_syst_PT[0] = 0.0598824; total_syst_PT[1] = 0.0400702; total_syst_PT[2] = 0.0418768; 
total_syst_PT[3] = 0.0394118; total_syst_PT[4] = 0.0391588; total_syst_PT[5] = 0.0388794; 
total_syst_PT[6] = 0.0300082; total_syst_PT[7] = 0.0304316; total_syst_PT[8] = 0.133682; 
total_syst_PT[9] = 0.12642; total_syst_PT[10] = 0.195756; 

// SSVHPT : S8
eff_SY[0]=0.0910; eff_SY[1]=0.1967; eff_SY[2]=0.2900; eff_SY[3]=0.3547; eff_SY[4]=0.3763; eff_SY[5]=0.4066; eff_SY[6]=0.4490; eff_SY[7]=0.3552; 
eff_stat_SY[0]=0.0055; eff_stat_SY[1]=0.0055; eff_stat_SY[2]=0.0076; eff_stat_SY[3]=0.0089; eff_stat_SY[4]=0.0096; eff_stat_SY[5]=0.0162; eff_stat_SY[6]=0.0281; eff_stat_SY[7]=0.0227; 
effMC_SY[0]=0.1080; effMC_SY[1]=0.2175; effMC_SY[2]=0.3198; effMC_SY[3]=0.3952; effMC_SY[4]=0.4541; effMC_SY[5]=0.4970; effMC_SY[6]=0.5075; effMC_SY[7]=0.5054; 
effMC_stat_SY[0]=0.0004; effMC_stat_SY[1]=0.0005; effMC_stat_SY[2]=0.0007; effMC_stat_SY[3]=0.0009; effMC_stat_SY[4]=0.0008; effMC_stat_SY[5]=0.0010; effMC_stat_SY[6]=0.0012; effMC_stat_SY[7]=0.0011; 
sf_SY[0]=0.8426; sf_SY[1]=0.9044; sf_SY[2]=0.9068; sf_SY[3]=0.8975; sf_SY[4]=0.8287; sf_SY[5]=0.8181; sf_SY[6]=0.8847; sf_SY[7]=0.7028; 
sf_stat_SY[0]=0.051; sf_stat_SY[1]=0.0254; sf_stat_SY[2]=0.0238; sf_stat_SY[3]=0.0226; sf_stat_SY[4]=0.0212; sf_stat_SY[5]=0.0326; sf_stat_SY[6]=0.0554; sf_stat_SY[7]=0.0449; 
away[0]=2.5813; away[1]=1.157; away[2]=0.6; away[3]=13.5501; 
mu[0]=1.6004; mu[1]=1.2121; mu[2]=2.3038; mu[3]=0.6233; 
ptrel[0]=1.1358; ptrel[1]=1.5702; ptrel[2]=2.6398; ptrel[3]=3.9295; 
gsplit[0]=0.9187; gsplit[1]=0.7741; gsplit[2]=0.3189; gsplit[3]=7.0663; 
closure[0]=0.1378; closure[1]=0.0704; closure[2]=0.1594; closure[3]=1.31; 
pu[0]=0.1033; pu[1]=0.2204; pu[2]=0.576; pu[3]=0.6504; 

// SSVHPT : IP3D
// float eff_IP[11];
eff_IP[0] = 0.0930049; eff_IP[1] = 0.212336; eff_IP[2] = 0.287841; 
eff_IP[3] = 0.344334; eff_IP[4] = 0.39556; eff_IP[5] = 0.428424; 
eff_IP[6] = 0.4338; eff_IP[7] = 0.446432; eff_IP[8] = 0.455589; 
eff_IP[9] = 0.421116; eff_IP[10] = 0.362453; 
// float eff_stat_IP[11];
eff_stat_IP[0] = 0.006668; eff_stat_IP[1] = 0.0524388; eff_stat_IP[2] = 0.0644098; 
eff_stat_IP[3] = 0.0716396; eff_stat_IP[4] = 0.0763235; eff_stat_IP[5] = 0.015052; 
eff_stat_IP[6] = 0.00872847; eff_stat_IP[7] = 0.0116345; eff_stat_IP[8] = 0.00942367; 
eff_stat_IP[9] = 0.0103904; eff_stat_IP[10] = 0.0182545; 
// float effMC_IP[11];
effMC_IP[0] = 0.129119; effMC_IP[1] = 0.247048; effMC_IP[2] = 0.342004; 
effMC_IP[3] = 0.418135; effMC_IP[4] = 0.447887; effMC_IP[5] = 0.478287; 
effMC_IP[6] = 0.503914; effMC_IP[7] = 0.497366; effMC_IP[8] = 0.494494; 
effMC_IP[9] = 0.456038; effMC_IP[10] = 0.388351; 
// float effMC_stat_IP[11];
effMC_stat_IP[0] = 0.00210543; effMC_stat_IP[1] = 0.00190416; effMC_stat_IP[2] = 0.00190586; 
effMC_stat_IP[3] = 0.00263696; effMC_stat_IP[4] = 0.00261107; effMC_stat_IP[5] = 0.00273334; 
effMC_stat_IP[6] = 0.00296065; effMC_stat_IP[7] = 0.00283996; effMC_stat_IP[8] = 0.00376257; 
effMC_stat_IP[9] = 0.00395054; effMC_stat_IP[10] = 0.00573518; 
// float sf_IP[11];
sf_IP[0] = 0.720303; sf_IP[1] = 0.859493; sf_IP[2] = 0.841631; 
sf_IP[3] = 0.8235; sf_IP[4] = 0.88317; sf_IP[5] = 0.895747; 
sf_IP[6] = 0.860862; sf_IP[7] = 0.897594; sf_IP[8] = 0.921324; 
sf_IP[9] = 0.923424; sf_IP[10] = 0.933312; 
// float sf_stat_IP[11];
sf_stat_IP[0] = 0.0529611; sf_stat_IP[1] = 0.212365; sf_stat_IP[2] = 0.188389; 
sf_stat_IP[3] = 0.17141; sf_stat_IP[4] = 0.170486; sf_stat_IP[5] = 0.0318843; 
sf_stat_IP[6] = 0.0180447; sf_stat_IP[7] = 0.0239471; sf_stat_IP[8] = 0.0203057; 
sf_stat_IP[9] = 0.0241476; sf_stat_IP[10] = 0.0489842; 
// float away_syst_IP[11];
away_syst_IP[0] = 0.0232093; away_syst_IP[1] = 0.0123634; away_syst_IP[2] = 0.0121064; 
away_syst_IP[3] = 0.0418689; away_syst_IP[4] = 0.0449027; away_syst_IP[5] = 0.0455421; 
away_syst_IP[6] = 0.0250467; away_syst_IP[7] = 0.0261154; away_syst_IP[8] = 0.0291712; 
away_syst_IP[9] = 0.0292377; away_syst_IP[10] = 0.132564; 
// float mupt_syst_IP[11];
mupt_syst_IP[0] = 0; mupt_syst_IP[1] = 0.00828638; mupt_syst_IP[2] = 0.00811417; 
mupt_syst_IP[3] = 0.028427; mupt_syst_IP[4] = 0.0304868; mupt_syst_IP[5] = 0.0309209; 
mupt_syst_IP[6] = 0.0116778; mupt_syst_IP[7] = 0.0121761; mupt_syst_IP[8] = 0.00815692; 
mupt_syst_IP[9] = 0.00817551; mupt_syst_IP[10] = 0.0401175; 
// float gluon_syst_IP[11];
gluon_syst_IP[0] = 0.30534; gluon_syst_IP[1] = 0; gluon_syst_IP[2] = 0; 
gluon_syst_IP[3] = 0.0249135; gluon_syst_IP[4] = 0.0267187; gluon_syst_IP[5] = 0.0270992; 
gluon_syst_IP[6] = 0.00281295; gluon_syst_IP[7] = 0.00293298; gluon_syst_IP[8] = 0.020293; 
gluon_syst_IP[9] = 0.0203392; gluon_syst_IP[10] = 0.101032; 
// float pu_syst_IP[11];
pu_syst_IP[0] = 0; pu_syst_IP[1] = 0.00781309; pu_syst_IP[2] = 0.00765072; 
pu_syst_IP[3] = 0.0112108; pu_syst_IP[4] = 0.0120231; pu_syst_IP[5] = 0.0121943; 
pu_syst_IP[6] = 0.00969307; pu_syst_IP[7] = 0.0101067; pu_syst_IP[8] = 0.00227816; 
pu_syst_IP[9] = 0.00228336; pu_syst_IP[10] = 0.0223098; 
// float total_syst_IP[11];
total_syst_IP[0] = 0.306221; total_syst_IP[1] = 0.0168096; total_syst_IP[2] = 0.0164602; 
total_syst_IP[3] = 0.0575105; total_syst_IP[4] = 0.0616777; total_syst_IP[5] = 0.062556; 
total_syst_IP[6] = 0.0294207; total_syst_IP[7] = 0.0306761; total_syst_IP[8] = 0.0365306; 
total_syst_IP[9] = 0.0366139; total_syst_IP[10] = 0.172881; 

// SSVHPT : JP
eff_JP[0] = 0.152745; eff_stat_JP[0] = 0.00517156; effMC_JP[0] = 0.168987; effMC_stat_JP[0] = 0.000868117; 
sf_JP[0] = 0.903884; sf_stat_JP[0] = 0.0306032; sf_syst_JP[0] = 0.15818; sf_eror_JP[0] = 0.161113; 
mupt_syst_JP[0] = 0.008; gluon_syst_JP[0] = 0.001; pu_syst_JP[0] = 0.003; total_syst_JP[0] = 0.175; 
cor_syst_JP[0] = 0.044; inc_syst_JP[0] = 0.169; bias_syst_JP[0] = 0; 
eff_JP[1] = 0.266863; eff_stat_JP[1] = 0.00618114; effMC_JP[1] = 0.285509; effMC_stat_JP[1] = 0.000791268; 
sf_JP[1] = 0.934691; sf_stat_JP[1] = 0.0216495; sf_syst_JP[1] = 0.0859916; sf_eror_JP[1] = 0.088675; 
mupt_syst_JP[1] = 0.008; gluon_syst_JP[1] = 0; pu_syst_JP[1] = 0.001; total_syst_JP[1] = 0.092; 
cor_syst_JP[1] = 0.022; inc_syst_JP[1] = 0.089; bias_syst_JP[1] = 0; 
eff_JP[2] = 0.344708; eff_stat_JP[2] = 0.00628496; effMC_JP[2] = 0.371306; effMC_stat_JP[2] = 0.0009104; 
sf_JP[2] = 0.928367; sf_stat_JP[2] = 0.0169267; sf_syst_JP[2] = 0.083553; sf_eror_JP[2] = 0.0852504; 
mupt_syst_JP[2] = 0.008; gluon_syst_JP[2] = 0.003; pu_syst_JP[2] = 0.001; total_syst_JP[2] = 0.09; 
cor_syst_JP[2] = 0.013; inc_syst_JP[2] = 0.089; bias_syst_JP[2] = 0; 
eff_JP[3] = 0.415525; eff_stat_JP[3] = 0.00730859; effMC_JP[3] = 0.434116; effMC_stat_JP[3] = 0.00129262; 
sf_JP[3] = 0.957176; sf_stat_JP[3] = 0.0168356; sf_syst_JP[3] = 0.0430729; sf_eror_JP[3] = 0.0462462; 
mupt_syst_JP[3] = 0.008; gluon_syst_JP[3] = 0.002; pu_syst_JP[3] = 0; total_syst_JP[3] = 0.045; 
cor_syst_JP[3] = 0.008; inc_syst_JP[3] = 0.043; bias_syst_JP[3] = 0; 
eff_JP[4] = 0.441712; eff_stat_JP[4] = 0.00784352; effMC_JP[4] = 0.466596; effMC_stat_JP[4] = 0.00130253; 
sf_JP[4] = 0.94667; sf_stat_JP[4] = 0.0168101; sf_syst_JP[4] = 0.0416535; sf_eror_JP[4] = 0.0449176; 
mupt_syst_JP[4] = 0.008; gluon_syst_JP[4] = 0.001; pu_syst_JP[4] = 0.001; total_syst_JP[4] = 0.044; 
cor_syst_JP[4] = 0.006; inc_syst_JP[4] = 0.043; bias_syst_JP[4] = 0; 
eff_JP[5] = 0.45935; eff_stat_JP[5] = 0.00798675; effMC_JP[5] = 0.484871; effMC_stat_JP[5] = 0.00135971; 
sf_JP[5] = 0.947366; sf_stat_JP[5] = 0.0164719; sf_syst_JP[5] = 0.0416841; sf_eror_JP[5] = 0.0448206; 
mupt_syst_JP[5] = 0.008; gluon_syst_JP[5] = 0.003; pu_syst_JP[5] = 0; total_syst_JP[5] = 0.044; 
cor_syst_JP[5] = 0.005; inc_syst_JP[5] = 0.043; bias_syst_JP[5] = 0; 
eff_JP[6] = 0.487694; eff_stat_JP[6] = 0.00693864; effMC_JP[6] = 0.505638; effMC_stat_JP[6] = 0.00130246; 
sf_JP[6] = 0.964511; sf_stat_JP[6] = 0.0137225; sf_syst_JP[6] = 0.0327934; sf_eror_JP[6] = 0.0355487; 
mupt_syst_JP[6] = 0.003; gluon_syst_JP[6] = 0.02; pu_syst_JP[6] = 0.001; total_syst_JP[6] = 0.034; 
cor_syst_JP[6] = 0.004; inc_syst_JP[6] = 0.027; bias_syst_JP[6] = 0; 
eff_JP[7] = 0.49138; eff_stat_JP[7] = 0.00800705; effMC_JP[7] = 0.510865; effMC_stat_JP[7] = 0.00130666; 
sf_JP[7] = 0.96186; sf_stat_JP[7] = 0.0156735; sf_syst_JP[7] = 0.0259702; sf_eror_JP[7] = 0.0303333; 
mupt_syst_JP[7] = 0.003; gluon_syst_JP[7] = 0.002; pu_syst_JP[7] = 0.001; total_syst_JP[7] = 0.027; 
cor_syst_JP[7] = 0.003; inc_syst_JP[7] = 0.027; bias_syst_JP[7] = 0; 
eff_JP[8] = 0.48933; eff_stat_JP[8] = 0.00842069; effMC_JP[8] = 0.505057; effMC_stat_JP[8] = 0.0012582; 
sf_JP[8] = 0.96886; sf_stat_JP[8] = 0.0166728; sf_syst_JP[8] = 0.0261592; sf_eror_JP[8] = 0.0310208; 
mupt_syst_JP[8] = 0.003; gluon_syst_JP[8] = 0.003; pu_syst_JP[8] = 0; total_syst_JP[8] = 0.027; 
cor_syst_JP[8] = 0.003; inc_syst_JP[8] = 0.027; bias_syst_JP[8] = 0; 
eff_JP[9] = 0.471904; eff_stat_JP[9] = 0.00782339; effMC_JP[9] = 0.485948; effMC_stat_JP[9] = 0.00153033; 
sf_JP[9] = 0.971099; sf_stat_JP[9] = 0.0160992; sf_syst_JP[9] = 0.0534104; sf_eror_JP[9] = 0.055784; 
mupt_syst_JP[9] = 0.003; gluon_syst_JP[9] = 0.016; pu_syst_JP[9] = 0; total_syst_JP[9] = 0.055; 
cor_syst_JP[9] = 0.002; inc_syst_JP[9] = 0.052; bias_syst_JP[9] = 0; 
eff_JP[10] = 0.433826; eff_stat_JP[10] = 0.0107112; effMC_JP[10] = 0.452853; effMC_stat_JP[10] = 0.00271564; 
sf_JP[10] = 0.957983; sf_stat_JP[10] = 0.0236527; sf_syst_JP[10] = 0.054605; sf_eror_JP[10] = 0.0595076; 
mupt_syst_JP[10] = 0.015; gluon_syst_JP[10] = 0.018; pu_syst_JP[10] = 0.001; total_syst_JP[10] = 0.057; 
cor_syst_JP[10] = 0.002; inc_syst_JP[10] = 0.052; bias_syst_JP[10] = 0; 
eff_JP[11] = 0.404835; eff_stat_JP[11] = 0.0147108; effMC_JP[11] = 0.432264; effMC_stat_JP[11] = 0.00454707; 
sf_JP[11] = 0.936547; sf_stat_JP[11] = 0.0340321; sf_syst_JP[11] = 0.0543197; sf_eror_JP[11] = 0.0641001; 
mupt_syst_JP[11] = 0.015; gluon_syst_JP[11] = 0.022; pu_syst_JP[11] = 0.001; total_syst_JP[11] = 0.058; 
cor_syst_JP[11] = 0.001; inc_syst_JP[11] = 0.052; bias_syst_JP[11] = 0; 
eff_JP[12] = 0.400656; eff_stat_JP[12] = 0.0255294; effMC_JP[12] = 0.408753; effMC_stat_JP[12] = 0.00766965; 
sf_JP[12] = 0.980189; sf_stat_JP[12] = 0.0624568; sf_syst_JP[12] = 0.0499896; sf_eror_JP[12] = 0.0799989; 
mupt_syst_JP[12] = 0.015; gluon_syst_JP[12] = 0.026; pu_syst_JP[12] = 0.019; total_syst_JP[12] = 0.051; 
cor_syst_JP[12] = 0.001; inc_syst_JP[12] = 0.037; bias_syst_JP[12] = 0; 
eff_JP[13] = 0.351274; eff_stat_JP[13] = 0.0254744; effMC_JP[13] = 0.382479; effMC_stat_JP[13] = 0.00665349; 
sf_JP[13] = 0.918413; sf_stat_JP[13] = 0.0666034; sf_syst_JP[13] = 0.0413286; sf_eror_JP[13] = 0.0783841; 
mupt_syst_JP[13] = 0.015; gluon_syst_JP[13] = 0.016; pu_syst_JP[13] = 0.013; total_syst_JP[13] = 0.045; 
cor_syst_JP[13] = 0.001; inc_syst_JP[13] = 0.037; bias_syst_JP[13] = 0; 
eff_JP[14] = 0.329342; eff_stat_JP[14] = 0.039901; effMC_JP[14] = 0.353174; effMC_stat_JP[14] = 0.00552762; 
sf_JP[14] = 0.932521; sf_stat_JP[14] = 0.112978; sf_syst_JP[14] = 0.042896; sf_eror_JP[14] = 0.120847; 
mupt_syst_JP[14] = 0.015; gluon_syst_JP[14] = 0.017; pu_syst_JP[14] = 0.014; total_syst_JP[14] = 0.046; 
cor_syst_JP[14] = 0.008; inc_syst_JP[14] = 0.037; bias_syst_JP[14] = 0; 
}

// #############################################################################
if ( title == "CSVL" ) { 
// CSVL : PTREL
// float eff_PT[11];
eff_PT[0] = 0.777917; eff_PT[1] = 0.862031; eff_PT[2] = 0.86328; 
eff_PT[3] = 0.860427; eff_PT[4] = 0.868143; eff_PT[5] = 0.860296; 
eff_PT[6] = 0.881937; eff_PT[7] = 0.860848; eff_PT[8] = 0.850732; 
eff_PT[9] = 0.82052; eff_PT[10] = 0.666204; 
// float eff_stat_PT[11];
eff_stat_PT[0] = 0.0106517; eff_stat_PT[1] = 0.00639494; eff_stat_PT[2] = 0.00712178; 
eff_stat_PT[3] = 0.00790591; eff_stat_PT[4] = 0.00987413; eff_stat_PT[5] = 0.0130904; 
eff_stat_PT[6] = 0.0105506; eff_stat_PT[7] = 0.0115368; eff_stat_PT[8] = 0.0158069; 
eff_stat_PT[9] = 0.0133697; eff_stat_PT[10] = 0.0243788; 
// float effMC_PT[11];
effMC_PT[0] = 0.804096; effMC_PT[1] = 0.865358; effMC_PT[2] = 0.859483; 
effMC_PT[3] = 0.871055; effMC_PT[4] = 0.870458; effMC_PT[5] = 0.871296; 
effMC_PT[6] = 0.879354; effMC_PT[7] = 0.86888; effMC_PT[8] = 0.878319; 
effMC_PT[9] = 0.865409; effMC_PT[10] = 0.831956; 
// float effMC_stat_PT[11];
effMC_stat_PT[0] = 0.0023706; effMC_stat_PT[1] = 0.00129307; effMC_stat_PT[2] = 0.00122742; 
effMC_stat_PT[3] = 0.00162192; effMC_stat_PT[4] = 0.00162303; effMC_stat_PT[5] = 0.00169486; 
effMC_stat_PT[6] = 0.00179348; effMC_stat_PT[7] = 0.00178649; effMC_stat_PT[8] = 0.00231828; 
effMC_stat_PT[9] = 0.00255049; effMC_stat_PT[10] = 0.00442785; 
// float sf_PT[11];
sf_PT[0] = 0.967444; sf_PT[1] = 0.996155; sf_PT[2] = 1.00442; 
sf_PT[3] = 0.987799; sf_PT[4] = 0.997341; sf_PT[5] = 0.987375; 
sf_PT[6] = 1.00294; sf_PT[7] = 0.990755; sf_PT[8] = 0.968591; 
sf_PT[9] = 0.948129; sf_PT[10] = 0.800768; 
// float sf_stat_PT[11];
sf_stat_PT[0] = 0.0135504; sf_stat_PT[1] = 0.00753836; sf_stat_PT[2] = 0.00840936; 
sf_stat_PT[3] = 0.00926074; sf_stat_PT[4] = 0.011495; sf_stat_PT[5] = 0.0151463; 
sf_stat_PT[6] = 0.0121713; sf_stat_PT[7] = 0.0134331; sf_stat_PT[8] = 0.0181774; 
sf_stat_PT[9] = 0.0156996; sf_stat_PT[10] = 0.0296114; 
// float away_syst_PT[11];
away_syst_PT[0] = 0.019948; away_syst_PT[1] = 0.0218029; away_syst_PT[2] = 0.0219838; 
away_syst_PT[3] = 0.0115416; away_syst_PT[4] = 0.0116531; away_syst_PT[5] = 0.0115366; 
away_syst_PT[6] = 0.00731702; away_syst_PT[7] = 0.00722812; away_syst_PT[8] = 0.0360105; 
away_syst_PT[9] = 0.0352498; away_syst_PT[10] = 0.145024; 
// float mupt_syst_PT[11];
mupt_syst_PT[0] = 0.028699; mupt_syst_PT[1] = 0.0061558; mupt_syst_PT[2] = 0.00620688; 
mupt_syst_PT[3] = 0.0117156; mupt_syst_PT[4] = 0.0118288; mupt_syst_PT[5] = 0.0117106; 
mupt_syst_PT[6] = 0.0230193; mupt_syst_PT[7] = 0.0227396; mupt_syst_PT[8] = 0.0619441; 
mupt_syst_PT[9] = 0.0606355; mupt_syst_PT[10] = 0.012438; 
// float lc_syst_PT[11];
lc_syst_PT[0] = 0.00348198; lc_syst_PT[1] = 0.00399694; lc_syst_PT[2] = 0.0040301; 
lc_syst_PT[3] = 0.00527208; lc_syst_PT[4] = 0.00532301; lc_syst_PT[5] = 0.00526982; 
lc_syst_PT[6] = 0.00322469; lc_syst_PT[7] = 0.00318551; lc_syst_PT[8] = 0.0182391; 
lc_syst_PT[9] = 0.0178538; lc_syst_PT[10] = 0.036974; 
// float gluon_syst_PT[11];
gluon_syst_PT[0] = 0.00105697; gluon_syst_PT[1] = 8.46775e-05; gluon_syst_PT[2] = 8.538e-05; 
gluon_syst_PT[3] = 0.000301122; gluon_syst_PT[4] = 0.000304031; gluon_syst_PT[5] = 0.000300993; 
gluon_syst_PT[6] = 0.0023812; gluon_syst_PT[7] = 0.00235227; gluon_syst_PT[8] = 0.0236203; 
gluon_syst_PT[9] = 0.0231214; gluon_syst_PT[10] = 0.035458; 
// float pu_syst_PT[11];
pu_syst_PT[0] = 0.00896603; pu_syst_PT[1] = 0.00637302; pu_syst_PT[2] = 0.00642589; 
pu_syst_PT[3] = 0.00325528; pu_syst_PT[4] = 0.00328672; pu_syst_PT[5] = 0.00325388; 
pu_syst_PT[6] = 0.00280652; pu_syst_PT[7] = 0.00277242; pu_syst_PT[8] = 0.00600419; 
pu_syst_PT[9] = 0.00587735; pu_syst_PT[10] = 0.012518; 
// float total_syst_PT[11];
total_syst_PT[0] = 0.0362655; total_syst_PT[1] = 0.0238717; total_syst_PT[2] = 0.0240698; 
total_syst_PT[3] = 0.0175768; total_syst_PT[4] = 0.0177466; total_syst_PT[5] = 0.0175693; 
total_syst_PT[6] = 0.0246449; total_syst_PT[7] = 0.0243455; total_syst_PT[8] = 0.077849; 
total_syst_PT[9] = 0.0762044; total_syst_PT[10] = 0.154815; 

// CSVL : S8
eff_SY[0]=0.7873; eff_SY[1]=0.8416; eff_SY[2]=0.8394; eff_SY[3]=0.8761; eff_SY[4]=0.8529; eff_SY[5]=0.9140; eff_SY[6]=0.8454; eff_SY[7]=0.8557; 
eff_stat_SY[0]=0.0089; eff_stat_SY[1]=0.0064; eff_stat_SY[2]=0.0082; eff_stat_SY[3]=0.0108; eff_stat_SY[4]=0.0129; eff_stat_SY[5]=0.0259; eff_stat_SY[6]=0.0316; eff_stat_SY[7]=0.0495; 
effMC_SY[0]=0.7722; effMC_SY[1]=0.8381; effMC_SY[2]=0.8465; effMC_SY[3]=0.8618; effMC_SY[4]=0.8698; effMC_SY[5]=0.8790; effMC_SY[6]=0.8745; effMC_SY[7]=0.8778; 
effMC_stat_SY[0]=0.0006; effMC_stat_SY[1]=0.0005; effMC_stat_SY[2]=0.0006; effMC_stat_SY[3]=0.0007; effMC_stat_SY[4]=0.0005; effMC_stat_SY[5]=0.0006; effMC_stat_SY[6]=0.0008; effMC_stat_SY[7]=0.0007; 
sf_SY[0]=1.0196; sf_SY[1]=1.0042; sf_SY[2]=0.9916; sf_SY[3]=1.0166; sf_SY[4]=0.9806; sf_SY[5]=1.0398; sf_SY[6]=0.9667; sf_SY[7]=0.9748; 
sf_stat_SY[0]=0.0116; sf_stat_SY[1]=0.0077; sf_stat_SY[2]=0.0097; sf_stat_SY[3]=0.0126; sf_stat_SY[4]=0.0148; sf_stat_SY[5]=0.0295; sf_stat_SY[6]=0.0361; sf_stat_SY[7]=0.0564; 
away[0]=0.8152; away[1]=1.6773; away[2]=0.618; away[3]=11.5292; 
mu[0]=0.1582; mu[1]=0.6989; mu[2]=0.9041; mu[3]=2.9574; 
ptrel[0]=0.4623; ptrel[1]=0.8154; ptrel[2]=4.097; ptrel[3]=6.3477; 
gsplit[0]=0.0487; gsplit[1]=0.0116; gsplit[2]=0.3766; gsplit[3]=2.8611; 
closure[0]=0.1339; closure[1]=0.0925; closure[2]=0.1255; closure[3]=0.7928; 
pu[0]=0.0365; pu[1]=0.1048; pu[2]=1.6136; pu[3]=1.7071; 

// CSVL : IP3D
// float eff_IP[11];
eff_IP[0] = 0.908644; eff_IP[1] = 0.908466; eff_IP[2] = 0.842129; 
eff_IP[3] = 0.790594; eff_IP[4] = 0.863286; eff_IP[5] = 0.815354; 
eff_IP[6] = 0.825874; eff_IP[7] = 0.771681; eff_IP[8] = 0.838201; 
eff_IP[9] = 0.84729; eff_IP[10] = 0.830375; 
// float eff_stat_IP[11];
eff_stat_IP[0] = 0.0625402; eff_stat_IP[1] = 0.0226021; eff_stat_IP[2] = 0.0368765; 
eff_stat_IP[3] = 0.0182494; eff_stat_IP[4] = 0.0217368; eff_stat_IP[5] = 0.0251022; 
eff_stat_IP[6] = 0.0116825; eff_stat_IP[7] = 0.0153187; eff_stat_IP[8] = 0.0109495; 
eff_stat_IP[9] = 0.0129985; eff_stat_IP[10] = 0.0225841; 
// float effMC_IP[11];
effMC_IP[0] = 0.808883; effMC_IP[1] = 0.865865; effMC_IP[2] = 0.858398; 
effMC_IP[3] = 0.874222; effMC_IP[4] = 0.871152; effMC_IP[5] = 0.870212; 
effMC_IP[6] = 0.878178; effMC_IP[7] = 0.867721; effMC_IP[8] = 0.876344; 
effMC_IP[9] = 0.864927; effMC_IP[10] = 0.840546; 
// float effMC_stat_IP[11];
effMC_stat_IP[0] = 0.00236328; effMC_stat_IP[1] = 0.00138263; effMC_stat_IP[2] = 0.00129665; 
effMC_stat_IP[3] = 0.00160871; effMC_stat_IP[4] = 0.00159845; effMC_stat_IP[5] = 0.00166619; 
effMC_stat_IP[6] = 0.00174458; effMC_stat_IP[7] = 0.00173313; effMC_stat_IP[8] = 0.00224939; 
effMC_stat_IP[9] = 0.00245578; effMC_stat_IP[10] = 0.0044294; 
// float sf_IP[11];
sf_IP[0] = 1.12333; sf_IP[1] = 1.0492; sf_IP[2] = 0.981047; 
sf_IP[3] = 0.90434; sf_IP[4] = 0.99097; sf_IP[5] = 0.936959; 
sf_IP[6] = 0.94044; sf_IP[7] = 0.88932; sf_IP[8] = 0.956475; 
sf_IP[9] = 0.979608; sf_IP[10] = 0.9879; 
// float sf_stat_IP[11];
sf_stat_IP[0] = 0.0773863; sf_stat_IP[1] = 0.0261572; sf_stat_IP[2] = 0.0429853; 
sf_stat_IP[3] = 0.0209413; sf_stat_IP[4] = 0.025018; sf_stat_IP[5] = 0.0289017; 
sf_stat_IP[6] = 0.0134337; sf_stat_IP[7] = 0.0177431; sf_stat_IP[8] = 0.0127334; 
sf_stat_IP[9] = 0.0152837; sf_stat_IP[10] = 0.0273681; 
// float away_syst_IP[11];
away_syst_IP[0] = 0.033603; away_syst_IP[1] = 0.110268; away_syst_IP[2] = 0.103105; 
away_syst_IP[3] = 0.0380171; away_syst_IP[4] = 0.0416589; away_syst_IP[5] = 0.0393883; 
away_syst_IP[6] = 0.0268172; away_syst_IP[7] = 0.0253595; away_syst_IP[8] = 0.00603025; 
away_syst_IP[9] = 0.0061761; away_syst_IP[10] = 0.0053363; 
// float mupt_syst_IP[11];
mupt_syst_IP[0] = 0.110314; mupt_syst_IP[1] = 0.0781879; mupt_syst_IP[2] = 0.073109; 
mupt_syst_IP[3] = 0.0213076; mupt_syst_IP[4] = 0.0233487; mupt_syst_IP[5] = 0.0220761; 
mupt_syst_IP[6] = 0.0153233; mupt_syst_IP[7] = 0.0144903; mupt_syst_IP[8] = 0.0110714; 
mupt_syst_IP[9] = 0.0113392; mupt_syst_IP[10] = 0.00530951; 
// float gluon_syst_IP[11];
gluon_syst_IP[0] = 0.0122367; gluon_syst_IP[1] = 0.00237171; gluon_syst_IP[2] = 0.00221765; 
gluon_syst_IP[3] = 0.00355744; gluon_syst_IP[4] = 0.00389822; gluon_syst_IP[5] = 0.00368576; 
gluon_syst_IP[6] = 0.00026031; gluon_syst_IP[7] = 0.000246161; gluon_syst_IP[8] = 0.00546305; 
gluon_syst_IP[9] = 0.00559518; gluon_syst_IP[10] = 0.0171463; 
// float pu_syst_IP[11];
pu_syst_IP[0] = 0.0228668; pu_syst_IP[1] = 0.0345334; pu_syst_IP[2] = 0.0322902; 
pu_syst_IP[3] = 0.0177589; pu_syst_IP[4] = 0.0194601; pu_syst_IP[5] = 0.0183995; 
pu_syst_IP[6] = 0.00312689; pu_syst_IP[7] = 0.00295692; pu_syst_IP[8] = 0.00243802; 
pu_syst_IP[9] = 0.00249699; pu_syst_IP[10] = 0.0150145; 
// float total_syst_IP[11];
total_syst_IP[0] = 0.118199; total_syst_IP[1] = 0.139537; total_syst_IP[2] = 0.130473; 
total_syst_IP[3] = 0.0471948; total_syst_IP[4] = 0.0517157; total_syst_IP[5] = 0.0488971; 
total_syst_IP[6] = 0.0310453; total_syst_IP[7] = 0.0293577; total_syst_IP[8] = 0.0139545; 
total_syst_IP[9] = 0.014292; total_syst_IP[10] = 0.024002; 

// CSVL : JP
eff_JP[0] = 0.791185; eff_stat_JP[0] = 0.0139995; effMC_JP[0] = 0.809235; effMC_stat_JP[0] = 0.000910193; 
sf_JP[0] = 0.977695; sf_stat_JP[0] = 0.0172996; sf_syst_JP[0] = 0.0821264; sf_eror_JP[0] = 0.0839287; 
mupt_syst_JP[0] = 0.006; gluon_syst_JP[0] = 0.003; pu_syst_JP[0] = 0.002; total_syst_JP[0] = 0.084; 
cor_syst_JP[0] = 0.044; inc_syst_JP[0] = 0.071; bias_syst_JP[0] = 0; 
eff_JP[1] = 0.857213; eff_stat_JP[1] = 0.02551; effMC_JP[1] = 0.864676; effMC_stat_JP[1] = 0.00059928; 
sf_JP[1] = 0.991369; sf_stat_JP[1] = 0.0295024; sf_syst_JP[1] = 0.0297411; sf_eror_JP[1] = 0.0418918; 
mupt_syst_JP[1] = 0.006; gluon_syst_JP[1] = 0.003; pu_syst_JP[1] = 0.001; total_syst_JP[1] = 0.03; 
cor_syst_JP[1] = 0.022; inc_syst_JP[1] = 0.02; bias_syst_JP[1] = 0; 
eff_JP[2] = 0.833568; eff_stat_JP[2] = 0.00808785; effMC_JP[2] = 0.84846; effMC_stat_JP[2] = 0.000675656; 
sf_JP[2] = 0.982448; sf_stat_JP[2] = 0.00953239; sf_syst_JP[2] = 0.0245612; sf_eror_JP[2] = 0.0263461; 
mupt_syst_JP[2] = 0.006; gluon_syst_JP[2] = 0.003; pu_syst_JP[2] = 0; total_syst_JP[2] = 0.025; 
cor_syst_JP[2] = 0.013; inc_syst_JP[2] = 0.02; bias_syst_JP[2] = 0; 
eff_JP[3] = 0.872241; eff_stat_JP[3] = 0.00765432; effMC_JP[3] = 0.882744; effMC_stat_JP[3] = 0.000839049; 
sf_JP[3] = 0.988102; sf_stat_JP[3] = 0.00867105; sf_syst_JP[3] = 0.0247026; sf_eror_JP[3] = 0.0261802; 
mupt_syst_JP[3] = 0.006; gluon_syst_JP[3] = 0.004; pu_syst_JP[3] = 0; total_syst_JP[3] = 0.025; 
cor_syst_JP[3] = 0.008; inc_syst_JP[3] = 0.023; bias_syst_JP[3] = 0; 
eff_JP[4] = 0.850999; eff_stat_JP[4] = 0.00999902; effMC_JP[4] = 0.870583; effMC_stat_JP[4] = 0.000876378; 
sf_JP[4] = 0.977506; sf_stat_JP[4] = 0.0114854; sf_syst_JP[4] = 0.0244376; sf_eror_JP[4] = 0.0270021; 
mupt_syst_JP[4] = 0.006; gluon_syst_JP[4] = 0.005; pu_syst_JP[4] = 0; total_syst_JP[4] = 0.025; 
cor_syst_JP[4] = 0.006; inc_syst_JP[4] = 0.023; bias_syst_JP[4] = 0; 
eff_JP[5] = 0.855536; eff_stat_JP[5] = 0.00958563; effMC_JP[5] = 0.875968; effMC_stat_JP[5] = 0.000896783; 
sf_JP[5] = 0.976675; sf_stat_JP[5] = 0.0109429; sf_syst_JP[5] = 0.0244169; sf_eror_JP[5] = 0.0267569; 
mupt_syst_JP[5] = 0.006; gluon_syst_JP[5] = 0.005; pu_syst_JP[5] = 0; total_syst_JP[5] = 0.025; 
cor_syst_JP[5] = 0.005; inc_syst_JP[5] = 0.023; bias_syst_JP[5] = 0; 
eff_JP[6] = 0.867532; eff_stat_JP[6] = 0.0075736; effMC_JP[6] = 0.885392; effMC_stat_JP[6] = 0.000829844; 
sf_JP[6] = 0.979829; sf_stat_JP[6] = 0.00855395; sf_syst_JP[6] = 0.0205764; sf_eror_JP[6] = 0.0222836; 
mupt_syst_JP[6] = 0.003; gluon_syst_JP[6] = 0.006; pu_syst_JP[6] = 0; total_syst_JP[6] = 0.021; 
cor_syst_JP[6] = 0.004; inc_syst_JP[6] = 0.02; bias_syst_JP[6] = 0; 
eff_JP[7] = 0.854256; eff_stat_JP[7] = 0.00995237; effMC_JP[7] = 0.876279; effMC_stat_JP[7] = 0.00086067; 
sf_JP[7] = 0.974868; sf_stat_JP[7] = 0.0113575; sf_syst_JP[7] = 0.0214471; sf_eror_JP[7] = 0.0242687; 
mupt_syst_JP[7] = 0.003; gluon_syst_JP[7] = 0.007; pu_syst_JP[7] = 0; total_syst_JP[7] = 0.022; 
cor_syst_JP[7] = 0.003; inc_syst_JP[7] = 0.02; bias_syst_JP[7] = 0; 
eff_JP[8] = 0.858915; eff_stat_JP[8] = 0.0107836; effMC_JP[8] = 0.881023; effMC_stat_JP[8] = 0.000814757; 
sf_JP[8] = 0.974907; sf_stat_JP[8] = 0.0122399; sf_syst_JP[8] = 0.021448; sf_eror_JP[8] = 0.0246947; 
mupt_syst_JP[8] = 0.003; gluon_syst_JP[8] = 0.009; pu_syst_JP[8] = 0; total_syst_JP[8] = 0.022; 
cor_syst_JP[8] = 0.003; inc_syst_JP[8] = 0.02; bias_syst_JP[8] = 0; 
eff_JP[9] = 0.851312; eff_stat_JP[9] = 0.0106489; effMC_JP[9] = 0.877352; effMC_stat_JP[9] = 0.00100439; 
sf_JP[9] = 0.97032; sf_stat_JP[9] = 0.0121376; sf_syst_JP[9] = 0.0436644; sf_eror_JP[9] = 0.04532; 
mupt_syst_JP[9] = 0.003; gluon_syst_JP[9] = 0.01; pu_syst_JP[9] = 0; total_syst_JP[9] = 0.045; 
cor_syst_JP[9] = 0.002; inc_syst_JP[9] = 0.044; bias_syst_JP[9] = 0; 
eff_JP[10] = 0.838567; eff_stat_JP[10] = 0.0137222; effMC_JP[10] = 0.8718; effMC_stat_JP[10] = 0.00182387; 
sf_JP[10] = 0.96188; sf_stat_JP[10] = 0.0157401; sf_syst_JP[10] = 0.0452084; sf_eror_JP[10] = 0.0478701; 
mupt_syst_JP[10] = 0.004; gluon_syst_JP[10] = 0.017; pu_syst_JP[10] = 0; total_syst_JP[10] = 0.047; 
cor_syst_JP[10] = 0.002; inc_syst_JP[10] = 0.044; bias_syst_JP[10] = 0; 
eff_JP[11] = 0.820111; eff_stat_JP[11] = 0.0130472; effMC_JP[11] = 0.859662; effMC_stat_JP[11] = 0.00318813; 
sf_JP[11] = 0.953992; sf_stat_JP[11] = 0.0151772; sf_syst_JP[11] = 0.0457916; sf_eror_JP[11] = 0.0482413; 
mupt_syst_JP[11] = 0.004; gluon_syst_JP[11] = 0.018; pu_syst_JP[11] = 0.003; total_syst_JP[11] = 0.048; 
cor_syst_JP[11] = 0.001; inc_syst_JP[11] = 0.044; bias_syst_JP[11] = 0; 
eff_JP[12] = 0.825898; eff_stat_JP[12] = 0.0167122; effMC_JP[12] = 0.843868; effMC_stat_JP[12] = 0.00566296; 
sf_JP[12] = 0.978705; sf_stat_JP[12] = 0.0198043; sf_syst_JP[12] = 0.0322973; sf_eror_JP[12] = 0.0378857; 
mupt_syst_JP[12] = 0.004; gluon_syst_JP[12] = 0.006; pu_syst_JP[12] = 0.003; total_syst_JP[12] = 0.033; 
cor_syst_JP[12] = 0.001; inc_syst_JP[12] = 0.032; bias_syst_JP[12] = 0; 
eff_JP[13] = 0.792619; eff_stat_JP[13] = 0.012835; effMC_JP[13] = 0.842621; effMC_stat_JP[13] = 0.00498551; 
sf_JP[13] = 0.940659; sf_stat_JP[13] = 0.0152323; sf_syst_JP[13] = 0.0310417; sf_eror_JP[13] = 0.0345776; 
mupt_syst_JP[13] = 0.004; gluon_syst_JP[13] = 0.005; pu_syst_JP[13] = 0.006; total_syst_JP[13] = 0.033; 
cor_syst_JP[13] = 0.001; inc_syst_JP[13] = 0.032; bias_syst_JP[13] = 0; 
eff_JP[14] = 0.731248; eff_stat_JP[14] = 0.0478664; effMC_JP[14] = 0.794927; effMC_stat_JP[14] = 0.00466948; 
sf_JP[14] = 0.919893; sf_stat_JP[14] = 0.0602148; sf_syst_JP[14] = 0.034036; sf_eror_JP[14] = 0.0691684; 
mupt_syst_JP[14] = 0.004; gluon_syst_JP[14] = 0.017; pu_syst_JP[14] = 0.001; total_syst_JP[14] = 0.037; 
cor_syst_JP[14] = 0.008; inc_syst_JP[14] = 0.032; bias_syst_JP[14] = 0; 
}

// #############################################################################
if ( title == "CSVM" ) { 
// CSVM : PTREL
// float eff_PT[11];
eff_PT[0] = 0.427749; eff_PT[1] = 0.577202; eff_PT[2] = 0.640131; 
eff_PT[3] = 0.660343; eff_PT[4] = 0.690354; eff_PT[5] = 0.699875; 
eff_PT[6] = 0.710442; eff_PT[7] = 0.708815; eff_PT[8] = 0.698329; 
eff_PT[9] = 0.649629; eff_PT[10] = 0.513398; 
// float eff_stat_PT[11];
eff_stat_PT[0] = 0.010402; eff_stat_PT[1] = 0.00774761; eff_stat_PT[2] = 0.00856539; 
eff_stat_PT[3] = 0.00877695; eff_stat_PT[4] = 0.0108105; eff_stat_PT[5] = 0.0141391; 
eff_stat_PT[6] = 0.01127; eff_stat_PT[7] = 0.0120509; eff_stat_PT[8] = 0.0157443; 
eff_stat_PT[9] = 0.0132134; eff_stat_PT[10] = 0.0238353; 
// float effMC_PT[11];
effMC_PT[0] = 0.483367; effMC_PT[1] = 0.62171; effMC_PT[2] = 0.66664; 
effMC_PT[3] = 0.715432; effMC_PT[4] = 0.728413; effMC_PT[5] = 0.749377; 
effMC_PT[6] = 0.769722; effMC_PT[7] = 0.758534; effMC_PT[8] = 0.768543; 
effMC_PT[9] = 0.737991; effMC_PT[10] = 0.682895; 
// float effMC_stat_PT[11];
effMC_stat_PT[0] = 0.00338359; effMC_stat_PT[1] = 0.002047; effMC_stat_PT[2] = 0.00180723; 
effMC_stat_PT[3] = 0.00239412; effMC_stat_PT[4] = 0.00231625; effMC_stat_PT[5] = 0.00237728; 
effMC_stat_PT[6] = 0.00250384; effMC_stat_PT[7] = 0.00242423; effMC_stat_PT[8] = 0.00316062; 
effMC_stat_PT[9] = 0.00349352; effMC_stat_PT[10] = 0.00574427; 
// float sf_PT[11];
sf_PT[0] = 0.884936; sf_PT[1] = 0.92841; sf_PT[2] = 0.960235; 
sf_PT[3] = 0.922999; sf_PT[4] = 0.94775; sf_PT[5] = 0.933943; 
sf_PT[6] = 0.922986; sf_PT[7] = 0.934454; sf_PT[8] = 0.90864; 
sf_PT[9] = 0.880267; sf_PT[10] = 0.751796; 
// float sf_stat_PT[11];
sf_stat_PT[0] = 0.0223937; sf_stat_PT[1] = 0.0128312; sf_stat_PT[2] = 0.0131096; 
sf_stat_PT[3] = 0.0126509; sf_stat_PT[4] = 0.015144; sf_stat_PT[5] = 0.0190991; 
sf_stat_PT[6] = 0.0149463; sf_stat_PT[7] = 0.0161653; sf_stat_PT[8] = 0.020824; 
sf_stat_PT[9] = 0.0183831; sf_stat_PT[10] = 0.0354716; 
// float away_syst_PT[11];
away_syst_PT[0] = 0.012944; away_syst_PT[1] = 0.0331379; away_syst_PT[2] = 0.0342739; 
away_syst_PT[3] = 0.0217786; away_syst_PT[4] = 0.0223626; away_syst_PT[5] = 0.0220368; 
away_syst_PT[6] = 0.0112504; away_syst_PT[7] = 0.0113902; away_syst_PT[8] = 0.0742381; 
away_syst_PT[9] = 0.0719199; away_syst_PT[10] = 0.070153; 
// float mupt_syst_PT[11];
mupt_syst_PT[0] = 0.037003; mupt_syst_PT[1] = 0.00312242; mupt_syst_PT[2] = 0.00322946; 
mupt_syst_PT[3] = 0.0114293; mupt_syst_PT[4] = 0.0117358; mupt_syst_PT[5] = 0.0115648; 
mupt_syst_PT[6] = 0.0120678; mupt_syst_PT[7] = 0.0122177; mupt_syst_PT[8] = 0.0692652; 
mupt_syst_PT[9] = 0.0671024; mupt_syst_PT[10] = 0.039248; 
// float lc_syst_PT[11];
lc_syst_PT[0] = 0.00539196; lc_syst_PT[1] = 0.00488161; lc_syst_PT[2] = 0.00504894; 
lc_syst_PT[3] = 0.00622145; lc_syst_PT[4] = 0.00638829; lc_syst_PT[5] = 0.00629522; 
lc_syst_PT[6] = 0.00411453; lc_syst_PT[7] = 0.00416566; lc_syst_PT[8] = 0.0215825; 
lc_syst_PT[9] = 0.0209086; lc_syst_PT[10] = 0.054694; 
// float gluon_syst_PT[11];
gluon_syst_PT[0] = 0.004417; gluon_syst_PT[1] = 0.000448191; gluon_syst_PT[2] = 0.000463554; 
gluon_syst_PT[3] = 0.00100571; gluon_syst_PT[4] = 0.00103268; gluon_syst_PT[5] = 0.00101763; 
gluon_syst_PT[6] = 0.0114995; gluon_syst_PT[7] = 0.0116424; gluon_syst_PT[8] = 0.0409549; 
gluon_syst_PT[9] = 0.0396761; gluon_syst_PT[10] = 0.043554; 
// float pu_syst_PT[11];
pu_syst_PT[0] = 0.020309; pu_syst_PT[1] = 0.0143572; pu_syst_PT[2] = 0.0148493; 
pu_syst_PT[3] = 0.00835967; pu_syst_PT[4] = 0.00858384; pu_syst_PT[5] = 0.00845879; 
pu_syst_PT[6] = 0.0124245; pu_syst_PT[7] = 0.0125789; pu_syst_PT[8] = 0.00805344; 
pu_syst_PT[9] = 0.00780197; pu_syst_PT[10] = 0.011188; 
// float total_syst_PT[11];
total_syst_PT[0] = 0.0446969; total_syst_PT[1] = 0.0365791; total_syst_PT[2] = 0.037833; 
total_syst_PT[3] = 0.0267308; total_syst_PT[4] = 0.0274476; total_syst_PT[5] = 0.0270478; 
total_syst_PT[6] = 0.0239946; total_syst_PT[7] = 0.0242927; total_syst_PT[8] = 0.111879; 
total_syst_PT[9] = 0.108386; total_syst_PT[10] = 0.107123; 

// CSVM : S8
eff_SY[0]=0.4316; eff_SY[1]=0.5630; eff_SY[2]=0.6231; eff_SY[3]=0.6703; eff_SY[4]=0.6720; eff_SY[5]=0.7428; eff_SY[6]=0.7390; eff_SY[7]=0.6214; 
eff_stat_SY[0]=0.0099; eff_stat_SY[1]=0.0074; eff_stat_SY[2]=0.0098; eff_stat_SY[3]=0.0109; eff_stat_SY[4]=0.0119; eff_stat_SY[5]=0.0244; eff_stat_SY[6]=0.0392; eff_stat_SY[7]=0.0308; 
effMC_SY[0]=0.4260; effMC_SY[1]=0.5729; effMC_SY[2]=0.6422; effMC_SY[3]=0.6945; effMC_SY[4]=0.7329; effMC_SY[5]=0.7631; effMC_SY[6]=0.7691; effMC_SY[7]=0.7743; 
effMC_stat_SY[0]=0.0007; effMC_stat_SY[1]=0.0007; effMC_stat_SY[2]=0.0008; effMC_stat_SY[3]=0.0008; effMC_stat_SY[4]=0.0007; effMC_stat_SY[5]=0.0008; effMC_stat_SY[6]=0.0010; effMC_stat_SY[7]=0.0009; 
sf_SY[0]=1.0131; sf_SY[1]=0.9827; sf_SY[2]=0.9703; sf_SY[3]=0.9652; sf_SY[4]=0.9169; sf_SY[5]=0.9734; sf_SY[6]=0.9609; sf_SY[7]=0.8025; 
sf_stat_SY[0]=0.0233; sf_stat_SY[1]=0.013; sf_stat_SY[2]=0.0153; sf_stat_SY[3]=0.0157; sf_stat_SY[4]=0.0163; sf_stat_SY[5]=0.032; sf_stat_SY[6]=0.051; sf_stat_SY[7]=0.0398; 
away[0]=2.5556; away[1]=0.1803; away[2]=0.7684; away[3]=14.8366; 
mu[0]=0.0741; mu[1]=0.4958; mu[2]=1.2898; mu[3]=0.1471; 
ptrel[0]=1.0556; ptrel[1]=1.6977; ptrel[2]=4.4457; ptrel[3]=0.9314; 
gsplit[0]=0.1811; gsplit[1]=0.1821; gsplit[2]=0.3525; gsplit[3]=3.7912; 
closure[0]=0.2536; closure[1]=0.084; closure[2]=0.0653; closure[3]=0.2856; 
pu[0]=0.0185; pu[1]=0.1803; pu[2]=1.5231; pu[3]=1.0948; 

// CSVM : IP3D
// float eff_IP[11];
eff_IP[0] = 0.405726; eff_IP[1] = 0.530795; eff_IP[2] = 0.572392; 
eff_IP[3] = 0.594873; eff_IP[4] = 0.674357; eff_IP[5] = 0.662509; 
eff_IP[6] = 0.70511; eff_IP[7] = 0.65881; eff_IP[8] = 0.721984; 
eff_IP[9] = 0.699425; eff_IP[10] = 0.719678; 
// float eff_stat_IP[11];
eff_stat_IP[0] = 0.080747; eff_stat_IP[1] = 0.0167383; eff_stat_IP[2] = 0.0140157; 
eff_stat_IP[3] = 0.0144631; eff_stat_IP[4] = 0.0198313; eff_stat_IP[5] = 0.023779; 
eff_stat_IP[6] = 0.0119878; eff_stat_IP[7] = 0.0137152; eff_stat_IP[8] = 0.0112713; 
eff_stat_IP[9] = 0.0130785; eff_stat_IP[10] = 0.0242813; 
// float effMC_IP[11];
effMC_IP[0] = 0.486491; effMC_IP[1] = 0.62092; effMC_IP[2] = 0.665676; 
effMC_IP[3] = 0.715524; effMC_IP[4] = 0.728406; effMC_IP[5] = 0.745225; 
effMC_IP[6] = 0.768805; effMC_IP[7] = 0.758541; effMC_IP[8] = 0.765092; 
effMC_IP[9] = 0.73642; effMC_IP[10] = 0.685582; 
// float effMC_stat_IP[11];
effMC_stat_IP[0] = 0.00336654; effMC_stat_IP[1] = 0.00219238; effMC_stat_IP[2] = 0.00190471; 
effMC_stat_IP[3] = 0.00238042; effMC_stat_IP[4] = 0.00228317; effMC_stat_IP[5] = 0.00233668; 
effMC_stat_IP[6] = 0.00243909; effMC_stat_IP[7] = 0.00235218; effMC_stat_IP[8] = 0.00306392; 
effMC_stat_IP[9] = 0.00336583; effMC_stat_IP[10] = 0.00574297; 
// float sf_IP[11];
sf_IP[0] = 0.833986; sf_IP[1] = 0.854853; sf_IP[2] = 0.859865; 
sf_IP[3] = 0.83138; sf_IP[4] = 0.925799; sf_IP[5] = 0.889005; 
sf_IP[6] = 0.917151; sf_IP[7] = 0.868522; sf_IP[8] = 0.943657; 
sf_IP[9] = 0.949764; sf_IP[10] = 1.04973; 
// float sf_stat_IP[11];
sf_stat_IP[0] = 0.166079; sf_stat_IP[1] = 0.0271258; sf_stat_IP[2] = 0.0211981; 
sf_stat_IP[3] = 0.0204016; sf_stat_IP[4] = 0.0273799; sf_stat_IP[5] = 0.0320301; 
sf_stat_IP[6] = 0.015862; sf_stat_IP[7] = 0.0182805; sf_stat_IP[8] = 0.015209; 
sf_stat_IP[9] = 0.0182824; sf_stat_IP[10] = 0.0364923; 
// float away_syst_IP[11];
away_syst_IP[0] = 0.0359282; away_syst_IP[1] = 0.0149812; away_syst_IP[2] = 0.015069; 
away_syst_IP[3] = 0.038133; away_syst_IP[4] = 0.0424637; away_syst_IP[5] = 0.0407761; 
away_syst_IP[6] = 0.0379435; away_syst_IP[7] = 0.0359317; away_syst_IP[8] = 0.0300011; 
away_syst_IP[9] = 0.0301952; away_syst_IP[10] = 0.0296933; 
// float mupt_syst_IP[11];
mupt_syst_IP[0] = 0.0085853; mupt_syst_IP[1] = 0.0181923; mupt_syst_IP[2] = 0.018299; 
mupt_syst_IP[3] = 0.0183783; mupt_syst_IP[4] = 0.0204655; mupt_syst_IP[5] = 0.0196521; 
mupt_syst_IP[6] = 0.0100962; mupt_syst_IP[7] = 0.00956084; mupt_syst_IP[8] = 0.000949506; 
mupt_syst_IP[9] = 0.000955651; mupt_syst_IP[10] = 0.0105754; 
// float gluon_syst_IP[11];
gluon_syst_IP[0] = 0.00707265; gluon_syst_IP[1] = 0.0277943; gluon_syst_IP[2] = 0.0279573; 
gluon_syst_IP[3] = 0.000363897; gluon_syst_IP[4] = 0.000405224; gluon_syst_IP[5] = 0.000389119; 
gluon_syst_IP[6] = 0.00194754; gluon_syst_IP[7] = 0.00184427; gluon_syst_IP[8] = 0.0218319; 
gluon_syst_IP[9] = 0.0219731; gluon_syst_IP[10] = 0.0741578; 
// float pu_syst_IP[11];
pu_syst_IP[0] = 0.0272092; pu_syst_IP[1] = 0.0138931; pu_syst_IP[2] = 0.0139746; 
pu_syst_IP[3] = 0.00639768; pu_syst_IP[4] = 0.00712426; pu_syst_IP[5] = 0.00684112; 
pu_syst_IP[6] = 0.0060448; pu_syst_IP[7] = 0.0057243; pu_syst_IP[8] = 0.00519886; 
pu_syst_IP[9] = 0.00523251; pu_syst_IP[10] = 0.0387592; 
// float total_syst_IP[11];
total_syst_IP[0] = 0.0464209; total_syst_IP[1] = 0.0389992; total_syst_IP[2] = 0.0392279; 
total_syst_IP[3] = 0.0428129; total_syst_IP[4] = 0.0476752; total_syst_IP[5] = 0.0457804; 
total_syst_IP[6] = 0.0397741; total_syst_IP[7] = 0.0376652; total_syst_IP[8] = 0.0374783; 
total_syst_IP[9] = 0.0377208; total_syst_IP[10] = 0.0894159; 

// CSVM : JP
eff_JP[0] = 0.50017; eff_stat_JP[0] = 0.00977812; effMC_JP[0] = 0.534872; effMC_stat_JP[0] = 0.00115547; 
sf_JP[0] = 0.93512; sf_stat_JP[0] = 0.0182812; sf_syst_JP[0] = 0.110344; sf_eror_JP[0] = 0.111848; 
mupt_syst_JP[0] = 0.002; gluon_syst_JP[0] = 0.004; pu_syst_JP[0] = 0.001; total_syst_JP[0] = 0.118; 
cor_syst_JP[0] = 0.044; inc_syst_JP[0] = 0.109; bias_syst_JP[0] = 0; 
eff_JP[1] = 0.620835; eff_stat_JP[1] = 0.00736126; effMC_JP[1] = 0.653073; effMC_stat_JP[1] = 0.000833902; 
sf_JP[1] = 0.950637; sf_stat_JP[1] = 0.0112717; sf_syst_JP[1] = 0.0617914; sf_eror_JP[1] = 0.0628111; 
mupt_syst_JP[1] = 0.002; gluon_syst_JP[1] = 0.002; pu_syst_JP[1] = 0; total_syst_JP[1] = 0.065; 
cor_syst_JP[1] = 0.022; inc_syst_JP[1] = 0.061; bias_syst_JP[1] = 0; 
eff_JP[2] = 0.641859; eff_stat_JP[2] = 0.00658748; effMC_JP[2] = 0.674577; effMC_stat_JP[2] = 0.000882849; 
sf_JP[2] = 0.951499; sf_stat_JP[2] = 0.00976535; sf_syst_JP[2] = 0.0589929; sf_eror_JP[2] = 0.0597957; 
mupt_syst_JP[2] = 0.002; gluon_syst_JP[2] = 0.002; pu_syst_JP[2] = 0; total_syst_JP[2] = 0.062; 
cor_syst_JP[2] = 0.013; inc_syst_JP[2] = 0.061; bias_syst_JP[2] = 0; 
eff_JP[3] = 0.704899; eff_stat_JP[3] = 0.0065762; effMC_JP[3] = 0.732992; effMC_stat_JP[3] = 0.00115376; 
sf_JP[3] = 0.961674; sf_stat_JP[3] = 0.00897173; sf_syst_JP[3] = 0.0326969; sf_eror_JP[3] = 0.0339055; 
mupt_syst_JP[3] = 0.002; gluon_syst_JP[3] = 0.004; pu_syst_JP[3] = 0; total_syst_JP[3] = 0.034; 
cor_syst_JP[3] = 0.008; inc_syst_JP[3] = 0.033; bias_syst_JP[3] = 0; 
eff_JP[4] = 0.71192; eff_stat_JP[4] = 0.00719148; effMC_JP[4] = 0.740993; effMC_stat_JP[4] = 0.00114381; 
sf_JP[4] = 0.960764; sf_stat_JP[4] = 0.00970519; sf_syst_JP[4] = 0.032666; sf_eror_JP[4] = 0.0340772; 
mupt_syst_JP[4] = 0.002; gluon_syst_JP[4] = 0.004; pu_syst_JP[4] = 0.001; total_syst_JP[4] = 0.034; 
cor_syst_JP[4] = 0.006; inc_syst_JP[4] = 0.033; bias_syst_JP[4] = 0; 
eff_JP[5] = 0.726684; eff_stat_JP[5] = 0.00706883; effMC_JP[5] = 0.754914; effMC_stat_JP[5] = 0.00117027; 
sf_JP[5] = 0.962605; sf_stat_JP[5] = 0.00936376; sf_syst_JP[5] = 0.0327286; sf_eror_JP[5] = 0.0340417; 
mupt_syst_JP[5] = 0.002; gluon_syst_JP[5] = 0.006; pu_syst_JP[5] = 0; total_syst_JP[5] = 0.034; 
cor_syst_JP[5] = 0.005; inc_syst_JP[5] = 0.033; bias_syst_JP[5] = 0; 
eff_JP[6] = 0.752098; eff_stat_JP[6] = 0.0051908; effMC_JP[6] = 0.775147; effMC_stat_JP[6] = 0.00108758; 
sf_JP[6] = 0.970265; sf_stat_JP[6] = 0.00669654; sf_syst_JP[6] = 0.0261972; sf_eror_JP[6] = 0.0270395; 
mupt_syst_JP[6] = 0.004; gluon_syst_JP[6] = 0.006; pu_syst_JP[6] = 0; total_syst_JP[6] = 0.027; 
cor_syst_JP[6] = 0.004; inc_syst_JP[6] = 0.026; bias_syst_JP[6] = 0; 
eff_JP[7] = 0.753491; eff_stat_JP[7] = 0.006892; effMC_JP[7] = 0.773253; effMC_stat_JP[7] = 0.00109452; 
sf_JP[7] = 0.974443; sf_stat_JP[7] = 0.00891299; sf_syst_JP[7] = 0.0272844; sf_eror_JP[7] = 0.0287033; 
mupt_syst_JP[7] = 0.004; gluon_syst_JP[7] = 0.009; pu_syst_JP[7] = 0.001; total_syst_JP[7] = 0.028; 
cor_syst_JP[7] = 0.003; inc_syst_JP[7] = 0.026; bias_syst_JP[7] = 0; 
eff_JP[8] = 0.755126; eff_stat_JP[8] = 0.00756375; effMC_JP[8] = 0.775292; effMC_stat_JP[8] = 0.00105038; 
sf_JP[8] = 0.973989; sf_stat_JP[8] = 0.009756; sf_syst_JP[8] = 0.0282457; sf_eror_JP[8] = 0.0298831; 
mupt_syst_JP[8] = 0.004; gluon_syst_JP[8] = 0.011; pu_syst_JP[8] = 0; total_syst_JP[8] = 0.029; 
cor_syst_JP[8] = 0.003; inc_syst_JP[8] = 0.026; bias_syst_JP[8] = 0; 
eff_JP[9] = 0.734344; eff_stat_JP[9] = 0.00668295; effMC_JP[9] = 0.759375; effMC_stat_JP[9] = 0.00130883; 
sf_JP[9] = 0.967038; sf_stat_JP[9] = 0.00880059; sf_syst_JP[9] = 0.0522201; sf_eror_JP[9] = 0.0529564; 
mupt_syst_JP[9] = 0.004; gluon_syst_JP[9] = 0.015; pu_syst_JP[9] = 0; total_syst_JP[9] = 0.054; 
cor_syst_JP[9] = 0.002; inc_syst_JP[9] = 0.052; bias_syst_JP[9] = 0; 
eff_JP[10] = 0.715155; eff_stat_JP[10] = 0.00899763; effMC_JP[10] = 0.742023; effMC_stat_JP[10] = 0.00238694; 
sf_JP[10] = 0.963791; sf_stat_JP[10] = 0.0121258; sf_syst_JP[10] = 0.0558999; sf_eror_JP[10] = 0.0571999; 
mupt_syst_JP[10] = 0.009; gluon_syst_JP[10] = 0.025; pu_syst_JP[10] = 0.001; total_syst_JP[10] = 0.058; 
cor_syst_JP[10] = 0.002; inc_syst_JP[10] = 0.052; bias_syst_JP[10] = 0; 
eff_JP[11] = 0.710254; eff_stat_JP[11] = 0.0119764; effMC_JP[11] = 0.729611; effMC_stat_JP[11] = 0.00407685; 
sf_JP[11] = 0.973469; sf_stat_JP[11] = 0.0164148; sf_syst_JP[11] = 0.0584081; sf_eror_JP[11] = 0.0606709; 
mupt_syst_JP[11] = 0.009; gluon_syst_JP[11] = 0.028; pu_syst_JP[11] = 0.005; total_syst_JP[11] = 0.06; 
cor_syst_JP[11] = 0.001; inc_syst_JP[11] = 0.052; bias_syst_JP[11] = 0; 
eff_JP[12] = 0.690667; eff_stat_JP[12] = 0.0159402; effMC_JP[12] = 0.706027; effMC_stat_JP[12] = 0.00710763; 
sf_JP[12] = 0.978244; sf_stat_JP[12] = 0.0225773; sf_syst_JP[12] = 0.0577164; sf_eror_JP[12] = 0.0619751; 
mupt_syst_JP[12] = 0.009; gluon_syst_JP[12] = 0.004; pu_syst_JP[12] = 0.009; total_syst_JP[12] = 0.059; 
cor_syst_JP[12] = 0.001; inc_syst_JP[12] = 0.057; bias_syst_JP[12] = 0; 
eff_JP[13] = 0.671027; eff_stat_JP[13] = 0.0167989; effMC_JP[13] = 0.687361; effMC_stat_JP[13] = 0.0063465; 
sf_JP[13] = 0.976237; sf_stat_JP[13] = 0.0244397; sf_syst_JP[13] = 0.0595505; sf_eror_JP[13] = 0.0643705; 
mupt_syst_JP[13] = 0.009; gluon_syst_JP[13] = 0.006; pu_syst_JP[13] = 0.02; total_syst_JP[13] = 0.061; 
cor_syst_JP[13] = 0.001; inc_syst_JP[13] = 0.057; bias_syst_JP[13] = 0; 
eff_JP[14] = 0.584277; eff_stat_JP[14] = 0.035261; effMC_JP[14] = 0.633593; effMC_stat_JP[14] = 0.00557234; 
sf_JP[14] = 0.922164; sf_stat_JP[14] = 0.0556524; sf_syst_JP[14] = 0.0544077; sf_eror_JP[14] = 0.0778292; 
mupt_syst_JP[14] = 0.009; gluon_syst_JP[14] = 0.007; pu_syst_JP[14] = 0.003; total_syst_JP[14] = 0.059; 
cor_syst_JP[14] = 0.008; inc_syst_JP[14] = 0.057; bias_syst_JP[14] = 0; 
}

// #############################################################################
if ( title == "CSVT" ) { 
// CSVT : PTREL
// float eff_PT[11];
eff_PT[0] = 0.294644; eff_PT[1] = 0.433698; eff_PT[2] = 0.490741; 
eff_PT[3] = 0.522193; eff_PT[4] = 0.512071; eff_PT[5] = 0.523333; 
eff_PT[6] = 0.527521; eff_PT[7] = 0.485651; eff_PT[8] = 0.455473; 
eff_PT[9] = 0.329993; eff_PT[10] = 0.2243; 
// float eff_stat_PT[11];
eff_stat_PT[0] = 0.0089674; eff_stat_PT[1] = 0.0074363; eff_stat_PT[2] = 0.00849531; 
eff_stat_PT[3] = 0.0076851; eff_stat_PT[4] = 0.0104497; eff_stat_PT[5] = 0.0129143; 
eff_stat_PT[6] = 0.00973132; eff_stat_PT[7] = 0.0110197; eff_stat_PT[8] = 0.0129372; 
eff_stat_PT[9] = 0.00913738; eff_stat_PT[10] = 0.00889217; 
// float effMC_PT[11];
effMC_PT[0] = 0.34475; effMC_PT[1] = 0.479638; effMC_PT[2] = 0.522866; 
effMC_PT[3] = 0.578587; effMC_PT[4] = 0.578296; effMC_PT[5] = 0.59453; 
effMC_PT[6] = 0.609894; effMC_PT[7] = 0.563087; effMC_PT[8] = 0.567889; 
effMC_PT[9] = 0.464457; effMC_PT[10] = 0.417102; 
// float effMC_stat_PT[11];
effMC_stat_PT[0] = 0.00321621; effMC_stat_PT[1] = 0.00214066; effMC_stat_PT[2] = 0.00194703; 
effMC_stat_PT[3] = 0.00237747; effMC_stat_PT[4] = 0.00264093; effMC_stat_PT[5] = 0.00275393; 
effMC_stat_PT[6] = 0.00298342; effMC_stat_PT[7] = 0.00290315; effMC_stat_PT[8] = 0.0038415; 
effMC_stat_PT[9] = 0.00409515; effMC_stat_PT[10] = 0.00586605; 
// float sf_PT[11];
sf_PT[0] = 0.85466; sf_PT[1] = 0.90422; sf_PT[2] = 0.938561; 
sf_PT[3] = 0.902531; sf_PT[4] = 0.885483; sf_PT[5] = 0.880246; 
sf_PT[6] = 0.864939; sf_PT[7] = 0.86248; sf_PT[8] = 0.802046; 
sf_PT[9] = 0.710493; sf_PT[10] = 0.537758; 
// float sf_stat_PT[11];
sf_stat_PT[0] = 0.0272059; sf_stat_PT[1] = 0.0160206; sf_stat_PT[2] = 0.0166193; 
sf_stat_PT[3] = 0.0137905; sf_stat_PT[4] = 0.0185168; sf_stat_PT[5] = 0.0221013; 
sf_stat_PT[6] = 0.0165072; sf_stat_PT[7] = 0.0200691; sf_stat_PT[8] = 0.0234183; 
sf_stat_PT[9] = 0.0206466; sf_stat_PT[10] = 0.0226207; 
// float away_syst_PT[11];
away_syst_PT[0] = 0.012025; away_syst_PT[1] = 0.0440312; away_syst_PT[2] = 0.0457034; 
away_syst_PT[3] = 0.0255597; away_syst_PT[4] = 0.0250769; away_syst_PT[5] = 0.0249286; 
away_syst_PT[6] = 0.0252903; away_syst_PT[7] = 0.0252184; away_syst_PT[8] = 0.0826287; 
away_syst_PT[9] = 0.0731967; away_syst_PT[10] = 0.028363; 
// float mupt_syst_PT[11];
mupt_syst_PT[0] = 0.03302; mupt_syst_PT[1] = 0.00783314; mupt_syst_PT[2] = 0.00813064; 
mupt_syst_PT[3] = 0.0111587; mupt_syst_PT[4] = 0.0109479; mupt_syst_PT[5] = 0.0108832; 
mupt_syst_PT[6] = 0.0261789; mupt_syst_PT[7] = 0.0261045; mupt_syst_PT[8] = 0.0897278; 
mupt_syst_PT[9] = 0.0794854; mupt_syst_PT[10] = 0.05236; 
// float lc_syst_PT[11];
lc_syst_PT[0] = 0.00541496; lc_syst_PT[1] = 0.00499899; lc_syst_PT[2] = 0.00518884; 
lc_syst_PT[3] = 0.00678567; lc_syst_PT[4] = 0.0066575; lc_syst_PT[5] = 0.00661812; 
lc_syst_PT[6] = 0.00516583; lc_syst_PT[7] = 0.00515115; lc_syst_PT[8] = 0.0174728; 
lc_syst_PT[9] = 0.0154783; lc_syst_PT[10] = 7.98702e-06; 
// float gluon_syst_PT[11];
gluon_syst_PT[0] = 0.00733697; gluon_syst_PT[1] = 0.00017988; gluon_syst_PT[2] = 0.000186712; 
gluon_syst_PT[3] = 0.000158253; gluon_syst_PT[4] = 0.000155264; gluon_syst_PT[5] = 0.000154346; 
gluon_syst_PT[6] = 0.0176365; gluon_syst_PT[7] = 0.0175864; gluon_syst_PT[8] = 0.0628114; 
gluon_syst_PT[9] = 0.0556415; gluon_syst_PT[10] = 0.011198; 
// float pu_syst_PT[11];
pu_syst_PT[0] = 0.025039; pu_syst_PT[1] = 0.0165684; pu_syst_PT[2] = 0.0171977; 
pu_syst_PT[3] = 0.00850065; pu_syst_PT[4] = 0.00834008; pu_syst_PT[5] = 0.00829076; 
pu_syst_PT[6] = 0.0181993; pu_syst_PT[7] = 0.0181476; pu_syst_PT[8] = 0.0122839; 
pu_syst_PT[9] = 0.0108817; pu_syst_PT[10] = 0.024948; 
// float total_syst_PT[11];
total_syst_PT[0] = 0.0441024; total_syst_PT[1] = 0.0479545; total_syst_PT[2] = 0.0497758; 
total_syst_PT[3] = 0.0299357; total_syst_PT[4] = 0.0293703; total_syst_PT[5] = 0.0291966; 
total_syst_PT[6] = 0.0446529; total_syst_PT[7] = 0.0445259; total_syst_PT[8] = 0.138853; 
total_syst_PT[9] = 0.123003; total_syst_PT[10] = 0.0655273; 

// CSVT : S8
eff_SY[0]=0.2967; eff_SY[1]=0.4198; eff_SY[2]=0.4815; eff_SY[3]=0.5345; eff_SY[4]=0.5103; eff_SY[5]=0.5530; eff_SY[6]=0.5367; eff_SY[7]=0.4728; eff_SY[9]=0.4684; 
eff_stat_SY[0]=0.0081; eff_stat_SY[1]=0.0067; eff_stat_SY[2]=0.0087; eff_stat_SY[3]=0.0102; eff_stat_SY[4]=0.0106; eff_stat_SY[5]=0.0194; eff_stat_SY[6]=0.0307; eff_stat_SY[7]=0.0332; eff_stat_SY[9]=0.0734; 
effMC_SY[0]=0.2960; effMC_SY[1]=0.4362; effMC_SY[2]=0.5026; effMC_SY[3]=0.5567; effMC_SY[4]=0.5820; effMC_SY[5]=0.6020; effMC_SY[6]=0.5805; effMC_SY[7]=0.5826; effMC_SY[9]=0.4972; 
effMC_stat_SY[0]=0.0006; effMC_stat_SY[1]=0.0006; effMC_stat_SY[2]=0.0008; effMC_stat_SY[3]=0.0009; effMC_stat_SY[4]=0.0008; effMC_stat_SY[5]=0.0010; effMC_stat_SY[6]=0.0012; effMC_stat_SY[7]=0.0010; effMC_stat_SY[9]=0.0014; 
sf_SY[0]=1.0024; sf_SY[1]=0.9624; sf_SY[2]=0.958; sf_SY[3]=0.9601; sf_SY[4]=0.8768; sf_SY[5]=0.9186; sf_SY[6]=0.9245; sf_SY[7]=0.8115; sf_SY[9]=0.9421; 
sf_stat_SY[0]=0.0274; sf_stat_SY[1]=0.0154; sf_stat_SY[2]=0.0174; sf_stat_SY[3]=0.0184; sf_stat_SY[4]=0.0183; sf_stat_SY[5]=0.0323; sf_stat_SY[6]=0.0529; sf_stat_SY[7]=0.057; sf_stat_SY[9]=0.1477; 
away[0]=2.8957; away[1]=0.4637; away[2]=0.5621; away[3]=21.2151; 
mu[0]=0.3245; mu[1]=0.3671; mu[2]=1.5552; mu[3]=1.8863; 
ptrel[0]=0.0999; ptrel[1]=2.0672; ptrel[2]=0.7307; ptrel[3]=3.2092; 
gsplit[0]=0.4559; gsplit[1]=0.2633; gsplit[2]=0.8413; gsplit[3]=4.205; 
closure[0]=0.2399; closure[1]=0.0351; closure[2]=0.0841; closure[3]=0.6136; 
pu[0]=0.025; pu[1]=0.1159; pu[2]=1.2929; pu[3]=0.4655; 

// CSVT : IP3D
// float eff_IP[11];
eff_IP[0] = 0.120326; eff_IP[1] = 0.400068; eff_IP[2] = 0.433583; 
eff_IP[3] = 0.457076; eff_IP[4] = 0.509317; eff_IP[5] = 0.287579; 
eff_IP[6] = 0.500155; eff_IP[7] = 0.493164; eff_IP[8] = 0.502225; 
eff_IP[9] = 0.416074; eff_IP[10] = 0.386578; 
// float eff_stat_IP[11];
eff_stat_IP[0] = 0.110003; eff_stat_IP[1] = 0.0113511; eff_stat_IP[2] = 0.0774332; 
eff_stat_IP[3] = 0.0786543; eff_stat_IP[4] = 0.0184767; eff_stat_IP[5] = 0.0622929; 
eff_stat_IP[6] = 0.0119416; eff_stat_IP[7] = 0.0115082; eff_stat_IP[8] = 0.0123466; 
eff_stat_IP[9] = 0.010811; eff_stat_IP[10] = 0.0208486; 
// float effMC_IP[11];
effMC_IP[0] = 0.345556; effMC_IP[1] = 0.478383; effMC_IP[2] = 0.52154; 
effMC_IP[3] = 0.579733; effMC_IP[4] = 0.57685; effMC_IP[5] = 0.59179; 
effMC_IP[6] = 0.607594; effMC_IP[7] = 0.563259; effMC_IP[8] = 0.565174; 
effMC_IP[9] = 0.46044; effMC_IP[10] = 0.407388; 
// float effMC_stat_IP[11];
effMC_stat_IP[0] = 0.00319646; effMC_stat_IP[1] = 0.00228962; effMC_stat_IP[2] = 0.00204763; 
effMC_stat_IP[3] = 0.00266847; effMC_stat_IP[4] = 0.00259999; effMC_stat_IP[5] = 0.00270341; 
effMC_stat_IP[6] = 0.00290582; effMC_stat_IP[7] = 0.00281874; effMC_stat_IP[8] = 0.00371929; 
effMC_stat_IP[9] = 0.0039517; effMC_stat_IP[10] = 0.00586605; 
// float sf_IP[11];
sf_IP[0] = 0.34821; sf_IP[1] = 0.836292; sf_IP[2] = 0.831352; 
sf_IP[3] = 0.788426; sf_IP[4] = 0.882928; sf_IP[5] = 0.485947; 
sf_IP[6] = 0.823174; sf_IP[7] = 0.875554; sf_IP[8] = 0.888621; 
sf_IP[9] = 0.903644; sf_IP[10] = 0.948919; 
// float sf_stat_IP[11];
sf_stat_IP[0] = 0.318354; sf_stat_IP[1] = 0.0240633; sf_stat_IP[2] = 0.148506; 
sf_stat_IP[3] = 0.135722; sf_stat_IP[4] = 0.0322766; sf_stat_IP[5] = 0.105285; 
sf_stat_IP[6] = 0.0200443; sf_stat_IP[7] = 0.020896; sf_stat_IP[8] = 0.0226148; 
sf_stat_IP[9] = 0.0247274; sf_stat_IP[10] = 0.0529691; 
// float away_syst_IP[11];
away_syst_IP[0] = 0.0266122; away_syst_IP[1] = 0.00621578; away_syst_IP[2] = 0.00617906; 
away_syst_IP[3] = 0.0373427; away_syst_IP[4] = 0.0418186; away_syst_IP[5] = 0.0230162; 
away_syst_IP[6] = 0.0237643; away_syst_IP[7] = 0.0252764; away_syst_IP[8] = 0.0489617; 
away_syst_IP[9] = 0.0497894; away_syst_IP[10] = 0.0477126; 
// float mupt_syst_IP[11];
mupt_syst_IP[0] = 0.0258705; mupt_syst_IP[1] = 0.0248721; mupt_syst_IP[2] = 0.0247251; 
mupt_syst_IP[3] = 0.0168075; mupt_syst_IP[4] = 0.0188221; mupt_syst_IP[5] = 0.0103593; 
mupt_syst_IP[6] = 0.0225614; mupt_syst_IP[7] = 0.0239971; mupt_syst_IP[8] = 0.0264111; 
mupt_syst_IP[9] = 0.0268576; mupt_syst_IP[10] = 0.0305917; 
// float gluon_syst_IP[11];
gluon_syst_IP[0] = 0.00766162; gluon_syst_IP[1] = 0.00116208; gluon_syst_IP[2] = 0.00115521; 
gluon_syst_IP[3] = 0.0032688; gluon_syst_IP[4] = 0.0036606; gluon_syst_IP[5] = 0.00201473; 
gluon_syst_IP[6] = 0.00520815; gluon_syst_IP[7] = 0.00553955; gluon_syst_IP[8] = 0.0247843; 
gluon_syst_IP[9] = 0.0252033; gluon_syst_IP[10] = 0.12064; 
// float pu_syst_IP[11];
pu_syst_IP[0] = 0.01227; pu_syst_IP[1] = 0.0285192; pu_syst_IP[2] = 0.0283508; 
pu_syst_IP[3] = 0.00830641; pu_syst_IP[4] = 0.00930203; pu_syst_IP[5] = 0.00511966; 
pu_syst_IP[6] = 0.0144097; pu_syst_IP[7] = 0.0153266; pu_syst_IP[8] = 0.0278932; 
pu_syst_IP[9] = 0.0283648; pu_syst_IP[10] = 0.0298972; 
// float total_syst_IP[11];
total_syst_IP[0] = 0.039834; total_syst_IP[1] = 0.038366; total_syst_IP[2] = 0.0381394; 
total_syst_IP[3] = 0.0419124; total_syst_IP[4] = 0.0469361; total_syst_IP[5] = 0.0258327; 
total_syst_IP[6] = 0.0361735; total_syst_IP[7] = 0.0384753; total_syst_IP[8] = 0.0669857; 
total_syst_IP[9] = 0.0681182; total_syst_IP[10] = 0.136602; 

// CSVT : JP
eff_JP[0] = 0.330923; eff_stat_JP[0] = 0.00903226; effMC_JP[0] = 0.396586; effMC_stat_JP[0] = 0.00113324; 
sf_JP[0] = 0.834429; sf_stat_JP[0] = 0.0227751; sf_syst_JP[0] = 0.144356; sf_eror_JP[0] = 0.146142; 
mupt_syst_JP[0] = 0.015; gluon_syst_JP[0] = 0.011; pu_syst_JP[0] = 0.005; total_syst_JP[0] = 0.173; 
cor_syst_JP[0] = 0.044; inc_syst_JP[0] = 0.166; bias_syst_JP[0] = 0; 
eff_JP[1] = 0.462666; eff_stat_JP[1] = 0.00681508; effMC_JP[1] = 0.516596; effMC_stat_JP[1] = 0.000875479; 
sf_JP[1] = 0.895605; sf_stat_JP[1] = 0.0131923; sf_syst_JP[1] = 0.072544; sf_eror_JP[1] = 0.0737338; 
mupt_syst_JP[1] = 0.015; gluon_syst_JP[1] = 0.002; pu_syst_JP[1] = 0.002; total_syst_JP[1] = 0.081; 
cor_syst_JP[1] = 0.022; inc_syst_JP[1] = 0.076; bias_syst_JP[1] = 0; 
eff_JP[2] = 0.476358; eff_stat_JP[2] = 0.00567016; effMC_JP[2] = 0.535874; effMC_stat_JP[2] = 0.000939714; 
sf_JP[2] = 0.888937; sf_stat_JP[2] = 0.0105811; sf_syst_JP[2] = 0.070226; sf_eror_JP[2] = 0.0710187; 
mupt_syst_JP[2] = 0.015; gluon_syst_JP[2] = 0.001; pu_syst_JP[2] = 0.002; total_syst_JP[2] = 0.079; 
cor_syst_JP[2] = 0.013; inc_syst_JP[2] = 0.076; bias_syst_JP[2] = 0; 
eff_JP[3] = 0.546433; eff_stat_JP[3] = 0.00565471; effMC_JP[3] = 0.597322; effMC_stat_JP[3] = 0.00127905; 
sf_JP[3] = 0.914805; sf_stat_JP[3] = 0.00946677; sf_syst_JP[3] = 0.0411662; sf_eror_JP[3] = 0.0422407; 
mupt_syst_JP[3] = 0.015; gluon_syst_JP[3] = 0.007; pu_syst_JP[3] = 0.001; total_syst_JP[3] = 0.045; 
cor_syst_JP[3] = 0.008; inc_syst_JP[3] = 0.041; bias_syst_JP[3] = 0; 
eff_JP[4] = 0.538992; eff_stat_JP[4] = 0.00608876; effMC_JP[4] = 0.586818; effMC_stat_JP[4] = 0.00128562; 
sf_JP[4] = 0.9185; sf_stat_JP[4] = 0.0103759; sf_syst_JP[4] = 0.0413325; sf_eror_JP[4] = 0.042615; 
mupt_syst_JP[4] = 0.015; gluon_syst_JP[4] = 0.009; pu_syst_JP[4] = 0.002; total_syst_JP[4] = 0.045; 
cor_syst_JP[4] = 0.006; inc_syst_JP[4] = 0.041; bias_syst_JP[4] = 0; 
eff_JP[5] = 0.560479; eff_stat_JP[5] = 0.00606078; effMC_JP[5] = 0.603878; effMC_stat_JP[5] = 0.00133066; 
sf_JP[5] = 0.928133; sf_stat_JP[5] = 0.0100364; sf_syst_JP[5] = 0.0436223; sf_eror_JP[5] = 0.0447619; 
mupt_syst_JP[5] = 0.015; gluon_syst_JP[5] = 0.016; pu_syst_JP[5] = 0.001; total_syst_JP[5] = 0.047; 
cor_syst_JP[5] = 0.005; inc_syst_JP[5] = 0.041; bias_syst_JP[5] = 0; 
eff_JP[6] = 0.571535; eff_stat_JP[6] = 0.00495701; effMC_JP[6] = 0.604916; effMC_stat_JP[6] = 0.00127354; 
sf_JP[6] = 0.944817; sf_stat_JP[6] = 0.00819454; sf_syst_JP[6] = 0.0406271; sf_eror_JP[6] = 0.0414453; 
mupt_syst_JP[6] = 0.003; gluon_syst_JP[6] = 0.016; pu_syst_JP[6] = 0.002; total_syst_JP[6] = 0.043; 
cor_syst_JP[6] = 0.004; inc_syst_JP[6] = 0.04; bias_syst_JP[6] = 0; 
eff_JP[7] = 0.551973; eff_stat_JP[7] = 0.0060102; effMC_JP[7] = 0.587329; effMC_stat_JP[7] = 0.00128688; 
sf_JP[7] = 0.939802; sf_stat_JP[7] = 0.0102331; sf_syst_JP[7] = 0.0422911; sf_eror_JP[7] = 0.0435115; 
mupt_syst_JP[7] = 0.003; gluon_syst_JP[7] = 0.021; pu_syst_JP[7] = 0.002; total_syst_JP[7] = 0.045; 
cor_syst_JP[7] = 0.003; inc_syst_JP[7] = 0.04; bias_syst_JP[7] = 0; 
eff_JP[8] = 0.543588; eff_stat_JP[8] = 0.00630713; effMC_JP[8] = 0.574628; effMC_stat_JP[8] = 0.00124417; 
sf_JP[8] = 0.945982; sf_stat_JP[8] = 0.010976; sf_syst_JP[8] = 0.050137; sf_eror_JP[8] = 0.0513244; 
mupt_syst_JP[8] = 0.003; gluon_syst_JP[8] = 0.035; pu_syst_JP[8] = 0.001; total_syst_JP[8] = 0.053; 
cor_syst_JP[8] = 0.003; inc_syst_JP[8] = 0.04; bias_syst_JP[8] = 0; 
eff_JP[9] = 0.469686; eff_stat_JP[9] = 0.00542137; effMC_JP[9] = 0.497408; effMC_stat_JP[9] = 0.00153091; 
sf_JP[9] = 0.944268; sf_stat_JP[9] = 0.0108993; sf_syst_JP[9] = 0.0934825; sf_eror_JP[9] = 0.0941158; 
mupt_syst_JP[9] = 0.003; gluon_syst_JP[9] = 0.083; pu_syst_JP[9] = 0.004; total_syst_JP[9] = 0.099; 
cor_syst_JP[9] = 0.002; inc_syst_JP[9] = 0.054; bias_syst_JP[9] = 0; 
eff_JP[10] = 0.444498; eff_stat_JP[10] = 0.00804122; effMC_JP[10] = 0.482056; effMC_stat_JP[10] = 0.00272604; 
sf_JP[10] = 0.922086; sf_stat_JP[10] = 0.0166811; sf_syst_JP[10] = 0.0792994; sf_eror_JP[10] = 0.0810349; 
mupt_syst_JP[10] = 0.021; gluon_syst_JP[10] = 0.063; pu_syst_JP[10] = 0.004; total_syst_JP[10] = 0.086; 
cor_syst_JP[10] = 0.002; inc_syst_JP[10] = 0.054; bias_syst_JP[10] = 0; 
eff_JP[11] = 0.429797; eff_stat_JP[11] = 0.0118747; effMC_JP[11] = 0.457922; effMC_stat_JP[11] = 0.0045731; 
sf_JP[11] = 0.93858; sf_stat_JP[11] = 0.0259318; sf_syst_JP[11] = 0.0957352; sf_eror_JP[11] = 0.0991851; 
mupt_syst_JP[11] = 0.021; gluon_syst_JP[11] = 0.084; pu_syst_JP[11] = 0.004; total_syst_JP[11] = 0.102; 
cor_syst_JP[11] = 0.001; inc_syst_JP[11] = 0.054; bias_syst_JP[11] = 0; 
eff_JP[12] = 0.455308; eff_stat_JP[12] = 0.0186278; effMC_JP[12] = 0.447982; effMC_stat_JP[12] = 0.00775832; 
sf_JP[12] = 1.01635; sf_stat_JP[12] = 0.0415816; sf_syst_JP[12] = 0.0965533; sf_eror_JP[12] = 0.105126; 
mupt_syst_JP[12] = 0.021; gluon_syst_JP[12] = 0.078; pu_syst_JP[12] = 0.019; total_syst_JP[12] = 0.095; 
cor_syst_JP[12] = 0.001; inc_syst_JP[12] = 0.047; bias_syst_JP[12] = 0; 
eff_JP[13] = 0.379615; eff_stat_JP[13] = 0.0185139; effMC_JP[13] = 0.412353; effMC_stat_JP[13] = 0.00673927; 
sf_JP[13] = 0.920606; sf_stat_JP[13] = 0.0448981; sf_syst_JP[13] = 0.0644424; sf_eror_JP[13] = 0.0785408; 
mupt_syst_JP[13] = 0.021; gluon_syst_JP[13] = 0.047; pu_syst_JP[13] = 0.007; total_syst_JP[13] = 0.07; 
cor_syst_JP[13] = 0.001; inc_syst_JP[13] = 0.047; bias_syst_JP[13] = 0; 
eff_JP[14] = 0.340054; eff_stat_JP[14] = 0.0314298; effMC_JP[14] = 0.380068; effMC_stat_JP[14] = 0.00561375; 
sf_JP[14] = 0.894719; sf_stat_JP[14] = 0.0826952; sf_syst_JP[14] = 0.0474201; sf_eror_JP[14] = 0.0953266; 
mupt_syst_JP[14] = 0.021; gluon_syst_JP[14] = 0.009; pu_syst_JP[14] = 0.007; total_syst_JP[14] = 0.053; 
cor_syst_JP[14] = 0.008; inc_syst_JP[14] = 0.047; bias_syst_JP[14] = 0; 
}

// #############################################################################
if ( title == "JPL" ) { 
// JPL : PTREL
// float eff_PT[11];
eff_PT[0] = 0.623004; eff_PT[1] = 0.74407; eff_PT[2] = 0.791315; 
eff_PT[3] = 0.815284; eff_PT[4] = 0.836422; eff_PT[5] = 0.838155; 
eff_PT[6] = 0.866501; eff_PT[7] = 0.862228; eff_PT[8] = 0.889257; 
eff_PT[9] = 0.832531; eff_PT[10] = 0.70446; 
// float eff_stat_PT[11];
eff_stat_PT[0] = 0.0115511; eff_stat_PT[1] = 0.00763368; eff_stat_PT[2] = 0.00801531; 
eff_stat_PT[3] = 0.00817475; eff_stat_PT[4] = 0.0101481; eff_stat_PT[5] = 0.0126561; 
eff_stat_PT[6] = 0.0101212; eff_stat_PT[7] = 0.0111749; eff_stat_PT[8] = 0.0113732; 
eff_stat_PT[9] = 0.011287; eff_stat_PT[10] = 0.0218298; 
// float effMC_PT[11];
effMC_PT[0] = 0.687791; effMC_PT[1] = 0.768155; effMC_PT[2] = 0.811886; 
effMC_PT[3] = 0.842518; effMC_PT[4] = 0.857748; effMC_PT[5] = 0.871163; 
effMC_PT[6] = 0.880951; effMC_PT[7] = 0.879172; effMC_PT[8] = 0.884625; 
effMC_PT[9] = 0.867664; effMC_PT[10] = 0.844274; 
// float effMC_stat_PT[11];
effMC_stat_PT[0] = 0.00299238; effMC_stat_PT[1] = 0.00169486; effMC_stat_PT[2] = 0.00141197; 
effMC_stat_PT[3] = 0.00179386; effMC_stat_PT[4] = 0.00169796; effMC_stat_PT[5] = 0.00174135; 
effMC_stat_PT[6] = 0.00182184; effMC_stat_PT[7] = 0.00175578; effMC_stat_PT[8] = 0.00228018; 
effMC_stat_PT[9] = 0.0025355; effMC_stat_PT[10] = 0.00437515; 
// float sf_PT[11];
sf_PT[0] = 0.905804; sf_PT[1] = 0.968645; sf_PT[2] = 0.974662; 
sf_PT[3] = 0.967675; sf_PT[4] = 0.975137; sf_PT[5] = 0.96211; 
sf_PT[6] = 0.983598; sf_PT[7] = 0.980728; sf_PT[8] = 1.00524; 
sf_PT[9] = 0.959509; sf_PT[10] = 0.834398; 
// float sf_stat_PT[11];
sf_stat_PT[0] = 0.0172506; sf_stat_PT[1] = 0.0101649; sf_stat_PT[2] = 0.0100169; 
sf_stat_PT[3] = 0.0099191; sf_stat_PT[4] = 0.0119876; sf_stat_PT[5] = 0.0146545; 
sf_stat_PT[6] = 0.0116676; sf_stat_PT[7] = 0.0128608; sf_stat_PT[8] = 0.0131151; 
sf_stat_PT[9] = 0.0133073; sf_stat_PT[10] = 0.0262153; 
// float away_syst_PT[11];
away_syst_PT[0] = 0.008753; away_syst_PT[1] = 0.0261315; away_syst_PT[2] = 0.0262938; 
away_syst_PT[3] = 0.0137802; away_syst_PT[4] = 0.0138865; away_syst_PT[5] = 0.0137009; 
away_syst_PT[6] = 0.023709; away_syst_PT[7] = 0.0236398; away_syst_PT[8] = 0.0420969; 
away_syst_PT[9] = 0.0401818; away_syst_PT[10] = 0.074864; 
// float mupt_syst_PT[11];
mupt_syst_PT[0] = 0.027855; mupt_syst_PT[1] = 0.00638951; mupt_syst_PT[2] = 0.0064292; 
mupt_syst_PT[3] = 0.0162615; mupt_syst_PT[4] = 0.0163869; mupt_syst_PT[5] = 0.016168; 
mupt_syst_PT[6] = 0.0236047; mupt_syst_PT[7] = 0.0235358; mupt_syst_PT[8] = 0.041987; 
mupt_syst_PT[9] = 0.0400769; mupt_syst_PT[10] = 0.012367; 
// float lc_syst_PT[11];
lc_syst_PT[0] = 0.005045; lc_syst_PT[1] = 0.00499811; lc_syst_PT[2] = 0.00502916; 
lc_syst_PT[3] = 0.00530097; lc_syst_PT[4] = 0.00534184; lc_syst_PT[5] = 0.00527048; 
lc_syst_PT[6] = 0.00171551; lc_syst_PT[7] = 0.0017105; lc_syst_PT[8] = 0.0155914; 
lc_syst_PT[9] = 0.0148821; lc_syst_PT[10] = 0.022383; 
// float gluon_syst_PT[11];
gluon_syst_PT[0] = 0.00193101; gluon_syst_PT[1] = 0.000803898; gluon_syst_PT[2] = 0.000808891; 
gluon_syst_PT[3] = 0.000740015; gluon_syst_PT[4] = 0.000745722; gluon_syst_PT[5] = 0.00073576; 
gluon_syst_PT[6] = 0.0041801; gluon_syst_PT[7] = 0.0041679; gluon_syst_PT[8] = 0.0125222; 
gluon_syst_PT[9] = 0.0119525; gluon_syst_PT[10] = 0.038577; 
// float pu_syst_PT[11];
pu_syst_PT[0] = 0.012915; pu_syst_PT[1] = 0.0108171; pu_syst_PT[2] = 0.0108843; 
pu_syst_PT[3] = 0.00916511; pu_syst_PT[4] = 0.00923579; pu_syst_PT[5] = 0.0091124; 
pu_syst_PT[6] = 0.00179075; pu_syst_PT[7] = 0.00178553; pu_syst_PT[8] = 0.00325823; 
pu_syst_PT[9] = 0.00311001; pu_syst_PT[10] = 0.018577; 
// float total_syst_PT[11];
total_syst_PT[0] = 0.0323805; total_syst_PT[1] = 0.0294333; total_syst_PT[2] = 0.0296161; 
total_syst_PT[3] = 0.0238113; total_syst_PT[4] = 0.0239949; total_syst_PT[5] = 0.0236743; 
total_syst_PT[6] = 0.0338071; total_syst_PT[7] = 0.0337085; total_syst_PT[8] = 0.0628138; 
total_syst_PT[9] = 0.0599563; total_syst_PT[10] = 0.0899548; 

// JPL : S8
eff_SY[0]=0.6503; eff_SY[1]=0.7228; eff_SY[2]=0.7855; eff_SY[3]=0.8172; eff_SY[4]=0.8056; eff_SY[5]=0.8874; eff_SY[6]=0.8546; eff_SY[7]=0.8356; eff_SY[9]=0.8324; 
eff_stat_SY[0]=0.0107; eff_stat_SY[1]=0.0074; eff_stat_SY[2]=0.0101; eff_stat_SY[3]=0.0105; eff_stat_SY[4]=0.0121; eff_stat_SY[5]=0.0270; eff_stat_SY[6]=0.0350; eff_stat_SY[7]=0.0459; eff_stat_SY[9]=0.0729; 
effMC_SY[0]=0.6295; effMC_SY[1]=0.7253; effMC_SY[2]=0.7894; effMC_SY[3]=0.8292; effMC_SY[4]=0.8561; effMC_SY[5]=0.8749; effMC_SY[6]=0.8793; effMC_SY[7]=0.8813; effMC_SY[9]=0.8785; 
effMC_stat_SY[0]=0.0007; effMC_stat_SY[1]=0.0006; effMC_stat_SY[2]=0.0007; effMC_stat_SY[3]=0.0007; effMC_stat_SY[4]=0.0006; effMC_stat_SY[5]=0.0007; effMC_stat_SY[6]=0.0008; effMC_stat_SY[7]=0.0007; effMC_stat_SY[9]=0.0009; 
sf_SY[0]=1.033; sf_SY[1]=0.9966; sf_SY[2]=0.9951; sf_SY[3]=0.9855; sf_SY[4]=0.9410; sf_SY[5]=1.0143; sf_SY[6]=0.9719; sf_SY[7]=0.9481; sf_SY[9]=0.9475; 
sf_stat_SY[0]=0.017; sf_stat_SY[1]=0.0102; sf_stat_SY[2]=0.0128; sf_stat_SY[3]=0.0127; sf_stat_SY[4]=0.0141; sf_stat_SY[5]=0.0309; sf_stat_SY[6]=0.0398; sf_stat_SY[7]=0.0521; sf_stat_SY[9]=0.083; 
away[0]=0.7684; away[1]=0.5339; away[2]=1.9747; away[3]=11.4817; 
mu[0]=0.0419; mu[1]=0.0621; mu[2]=0.2323; mu[3]=0.8852; 
ptrel[0]=0.1676; ptrel[1]=1.1547; ptrel[2]=2.9504; ptrel[3]=1.1674; 
gsplit[0]=0.0559; gsplit[1]=0.1424; gsplit[2]=0.0571; gsplit[3]=4.6467; 
closure[0]=0.2932; closure[1]=0.0831; closure[2]=0.0228; closure[3]=0.4862; 
pu[0]=0.0559; pu[1]=0.2235; pu[2]=1.8818; pu[3]=0.6671; 

// JPL : IP3D
// float eff_IP[11];
eff_IP[0] = 0.589981; eff_IP[1] = 0.718782; eff_IP[2] = 0.786051; 
eff_IP[3] = 0.752099; eff_IP[4] = 0.841184; eff_IP[5] = 0.800259; 
eff_IP[6] = 0.838722; eff_IP[7] = 0.831645; eff_IP[8] = 0.861657; 
eff_IP[9] = 0.856039; eff_IP[10] = 0.870792; 
// float eff_stat_IP[11];
eff_stat_IP[0] = 0.0512184; eff_stat_IP[1] = 0.0164999; eff_stat_IP[2] = 0.0154453; 
eff_stat_IP[3] = 0.0135656; eff_stat_IP[4] = 0.015136; eff_stat_IP[5] = 0.018249; 
eff_stat_IP[6] = 0.00879592; eff_stat_IP[7] = 0.00999546; eff_stat_IP[8] = 0.00840028; 
eff_stat_IP[9] = 0.00978528; eff_stat_IP[10] = 0.0172941; 
// float effMC_IP[11];
effMC_IP[0] = 0.691655; effMC_IP[1] = 0.769404; effMC_IP[2] = 0.811234; 
effMC_IP[3] = 0.845908; effMC_IP[4] = 0.856259; effMC_IP[5] = 0.867602; 
effMC_IP[6] = 0.878612; effMC_IP[7] = 0.875684; effMC_IP[8] = 0.880822; 
effMC_IP[9] = 0.865817; effMC_IP[10] = 0.841262; 
// float effMC_stat_IP[11];
effMC_stat_IP[0] = 0.00298121; effMC_stat_IP[1] = 0.00181579; effMC_stat_IP[2] = 0.00149186; 
effMC_stat_IP[3] = 0.00178419; effMC_stat_IP[4] = 0.00167422; effMC_stat_IP[5] = 0.00171625; 
effMC_stat_IP[6] = 0.00177446; effMC_stat_IP[7] = 0.00170549; effMC_stat_IP[8] = 0.00220783; 
effMC_stat_IP[9] = 0.00244198; effMC_stat_IP[10] = 0.00437724; 
// float sf_IP[11];
sf_IP[0] = 0.852999; sf_IP[1] = 0.934207; sf_IP[2] = 0.968958; 
sf_IP[3] = 0.889102; sf_IP[4] = 0.982394; sf_IP[5] = 0.92238; 
sf_IP[6] = 0.9546; sf_IP[7] = 0.949708; sf_IP[8] = 0.978241; 
sf_IP[9] = 0.988706; sf_IP[10] = 1.0351; 
// float sf_stat_IP[11];
sf_stat_IP[0] = 0.0741432; sf_stat_IP[1] = 0.0215581; sf_stat_IP[2] = 0.0191224; 
sf_stat_IP[3] = 0.016146; sf_stat_IP[4] = 0.0177809; sf_stat_IP[5] = 0.0211129; 
sf_stat_IP[6] = 0.0101951; sf_stat_IP[7] = 0.0115633; sf_stat_IP[8] = 0.00984703; 
sf_stat_IP[9] = 0.0116407; sf_stat_IP[10] = 0.0212512; 
// float away_syst_IP[11];
away_syst_IP[0] = 0.0848424; away_syst_IP[1] = 0.0447933; away_syst_IP[2] = 0.0464595; 
away_syst_IP[3] = 0.0211747; away_syst_IP[4] = 0.0233965; away_syst_IP[5] = 0.0219672; 
away_syst_IP[6] = 0.0121732; away_syst_IP[7] = 0.0121108; away_syst_IP[8] = 0.0101289; 
away_syst_IP[9] = 0.0102372; away_syst_IP[10] = 0.0218734; 
// float mupt_syst_IP[11];
mupt_syst_IP[0] = 0.0756571; mupt_syst_IP[1] = 0.0211931; mupt_syst_IP[2] = 0.0219815; 
mupt_syst_IP[3] = 0.0264245; mupt_syst_IP[4] = 0.0291972; mupt_syst_IP[5] = 0.0274135; 
mupt_syst_IP[6] = 0.0113098; mupt_syst_IP[7] = 0.0112518; mupt_syst_IP[8] = 0.00685593; 
mupt_syst_IP[9] = 0.00692927; mupt_syst_IP[10] = 0.00623106; 
// float gluon_syst_IP[11];
gluon_syst_IP[0] = 0.050109; gluon_syst_IP[1] = 0.00364101; gluon_syst_IP[2] = 0.00377645; 
gluon_syst_IP[3] = 0.0016978; gluon_syst_IP[4] = 0.00187595; gluon_syst_IP[5] = 0.00176135; 
gluon_syst_IP[6] = 0.000340319; gluon_syst_IP[7] = 0.000338575; gluon_syst_IP[8] = 0.00350995; 
gluon_syst_IP[9] = 0.0035475; gluon_syst_IP[10] = 0.0132597; 
// float pu_syst_IP[11];
pu_syst_IP[0] = 0.0612626; pu_syst_IP[1] = 0.00364592; pu_syst_IP[2] = 0.00378155; 
pu_syst_IP[3] = 0.00930651; pu_syst_IP[4] = 0.010283; pu_syst_IP[5] = 0.00965485; 
pu_syst_IP[6] = 0.00331618; pu_syst_IP[7] = 0.00329919; pu_syst_IP[8] = 0.001559; 
pu_syst_IP[9] = 0.00157567; pu_syst_IP[10] = 0.0170581; 
// float total_syst_IP[11];
total_syst_IP[0] = 0.138514; total_syst_IP[1] = 0.049821; total_syst_IP[2] = 0.0516743; 
total_syst_IP[3] = 0.0351584; total_syst_IP[4] = 0.0388475; total_syst_IP[5] = 0.0364744; 
total_syst_IP[6] = 0.0169473; total_syst_IP[7] = 0.0168604; total_syst_IP[8] = 0.0128198; 
total_syst_IP[9] = 0.012957; total_syst_IP[10] = 0.0313699; 

// JPL : JP
eff_JP[0] = 0.707891; eff_stat_JP[0] = 0.0123095; effMC_JP[0] = 0.707979; effMC_stat_JP[0] = 0.00105658; 
sf_JP[0] = 0.999874; sf_stat_JP[0] = 0.0173868; sf_syst_JP[0] = 0.123984; sf_eror_JP[0] = 0.125198; 
mupt_syst_JP[0] = 0.004; gluon_syst_JP[0] = 0.001; pu_syst_JP[0] = 0; total_syst_JP[0] = 0.124; 
cor_syst_JP[0] = 0.028; inc_syst_JP[0] = 0.119; bias_syst_JP[0] = 0.02; 
eff_JP[1] = 0.771977; eff_stat_JP[1] = 0.00823036; effMC_JP[1] = 0.779101; effMC_stat_JP[1] = 0.000727871; 
sf_JP[1] = 0.990858; sf_stat_JP[1] = 0.0105639; sf_syst_JP[1] = 0.0653967; sf_eror_JP[1] = 0.0662444; 
mupt_syst_JP[1] = 0.004; gluon_syst_JP[1] = 0.001; pu_syst_JP[1] = 0.001; total_syst_JP[1] = 0.066; 
cor_syst_JP[1] = 0.013; inc_syst_JP[1] = 0.061; bias_syst_JP[1] = 0.02; 
eff_JP[2] = 0.81533; eff_stat_JP[2] = 0.0050997; effMC_JP[2] = 0.817726; effMC_stat_JP[2] = 0.000727765; 
sf_JP[2] = 0.997072; sf_stat_JP[2] = 0.00623646; sf_syst_JP[2] = 0.0648097; sf_eror_JP[2] = 0.065109; 
mupt_syst_JP[2] = 0.004; gluon_syst_JP[2] = 0.001; pu_syst_JP[2] = 0; total_syst_JP[2] = 0.065; 
cor_syst_JP[2] = 0.007; inc_syst_JP[2] = 0.061; bias_syst_JP[2] = 0.02; 
eff_JP[3] = 0.852896; eff_stat_JP[3] = 0.00637538; effMC_JP[3] = 0.857264; effMC_stat_JP[3] = 0.000911459; 
sf_JP[3] = 0.994906; sf_stat_JP[3] = 0.0074369; sf_syst_JP[3] = 0.063674; sf_eror_JP[3] = 0.0641068; 
mupt_syst_JP[3] = 0.004; gluon_syst_JP[3] = 0; pu_syst_JP[3] = 0; total_syst_JP[3] = 0.064; 
cor_syst_JP[3] = 0.004; inc_syst_JP[3] = 0.061; bias_syst_JP[3] = 0.02; 
eff_JP[4] = 0.862954; eff_stat_JP[4] = 0.00606041; effMC_JP[4] = 0.866289; effMC_stat_JP[4] = 0.000890359; 
sf_JP[4] = 0.99615; sf_stat_JP[4] = 0.00699582; sf_syst_JP[4] = 0.0637536; sf_eror_JP[4] = 0.0641363; 
mupt_syst_JP[4] = 0.004; gluon_syst_JP[4] = 0.001; pu_syst_JP[4] = 0; total_syst_JP[4] = 0.064; 
cor_syst_JP[4] = 0.002; inc_syst_JP[4] = 0.061; bias_syst_JP[4] = 0.02; 
eff_JP[5] = 0.863522; eff_stat_JP[5] = 0.0208593; effMC_JP[5] = 0.870284; effMC_stat_JP[5] = 0.000912866; 
sf_JP[5] = 0.99223; sf_stat_JP[5] = 0.0239684; sf_syst_JP[5] = 0.0635028; sf_eror_JP[5] = 0.0678755; 
mupt_syst_JP[5] = 0.004; gluon_syst_JP[5] = 0.001; pu_syst_JP[5] = 0; total_syst_JP[5] = 0.064; 
cor_syst_JP[5] = 0.002; inc_syst_JP[5] = 0.061; bias_syst_JP[5] = 0.02; 
eff_JP[6] = 0.890167; eff_stat_JP[6] = 0.00910462; effMC_JP[6] = 0.885659; effMC_stat_JP[6] = 0.000828831; 
sf_JP[6] = 1.00509; sf_stat_JP[6] = 0.0102801; sf_syst_JP[6] = 0.0432188; sf_eror_JP[6] = 0.0444246; 
mupt_syst_JP[6] = 0.005; gluon_syst_JP[6] = 0.001; pu_syst_JP[6] = 0; total_syst_JP[6] = 0.043; 
cor_syst_JP[6] = 0.001; inc_syst_JP[6] = 0.042; bias_syst_JP[6] = 0.01; 
eff_JP[7] = 0.885543; eff_stat_JP[7] = 0.00583888; effMC_JP[7] = 0.883417; effMC_stat_JP[7] = 0.000837991; 
sf_JP[7] = 1.0024; sf_stat_JP[7] = 0.00660943; sf_syst_JP[7] = 0.0441058; sf_eror_JP[7] = 0.0445983; 
mupt_syst_JP[7] = 0.005; gluon_syst_JP[7] = 0.002; pu_syst_JP[7] = 0; total_syst_JP[7] = 0.044; 
cor_syst_JP[7] = 0.001; inc_syst_JP[7] = 0.042; bias_syst_JP[7] = 0.01; 
eff_JP[8] = 0.886085; eff_stat_JP[8] = 0.00703088; effMC_JP[8] = 0.885334; effMC_stat_JP[8] = 0.000801234; 
sf_JP[8] = 1.00085; sf_stat_JP[8] = 0.0079415; sf_syst_JP[8] = 0.0440374; sf_eror_JP[8] = 0.0447478; 
mupt_syst_JP[8] = 0.005; gluon_syst_JP[8] = 0.003; pu_syst_JP[8] = 0; total_syst_JP[8] = 0.044; 
cor_syst_JP[8] = 0.001; inc_syst_JP[8] = 0.042; bias_syst_JP[8] = 0.01; 
eff_JP[9] = 0.877209; eff_stat_JP[9] = 0.0102112; effMC_JP[9] = 0.879117; effMC_stat_JP[9] = 0.00099797; 
sf_JP[9] = 0.997831; sf_stat_JP[9] = 0.0116153; sf_syst_JP[9] = 0.0548807; sf_eror_JP[9] = 0.0560964; 
mupt_syst_JP[9] = 0.005; gluon_syst_JP[9] = 0.004; pu_syst_JP[9] = 0; total_syst_JP[9] = 0.055; 
cor_syst_JP[9] = 0.001; inc_syst_JP[9] = 0.054; bias_syst_JP[9] = 0.01; 
eff_JP[10] = 0.876891; eff_stat_JP[10] = 0.0122444; effMC_JP[10] = 0.87036; effMC_stat_JP[10] = 0.00183257; 
sf_JP[10] = 1.0075; sf_stat_JP[10] = 0.0140683; sf_syst_JP[10] = 0.0554125; sf_eror_JP[10] = 0.0571705; 
mupt_syst_JP[10] = 0; gluon_syst_JP[10] = 0.008; pu_syst_JP[10] = 0; total_syst_JP[10] = 0.055; 
cor_syst_JP[10] = 0; inc_syst_JP[10] = 0.054; bias_syst_JP[10] = 0.01; 
eff_JP[11] = 0.879428; eff_stat_JP[11] = 0.0140546; effMC_JP[11] = 0.866434; effMC_stat_JP[11] = 0.00312239; 
sf_JP[11] = 1.015; sf_stat_JP[11] = 0.0162212; sf_syst_JP[11] = 0.05684; sf_eror_JP[11] = 0.0591093; 
mupt_syst_JP[11] = 0; gluon_syst_JP[11] = 0.01; pu_syst_JP[11] = 0.002; total_syst_JP[11] = 0.056; 
cor_syst_JP[11] = 0; inc_syst_JP[11] = 0.054; bias_syst_JP[11] = 0.01; 
eff_JP[12] = 0.848048; eff_stat_JP[12] = 0.0192939; effMC_JP[12] = 0.85857; effMC_stat_JP[12] = 0.00532296; 
sf_JP[12] = 0.987745; sf_stat_JP[12] = 0.0224721; sf_syst_JP[12] = 0.0671667; sf_eror_JP[12] = 0.0708262; 
mupt_syst_JP[12] = 0; gluon_syst_JP[12] = 0.013; pu_syst_JP[12] = 0.004; total_syst_JP[12] = 0.068; 
cor_syst_JP[12] = 0; inc_syst_JP[12] = 0.066; bias_syst_JP[12] = 0.01; 
eff_JP[13] = 0.857379; eff_stat_JP[13] = 0.0223778; effMC_JP[13] = 0.845517; effMC_stat_JP[13] = 0.00476504; 
sf_JP[13] = 1.01403; sf_stat_JP[13] = 0.0264664; sf_syst_JP[13] = 0.06794; sf_eror_JP[13] = 0.0729131; 
mupt_syst_JP[13] = 0; gluon_syst_JP[13] = 0.001; pu_syst_JP[13] = 0.008; total_syst_JP[13] = 0.067; 
cor_syst_JP[13] = 0; inc_syst_JP[13] = 0.066; bias_syst_JP[13] = 0.01; 
eff_JP[14] = 0.803299; eff_stat_JP[14] = 0.0287345; effMC_JP[14] = 0.815536; effMC_stat_JP[14] = 0.00485472; 
sf_JP[14] = 0.984994; sf_stat_JP[14] = 0.0352339; sf_syst_JP[14] = 0.0659946; sf_eror_JP[14] = 0.0748112; 
mupt_syst_JP[14] = 0; gluon_syst_JP[14] = 0.002; pu_syst_JP[14] = 0.005; total_syst_JP[14] = 0.067; 
cor_syst_JP[14] = 0.007; inc_syst_JP[14] = 0.066; bias_syst_JP[14] = 0.01; 
}

// #############################################################################
if ( title == "JPM" ) { 
// JPM : PTREL
// float eff_PT[11];
eff_PT[0] = 0.359386; eff_PT[1] = 0.465029; eff_PT[2] = 0.533829; 
eff_PT[3] = 0.567208; eff_PT[4] = 0.591411; eff_PT[5] = 0.606508; 
eff_PT[6] = 0.630404; eff_PT[7] = 0.620887; eff_PT[8] = 0.618303; 
eff_PT[9] = 0.563094; eff_PT[10] = 0.509885; 
// float eff_stat_PT[11];
eff_stat_PT[0] = 0.00980683; eff_stat_PT[1] = 0.00752499; eff_stat_PT[2] = 0.00852566; 
eff_stat_PT[3] = 0.00879049; eff_stat_PT[4] = 0.0107238; eff_stat_PT[5] = 0.0137759; 
eff_stat_PT[6] = 0.0106183; eff_stat_PT[7] = 0.0115145; eff_stat_PT[8] = 0.0133963; 
eff_stat_PT[9] = 0.0113811; eff_stat_PT[10] = 0.0221644; 
// float effMC_PT[11];
effMC_PT[0] = 0.405123; effMC_PT[1] = 0.509378; effMC_PT[2] = 0.578592; 
effMC_PT[3] = 0.637827; effMC_PT[4] = 0.657585; effMC_PT[5] = 0.677864; 
effMC_PT[6] = 0.704275; effMC_PT[7] = 0.686766; effMC_PT[8] = 0.698757; 
effMC_PT[9] = 0.675192; effMC_PT[10] = 0.621738; 
// float effMC_stat_PT[11];
effMC_stat_PT[0] = 0.00334186; effMC_stat_PT[1] = 0.00213701; effMC_stat_PT[2] = 0.00191117; 
effMC_stat_PT[3] = 0.00258042; effMC_stat_PT[4] = 0.00250174; effMC_stat_PT[5] = 0.00258988; 
effMC_stat_PT[6] = 0.00276148; effMC_stat_PT[7] = 0.00267711; effMC_stat_PT[8] = 0.00350254; 
effMC_stat_PT[9] = 0.00380284; effMC_stat_PT[10] = 0.00595454; 
// float sf_PT[11];
sf_PT[0] = 0.887105; sf_PT[1] = 0.912934; sf_PT[2] = 0.922635; 
sf_PT[3] = 0.889282; sf_PT[4] = 0.899368; sf_PT[5] = 0.894734; 
sf_PT[6] = 0.89511; sf_PT[7] = 0.904073; sf_PT[8] = 0.884862; 
sf_PT[9] = 0.833977; sf_PT[10] = 0.820097; 
// float sf_stat_PT[11];
sf_stat_PT[0] = 0.025289; sf_stat_PT[1] = 0.0152613; sf_stat_PT[2] = 0.015047; 
sf_stat_PT[3] = 0.0142438; sf_stat_PT[4] = 0.0166629; sf_stat_PT[5] = 0.020608; 
sf_stat_PT[6] = 0.01548; sf_stat_PT[7] = 0.0171327; sf_stat_PT[8] = 0.019678; 
sf_stat_PT[9] = 0.0174984; sf_stat_PT[10] = 0.0365042; 
// float away_syst_PT[11];
away_syst_PT[0] = 0.027883; away_syst_PT[1] = 0.0426228; away_syst_PT[2] = 0.0430757; 
away_syst_PT[3] = 0.020344; away_syst_PT[4] = 0.0205748; away_syst_PT[5] = 0.0204688; 
away_syst_PT[6] = 0.0170476; away_syst_PT[7] = 0.0172183; away_syst_PT[8] = 0.0791308; 
away_syst_PT[9] = 0.0745802; away_syst_PT[10] = 0.113708; 
// float mupt_syst_PT[11];
mupt_syst_PT[0] = 0.012772; mupt_syst_PT[1] = 0.00626271; mupt_syst_PT[2] = 0.00632926; 
mupt_syst_PT[3] = 0.0139745; mupt_syst_PT[4] = 0.014133; mupt_syst_PT[5] = 0.0140601; 
mupt_syst_PT[6] = 0.00956658; mupt_syst_PT[7] = 0.00966238; mupt_syst_PT[8] = 0.0618539; 
mupt_syst_PT[9] = 0.0582969; mupt_syst_PT[10] = 0.038637; 
// float lc_syst_PT[11];
lc_syst_PT[0] = 0.00586498; lc_syst_PT[1] = 0.00498267; lc_syst_PT[2] = 0.00503561; 
lc_syst_PT[3] = 0.00615427; lc_syst_PT[4] = 0.00622407; lc_syst_PT[5] = 0.006192; 
lc_syst_PT[6] = 0.00345575; lc_syst_PT[7] = 0.00349036; lc_syst_PT[8] = 0.0203962; 
lc_syst_PT[9] = 0.0192233; lc_syst_PT[10] = 0.0513099; 
// float gluon_syst_PT[11];
gluon_syst_PT[0] = 0.00266898; gluon_syst_PT[1] = 0.000223961; gluon_syst_PT[2] = 0.000226341; 
gluon_syst_PT[3] = 0.000275131; gluon_syst_PT[4] = 0.000278252; gluon_syst_PT[5] = 0.000276818; 
gluon_syst_PT[6] = 0.0141762; gluon_syst_PT[7] = 0.0143182; gluon_syst_PT[8] = 0.0375659; 
gluon_syst_PT[9] = 0.0354056; gluon_syst_PT[10] = 0.0331711; 
// float pu_syst_PT[11];
pu_syst_PT[0] = 0.02313; pu_syst_PT[1] = 0.0138194; pu_syst_PT[2] = 0.0139663; 
pu_syst_PT[3] = 0.00746028; pu_syst_PT[4] = 0.00754489; pu_syst_PT[5] = 0.00750602; 
pu_syst_PT[6] = 0.0145076; pu_syst_PT[7] = 0.0146529; pu_syst_PT[8] = 0.00814101; 
pu_syst_PT[9] = 0.00767285; pu_syst_PT[10] = 0.023144; 
// float total_syst_PT[11];
total_syst_PT[0] = 0.03895; total_syst_PT[1] = 0.0455168; total_syst_PT[2] = 0.0460004; 
total_syst_PT[3] = 0.0265098; total_syst_PT[4] = 0.0268105; total_syst_PT[5] = 0.0266724; 
total_syst_PT[6] = 0.0283817; total_syst_PT[7] = 0.0286658; total_syst_PT[8] = 0.109458; 
total_syst_PT[9] = 0.103164; total_syst_PT[10] = 0.136715; 

// JPM : S8
eff_SY[0]=0.3433; eff_SY[1]=0.4476; eff_SY[2]=0.5183; eff_SY[3]=0.5920; eff_SY[4]=0.5648; eff_SY[5]=0.6468; eff_SY[6]=0.6373; eff_SY[7]=0.5729; eff_SY[9]=0.6104; 
eff_stat_SY[0]=0.0087; eff_stat_SY[1]=0.0071; eff_stat_SY[2]=0.0091; eff_stat_SY[3]=0.0116; eff_stat_SY[4]=0.0124; eff_stat_SY[5]=0.0236; eff_stat_SY[6]=0.0355; eff_stat_SY[7]=0.0410; eff_stat_SY[9]=0.0763; 
effMC_SY[0]=0.3463; effMC_SY[1]=0.4656; effMC_SY[2]=0.5547; effMC_SY[3]=0.6169; effMC_SY[4]=0.6510; effMC_SY[5]=0.6917; effMC_SY[6]=0.6997; effMC_SY[7]=0.7023; effMC_SY[9]=0.6929; 
effMC_stat_SY[0]=0.0006; effMC_stat_SY[1]=0.0007; effMC_stat_SY[2]=0.0008; effMC_stat_SY[3]=0.0009; effMC_stat_SY[4]=0.0010; effMC_stat_SY[5]=0.0009; effMC_stat_SY[6]=0.0011; effMC_stat_SY[7]=0.0010; effMC_stat_SY[9]=0.0012; 
sf_SY[0]=0.9913; sf_SY[1]=0.9613; sf_SY[2]=0.9344; sf_SY[3]=0.9596; sf_SY[4]=0.8676; sf_SY[5]=0.9351; sf_SY[6]=0.9108; sf_SY[7]=0.8157; sf_SY[9]=0.8809; 
sf_stat_SY[0]=0.0252; sf_stat_SY[1]=0.0153; sf_stat_SY[2]=0.0165; sf_stat_SY[3]=0.0189; sf_stat_SY[4]=0.0191; sf_stat_SY[5]=0.0341; sf_stat_SY[6]=0.0508; sf_stat_SY[7]=0.0584; sf_stat_SY[9]=0.1101; 
away[0]=2.9851; away[1]=0.1376; away[2]=0.8077; away[3]=12.0313; 
mu[0]=0.4133; mu[1]=0.2751; mu[2]=0.3168; mu[3]=1.2311; 
ptrel[0]=1.4007; ptrel[1]=2.7682; ptrel[2]=1.5996; ptrel[3]=2.7607; 
gsplit[0]=0.2399; gsplit[1]=0.2662; gsplit[2]=0.5606; gsplit[3]=4.8474; 
closure[0]=0.2181; closure[1]=0.047; closure[2]=0.1294; closure[3]=0.3707; 
pu[0]=0.1607; pu[1]=0.2063; pu[2]=0.8711; pu[3]=0.6715; 

// JPM : IP3D
// float eff_IP[11];
eff_IP[0] = 0.311417; eff_IP[1] = 0.380209; eff_IP[2] = 0.479005; 
eff_IP[3] = 0.50646; eff_IP[4] = 0.605927; eff_IP[5] = 0.630877; 
eff_IP[6] = 0.633477; eff_IP[7] = 0.621776; eff_IP[8] = 0.644164; 
eff_IP[9] = 0.621978; eff_IP[10] = 0.608524; 
// float eff_stat_IP[11];
eff_stat_IP[0] = 0.0702541; eff_stat_IP[1] = 0.0231153; eff_stat_IP[2] = 0.0813767; 
eff_stat_IP[3] = 0.0826222; eff_stat_IP[4] = 0.016869; eff_stat_IP[5] = 0.0194923; 
eff_stat_IP[6] = 0.0102988; eff_stat_IP[7] = 0.0124926; eff_stat_IP[8] = 0.0106864; 
eff_stat_IP[9] = 0.0120071; eff_stat_IP[10] = 0.0225995; 
// float effMC_IP[11];
effMC_IP[0] = 0.404179; effMC_IP[1] = 0.510247; effMC_IP[2] = 0.577151; 
effMC_IP[3] = 0.640747; effMC_IP[4] = 0.656815; effMC_IP[5] = 0.678769; 
effMC_IP[6] = 0.700908; effMC_IP[7] = 0.685144; effMC_IP[8] = 0.69507; 
effMC_IP[9] = 0.668773; effMC_IP[10] = 0.61219; 
// float effMC_stat_IP[11];
effMC_stat_IP[0] = 0.00332357; effMC_stat_IP[1] = 0.00228643; effMC_stat_IP[2] = 0.0020108; 
effMC_stat_IP[3] = 0.00256301; effMC_stat_IP[4] = 0.00246496; effMC_stat_IP[5] = 0.00254452; 
effMC_stat_IP[6] = 0.00268923; effMC_stat_IP[7] = 0.002597; effMC_stat_IP[8] = 0.00338947; 
effMC_stat_IP[9] = 0.00366722; effMC_stat_IP[10] = 0.00595363; 
// float sf_IP[11];
sf_IP[0] = 0.770493; sf_IP[1] = 0.745147; sf_IP[2] = 0.829948; 
sf_IP[3] = 0.790421; sf_IP[4] = 0.922523; sf_IP[5] = 0.929443; 
sf_IP[6] = 0.903794; sf_IP[7] = 0.907511; sf_IP[8] = 0.926762; 
sf_IP[9] = 0.930029; sf_IP[10] = 0.994012; 
// float sf_stat_IP[11];
sf_stat_IP[0] = 0.173935; sf_stat_IP[1] = 0.0454251; sf_stat_IP[2] = 0.141027; 
sf_stat_IP[3] = 0.128985; sf_stat_IP[4] = 0.0259153; sf_stat_IP[5] = 0.0289277; 
sf_stat_IP[6] = 0.0150972; sf_stat_IP[7] = 0.0185552; sf_stat_IP[8] = 0.016025; 
sf_stat_IP[9] = 0.0186642; sf_stat_IP[10] = 0.0381606; 
// float away_syst_IP[11];
away_syst_IP[0] = 0.212314; away_syst_IP[1] = 0.0243688; away_syst_IP[2] = 0.027142; 
away_syst_IP[3] = 0.00785621; away_syst_IP[4] = 0.0091692; away_syst_IP[5] = 0.00923798; 
away_syst_IP[6] = 0.0274535; away_syst_IP[7] = 0.0275664; away_syst_IP[8] = 0.02927; 
away_syst_IP[9] = 0.0293731; away_syst_IP[10] = 0.0138647; 
// float mupt_syst_IP[11];
mupt_syst_IP[0] = 0.0947668; mupt_syst_IP[1] = 0.0235424; mupt_syst_IP[2] = 0.0262216; 
mupt_syst_IP[3] = 0.0433777; mupt_syst_IP[4] = 0.0506274; mupt_syst_IP[5] = 0.0510071; 
mupt_syst_IP[6] = 0.00433743; mupt_syst_IP[7] = 0.00435527; mupt_syst_IP[8] = 0.000865509; 
mupt_syst_IP[9] = 0.00086856; mupt_syst_IP[10] = 0.0167754; 
// float gluon_syst_IP[11];
gluon_syst_IP[0] = 0.00533425; gluon_syst_IP[1] = 0.00236409; gluon_syst_IP[2] = 0.00263313; 
gluon_syst_IP[3] = 0.03187; gluon_syst_IP[4] = 0.0371963; gluon_syst_IP[5] = 0.0374754; 
gluon_syst_IP[6] = 0.00244272; gluon_syst_IP[7] = 0.00245277; gluon_syst_IP[8] = 0.00611104; 
gluon_syst_IP[9] = 0.00613258; gluon_syst_IP[10] = 0.0620372; 
// float pu_syst_IP[11];
pu_syst_IP[0] = 0.0814644; pu_syst_IP[1] = 0.00577193; pu_syst_IP[2] = 0.0064288; 
pu_syst_IP[3] = 0.0328044; pu_syst_IP[4] = 0.0382869; pu_syst_IP[5] = 0.0385741; 
pu_syst_IP[6] = 0.00992641; pu_syst_IP[7] = 0.00996724; pu_syst_IP[8] = 0.00114936; 
pu_syst_IP[9] = 0.00115341; pu_syst_IP[10] = 0.012139; 
// float total_syst_IP[11];
total_syst_IP[0] = 0.24642; total_syst_IP[1] = 0.0344526; total_syst_IP[2] = 0.0383735; 
total_syst_IP[3] = 0.063523; total_syst_IP[4] = 0.0741395; total_syst_IP[5] = 0.0746956; 
total_syst_IP[6] = 0.0296143; total_syst_IP[7] = 0.0297361; total_syst_IP[8] = 0.0299357; 
total_syst_IP[9] = 0.0300412; total_syst_IP[10] = 0.0668552; 

// JPM : JP
eff_JP[0] = 0.442165; eff_stat_JP[0] = 0.0165721; effMC_JP[0] = 0.448076; effMC_stat_JP[0] = 0.00115559; 
sf_JP[0] = 0.986803; sf_stat_JP[0] = 0.0369851; sf_syst_JP[0] = 0.238806; sf_eror_JP[0] = 0.241653; 
mupt_syst_JP[0] = 0.015; gluon_syst_JP[0] = 0; pu_syst_JP[0] = 0; total_syst_JP[0] = 0.242; 
cor_syst_JP[0] = 0.028; inc_syst_JP[0] = 0.232; bias_syst_JP[0] = 0.06; 
eff_JP[1] = 0.53095; eff_stat_JP[1] = 0.0138112; effMC_JP[1] = 0.544231; effMC_stat_JP[1] = 0.000873825; 
sf_JP[1] = 0.975598; sf_stat_JP[1] = 0.0253776; sf_syst_JP[1] = 0.0829258; sf_eror_JP[1] = 0.086722; 
mupt_syst_JP[1] = 0.015; gluon_syst_JP[1] = 0.003; pu_syst_JP[1] = 0.001; total_syst_JP[1] = 0.085; 
cor_syst_JP[1] = 0.013; inc_syst_JP[1] = 0.057; bias_syst_JP[1] = 0.06; 
eff_JP[2] = 0.586153; eff_stat_JP[2] = 0.0151536; effMC_JP[2] = 0.602056; effMC_stat_JP[2] = 0.000922686; 
sf_JP[2] = 0.973586; sf_stat_JP[2] = 0.0251696; sf_syst_JP[2] = 0.0817812; sf_eror_JP[2] = 0.0855668; 
mupt_syst_JP[2] = 0.015; gluon_syst_JP[2] = 0.003; pu_syst_JP[2] = 0.001; total_syst_JP[2] = 0.084; 
cor_syst_JP[2] = 0.007; inc_syst_JP[2] = 0.057; bias_syst_JP[2] = 0.06; 
eff_JP[3] = 0.64293; eff_stat_JP[3] = 0.0157556; effMC_JP[3] = 0.659442; effMC_stat_JP[3] = 0.0012348; 
sf_JP[3] = 0.974959; sf_stat_JP[3] = 0.0238924; sf_syst_JP[3] = 0.0750718; sf_eror_JP[3] = 0.0787821; 
mupt_syst_JP[3] = 0.015; gluon_syst_JP[3] = 0.004; pu_syst_JP[3] = 0.001; total_syst_JP[3] = 0.077; 
cor_syst_JP[3] = 0.004; inc_syst_JP[3] = 0.046; bias_syst_JP[3] = 0.06; 
eff_JP[4] = 0.65807; eff_stat_JP[4] = 0.0156281; effMC_JP[4] = 0.67451; effMC_stat_JP[4] = 0.00122578; 
sf_JP[4] = 0.975626; sf_stat_JP[4] = 0.0231695; sf_syst_JP[4] = 0.0751232; sf_eror_JP[4] = 0.078615; 
mupt_syst_JP[4] = 0.015; gluon_syst_JP[4] = 0.004; pu_syst_JP[4] = 0; total_syst_JP[4] = 0.077; 
cor_syst_JP[4] = 0.002; inc_syst_JP[4] = 0.046; bias_syst_JP[4] = 0.06; 
eff_JP[5] = 0.664011; eff_stat_JP[5] = 0.0137116; effMC_JP[5] = 0.685679; effMC_stat_JP[5] = 0.00126133; 
sf_JP[5] = 0.968397; sf_stat_JP[5] = 0.0199971; sf_syst_JP[5] = 0.0745666; sf_eror_JP[5] = 0.0772014; 
mupt_syst_JP[5] = 0.015; gluon_syst_JP[5] = 0.005; pu_syst_JP[5] = 0; total_syst_JP[5] = 0.077; 
cor_syst_JP[5] = 0.002; inc_syst_JP[5] = 0.046; bias_syst_JP[5] = 0.06; 
eff_JP[6] = 0.705411; eff_stat_JP[6] = 0.0158979; effMC_JP[6] = 0.709291; effMC_stat_JP[6] = 0.0011827; 
sf_JP[6] = 0.994531; sf_stat_JP[6] = 0.0224138; sf_syst_JP[6] = 0.0546992; sf_eror_JP[6] = 0.0591133; 
mupt_syst_JP[6] = 0.008; gluon_syst_JP[6] = 0.005; pu_syst_JP[6] = 0.001; total_syst_JP[6] = 0.055; 
cor_syst_JP[6] = 0.001; inc_syst_JP[6] = 0.037; bias_syst_JP[6] = 0.04; 
eff_JP[7] = 0.697647; eff_stat_JP[7] = 0.0140772; effMC_JP[7] = 0.706105; effMC_stat_JP[7] = 0.00118951; 
sf_JP[7] = 0.988022; sf_stat_JP[7] = 0.0199363; sf_syst_JP[7] = 0.0543412; sf_eror_JP[7] = 0.0578829; 
mupt_syst_JP[7] = 0.008; gluon_syst_JP[7] = 0.004; pu_syst_JP[7] = 0; total_syst_JP[7] = 0.055; 
cor_syst_JP[7] = 0.001; inc_syst_JP[7] = 0.037; bias_syst_JP[7] = 0.04; 
eff_JP[8] = 0.702727; eff_stat_JP[8] = 0.0140557; effMC_JP[8] = 0.711045; effMC_stat_JP[8] = 0.00113986; 
sf_JP[8] = 0.988301; sf_stat_JP[8] = 0.0197677; sf_syst_JP[8] = 0.0543565; sf_eror_JP[8] = 0.0578394; 
mupt_syst_JP[8] = 0.008; gluon_syst_JP[8] = 0.003; pu_syst_JP[8] = 0; total_syst_JP[8] = 0.055; 
cor_syst_JP[8] = 0.001; inc_syst_JP[8] = 0.037; bias_syst_JP[8] = 0.04; 
eff_JP[9] = 0.685134; eff_stat_JP[9] = 0.0145318; effMC_JP[9] = 0.694155; effMC_stat_JP[9] = 0.00141056; 
sf_JP[9] = 0.987005; sf_stat_JP[9] = 0.0209345; sf_syst_JP[9] = 0.0592203; sf_eror_JP[9] = 0.0628116; 
mupt_syst_JP[9] = 0.008; gluon_syst_JP[9] = 0.003; pu_syst_JP[9] = 0.001; total_syst_JP[9] = 0.06; 
cor_syst_JP[9] = 0.001; inc_syst_JP[9] = 0.044; bias_syst_JP[9] = 0.04; 
eff_JP[10] = 0.661152; eff_stat_JP[10] = 0.0168804; effMC_JP[10] = 0.675449; effMC_stat_JP[10] = 0.00255434; 
sf_JP[10] = 0.978835; sf_stat_JP[10] = 0.0249914; sf_syst_JP[10] = 0.0587301; sf_eror_JP[10] = 0.0638263; 
mupt_syst_JP[10] = 0.003; gluon_syst_JP[10] = 0.01; pu_syst_JP[10] = 0.002; total_syst_JP[10] = 0.06; 
cor_syst_JP[10] = 0; inc_syst_JP[10] = 0.044; bias_syst_JP[10] = 0.04; 
eff_JP[11] = 0.655055; eff_stat_JP[11] = 0.0199992; effMC_JP[11] = 0.661773; effMC_stat_JP[11] = 0.0043424; 
sf_JP[11] = 0.989846; sf_stat_JP[11] = 0.0302206; sf_syst_JP[11] = 0.0593908; sf_eror_JP[11] = 0.0666375; 
mupt_syst_JP[11] = 0.003; gluon_syst_JP[11] = 0.01; pu_syst_JP[11] = 0.001; total_syst_JP[11] = 0.06; 
cor_syst_JP[11] = 0; inc_syst_JP[11] = 0.044; bias_syst_JP[11] = 0.04; 
eff_JP[12] = 0.634592; eff_stat_JP[12] = 0.026733; effMC_JP[12] = 0.652621; effMC_stat_JP[12] = 0.00727321; 
sf_JP[12] = 0.972374; sf_stat_JP[12] = 0.0409624; sf_syst_JP[12] = 0.0661215; sf_eror_JP[12] = 0.0777815; 
mupt_syst_JP[12] = 0.003; gluon_syst_JP[12] = 0.002; pu_syst_JP[12] = 0; total_syst_JP[12] = 0.068; 
cor_syst_JP[12] = 0; inc_syst_JP[12] = 0.055; bias_syst_JP[12] = 0.04; 
eff_JP[13] = 0.59319; eff_stat_JP[13] = 0.0323508; effMC_JP[13] = 0.616988; effMC_stat_JP[13] = 0.00640929; 
sf_JP[13] = 0.96143; sf_stat_JP[13] = 0.0524336; sf_syst_JP[13] = 0.0673001; sf_eror_JP[13] = 0.0853146; 
mupt_syst_JP[13] = 0.003; gluon_syst_JP[13] = 0.003; pu_syst_JP[13] = 0.018; total_syst_JP[13] = 0.07; 
cor_syst_JP[13] = 0; inc_syst_JP[13] = 0.055; bias_syst_JP[13] = 0.04; 
eff_JP[14] = 0.583035; eff_stat_JP[14] = 0.0444833; effMC_JP[14] = 0.57958; effMC_stat_JP[14] = 0.00617854; 
sf_JP[14] = 1.00596; sf_stat_JP[14] = 0.076751; sf_syst_JP[14] = 0.0714235; sf_eror_JP[14] = 0.104843; 
mupt_syst_JP[14] = 0.003; gluon_syst_JP[14] = 0.018; pu_syst_JP[14] = 0.002; total_syst_JP[14] = 0.071; 
cor_syst_JP[14] = 0.007; inc_syst_JP[14] = 0.055; bias_syst_JP[14] = 0.04; 
}

// #############################################################################
if ( title == "JPT" ) { 
// JPT : PTREL
// float eff_PT[11];
eff_PT[0] = 0.143378; eff_PT[1] = 0.235408; eff_PT[2] = 0.313288; 
eff_PT[3] = 0.346904; eff_PT[4] = 0.372452; eff_PT[5] = 0.385939; 
eff_PT[6] = 0.40974; eff_PT[7] = 0.394103; eff_PT[8] = 0.389766; 
eff_PT[9] = 0.337588; eff_PT[10] = 0.322982; 
// float eff_stat_PT[11];
eff_stat_PT[0] = 0.00636435; eff_stat_PT[1] = 0.00575911; eff_stat_PT[2] = 0.00658724; 
eff_stat_PT[3] = 0.00761191; eff_stat_PT[4] = 0.00911606; eff_stat_PT[5] = 0.0107991; 
eff_stat_PT[6] = 0.00787072; eff_stat_PT[7] = 0.00911358; eff_stat_PT[8] = 0.0104803; 
eff_stat_PT[9] = 0.00847121; eff_stat_PT[10] = 0.0142244; 
// float effMC_PT[11];
effMC_PT[0] = 0.178837; effMC_PT[1] = 0.280125; effMC_PT[2] = 0.348681; 
effMC_PT[3] = 0.416489; effMC_PT[4] = 0.443961; effMC_PT[5] = 0.472705; 
effMC_PT[6] = 0.493821; effMC_PT[7] = 0.482667; effMC_PT[8] = 0.494998; 
effMC_PT[9] = 0.4715; effMC_PT[10] = 0.427451; 
// float effMC_stat_PT[11];
effMC_stat_PT[0] = 0.00254297; effMC_stat_PT[1] = 0.00190625; effMC_stat_PT[2] = 0.0018503; 
effMC_stat_PT[3] = 0.00267633; effMC_stat_PT[4] = 0.00265938; effMC_stat_PT[5] = 0.00278645; 
effMC_stat_PT[6] = 0.00303876; effMC_stat_PT[7] = 0.0029229; effMC_stat_PT[8] = 0.00388524; 
effMC_stat_PT[9] = 0.00409852; effMC_stat_PT[10] = 0.00588951; 
// float sf_PT[11];
sf_PT[0] = 0.801723; sf_PT[1] = 0.840367; sf_PT[2] = 0.898495; 
sf_PT[3] = 0.832925; sf_PT[4] = 0.83893; sf_PT[5] = 0.816448; 
sf_PT[6] = 0.829733; sf_PT[7] = 0.816511; sf_PT[8] = 0.787409; 
sf_PT[9] = 0.715988; sf_PT[10] = 0.755599; 
// float sf_stat_PT[11];
sf_stat_PT[0] = 0.0373688; sf_stat_PT[1] = 0.0213396; sf_stat_PT[2] = 0.0194843; 
sf_stat_PT[3] = 0.019044; sf_stat_PT[4] = 0.0211395; sf_stat_PT[5] = 0.0233467; 
sf_stat_PT[6] = 0.0167363; sf_stat_PT[7] = 0.0195184; sf_stat_PT[8] = 0.022056; 
sf_stat_PT[9] = 0.019014; sf_stat_PT[10] = 0.0348677; 
// float away_syst_PT[11];
away_syst_PT[0] = 0.015053; away_syst_PT[1] = 0.0582126; away_syst_PT[2] = 0.0622391; 
away_syst_PT[3] = 0.0172324; away_syst_PT[4] = 0.0173567; away_syst_PT[5] = 0.0168915; 
away_syst_PT[6] = 0.0358913; away_syst_PT[7] = 0.0353193; away_syst_PT[8] = 0.0695693; 
away_syst_PT[9] = 0.0632591; away_syst_PT[10] = 0.114568; 
// float mupt_syst_PT[11];
mupt_syst_PT[0] = 0.025016; mupt_syst_PT[1] = 0.0100613; mupt_syst_PT[2] = 0.0107572; 
mupt_syst_PT[3] = 0.018847; mupt_syst_PT[4] = 0.0189828; mupt_syst_PT[5] = 0.0184741; 
mupt_syst_PT[6] = 0.0322146; mupt_syst_PT[7] = 0.0317012; mupt_syst_PT[8] = 0.0790492; 
mupt_syst_PT[9] = 0.0718791; mupt_syst_PT[10] = 0.107235; 
// float lc_syst_PT[11];
lc_syst_PT[0] = 0.00497001; lc_syst_PT[1] = 0.00447694; lc_syst_PT[2] = 0.0047866; 
lc_syst_PT[3] = 0.00603699; lc_syst_PT[4] = 0.00608051; lc_syst_PT[5] = 0.00591757; 
lc_syst_PT[6] = 0.00454269; lc_syst_PT[7] = 0.0044703; lc_syst_PT[8] = 0.0172608; 
lc_syst_PT[9] = 0.0156952; lc_syst_PT[10] = 0.056356; 
// float gluon_syst_PT[11];
gluon_syst_PT[0] = 0.000333011; gluon_syst_PT[1] = 0.00170795; gluon_syst_PT[2] = 0.00182609; 
gluon_syst_PT[3] = 0.00216641; gluon_syst_PT[4] = 0.00218203; gluon_syst_PT[5] = 0.00212355; 
gluon_syst_PT[6] = 0.0161203; gluon_syst_PT[7] = 0.0158635; gluon_syst_PT[8] = 0.0452712; 
gluon_syst_PT[9] = 0.0411649; gluon_syst_PT[10] = 0.027246; 
// float pu_syst_PT[11];
pu_syst_PT[0] = 0.026725; pu_syst_PT[1] = 0.0186911; pu_syst_PT[2] = 0.019984; 
pu_syst_PT[3] = 0.0141208; pu_syst_PT[4] = 0.0142226; pu_syst_PT[5] = 0.0138414; 
pu_syst_PT[6] = 0.00422121; pu_syst_PT[7] = 0.00415395; pu_syst_PT[8] = 0.00760333; 
pu_syst_PT[9] = 0.00691368; pu_syst_PT[10] = 0.044388; 
// float total_syst_PT[11];
total_syst_PT[0] = 0.0398927; total_syst_PT[1] = 0.062147; total_syst_PT[2] = 0.0664457; 
total_syst_PT[3] = 0.0298781; total_syst_PT[4] = 0.0300935; total_syst_PT[5] = 0.029287; 
total_syst_PT[6] = 0.0512278; total_syst_PT[7] = 0.0504114; total_syst_PT[8] = 0.116163; 
total_syst_PT[9] = 0.105627; total_syst_PT[10] = 0.174682; 

// JPT : S8
eff_SY[0]=0.1348; eff_SY[1]=0.2233; eff_SY[2]=0.3044; eff_SY[3]=0.3604; eff_SY[4]=0.3595; eff_SY[5]=0.4138; eff_SY[6]=0.4291; eff_SY[7]=0.3697; 
eff_stat_SY[0]=0.0060; eff_stat_SY[1]=0.0057; eff_stat_SY[2]=0.0077; eff_stat_SY[3]=0.0086; eff_stat_SY[4]=0.0092; eff_stat_SY[5]=0.0177; eff_stat_SY[6]=0.0252; eff_stat_SY[7]=0.0334; 
effMC_SY[0]=0.1472; effMC_SY[1]=0.2489; effMC_SY[2]=0.3346; effMC_SY[3]=0.3980; effMC_SY[4]=0.4474; effMC_SY[5]=0.4875; effMC_SY[6]=0.4991; effMC_SY[7]=0.5051; 
effMC_stat_SY[0]=0.0004; effMC_stat_SY[1]=0.0005; effMC_stat_SY[2]=0.0007; effMC_stat_SY[3]=0.0009; effMC_stat_SY[4]=0.0008; effMC_stat_SY[5]=0.0010; effMC_stat_SY[6]=0.0012; effMC_stat_SY[7]=0.0011; 
sf_SY[0]=0.9158; sf_SY[1]=0.8971; sf_SY[2]=0.9097; sf_SY[3]=0.9055; sf_SY[4]=0.8035; sf_SY[5]=0.8488; sf_SY[6]=0.8597; sf_SY[7]=0.7319; 
sf_stat_SY[0]=0.0408; sf_stat_SY[1]=0.023; sf_stat_SY[2]=0.0231; sf_stat_SY[3]=0.0217; sf_stat_SY[4]=0.0206; sf_stat_SY[5]=0.0363; sf_stat_SY[6]=0.0505; sf_stat_SY[7]=0.0661; 
away[0]=4.6691; away[1]=0.8128; away[2]=1.0613; away[3]=12.109; 
mu[0]=0.136; mu[1]=1.3173; mu[2]=0.5065; mu[3]=2.1036; 
ptrel[0]=0.9973; ptrel[1]=0.7287; ptrel[2]=2.6773; ptrel[3]=6.5264; 
gsplit[0]=0.3656; gsplit[1]=0.1888; gsplit[2]=0.9137; gsplit[3]=5.9524; 
closure[0]=0.1625; closure[1]=0.0236; closure[2]=0.1218; closure[3]=1.9516; 
pu[0]=0.0907; pu[1]=0.3083; pu[2]=1.2783; pu[3]=0.2967; 

// JPT : IP3D
// float eff_IP[11];
eff_IP[0] = 0.137104; eff_IP[1] = 0.217602; eff_IP[2] = 0.283445; 
eff_IP[3] = 0.0187397; eff_IP[4] = 0.369722; eff_IP[5] = 0.39519; 
eff_IP[6] = 0.392183; eff_IP[7] = 0.426907; eff_IP[8] = 0.433683; 
eff_IP[9] = 0.415515; eff_IP[10] = 0.394583; 
// float eff_stat_IP[11];
eff_stat_IP[0] = 0.00617352; eff_stat_IP[1] = 0.0529977; eff_stat_IP[2] = 0.06334; 
eff_stat_IP[3] = 0.164028; eff_stat_IP[4] = 0.0112156; eff_stat_IP[5] = 0.075677; 
eff_stat_IP[6] = 0.00952218; eff_stat_IP[7] = 0.00959186; eff_stat_IP[8] = 0.00972168; 
eff_stat_IP[9] = 0.0109914; eff_stat_IP[10] = 0.0197978; 
// float effMC_IP[11];
effMC_IP[0] = 0.181394; effMC_IP[1] = 0.282163; effMC_IP[2] = 0.351522; 
effMC_IP[3] = 0.420792; effMC_IP[4] = 0.444959; effMC_IP[5] = 0.470917; 
effMC_IP[6] = 0.489625; effMC_IP[7] = 0.48326; effMC_IP[8] = 0.492378; 
effMC_IP[9] = 0.464578; effMC_IP[10] = 0.417994; 
// float effMC_stat_IP[11];
effMC_stat_IP[0] = 0.00252407; effMC_stat_IP[1] = 0.00203604; effMC_stat_IP[2] = 0.00194171; 
effMC_stat_IP[3] = 0.00265552; effMC_stat_IP[4] = 0.00261451; effMC_stat_IP[5] = 0.00273134; 
effMC_stat_IP[6] = 0.00295885; effMC_stat_IP[7] = 0.00283735; effMC_stat_IP[8] = 0.00376213; 
effMC_stat_IP[9] = 0.0039545; effMC_stat_IP[10] = 0.0058885; 
// float sf_IP[11];
sf_IP[0] = 0.755835; sf_IP[1] = 0.771193; sf_IP[2] = 0.806336; 
sf_IP[3] = 0.0445344; sf_IP[4] = 0.830911; sf_IP[5] = 0.839192; 
sf_IP[6] = 0.800986; sf_IP[7] = 0.88339; sf_IP[8] = 0.880793; 
sf_IP[9] = 0.894391; sf_IP[10] = 0.943994; 
// float sf_stat_IP[11];
sf_stat_IP[0] = 0.0356219; sf_stat_IP[1] = 0.187909; sf_stat_IP[2] = 0.180243; 
sf_stat_IP[3] = 0.389808; sf_stat_IP[4] = 0.0256743; sf_stat_IP[5] = 0.160775; 
sf_stat_IP[6] = 0.0200412; sf_stat_IP[7] = 0.0205147; sf_stat_IP[8] = 0.0208598; 
sf_stat_IP[9] = 0.0248537; sf_stat_IP[10] = 0.0491954; 
// float away_syst_IP[11];
away_syst_IP[0] = 0.0425091; away_syst_IP[1] = 0.0272262; away_syst_IP[2] = 0.0284669; 
away_syst_IP[3] = 0.000400076; away_syst_IP[4] = 0.00746452; away_syst_IP[5] = 0.00753891; 
away_syst_IP[6] = 0.0161113; away_syst_IP[7] = 0.0177688; away_syst_IP[8] = 0.0373499; 
away_syst_IP[9] = 0.0379266; away_syst_IP[10] = 0.111587; 
// float mupt_syst_IP[11];
mupt_syst_IP[0] = 0.0487599; mupt_syst_IP[1] = 0.00177632; mupt_syst_IP[2] = 0.00185727; 
mupt_syst_IP[3] = 0.000159282; mupt_syst_IP[4] = 0.00297185; mupt_syst_IP[5] = 0.00300147; 
mupt_syst_IP[6] = 0.0146773; mupt_syst_IP[7] = 0.0161873; mupt_syst_IP[8] = 0.00548109; 
mupt_syst_IP[9] = 0.00556571; mupt_syst_IP[10] = 0.0343767; 
// float gluon_syst_IP[11];
gluon_syst_IP[0] = 0.164109; gluon_syst_IP[1] = 0.00424586; gluon_syst_IP[2] = 0.00443934; 
gluon_syst_IP[3] = 0.000282079; gluon_syst_IP[4] = 0.00526295; gluon_syst_IP[5] = 0.0053154; 
gluon_syst_IP[6] = 0.00682301; gluon_syst_IP[7] = 0.00752495; gluon_syst_IP[8] = 0.00236658; 
gluon_syst_IP[9] = 0.00240312; gluon_syst_IP[10] = 0.108061; 
// float pu_syst_IP[11];
pu_syst_IP[0] = 0.026078; pu_syst_IP[1] = 0.00702166; pu_syst_IP[2] = 0.00734164; 
pu_syst_IP[3] = 0.000402028; pu_syst_IP[4] = 0.00750094; pu_syst_IP[5] = 0.00757569; 
pu_syst_IP[6] = 0.00796224; pu_syst_IP[7] = 0.00878138; pu_syst_IP[8] = 0.0185095; 
pu_syst_IP[9] = 0.0187952; pu_syst_IP[10] = 0.0244211; 
// float total_syst_IP[11];
total_syst_IP[0] = 0.178315; total_syst_IP[1] = 0.0284913; total_syst_IP[2] = 0.0297896; 
total_syst_IP[3] = 0.000653167; total_syst_IP[4] = 0.0121866; total_syst_IP[5] = 0.0123081; 
total_syst_IP[6] = 0.0241857; total_syst_IP[7] = 0.0266739; total_syst_IP[8] = 0.0421101; 
total_syst_IP[9] = 0.0427602; total_syst_IP[10] = 0.160956; 

// JPT : JP
eff_JP[0] = 0.201749; eff_stat_JP[0] = 0.00558609; effMC_JP[0] = 0.221567; effMC_stat_JP[0] = 0.000965051; 
sf_JP[0] = 0.910555; sf_stat_JP[0] = 0.0252117; sf_syst_JP[0] = 0.200322; sf_eror_JP[0] = 0.201902; 
mupt_syst_JP[0] = 0.018; gluon_syst_JP[0] = 0; pu_syst_JP[0] = 0.001; total_syst_JP[0] = 0.22; 
cor_syst_JP[0] = 0.028; inc_syst_JP[0] = 0.193; bias_syst_JP[0] = 0.1; 
eff_JP[1] = 0.287674; eff_stat_JP[1] = 0.00571061; effMC_JP[1] = 0.3182; effMC_stat_JP[1] = 0.00081722; 
sf_JP[1] = 0.904067; sf_stat_JP[1] = 0.0179467; sf_syst_JP[1] = 0.112104; sf_eror_JP[1] = 0.113532; 
mupt_syst_JP[1] = 0.018; gluon_syst_JP[1] = 0.004; pu_syst_JP[1] = 0; total_syst_JP[1] = 0.124; 
cor_syst_JP[1] = 0.013; inc_syst_JP[1] = 0.069; bias_syst_JP[1] = 0.1; 
eff_JP[2] = 0.349201; eff_stat_JP[2] = 0.00670989; effMC_JP[2] = 0.38262; effMC_stat_JP[2] = 0.000916189; 
sf_JP[2] = 0.912657; sf_stat_JP[2] = 0.0175367; sf_syst_JP[2] = 0.112257; sf_eror_JP[2] = 0.113618; 
mupt_syst_JP[2] = 0.018; gluon_syst_JP[2] = 0.005; pu_syst_JP[2] = 0; total_syst_JP[2] = 0.123; 
cor_syst_JP[2] = 0.007; inc_syst_JP[2] = 0.069; bias_syst_JP[2] = 0.1; 
eff_JP[3] = 0.408399; eff_stat_JP[3] = 0.00758946; effMC_JP[3] = 0.441922; effMC_stat_JP[3] = 0.001294; 
sf_JP[3] = 0.924141; sf_stat_JP[3] = 0.0171738; sf_syst_JP[3] = 0.1072; sf_eror_JP[3] = 0.108567; 
mupt_syst_JP[3] = 0.018; gluon_syst_JP[3] = 0.006; pu_syst_JP[3] = 0; total_syst_JP[3] = 0.116; 
cor_syst_JP[3] = 0.004; inc_syst_JP[3] = 0.055; bias_syst_JP[3] = 0.1; 
eff_JP[4] = 0.431745; eff_stat_JP[4] = 0.00781069; effMC_JP[4] = 0.464415; effMC_stat_JP[4] = 0.00130472; 
sf_JP[4] = 0.929656; sf_stat_JP[4] = 0.0168183; sf_syst_JP[4] = 0.10784; sf_eror_JP[4] = 0.109144; 
mupt_syst_JP[4] = 0.018; gluon_syst_JP[4] = 0.007; pu_syst_JP[4] = 0.001; total_syst_JP[4] = 0.116; 
cor_syst_JP[4] = 0.002; inc_syst_JP[4] = 0.055; bias_syst_JP[4] = 0.1; 
eff_JP[5] = 0.442844; eff_stat_JP[5] = 0.00708612; effMC_JP[5] = 0.479687; effMC_stat_JP[5] = 0.00135735; 
sf_JP[5] = 0.923195; sf_stat_JP[5] = 0.0147724; sf_syst_JP[5] = 0.107091; sf_eror_JP[5] = 0.108105; 
mupt_syst_JP[5] = 0.018; gluon_syst_JP[5] = 0.009; pu_syst_JP[5] = 0; total_syst_JP[5] = 0.116; 
cor_syst_JP[5] = 0.002; inc_syst_JP[5] = 0.055; bias_syst_JP[5] = 0.1; 
eff_JP[6] = 0.487008; eff_stat_JP[6] = 0.00851098; effMC_JP[6] = 0.5074; effMC_stat_JP[6] = 0.00130213; 
sf_JP[6] = 0.959815; sf_stat_JP[6] = 0.0167737; sf_syst_JP[6] = 0.0710263; sf_eror_JP[6] = 0.0729801; 
mupt_syst_JP[6] = 0.012; gluon_syst_JP[6] = 0.009; pu_syst_JP[6] = 0; total_syst_JP[6] = 0.074; 
cor_syst_JP[6] = 0.001; inc_syst_JP[6] = 0.04; bias_syst_JP[6] = 0.06; 
eff_JP[7] = 0.485079; eff_stat_JP[7] = 0.00746136; effMC_JP[7] = 0.50636; effMC_stat_JP[7] = 0.00130549; 
sf_JP[7] = 0.957973; sf_stat_JP[7] = 0.0147353; sf_syst_JP[7] = 0.069932; sf_eror_JP[7] = 0.0714676; 
mupt_syst_JP[7] = 0.012; gluon_syst_JP[7] = 0.007; pu_syst_JP[7] = 0.002; total_syst_JP[7] = 0.073; 
cor_syst_JP[7] = 0.001; inc_syst_JP[7] = 0.04; bias_syst_JP[7] = 0.06; 
eff_JP[8] = 0.494893; eff_stat_JP[8] = 0.00743646; effMC_JP[8] = 0.51648; effMC_stat_JP[8] = 0.00125667; 
sf_JP[8] = 0.958208; sf_stat_JP[8] = 0.0143984; sf_syst_JP[8] = 0.0699492; sf_eror_JP[8] = 0.0714157; 
mupt_syst_JP[8] = 0.012; gluon_syst_JP[8] = 0.001; pu_syst_JP[8] = 0.001; total_syst_JP[8] = 0.073; 
cor_syst_JP[8] = 0.001; inc_syst_JP[8] = 0.04; bias_syst_JP[8] = 0.06; 
eff_JP[9] = 0.479115; eff_stat_JP[9] = 0.00760533; effMC_JP[9] = 0.498279; effMC_stat_JP[9] = 0.00153066; 
sf_JP[9] = 0.961545; sf_stat_JP[9] = 0.0152632; sf_syst_JP[9] = 0.0730774; sf_eror_JP[9] = 0.0746544; 
mupt_syst_JP[9] = 0.012; gluon_syst_JP[9] = 0.005; pu_syst_JP[9] = 0.001; total_syst_JP[9] = 0.076; 
cor_syst_JP[9] = 0.001; inc_syst_JP[9] = 0.045; bias_syst_JP[9] = 0.06; 
eff_JP[10] = 0.446049; eff_stat_JP[10] = 0.00873969; effMC_JP[10] = 0.484511; effMC_stat_JP[10] = 0.00272648; 
sf_JP[10] = 0.92062; sf_stat_JP[10] = 0.0180382; sf_syst_JP[10] = 0.103109; sf_eror_JP[10] = 0.104675; 
mupt_syst_JP[10] = 0.007; gluon_syst_JP[10] = 0.022; pu_syst_JP[10] = 0.003; total_syst_JP[10] = 0.112; 
cor_syst_JP[10] = 0; inc_syst_JP[10] = 0.045; bias_syst_JP[10] = 0.1; 
eff_JP[11] = 0.435916; eff_stat_JP[11] = 0.0111323; effMC_JP[11] = 0.459159; effMC_stat_JP[11] = 0.00457391; 
sf_JP[11] = 0.949385; sf_stat_JP[11] = 0.024245; sf_syst_JP[11] = 0.10823; sf_eror_JP[11] = 0.110912; 
mupt_syst_JP[11] = 0.007; gluon_syst_JP[11] = 0.029; pu_syst_JP[11] = 0.003; total_syst_JP[11] = 0.114; 
cor_syst_JP[11] = 0; inc_syst_JP[11] = 0.045; bias_syst_JP[11] = 0.1; 
eff_JP[12] = 0.445118; eff_stat_JP[12] = 0.0158891; effMC_JP[12] = 0.45759; effMC_stat_JP[12] = 0.0076102; 
sf_JP[12] = 0.972745; sf_stat_JP[12] = 0.0347235; sf_syst_JP[12] = 0.114784; sf_eror_JP[12] = 0.119921; 
mupt_syst_JP[12] = 0.007; gluon_syst_JP[12] = 0.023; pu_syst_JP[12] = 0.003; total_syst_JP[12] = 0.118; 
cor_syst_JP[12] = 0; inc_syst_JP[12] = 0.058; bias_syst_JP[12] = 0.1; 
eff_JP[13] = 0.393381; eff_stat_JP[13] = 0.017821; effMC_JP[13] = 0.420733; effMC_stat_JP[13] = 0.0065089; 
sf_JP[13] = 0.934989; sf_stat_JP[13] = 0.042357; sf_syst_JP[13] = 0.115939; sf_eror_JP[13] = 0.123434; 
mupt_syst_JP[13] = 0.007; gluon_syst_JP[13] = 0.043; pu_syst_JP[13] = 0.01; total_syst_JP[13] = 0.124; 
cor_syst_JP[13] = 0; inc_syst_JP[13] = 0.058; bias_syst_JP[13] = 0.1; 
eff_JP[14] = 0.387496; eff_stat_JP[14] = 0.0237165; effMC_JP[14] = 0.394128; effMC_stat_JP[14] = 0.00611641; 
sf_JP[14] = 0.983173; sf_stat_JP[14] = 0.0601746; sf_syst_JP[14] = 0.144526; sf_eror_JP[14] = 0.156553; 
mupt_syst_JP[14] = 0.007; gluon_syst_JP[14] = 0.089; pu_syst_JP[14] = 0.013; total_syst_JP[14] = 0.147; 
cor_syst_JP[14] = 0.007; inc_syst_JP[14] = 0.058; bias_syst_JP[14] = 0.1; 
}

// #############################################################################
if ( title == "JBPL" ) { 
// JBPL : PTREL
// float eff_PT[11];
eff_PT[0] = 0.625225; eff_PT[1] = 0.742144; eff_PT[2] = 0.799421; 
eff_PT[3] = 0.814032; eff_PT[4] = 0.844176; eff_PT[5] = 0.841144; 
eff_PT[6] = 0.87322; eff_PT[7] = 0.884921; eff_PT[8] = 0.910747; 
eff_PT[9] = 0.904489; eff_PT[10] = 0.864732; 
// float eff_stat_PT[11];
eff_stat_PT[0] = 0.0114512; eff_stat_PT[1] = 0.00757409; eff_stat_PT[2] = 0.007912; 
eff_stat_PT[3] = 0.00818892; eff_stat_PT[4] = 0.00979731; eff_stat_PT[5] = 0.0127586; 
eff_stat_PT[6] = 0.0100035; eff_stat_PT[7] = 0.0118992; eff_stat_PT[8] = 0.013124; 
eff_stat_PT[9] = 0.0134864; eff_stat_PT[10] = 0.0290183; 
// float effMC_PT[11];
effMC_PT[0] = 0.700696; effMC_PT[1] = 0.77446; effMC_PT[2] = 0.813454; 
effMC_PT[3] = 0.843395; effMC_PT[4] = 0.86014; effMC_PT[5] = 0.872641; 
effMC_PT[6] = 0.88673; effMC_PT[7] = 0.889833; effMC_PT[8] = 0.902818; 
effMC_PT[9] = 0.907414; effMC_PT[10] = 0.91281; 
// float effMC_stat_PT[11];
effMC_stat_PT[0] = 0.00293031; effMC_stat_PT[1] = 0.00166865; effMC_stat_PT[2] = 0.00139984; 
effMC_stat_PT[3] = 0.00178491; effMC_stat_PT[4] = 0.00167209; effMC_stat_PT[5] = 0.00171082; 
effMC_stat_PT[6] = 0.00175561; effMC_stat_PT[7] = 0.00165617; effMC_stat_PT[8] = 0.00206615; 
effMC_stat_PT[9] = 0.0021098; effMC_stat_PT[10] = 0.00307741; 
// float sf_PT[11];
sf_PT[0] = 0.892292; sf_PT[1] = 0.958272; sf_PT[2] = 0.98275; 
sf_PT[3] = 0.965185; sf_PT[4] = 0.98144; sf_PT[5] = 0.963906; 
sf_PT[6] = 0.984765; sf_PT[7] = 0.994479; sf_PT[8] = 1.00878; 
sf_PT[9] = 0.996776; sf_PT[10] = 0.94733; 
// float sf_stat_PT[11];
sf_stat_PT[0] = 0.0167633; sf_stat_PT[1] = 0.0099954; sf_stat_PT[2] = 0.00987236; 
sf_stat_PT[3] = 0.00992201; sf_stat_PT[4] = 0.011549; sf_stat_PT[5] = 0.0147423; 
sf_stat_PT[6] = 0.0114486; sf_stat_PT[7] = 0.0134999; sf_stat_PT[8] = 0.0147189; 
sf_stat_PT[9] = 0.0150421; sf_stat_PT[10] = 0.0319501; 
// float away_syst_PT[11];
away_syst_PT[0] = 0.00313401; away_syst_PT[1] = 0.0348468; away_syst_PT[2] = 0.035737; 
away_syst_PT[3] = 0.0227973; away_syst_PT[4] = 0.0231813; away_syst_PT[5] = 0.0227671; 
away_syst_PT[6] = 0.0049271; away_syst_PT[7] = 0.0049757; away_syst_PT[8] = 0.0408379; 
away_syst_PT[9] = 0.0403519; away_syst_PT[10] = 0.14913; 
// float mupt_syst_PT[11];
mupt_syst_PT[0] = 0.023902; mupt_syst_PT[1] = 0.00493491; mupt_syst_PT[2] = 0.00506097; 
mupt_syst_PT[3] = 0.0115436; mupt_syst_PT[4] = 0.011738; mupt_syst_PT[5] = 0.0115283; 
mupt_syst_PT[6] = 0.0217691; mupt_syst_PT[7] = 0.0219839; mupt_syst_PT[8] = 0.0520299; 
mupt_syst_PT[9] = 0.0514108; mupt_syst_PT[10] = 0.027427; 
// float lc_syst_PT[11];
lc_syst_PT[0] = 0.00482303; lc_syst_PT[1] = 0.00514544; lc_syst_PT[2] = 0.00527688; 
lc_syst_PT[3] = 0.00617863; lc_syst_PT[4] = 0.00628268; lc_syst_PT[5] = 0.00617044; 
lc_syst_PT[6] = 0.00482266; lc_syst_PT[7] = 0.00487023; lc_syst_PT[8] = 0.0104099; 
lc_syst_PT[9] = 0.0102861; lc_syst_PT[10] = 0.017942; 
// float gluon_syst_PT[11];
gluon_syst_PT[0] = 0.00190097; gluon_syst_PT[1] = 0.00103704; gluon_syst_PT[2] = 0.00106353; 
gluon_syst_PT[3] = 0.000253945; gluon_syst_PT[4] = 0.000258222; gluon_syst_PT[5] = 0.000253609; 
gluon_syst_PT[6] = 0.00261982; gluon_syst_PT[7] = 0.00264566; gluon_syst_PT[8] = 0.00554284; 
gluon_syst_PT[9] = 0.00547688; gluon_syst_PT[10] = 0.049826; 
// float pu_syst_PT[11];
pu_syst_PT[0] = 0.00908804; pu_syst_PT[1] = 0.00873691; pu_syst_PT[2] = 0.00896009; 
pu_syst_PT[3] = 0.0026989; pu_syst_PT[4] = 0.00274435; pu_syst_PT[5] = 0.00269532; 
pu_syst_PT[6] = 0.00273624; pu_syst_PT[7] = 0.00276323; pu_syst_PT[8] = 0.00819845; 
pu_syst_PT[9] = 0.00810089; pu_syst_PT[10] = 0.035814; 
// float total_syst_PT[11];
total_syst_PT[0] = 0.0262792; total_syst_PT[1] = 0.0366407; total_syst_PT[2] = 0.0375766; 
total_syst_PT[3] = 0.0264291; total_syst_PT[4] = 0.0268742; total_syst_PT[5] = 0.0263941; 
total_syst_PT[6] = 0.0231469; total_syst_PT[7] = 0.0233753; total_syst_PT[8] = 0.0676842; 
total_syst_PT[9] = 0.0668788; total_syst_PT[10] = 0.164558; 

// JBPL : S8
eff_SY[0]=0.6677; eff_SY[1]=0.7239; eff_SY[2]=0.7821; eff_SY[3]=0.8140; eff_SY[4]=0.8226; eff_SY[5]=0.8720; eff_SY[6]=0.8969; eff_SY[7]=0.8542; eff_SY[9]=0.9099; 
eff_stat_SY[0]=0.0103; eff_stat_SY[1]=0.0073; eff_stat_SY[2]=0.0091; eff_stat_SY[3]=0.0095; eff_stat_SY[4]=0.0119; eff_stat_SY[5]=0.0198; eff_stat_SY[6]=0.0310; eff_stat_SY[7]=0.0265; eff_stat_SY[9]=0.0482; 
effMC_SY[0]=0.6591; effMC_SY[1]=0.7458; effMC_SY[2]=0.8034; effMC_SY[3]=0.8410; effMC_SY[4]=0.8691; effMC_SY[5]=0.8940; effMC_SY[6]=0.9038; effMC_SY[7]=0.9138; effMC_SY[9]=0.9203; 
effMC_stat_SY[0]=0.0006; effMC_stat_SY[1]=0.0006; effMC_stat_SY[2]=0.0006; effMC_stat_SY[3]=0.0007; effMC_stat_SY[4]=0.0005; effMC_stat_SY[5]=0.0006; effMC_stat_SY[6]=0.0007; effMC_stat_SY[7]=0.0006; effMC_stat_SY[9]=0.0007; 
sf_SY[0]=1.013; sf_SY[1]=0.9706; sf_SY[2]=0.9735; sf_SY[3]=0.9679; sf_SY[4]=0.9465; sf_SY[5]=0.9754; sf_SY[6]=0.9924; sf_SY[7]=0.9348; sf_SY[9]=0.9887; 
sf_stat_SY[0]=0.0157; sf_stat_SY[1]=0.0098; sf_stat_SY[2]=0.0114; sf_stat_SY[3]=0.0113; sf_stat_SY[4]=0.0137; sf_stat_SY[5]=0.0222; sf_stat_SY[6]=0.0343; sf_stat_SY[7]=0.029; sf_stat_SY[9]=0.0524; 
away[0]=1.4304; away[1]=0.2708; away[2]=0.9066; away[3]=6.1861; 
mu[0]=0.1528; mu[1]=0.1477; mu[2]=0.3328; mu[3]=0.5295; 
ptrel[0]=0.5277; ptrel[1]=0.2954; ptrel[2]=1.8017; ptrel[3]=1.5646; 
gsplit[0]=0.1898; gsplit[1]=0.3626; gsplit[2]=0.7123; gsplit[3]=0.6969; 
closure[0]=0.2033; closure[1]=0.0936; closure[2]=0.089; closure[3]=0.2287; 
pu[0]=0.0972; pu[1]=0.1477; pu[2]=1.0328; pu[3]=0.1805; 

// JBPL : IP3D
// float eff_IP[11];
eff_IP[0] = 0.602067; eff_IP[1] = 0.714728; eff_IP[2] = 0.790594; 
eff_IP[3] = 0.765914; eff_IP[4] = 0.854815; eff_IP[5] = 0.829299; 
eff_IP[6] = 0.848216; eff_IP[7] = 0.857053; eff_IP[8] = 0.882763; 
eff_IP[9] = 0.904504; eff_IP[10] = 0.930323; 
// float eff_stat_IP[11];
eff_stat_IP[0] = 0.0949281; eff_stat_IP[1] = 0.0181849; eff_stat_IP[2] = 0.0176732; 
eff_stat_IP[3] = 0.0148501; eff_stat_IP[4] = 0.0162948; eff_stat_IP[5] = 0.0198446; 
eff_stat_IP[6] = 0.00866533; eff_stat_IP[7] = 0.00927683; eff_stat_IP[8] = 0.00737873; 
eff_stat_IP[9] = 0.00786534; eff_stat_IP[10] = 0.0135923; 
// float effMC_IP[11];
effMC_IP[0] = 0.703702; effMC_IP[1] = 0.774712; effMC_IP[2] = 0.812828; 
effMC_IP[3] = 0.848288; effMC_IP[4] = 0.859782; effMC_IP[5] = 0.870824; 
effMC_IP[6] = 0.886115; effMC_IP[7] = 0.888829; effMC_IP[8] = 0.901168; 
effMC_IP[9] = 0.906641; effMC_IP[10] = 0.920208; 
// float effMC_stat_IP[11];
effMC_stat_IP[0] = 0.00291961; effMC_stat_IP[1] = 0.00178762; effMC_stat_IP[2] = 0.00147706; 
effMC_stat_IP[3] = 0.0017707; effMC_stat_IP[4] = 0.00164784; effMC_stat_IP[5] = 0.00168252; 
effMC_stat_IP[6] = 0.00170918; effMC_stat_IP[7] = 0.00160554; effMC_stat_IP[8] = 0.00199808; 
effMC_stat_IP[9] = 0.002028; effMC_stat_IP[10] = 0.00308273; 
// float sf_IP[11];
sf_IP[0] = 0.855571; sf_IP[1] = 0.922573; sf_IP[2] = 0.972646; 
sf_IP[3] = 0.902894; sf_IP[4] = 0.994223; sf_IP[5] = 0.952315; 
sf_IP[6] = 0.95723; sf_IP[7] = 0.964249; sf_IP[8] = 0.979577; 
sf_IP[9] = 0.997643; sf_IP[10] = 1.01099; 
// float sf_stat_IP[11];
sf_stat_IP[0] = 0.134945; sf_stat_IP[1] = 0.0235695; sf_stat_IP[2] = 0.0218145; 
sf_stat_IP[3] = 0.0176072; sf_stat_IP[4] = 0.0190478; sf_stat_IP[5] = 0.0228625; 
sf_stat_IP[6] = 0.00995179; sf_stat_IP[7] = 0.0105815; sf_stat_IP[8] = 0.00847113; 
sf_stat_IP[9] = 0.00895767; sf_stat_IP[10] = 0.0151542; 
// float away_syst_IP[11];
away_syst_IP[0] = 0.143189; away_syst_IP[1] = 0.0331337; away_syst_IP[2] = 0.0349321; 
away_syst_IP[3] = 0.030602; away_syst_IP[4] = 0.0336975; away_syst_IP[5] = 0.0322771; 
away_syst_IP[6] = 0.00813208; away_syst_IP[7] = 0.00819171; away_syst_IP[8] = 0.00588617; 
away_syst_IP[9] = 0.00599472; away_syst_IP[10] = 0.0126832; 
// float mupt_syst_IP[11];
mupt_syst_IP[0] = 0.0457331; mupt_syst_IP[1] = 0.0260302; mupt_syst_IP[2] = 0.0274431; 
mupt_syst_IP[3] = 0.00895619; mupt_syst_IP[4] = 0.00986213; mupt_syst_IP[5] = 0.00944642; 
mupt_syst_IP[6] = 0.0114235; mupt_syst_IP[7] = 0.0115073; mupt_syst_IP[8] = 0.00745506; 
mupt_syst_IP[9] = 0.00759255; mupt_syst_IP[10] = 0.00798277; 
// float gluon_syst_IP[11];
gluon_syst_IP[0] = 0.0352916; gluon_syst_IP[1] = 0.00101582; gluon_syst_IP[2] = 0.00107095; 
gluon_syst_IP[3] = 0.000653448; gluon_syst_IP[4] = 0.000719545; gluon_syst_IP[5] = 0.000689216; 
gluon_syst_IP[6] = 0.00228726; gluon_syst_IP[7] = 0.00230404; gluon_syst_IP[8] = 0.00665122; 
gluon_syst_IP[9] = 0.00677388; gluon_syst_IP[10] = 0.0188513; 
// float pu_syst_IP[11];
pu_syst_IP[0] = 0.0392405; pu_syst_IP[1] = 0.0091133; pu_syst_IP[2] = 0.00960792; 
pu_syst_IP[3] = 0.00432103; pu_syst_IP[4] = 0.00475811; pu_syst_IP[5] = 0.00455755; 
pu_syst_IP[6] = 0.00453575; pu_syst_IP[7] = 0.00456901; pu_syst_IP[8] = 0.00165038; 
pu_syst_IP[9] = 0.00168082; pu_syst_IP[10] = 0.0132584; 
// float total_syst_IP[11];
total_syst_IP[0] = 0.159311; total_syst_IP[1] = 0.0431219; total_syst_IP[2] = 0.0454624; 
total_syst_IP[3] = 0.0321838; total_syst_IP[4] = 0.0354392; total_syst_IP[5] = 0.0339454; 
total_syst_IP[6] = 0.0149141; total_syst_IP[7] = 0.0150235; total_syst_IP[8] = 0.0117127; 
total_syst_IP[9] = 0.0119287; total_syst_IP[10] = 0.0274908; 

// JBPL : JP
eff_JP[0] = 0.712799; eff_stat_JP[0] = 0.0124288; effMC_JP[0] = 0.718862; effMC_stat_JP[0] = 0.00104465; 
sf_JP[0] = 0.991564; sf_stat_JP[0] = 0.0172895; sf_syst_JP[0] = 0.122954; sf_eror_JP[0] = 0.124164; 
mupt_syst_JP[0] = 0.006; gluon_syst_JP[0] = 0.002; pu_syst_JP[0] = 0.001; total_syst_JP[0] = 0.124; 
cor_syst_JP[0] = 0.028; inc_syst_JP[0] = 0.119; bias_syst_JP[0] = 0.02; 
eff_JP[1] = 0.78106; eff_stat_JP[1] = 0.00898451; effMC_JP[1] = 0.785668; effMC_stat_JP[1] = 0.000719984; 
sf_JP[1] = 0.994132; sf_stat_JP[1] = 0.0114355; sf_syst_JP[1] = 0.060642; sf_eror_JP[1] = 0.0617108; 
mupt_syst_JP[1] = 0.006; gluon_syst_JP[1] = 0.001; pu_syst_JP[1] = 0.001; total_syst_JP[1] = 0.061; 
cor_syst_JP[1] = 0.013; inc_syst_JP[1] = 0.056; bias_syst_JP[1] = 0.02; 
eff_JP[2] = 0.824662; eff_stat_JP[2] = 0.00844805; effMC_JP[2] = 0.822938; effMC_stat_JP[2] = 0.000719566; 
sf_JP[2] = 1.0021; sf_stat_JP[2] = 0.0102657; sf_syst_JP[2] = 0.0601259; sf_eror_JP[2] = 0.060996; 
mupt_syst_JP[2] = 0.006; gluon_syst_JP[2] = 0.001; pu_syst_JP[2] = 0.001; total_syst_JP[2] = 0.06; 
cor_syst_JP[2] = 0.007; inc_syst_JP[2] = 0.056; bias_syst_JP[2] = 0.02; 
eff_JP[3] = 0.864478; eff_stat_JP[3] = 0.00813698; effMC_JP[3] = 0.86368; effMC_stat_JP[3] = 0.000894064; 
sf_JP[3] = 1.00092; sf_stat_JP[3] = 0.00942129; sf_syst_JP[3] = 0.0530489; sf_eror_JP[3] = 0.053879; 
mupt_syst_JP[3] = 0.006; gluon_syst_JP[3] = 0.001; pu_syst_JP[3] = 0.001; total_syst_JP[3] = 0.053; 
cor_syst_JP[3] = 0.004; inc_syst_JP[3] = 0.049; bias_syst_JP[3] = 0.02; 
eff_JP[4] = 0.875248; eff_stat_JP[4] = 0.00631372; effMC_JP[4] = 0.876501; effMC_stat_JP[4] = 0.000860713; 
sf_JP[4] = 0.998571; sf_stat_JP[4] = 0.00720331; sf_syst_JP[4] = 0.0529243; sf_eror_JP[4] = 0.0534122; 
mupt_syst_JP[4] = 0.006; gluon_syst_JP[4] = 0.001; pu_syst_JP[4] = 0; total_syst_JP[4] = 0.053; 
cor_syst_JP[4] = 0.002; inc_syst_JP[4] = 0.049; bias_syst_JP[4] = 0.02; 
eff_JP[5] = 0.879474; eff_stat_JP[5] = 0.00509964; effMC_JP[5] = 0.882735; effMC_stat_JP[5] = 0.000874136; 
sf_JP[5] = 0.996307; sf_stat_JP[5] = 0.00577708; sf_syst_JP[5] = 0.0528043; sf_eror_JP[5] = 0.0531194; 
mupt_syst_JP[5] = 0.006; gluon_syst_JP[5] = 0.001; pu_syst_JP[5] = 0; total_syst_JP[5] = 0.053; 
cor_syst_JP[5] = 0.002; inc_syst_JP[5] = 0.049; bias_syst_JP[5] = 0.02; 
eff_JP[6] = 0.909563; eff_stat_JP[6] = 0.00881127; effMC_JP[6] = 0.9023; effMC_stat_JP[6] = 0.000773311; 
sf_JP[6] = 1.00805; sf_stat_JP[6] = 0.00976533; sf_syst_JP[6] = 0.0443541; sf_eror_JP[6] = 0.0454164; 
mupt_syst_JP[6] = 0.004; gluon_syst_JP[6] = 0.001; pu_syst_JP[6] = 0.001; total_syst_JP[6] = 0.044; 
cor_syst_JP[6] = 0.001; inc_syst_JP[6] = 0.043; bias_syst_JP[6] = 0.01; 
eff_JP[7] = 0.909859; eff_stat_JP[7] = 0.00945682; effMC_JP[7] = 0.906955; effMC_stat_JP[7] = 0.00075854; 
sf_JP[7] = 1.00321; sf_stat_JP[7] = 0.010427; sf_syst_JP[7] = 0.0441411; sf_eror_JP[7] = 0.0453559; 
mupt_syst_JP[7] = 0.004; gluon_syst_JP[7] = 0.001; pu_syst_JP[7] = 0; total_syst_JP[7] = 0.044; 
cor_syst_JP[7] = 0.001; inc_syst_JP[7] = 0.043; bias_syst_JP[7] = 0.01; 
eff_JP[8] = 0.915801; eff_stat_JP[8] = 0.00685511; effMC_JP[8] = 0.916611; effMC_stat_JP[8] = 0.000695241; 
sf_JP[8] = 0.999118; sf_stat_JP[8] = 0.00747876; sf_syst_JP[8] = 0.0439612; sf_eror_JP[8] = 0.0445928; 
mupt_syst_JP[8] = 0.004; gluon_syst_JP[8] = 0.002; pu_syst_JP[8] = 0; total_syst_JP[8] = 0.044; 
cor_syst_JP[8] = 0.001; inc_syst_JP[8] = 0.043; bias_syst_JP[8] = 0.01; 
eff_JP[9] = 0.918183; eff_stat_JP[9] = 0.010487; effMC_JP[9] = 0.921375; effMC_stat_JP[9] = 0.000823968; 
sf_JP[9] = 0.996534; sf_stat_JP[9] = 0.0113819; sf_syst_JP[9] = 0.0538128; sf_eror_JP[9] = 0.0550034; 
mupt_syst_JP[9] = 0.004; gluon_syst_JP[9] = 0.002; pu_syst_JP[9] = 0.001; total_syst_JP[9] = 0.054; 
cor_syst_JP[9] = 0.001; inc_syst_JP[9] = 0.053; bias_syst_JP[9] = 0.01; 
eff_JP[10] = 0.930344; eff_stat_JP[10] = 0.0120334; effMC_JP[10] = 0.925099; effMC_stat_JP[10] = 0.00143608; 
sf_JP[10] = 1.00567; sf_stat_JP[10] = 0.0130077; sf_syst_JP[10] = 0.0543062; sf_eror_JP[10] = 0.0558423; 
mupt_syst_JP[10] = 0.002; gluon_syst_JP[10] = 0.003; pu_syst_JP[10] = 0; total_syst_JP[10] = 0.054; 
cor_syst_JP[10] = 0; inc_syst_JP[10] = 0.053; bias_syst_JP[10] = 0.01; 
eff_JP[11] = 0.93139; eff_stat_JP[11] = 0.0116994; effMC_JP[11] = 0.929691; effMC_stat_JP[11] = 0.00234663; 
sf_JP[11] = 1.00183; sf_stat_JP[11] = 0.0125842; sf_syst_JP[11] = 0.0540988; sf_eror_JP[11] = 0.0555432; 
mupt_syst_JP[11] = 0.002; gluon_syst_JP[11] = 0.004; pu_syst_JP[11] = 0.001; total_syst_JP[11] = 0.054; 
cor_syst_JP[11] = 0; inc_syst_JP[11] = 0.053; bias_syst_JP[11] = 0.01; 
eff_JP[12] = 0.927758; eff_stat_JP[12] = 0.0129392; effMC_JP[12] = 0.935072; effMC_stat_JP[12] = 0.00376386; 
sf_JP[12] = 0.992178; sf_stat_JP[12] = 0.0138377; sf_syst_JP[12] = 0.083343; sf_eror_JP[12] = 0.0844839; 
mupt_syst_JP[12] = 0.002; gluon_syst_JP[12] = 0.004; pu_syst_JP[12] = 0.001; total_syst_JP[12] = 0.084; 
cor_syst_JP[12] = 0; inc_syst_JP[12] = 0.083; bias_syst_JP[12] = 0.01; 
eff_JP[13] = 0.940462; eff_stat_JP[13] = 0.0224368; effMC_JP[13] = 0.936155; effMC_stat_JP[13] = 0.00322332; 
sf_JP[13] = 1.0046; sf_stat_JP[13] = 0.023967; sf_syst_JP[13] = 0.0843864; sf_eror_JP[13] = 0.0877239; 
mupt_syst_JP[13] = 0.002; gluon_syst_JP[13] = 0.006; pu_syst_JP[13] = 0.006; total_syst_JP[13] = 0.084; 
cor_syst_JP[13] = 0; inc_syst_JP[13] = 0.083; bias_syst_JP[13] = 0.01; 
eff_JP[14] = 0.944294; eff_stat_JP[14] = 0.0253142; effMC_JP[14] = 0.934751; effMC_stat_JP[14] = 0.00309116; 
sf_JP[14] = 1.01021; sf_stat_JP[14] = 0.0270812; sf_syst_JP[14] = 0.0858679; sf_eror_JP[14] = 0.0900371; 
mupt_syst_JP[14] = 0.002; gluon_syst_JP[14] = 0.01; pu_syst_JP[14] = 0.008; total_syst_JP[14] = 0.085; 
cor_syst_JP[14] = 0.007; inc_syst_JP[14] = 0.083; bias_syst_JP[14] = 0.01; 
}

// #############################################################################
if ( title == "JBPM" ) { 
// JBPM : PTREL
// float eff_PT[11];
eff_PT[0] = 0.359006; eff_PT[1] = 0.465528; eff_PT[2] = 0.53428; 
eff_PT[3] = 0.566854; eff_PT[4] = 0.592257; eff_PT[5] = 0.609492; 
eff_PT[6] = 0.643292; eff_PT[7] = 0.650207; eff_PT[8] = 0.641303; 
eff_PT[9] = 0.623755; eff_PT[10] = 0.62268; 
// float eff_stat_PT[11];
eff_stat_PT[0] = 0.00978639; eff_stat_PT[1] = 0.00748801; eff_stat_PT[2] = 0.0085123; 
eff_stat_PT[3] = 0.00870705; eff_stat_PT[4] = 0.0106063; eff_stat_PT[5] = 0.0129643; 
eff_stat_PT[6] = 0.0107143; eff_stat_PT[7] = 0.0117344; eff_stat_PT[8] = 0.0146605; 
eff_stat_PT[9] = 0.0134934; eff_stat_PT[10] = 0.027125; 
// float effMC_PT[11];
effMC_PT[0] = 0.420545; effMC_PT[1] = 0.518706; effMC_PT[2] = 0.579137; 
effMC_PT[3] = 0.633271; effMC_PT[4] = 0.653868; effMC_PT[5] = 0.678984; 
effMC_PT[6] = 0.710783; effMC_PT[7] = 0.701237; effMC_PT[8] = 0.734198; 
effMC_PT[9] = 0.737901; effMC_PT[10] = 0.750494; 
// float effMC_stat_PT[11];
effMC_stat_PT[0] = 0.00335965; effMC_stat_PT[1] = 0.00213464; effMC_stat_PT[2] = 0.00191198; 
effMC_stat_PT[3] = 0.00259212; effMC_stat_PT[4] = 0.00250695; effMC_stat_PT[5] = 0.00258624; 
effMC_stat_PT[6] = 0.00273972; effMC_stat_PT[7] = 0.00262084; effMC_stat_PT[8] = 0.00335864; 
effMC_stat_PT[9] = 0.00350992; effMC_stat_PT[10] = 0.00514303; 
// float sf_PT[11];
sf_PT[0] = 0.853667; sf_PT[1] = 0.897478; sf_PT[2] = 0.922546; 
sf_PT[3] = 0.89512; sf_PT[4] = 0.905774; sf_PT[5] = 0.897653; 
sf_PT[6] = 0.905047; sf_PT[7] = 0.927229; sf_PT[8] = 0.873474; 
sf_PT[9] = 0.845309; sf_PT[10] = 0.829694; 
// float sf_stat_PT[11];
sf_stat_PT[0] = 0.0242495; sf_stat_PT[1] = 0.0149009; sf_stat_PT[2] = 0.0150105; 
sf_stat_PT[3] = 0.0142291; sf_stat_PT[4] = 0.0165884; sf_stat_PT[5] = 0.0193974; 
sf_stat_PT[6] = 0.0154723; sf_stat_PT[7] = 0.0170889; sf_stat_PT[8] = 0.0203639; 
sf_stat_PT[9] = 0.018723; sf_stat_PT[10] = 0.0365874; 
// float away_syst_PT[11];
away_syst_PT[0] = 0.026753; away_syst_PT[1] = 0.0372876; away_syst_PT[2] = 0.0383291; 
away_syst_PT[3] = 0.0217542; away_syst_PT[4] = 0.0220132; away_syst_PT[5] = 0.0218158; 
away_syst_PT[6] = 0.0104402; away_syst_PT[7] = 0.0106961; away_syst_PT[8] = 0.0755603; 
away_syst_PT[9] = 0.0731239; away_syst_PT[10] = 0.016013; 
// float mupt_syst_PT[11];
mupt_syst_PT[0] = 0.013924; mupt_syst_PT[1] = 0.00476967; mupt_syst_PT[2] = 0.00490289; 
mupt_syst_PT[3] = 0.010326; mupt_syst_PT[4] = 0.0104489; mupt_syst_PT[5] = 0.0103552; 
mupt_syst_PT[6] = 0.0166145; mupt_syst_PT[7] = 0.0170217; mupt_syst_PT[8] = 0.0662552; 
mupt_syst_PT[9] = 0.0641188; mupt_syst_PT[10] = 0.032049; 
// float lc_syst_PT[11];
lc_syst_PT[0] = 0.00541002; lc_syst_PT[1] = 0.00478545; lc_syst_PT[2] = 0.00491912; 
lc_syst_PT[3] = 0.00638232; lc_syst_PT[4] = 0.00645828; lc_syst_PT[5] = 0.00640038; 
lc_syst_PT[6] = 0.00466093; lc_syst_PT[7] = 0.00477516; lc_syst_PT[8] = 0.0157563; 
lc_syst_PT[9] = 0.0152483; lc_syst_PT[10] = 0.042382; 
// float gluon_syst_PT[11];
gluon_syst_PT[0] = 0.000160992; gluon_syst_PT[1] = 0.00137346; gluon_syst_PT[2] = 0.00141183; 
gluon_syst_PT[3] = 0.00363849; gluon_syst_PT[4] = 0.00368179; gluon_syst_PT[5] = 0.00364878; 
gluon_syst_PT[6] = 0.00416996; gluon_syst_PT[7] = 0.00427216; gluon_syst_PT[8] = 0.0091773; 
gluon_syst_PT[9] = 0.00888138; gluon_syst_PT[10] = 0.081834; 
// float pu_syst_PT[11];
pu_syst_PT[0] = 0.018625; pu_syst_PT[1] = 0.0112507; pu_syst_PT[2] = 0.0115649; 
pu_syst_PT[3] = 0.00507075; pu_syst_PT[4] = 0.0051311; pu_syst_PT[5] = 0.00508509; 
pu_syst_PT[6] = 0.00843272; pu_syst_PT[7] = 0.0086394; pu_syst_PT[8] = 0.0135715; 
pu_syst_PT[9] = 0.0131339; pu_syst_PT[10] = 0.036163; 
// float total_syst_PT[11];
total_syst_PT[0] = 0.0358579; total_syst_PT[1] = 0.0395535; total_syst_PT[2] = 0.0406583; 
total_syst_PT[3] = 0.0256819; total_syst_PT[4] = 0.0259875; total_syst_PT[5] = 0.0257545; 
total_syst_PT[6] = 0.0222545; total_syst_PT[7] = 0.0228; total_syst_PT[8] = 0.103033; 
total_syst_PT[9] = 0.0997107; total_syst_PT[10] = 0.105282; 

// JBPM : S8
eff_SY[0]=0.3698; eff_SY[1]=0.4567; eff_SY[2]=0.5250; eff_SY[3]=0.5884; eff_SY[4]=0.5987; eff_SY[5]=0.6798; eff_SY[6]=0.7340; eff_SY[7]=0.6331; eff_SY[9]=0.7006; 
eff_stat_SY[0]=0.0092; eff_stat_SY[1]=0.0069; eff_stat_SY[2]=0.0089; eff_stat_SY[3]=0.0105; eff_stat_SY[4]=0.0116; eff_stat_SY[5]=0.0219; eff_stat_SY[6]=0.0388; eff_stat_SY[7]=0.0310; eff_stat_SY[9]=0.0620; 
effMC_SY[0]=0.3850; effMC_SY[1]=0.4928; effMC_SY[2]=0.5717; effMC_SY[3]=0.6287; effMC_SY[4]=0.6757; effMC_SY[5]=0.7186; effMC_SY[6]=0.7368; effMC_SY[7]=0.7575; effMC_SY[9]=0.7685; 
effMC_stat_SY[0]=0.0006; effMC_stat_SY[1]=0.0007; effMC_stat_SY[2]=0.0008; effMC_stat_SY[3]=0.0009; effMC_stat_SY[4]=0.0007; effMC_stat_SY[5]=0.0009; effMC_stat_SY[6]=0.0010; effMC_stat_SY[7]=0.0009; effMC_stat_SY[9]=0.0011; 
sf_SY[0]=0.9605; sf_SY[1]=0.9267; sf_SY[2]=0.9183; sf_SY[3]=0.9359; sf_SY[4]=0.8860; sf_SY[5]=0.946; sf_SY[6]=0.9962; sf_SY[7]=0.8358; sf_SY[9]=0.9116; 
sf_stat_SY[0]=0.0239; sf_stat_SY[1]=0.0141; sf_stat_SY[2]=0.0156; sf_stat_SY[3]=0.0168; sf_stat_SY[4]=0.0172; sf_stat_SY[5]=0.0305; sf_stat_SY[6]=0.0527; sf_stat_SY[7]=0.0409; sf_stat_SY[9]=0.0807; 
away[0]=2.0036; away[1]=0.9687; away[2]=1.5036; away[3]=16.5653; 
mu[0]=0.8014; mu[1]=0.6288; mu[2]=1.1679; mu[3]=1.0025; 
ptrel[0]=0.2894; ptrel[1]=2.5323; ptrel[2]=3.5766; ptrel[3]=1.1175; 
gsplit[0]=0.4115; gsplit[1]=0.9194; gsplit[2]=1.8748; gsplit[3]=0.837; 
closure[0]=0.1852; closure[1]=0.0919; closure[2]=0.0276; closure[3]=0.3924; 
pu[0]=0.1336; pu[1]=0.136; pu[2]=1.4453; pu[3]=0.5588; 

// JBPM : IP3D
// float eff_IP[11];
eff_IP[0] = 0.198069; eff_IP[1] = 0.405917; eff_IP[2] = 0.478958; 
eff_IP[3] = 0.522153; eff_IP[4] = 0.603872; eff_IP[5] = 0.620561; 
eff_IP[6] = 0.646192; eff_IP[7] = 0.626955; eff_IP[8] = 0.690756; 
eff_IP[9] = 0.71078; eff_IP[10] = 0.775837; 
// float eff_stat_IP[11];
eff_stat_IP[0] = 0.0684843; eff_stat_IP[1] = 0.0775523; eff_stat_IP[2] = 0.0808511; 
eff_stat_IP[3] = 0.0117437; eff_stat_IP[4] = 0.0160494; eff_stat_IP[5] = 0.0196242; 
eff_stat_IP[6] = 0.0102369; eff_stat_IP[7] = 0.0123082; eff_stat_IP[8] = 0.0102882; 
eff_stat_IP[9] = 0.0114164; eff_stat_IP[10] = 0.0212283; 
// float effMC_IP[11];
effMC_IP[0] = 0.418808; effMC_IP[1] = 0.519792; effMC_IP[2] = 0.577078; 
effMC_IP[3] = 0.635105; effMC_IP[4] = 0.653273; effMC_IP[5] = 0.678594; 
effMC_IP[6] = 0.708263; effMC_IP[7] = 0.701916; effMC_IP[8] = 0.7297; 
effMC_IP[9] = 0.734697; effMC_IP[10] = 0.753169; 
// float effMC_stat_IP[11];
effMC_stat_IP[0] = 0.00334178; effMC_stat_IP[1] = 0.00228369; effMC_stat_IP[2] = 0.00201138; 
effMC_stat_IP[3] = 0.00257356; effMC_stat_IP[4] = 0.00246792; effMC_stat_IP[5] = 0.00254023; 
effMC_stat_IP[6] = 0.00266916; effMC_stat_IP[7] = 0.00254278; effMC_stat_IP[8] = 0.00325197; 
effMC_stat_IP[9] = 0.0033829; effMC_stat_IP[10] = 0.00514241; 
// float sf_IP[11];
sf_IP[0] = 0.472934; sf_IP[1] = 0.780922; sf_IP[2] = 0.829971; 
sf_IP[3] = 0.822152; sf_IP[4] = 0.924379; sf_IP[5] = 0.914481; 
sf_IP[6] = 0.912363; sf_IP[7] = 0.893205; sf_IP[8] = 0.94663; 
sf_IP[9] = 0.967447; sf_IP[10] = 1.0301; 
// float sf_stat_IP[11];
sf_stat_IP[0] = 0.163565; sf_stat_IP[1] = 0.149238; sf_stat_IP[2] = 0.140134; 
sf_stat_IP[3] = 0.0187886; sf_stat_IP[4] = 0.0248146; sf_stat_IP[5] = 0.0291208; 
sf_stat_IP[6] = 0.0148569; sf_stat_IP[7] = 0.0178311; sf_stat_IP[8] = 0.0147169; 
sf_stat_IP[9] = 0.0161649; sf_stat_IP[10] = 0.0290495; 
// float away_syst_IP[11];
away_syst_IP[0] = 0.0447266; away_syst_IP[1] = 0.0403461; away_syst_IP[2] = 0.0428802; 
away_syst_IP[3] = 0.0079931; away_syst_IP[4] = 0.00898696; away_syst_IP[5] = 0.00889073; 
away_syst_IP[6] = 0.0159624; away_syst_IP[7] = 0.0156272; away_syst_IP[8] = 0.0264368; 
away_syst_IP[9] = 0.0270182; away_syst_IP[10] = 0.0617519; 
// float mupt_syst_IP[11];
mupt_syst_IP[0] = 0.0269484; mupt_syst_IP[1] = 0.0320877; mupt_syst_IP[2] = 0.0341031; 
mupt_syst_IP[3] = 0.0463739; mupt_syst_IP[4] = 0.05214; mupt_syst_IP[5] = 0.0515817; 
mupt_syst_IP[6] = 0.0073377; mupt_syst_IP[7] = 0.00718362; mupt_syst_IP[8] = 0.00397467; 
mupt_syst_IP[9] = 0.00406208; mupt_syst_IP[10] = 0.00391113; 
// float gluon_syst_IP[11];
gluon_syst_IP[0] = 0.00368183; gluon_syst_IP[1] = 0.00107853; gluon_syst_IP[2] = 0.00114627; 
gluon_syst_IP[3] = 0.0268446; gluon_syst_IP[4] = 0.0301825; gluon_syst_IP[5] = 0.0298593; 
gluon_syst_IP[6] = 0.00433544; gluon_syst_IP[7] = 0.0042444; gluon_syst_IP[8] = 0.0159772; 
gluon_syst_IP[9] = 0.0163286; gluon_syst_IP[10] = 0.0652504; 
// float pu_syst_IP[11];
pu_syst_IP[0] = 0.0113479; pu_syst_IP[1] = 0.00747931; pu_syst_IP[2] = 0.00794907; 
pu_syst_IP[3] = 0.0336469; pu_syst_IP[4] = 0.0378306; pu_syst_IP[5] = 0.0374255; 
pu_syst_IP[6] = 0.0108236; pu_syst_IP[7] = 0.0105963; pu_syst_IP[8] = 0.00212729; 
pu_syst_IP[9] = 0.00217407; pu_syst_IP[10] = 0.0230046; 
// float total_syst_IP[11];
total_syst_IP[0] = 0.0535631; total_syst_IP[1] = 0.0521012; total_syst_IP[2] = 0.0553736; 
total_syst_IP[3] = 0.0637744; total_syst_IP[4] = 0.0717041; total_syst_IP[5] = 0.0709363; 
total_syst_IP[6] = 0.0210852; total_syst_IP[7] = 0.0206424; total_syst_IP[8] = 0.031217; 
total_syst_IP[9] = 0.0319035; total_syst_IP[10] = 0.0928193; 

// JBPM : JP
eff_JP[0] = 0.449797; eff_stat_JP[0] = 0.0167349; effMC_JP[0] = 0.463249; effMC_stat_JP[0] = 0.00115873; 
sf_JP[0] = 0.970964; sf_stat_JP[0] = 0.036125; sf_syst_JP[0] = 0.278667; sf_eror_JP[0] = 0.280998; 
mupt_syst_JP[0] = 0.013; gluon_syst_JP[0] = 0; pu_syst_JP[0] = 0; total_syst_JP[0] = 0.287; 
cor_syst_JP[0] = 0.028; inc_syst_JP[0] = 0.279; bias_syst_JP[0] = 0.06; 
eff_JP[1] = 0.534605; eff_stat_JP[1] = 0.014156; effMC_JP[1] = 0.551262; effMC_stat_JP[1] = 0.000872641; 
sf_JP[1] = 0.969789; sf_stat_JP[1] = 0.0256792; sf_syst_JP[1] = 0.082432; sf_eror_JP[1] = 0.0863392; 
mupt_syst_JP[1] = 0.013; gluon_syst_JP[1] = 0.002; pu_syst_JP[1] = 0.001; total_syst_JP[1] = 0.085; 
cor_syst_JP[1] = 0.013; inc_syst_JP[1] = 0.057; bias_syst_JP[1] = 0.06; 
eff_JP[2] = 0.591183; eff_stat_JP[2] = 0.0155203; effMC_JP[2] = 0.607315; effMC_stat_JP[2] = 0.000920564; 
sf_JP[2] = 0.973436; sf_stat_JP[2] = 0.0255556; sf_syst_JP[2] = 0.0817686; sf_eror_JP[2] = 0.0856691; 
mupt_syst_JP[2] = 0.013; gluon_syst_JP[2] = 0.003; pu_syst_JP[2] = 0.001; total_syst_JP[2] = 0.084; 
cor_syst_JP[2] = 0.007; inc_syst_JP[2] = 0.057; bias_syst_JP[2] = 0.06; 
eff_JP[3] = 0.646734; eff_stat_JP[3] = 0.0162871; effMC_JP[3] = 0.664714; effMC_stat_JP[3] = 0.00123009; 
sf_JP[3] = 0.972947; sf_stat_JP[3] = 0.0245025; sf_syst_JP[3] = 0.0749169; sf_eror_JP[3] = 0.0788221; 
mupt_syst_JP[3] = 0.013; gluon_syst_JP[3] = 0.004; pu_syst_JP[3] = 0; total_syst_JP[3] = 0.077; 
cor_syst_JP[3] = 0.004; inc_syst_JP[3] = 0.046; bias_syst_JP[3] = 0.06; 
eff_JP[4] = 0.668901; eff_stat_JP[4] = 0.0161092; effMC_JP[4] = 0.686245; effMC_stat_JP[4] = 0.00121391; 
sf_JP[4] = 0.974724; sf_stat_JP[4] = 0.0234744; sf_syst_JP[4] = 0.0750537; sf_eror_JP[4] = 0.0786391; 
mupt_syst_JP[4] = 0.013; gluon_syst_JP[4] = 0.005; pu_syst_JP[4] = 0; total_syst_JP[4] = 0.077; 
cor_syst_JP[4] = 0.002; inc_syst_JP[4] = 0.046; bias_syst_JP[4] = 0.06; 
eff_JP[5] = 0.681444; eff_stat_JP[5] = 0.0138893; effMC_JP[5] = 0.701431; effMC_stat_JP[5] = 0.00124335; 
sf_JP[5] = 0.971509; sf_stat_JP[5] = 0.0198014; sf_syst_JP[5] = 0.0748062; sf_eror_JP[5] = 0.0773825; 
mupt_syst_JP[5] = 0.013; gluon_syst_JP[5] = 0.007; pu_syst_JP[5] = 0; total_syst_JP[5] = 0.077; 
cor_syst_JP[5] = 0.002; inc_syst_JP[5] = 0.046; bias_syst_JP[5] = 0.06; 
eff_JP[6] = 0.732045; eff_stat_JP[6] = 0.0164967; effMC_JP[6] = 0.734595; effMC_stat_JP[6] = 0.00115003; 
sf_JP[6] = 0.996528; sf_stat_JP[6] = 0.0224569; sf_syst_JP[6] = 0.054809; sf_eror_JP[6] = 0.0592313; 
mupt_syst_JP[6] = 0.009; gluon_syst_JP[6] = 0.006; pu_syst_JP[6] = 0.001; total_syst_JP[6] = 0.055; 
cor_syst_JP[6] = 0.001; inc_syst_JP[6] = 0.036; bias_syst_JP[6] = 0.04; 
eff_JP[7] = 0.736179; eff_stat_JP[7] = 0.0146934; effMC_JP[7] = 0.744059; effMC_stat_JP[7] = 0.0011395; 
sf_JP[7] = 0.989414; sf_stat_JP[7] = 0.0197477; sf_syst_JP[7] = 0.0544178; sf_eror_JP[7] = 0.0578901; 
mupt_syst_JP[7] = 0.009; gluon_syst_JP[7] = 0.006; pu_syst_JP[7] = 0.001; total_syst_JP[7] = 0.055; 
cor_syst_JP[7] = 0.001; inc_syst_JP[7] = 0.036; bias_syst_JP[7] = 0.04; 
eff_JP[8] = 0.755536; eff_stat_JP[8] = 0.0147141; effMC_JP[8] = 0.764007; effMC_stat_JP[8] = 0.00106779; 
sf_JP[8] = 0.988915; sf_stat_JP[8] = 0.019259; sf_syst_JP[8] = 0.0543903; sf_eror_JP[8] = 0.0576994; 
mupt_syst_JP[8] = 0.009; gluon_syst_JP[8] = 0.005; pu_syst_JP[8] = 0; total_syst_JP[8] = 0.055; 
cor_syst_JP[8] = 0.001; inc_syst_JP[8] = 0.036; bias_syst_JP[8] = 0.04; 
eff_JP[9] = 0.75873; eff_stat_JP[9] = 0.0156463; effMC_JP[9] = 0.769813; effMC_stat_JP[9] = 0.00128868; 
sf_JP[9] = 0.985603; sf_stat_JP[9] = 0.0203247; sf_syst_JP[9] = 0.0581506; sf_eror_JP[9] = 0.0616002; 
mupt_syst_JP[9] = 0.009; gluon_syst_JP[9] = 0.003; pu_syst_JP[9] = 0; total_syst_JP[9] = 0.059; 
cor_syst_JP[9] = 0.001; inc_syst_JP[9] = 0.043; bias_syst_JP[9] = 0.04; 
eff_JP[10] = 0.759769; eff_stat_JP[10] = 0.0187045; effMC_JP[10] = 0.772144; effMC_stat_JP[10] = 0.00228834; 
sf_JP[10] = 0.983971; sf_stat_JP[10] = 0.0242242; sf_syst_JP[10] = 0.0580543; sf_eror_JP[10] = 0.0629056; 
mupt_syst_JP[10] = 0.007; gluon_syst_JP[10] = 0.003; pu_syst_JP[10] = 0.001; total_syst_JP[10] = 0.059; 
cor_syst_JP[10] = 0; inc_syst_JP[10] = 0.043; bias_syst_JP[10] = 0.04; 
eff_JP[11] = 0.76926; eff_stat_JP[11] = 0.0175641; effMC_JP[11] = 0.778353; effMC_stat_JP[11] = 0.00381233; 
sf_JP[11] = 0.98832; sf_stat_JP[11] = 0.0225657; sf_syst_JP[11] = 0.0583109; sf_eror_JP[11] = 0.0625249; 
mupt_syst_JP[11] = 0.007; gluon_syst_JP[11] = 0.005; pu_syst_JP[11] = 0.002; total_syst_JP[11] = 0.059; 
cor_syst_JP[11] = 0; inc_syst_JP[11] = 0.043; bias_syst_JP[11] = 0.04; 
eff_JP[12] = 0.766537; eff_stat_JP[12] = 0.0282228; effMC_JP[12] = 0.787859; effMC_stat_JP[12] = 0.00624498; 
sf_JP[12] = 0.972941; sf_stat_JP[12] = 0.0358221; sf_syst_JP[12] = 0.065187; sf_eror_JP[12] = 0.0743813; 
mupt_syst_JP[12] = 0.007; gluon_syst_JP[12] = 0.002; pu_syst_JP[12] = 0.011; total_syst_JP[12] = 0.067; 
cor_syst_JP[12] = 0; inc_syst_JP[12] = 0.052; bias_syst_JP[12] = 0.04; 
eff_JP[13] = 0.782338; eff_stat_JP[13] = 0.0306797; effMC_JP[13] = 0.783055; effMC_stat_JP[13] = 0.00543421; 
sf_JP[13] = 0.999082; sf_stat_JP[13] = 0.0391794; sf_syst_JP[13] = 0.0669385; sf_eror_JP[13] = 0.0775615; 
mupt_syst_JP[13] = 0.007; gluon_syst_JP[13] = 0.011; pu_syst_JP[13] = 0.003; total_syst_JP[13] = 0.067; 
cor_syst_JP[13] = 0; inc_syst_JP[13] = 0.052; bias_syst_JP[13] = 0.04; 
eff_JP[14] = 0.81942; eff_stat_JP[14] = 0.0346677; effMC_JP[14] = 0.773639; effMC_stat_JP[14] = 0.0052379; 
sf_JP[14] = 1.05918; sf_stat_JP[14] = 0.0448112; sf_syst_JP[14] = 0.0730832; sf_eror_JP[14] = 0.0857275; 
mupt_syst_JP[14] = 0.007; gluon_syst_JP[14] = 0.01; pu_syst_JP[14] = 0.016; total_syst_JP[14] = 0.069; 
cor_syst_JP[14] = 0.007; inc_syst_JP[14] = 0.052; bias_syst_JP[14] = 0.04; 
}

// #############################################################################
if ( title == "JBPT" ) { 
// JBPT : PTREL
// float eff_PT[11];
eff_PT[0] = 0.148079; eff_PT[1] = 0.227749; eff_PT[2] = 0.303963; 
eff_PT[3] = 0.342541; eff_PT[4] = 0.363873; eff_PT[5] = 0.379256; 
eff_PT[6] = 0.410434; eff_PT[7] = 0.400071; eff_PT[8] = 0.414814; 
eff_PT[9] = 0.383199; eff_PT[10] = 0.430362; 
// float eff_stat_PT[11];
eff_stat_PT[0] = 0.00646738; eff_stat_PT[1] = 0.0057655; eff_stat_PT[2] = 0.00647746; 
eff_stat_PT[3] = 0.00784235; eff_stat_PT[4] = 0.0091602; eff_stat_PT[5] = 0.0107119; 
eff_stat_PT[6] = 0.00792703; eff_stat_PT[7] = 0.00936377; eff_stat_PT[8] = 0.0114147; 
eff_stat_PT[9] = 0.00992931; eff_stat_PT[10] = 0.018871; 
// float effMC_PT[11];
effMC_PT[0] = 0.187365; effMC_PT[1] = 0.277345; effMC_PT[2] = 0.338719; 
effMC_PT[3] = 0.400267; effMC_PT[4] = 0.4298; effMC_PT[5] = 0.464243; 
effMC_PT[6] = 0.491301; effMC_PT[7] = 0.492592; effMC_PT[8] = 0.52922; 
effMC_PT[9] = 0.529695; effMC_PT[10] = 0.547142; 
// float effMC_stat_PT[11];
effMC_stat_PT[0] = 0.00259759; effMC_stat_PT[1] = 0.00190019; effMC_stat_PT[2] = 0.00183456; 
effMC_stat_PT[3] = 0.00265678; effMC_stat_PT[4] = 0.00264728; effMC_stat_PT[5] = 0.00278184; 
effMC_stat_PT[6] = 0.00303831; effMC_stat_PT[7] = 0.00292628; effMC_stat_PT[8] = 0.0038819; 
effMC_stat_PT[9] = 0.00410244; effMC_stat_PT[10] = 0.00604523; 
// float sf_PT[11];
sf_PT[0] = 0.790322; sf_PT[1] = 0.821174; sf_PT[2] = 0.89739; 
sf_PT[3] = 0.855781; sf_PT[4] = 0.846611; sf_PT[5] = 0.816935; 
sf_PT[6] = 0.835403; sf_PT[7] = 0.812176; sf_PT[8] = 0.783822; 
sf_PT[9] = 0.723433; sf_PT[10] = 0.786563; 
// float sf_stat_PT[11];
sf_stat_PT[0] = 0.0362148; sf_stat_PT[1] = 0.0215361; sf_stat_PT[2] = 0.0197314; 
sf_stat_PT[3] = 0.0203996; sf_stat_PT[4] = 0.0219414; sf_stat_PT[5] = 0.0235875; 
sf_stat_PT[6] = 0.0169417; sf_stat_PT[7] = 0.0196119; sf_stat_PT[8] = 0.0223221; 
sf_stat_PT[9] = 0.0195648; sf_stat_PT[10] = 0.0355681; 
// float away_syst_PT[11];
away_syst_PT[0] = 0.030183; away_syst_PT[1] = 0.0539373; away_syst_PT[2] = 0.0589434; 
away_syst_PT[3] = 0.017481; away_syst_PT[4] = 0.0172937; away_syst_PT[5] = 0.0166875; 
away_syst_PT[6] = 0.0418964; away_syst_PT[7] = 0.0407315; away_syst_PT[8] = 0.0702237; 
away_syst_PT[9] = 0.0648134; away_syst_PT[10] = 0.119389; 
// float mupt_syst_PT[11];
mupt_syst_PT[0] = 0.030496; mupt_syst_PT[1] = 0.0133837; mupt_syst_PT[2] = 0.0146259; 
mupt_syst_PT[3] = 0.0189736; mupt_syst_PT[4] = 0.0187702; mupt_syst_PT[5] = 0.0181123; 
mupt_syst_PT[6] = 0.0373401; mupt_syst_PT[7] = 0.0363019; mupt_syst_PT[8] = 0.0837653; 
mupt_syst_PT[9] = 0.0773117; mupt_syst_PT[10] = 0.10844; 
// float lc_syst_PT[11];
lc_syst_PT[0] = 0.004731; lc_syst_PT[1] = 0.00424584; lc_syst_PT[2] = 0.00463991; 
lc_syst_PT[3] = 0.00643389; lc_syst_PT[4] = 0.00636495; lc_syst_PT[5] = 0.00614184; 
lc_syst_PT[6] = 0.00490863; lc_syst_PT[7] = 0.00477215; lc_syst_PT[8] = 0.0164237; 
lc_syst_PT[9] = 0.0151583; lc_syst_PT[10] = 0.059725; 
// float gluon_syst_PT[11];
gluon_syst_PT[0] = 0.00256097; gluon_syst_PT[1] = 0.00392474; gluon_syst_PT[2] = 0.00428901; 
gluon_syst_PT[3] = 0.00160066; gluon_syst_PT[4] = 0.00158351; gluon_syst_PT[5] = 0.001528; 
gluon_syst_PT[6] = 0.00412928; gluon_syst_PT[7] = 0.00401447; gluon_syst_PT[8] = 0.013033; 
gluon_syst_PT[9] = 0.0120289; gluon_syst_PT[10] = 0.036967; 
// float pu_syst_PT[11];
pu_syst_PT[0] = 0.010665; pu_syst_PT[1] = 0.0123667; pu_syst_PT[2] = 0.0135145; 
pu_syst_PT[3] = 0.0185595; pu_syst_PT[4] = 0.0183606; pu_syst_PT[5] = 0.017717; 
pu_syst_PT[6] = 0.00882187; pu_syst_PT[7] = 0.00857659; pu_syst_PT[8] = 0.0113891; 
pu_syst_PT[9] = 0.0105117; pu_syst_PT[10] = 0.055474; 
// float total_syst_PT[11];
total_syst_PT[0] = 0.0445388; total_syst_PT[1] = 0.0572252; total_syst_PT[2] = 0.0625364; 
total_syst_PT[3] = 0.0324653; total_syst_PT[4] = 0.0321174; total_syst_PT[5] = 0.0309916; 
total_syst_PT[6] = 0.0571713; total_syst_PT[7] = 0.0555818; total_syst_PT[8] = 0.111881; 
total_syst_PT[9] = 0.103261; total_syst_PT[10] = 0.184456; 

// JBPT : S8
eff_SY[0]=0.1425; eff_SY[1]=0.2155; eff_SY[2]=0.2977; eff_SY[3]=0.3580; eff_SY[4]=0.3785; eff_SY[5]=0.4411; eff_SY[6]=0.5232; eff_SY[7]=0.4655; 
eff_stat_SY[0]=0.0062; eff_stat_SY[1]=0.0054; eff_stat_SY[2]=0.0075; eff_stat_SY[3]=0.0087; eff_stat_SY[4]=0.0097; eff_stat_SY[5]=0.0166; eff_stat_SY[6]=0.0309; eff_stat_SY[7]=0.0329; 
effMC_SY[0]=0.1726; effMC_SY[1]=0.2634; effMC_SY[2]=0.3390; effMC_SY[3]=0.3991; effMC_SY[4]=0.4546; effMC_SY[5]=0.5091; effMC_SY[6]=0.5350; effMC_SY[7]=0.5616; 
effMC_stat_SY[0]=0.0005; effMC_stat_SY[1]=0.0006; effMC_stat_SY[2]=0.0007; effMC_stat_SY[3]=0.0009; effMC_stat_SY[4]=0.0008; effMC_stat_SY[5]=0.0010; effMC_stat_SY[6]=0.0012; effMC_stat_SY[7]=0.0011; 
sf_SY[0]=0.8256; sf_SY[1]=0.8181; sf_SY[2]=0.8782; sf_SY[3]=0.897; sf_SY[4]=0.8326; sf_SY[5]=0.8664; sf_SY[6]=0.9779; sf_SY[7]=0.8289; 
sf_stat_SY[0]=0.036; sf_stat_SY[1]=0.0206; sf_stat_SY[2]=0.0222; sf_stat_SY[3]=0.0219; sf_stat_SY[4]=0.0214; sf_stat_SY[5]=0.0327; sf_stat_SY[6]=0.0578; sf_stat_SY[7]=0.0586; 
away[0]=4.0714; away[1]=2.2447; away[2]=1.6164; away[3]=22.0187; 
mu[0]=0.9607; mu[1]=1.6699; mu[2]=2.5431; mu[3]=0.724; 
ptrel[0]=0.7777; ptrel[1]=0.7117; ptrel[2]=2.1552; ptrel[3]=0.8518; 
gsplit[0]=1.0341; gsplit[1]=1.9159; gsplit[2]=3.4655; gsplit[3]=0.8004; 
closure[0]=0.0; closure[1]=0.0234; closure[2]=0.077; closure[3]=1.305; 
pu[0]=0.1372; pu[1]=0.1916; pu[2]=1.7026; pu[3]=0.0426; 

// JBPT : IP3D
// float eff_IP[11];
eff_IP[0] = 0.142614; eff_IP[1] = 0.214761; eff_IP[2] = 0.276254; 
eff_IP[3] = 0.319317; eff_IP[4] = 0.363301; eff_IP[5] = 0.387019; 
eff_IP[6] = 0.409068; eff_IP[7] = 0.424023; eff_IP[8] = 0.461858; 
eff_IP[9] = 0.479511; eff_IP[10] = 0.54012; 
// float eff_stat_IP[11];
eff_stat_IP[0] = 0.00620448; eff_stat_IP[1] = 0.00677319; eff_stat_IP[2] = 0.0621496; 
eff_stat_IP[3] = 0.00798231; eff_stat_IP[4] = 0.0106629; eff_stat_IP[5] = 0.075016; 
eff_stat_IP[6] = 0.00634851; eff_stat_IP[7] = 0.0130584; eff_stat_IP[8] = 0.00990085; 
eff_stat_IP[9] = 0.011765; eff_stat_IP[10] = 0.0223062; 
// float effMC_IP[11];
effMC_IP[0] = 0.19057; effMC_IP[1] = 0.281215; effMC_IP[2] = 0.3405; 
effMC_IP[3] = 0.402025; effMC_IP[4] = 0.429765; effMC_IP[5] = 0.463133; 
effMC_IP[6] = 0.489063; effMC_IP[7] = 0.493108; effMC_IP[8] = 0.524437; 
effMC_IP[9] = 0.526817; effMC_IP[10] = 0.543327; 
// float effMC_stat_IP[11];
effMC_stat_IP[0] = 0.00257865; effMC_stat_IP[1] = 0.00202938; effMC_stat_IP[2] = 0.00192427; 
effMC_stat_IP[3] = 0.00263606; effMC_stat_IP[4] = 0.00260231; effMC_stat_IP[5] = 0.00272693; 
effMC_stat_IP[6] = 0.00295843; effMC_stat_IP[7] = 0.00284063; effMC_stat_IP[8] = 0.00375883; 
effMC_stat_IP[9] = 0.00395807; effMC_stat_IP[10] = 0.00604355; 
// float sf_IP[11];
sf_IP[0] = 0.748359; sf_IP[1] = 0.763691; sf_IP[2] = 0.811318; 
sf_IP[3] = 0.794272; sf_IP[4] = 0.845348; sf_IP[5] = 0.835654; 
sf_IP[6] = 0.836433; sf_IP[7] = 0.8599; sf_IP[8] = 0.880674; 
sf_IP[9] = 0.910203; sf_IP[10] = 0.994097; 
// float sf_stat_IP[11];
sf_stat_IP[0] = 0.034096; sf_stat_IP[1] = 0.024708; sf_stat_IP[2] = 0.182582; 
sf_stat_IP[3] = 0.020527; sf_stat_IP[4] = 0.0253334; sf_stat_IP[5] = 0.16205; 
sf_stat_IP[6] = 0.0139322; sf_stat_IP[7] = 0.0269411; sf_stat_IP[8] = 0.0199063; 
sf_stat_IP[9] = 0.0233558; sf_stat_IP[10] = 0.0425178; 
// float away_syst_IP[11];
away_syst_IP[0] = 0.0491698; away_syst_IP[1] = 0.00437268; away_syst_IP[2] = 0.00464538; 
away_syst_IP[3] = 0.0142431; away_syst_IP[4] = 0.015159; away_syst_IP[5] = 0.0149852; 
away_syst_IP[6] = 0.0258814; away_syst_IP[7] = 0.0266076; away_syst_IP[8] = 0.0181514; 
away_syst_IP[9] = 0.01876; away_syst_IP[10] = 0.0750535; 
// float mupt_syst_IP[11];
mupt_syst_IP[0] = 0.0873047; mupt_syst_IP[1] = 0.0228176; mupt_syst_IP[2] = 0.0242406; 
mupt_syst_IP[3] = 0.00507187; mupt_syst_IP[4] = 0.00539801; mupt_syst_IP[5] = 0.00533611; 
mupt_syst_IP[6] = 0.00719952; mupt_syst_IP[7] = 0.00740151; mupt_syst_IP[8] = 0.0105953; 
mupt_syst_IP[9] = 0.0109505; mupt_syst_IP[10] = 0.0434275; 
// float gluon_syst_IP[11];
gluon_syst_IP[0] = 0.00415621; gluon_syst_IP[1] = 0.0197719; gluon_syst_IP[2] = 0.021005; 
gluon_syst_IP[3] = 0.0240427; gluon_syst_IP[4] = 0.0255888; gluon_syst_IP[5] = 0.0252953; 
gluon_syst_IP[6] = 0.0248196; gluon_syst_IP[7] = 0.0255159; gluon_syst_IP[8] = 0.0210735; 
gluon_syst_IP[9] = 0.0217801; gluon_syst_IP[10] = 0.127424; 
// float pu_syst_IP[11];
pu_syst_IP[0] = 0.013251; pu_syst_IP[1] = 0.0120261; pu_syst_IP[2] = 0.0127761; 
pu_syst_IP[3] = 0.0177339; pu_syst_IP[4] = 0.0188743; pu_syst_IP[5] = 0.0186578; 
pu_syst_IP[6] = 0.00384673; pu_syst_IP[7] = 0.00395465; pu_syst_IP[8] = 0.00334768; 
pu_syst_IP[9] = 0.00345992; pu_syst_IP[10] = 0.0479588; 
// float total_syst_IP[11];
total_syst_IP[0] = 0.101157; total_syst_IP[1] = 0.032792; total_syst_IP[2] = 0.0348371; 
total_syst_IP[3] = 0.0334833; total_syst_IP[4] = 0.0356365; total_syst_IP[5] = 0.0352278; 
total_syst_IP[6] = 0.0367762; total_syst_IP[7] = 0.037808; total_syst_IP[8] = 0.0299505; 
total_syst_IP[9] = 0.0309547; total_syst_IP[10] = 0.161419; 

// JBPT : JP
eff_JP[0] = 0.209555; eff_stat_JP[0] = 0.00581037; effMC_JP[0] = 0.230027; effMC_stat_JP[0] = 0.000977944; 
sf_JP[0] = 0.911001; sf_stat_JP[0] = 0.0252595; sf_syst_JP[0] = 0.206797; sf_eror_JP[0] = 0.208334; 
mupt_syst_JP[0] = 0.017; gluon_syst_JP[0] = 0.001; pu_syst_JP[0] = 0.001; total_syst_JP[0] = 0.227; 
cor_syst_JP[0] = 0.028; inc_syst_JP[0] = 0.201; bias_syst_JP[0] = 0.1; 
eff_JP[1] = 0.284664; eff_stat_JP[1] = 0.00570982; effMC_JP[1] = 0.316202; effMC_stat_JP[1] = 0.000815843; 
sf_JP[1] = 0.900262; sf_stat_JP[1] = 0.0180575; sf_syst_JP[1] = 0.128737; sf_eror_JP[1] = 0.129998; 
mupt_syst_JP[1] = 0.017; gluon_syst_JP[1] = 0.007; pu_syst_JP[1] = 0; total_syst_JP[1] = 0.143; 
cor_syst_JP[1] = 0.013; inc_syst_JP[1] = 0.099; bias_syst_JP[1] = 0.1; 
eff_JP[2] = 0.341863; eff_stat_JP[2] = 0.00663471; effMC_JP[2] = 0.376988; effMC_stat_JP[2] = 0.000913559; 
sf_JP[2] = 0.906829; sf_stat_JP[2] = 0.0175993; sf_syst_JP[2] = 0.12877; sf_eror_JP[2] = 0.129967; 
mupt_syst_JP[2] = 0.017; gluon_syst_JP[2] = 0.008; pu_syst_JP[2] = 0; total_syst_JP[2] = 0.142; 
cor_syst_JP[2] = 0.007; inc_syst_JP[2] = 0.099; bias_syst_JP[2] = 0.1; 
eff_JP[3] = 0.403122; eff_stat_JP[3] = 0.0075073; effMC_JP[3] = 0.437618; effMC_stat_JP[3] = 0.00129264; 
sf_JP[3] = 0.921175; sf_stat_JP[3] = 0.017155; sf_syst_JP[3] = 0.116068; sf_eror_JP[3] = 0.117329; 
mupt_syst_JP[3] = 0.017; gluon_syst_JP[3] = 0.005; pu_syst_JP[3] = 0; total_syst_JP[3] = 0.126; 
cor_syst_JP[3] = 0.004; inc_syst_JP[3] = 0.075; bias_syst_JP[3] = 0.1; 
eff_JP[4] = 0.432648; eff_stat_JP[4] = 0.00783207; effMC_JP[4] = 0.467563; effMC_stat_JP[4] = 0.00130528; 
sf_JP[4] = 0.925324; sf_stat_JP[4] = 0.0167508; sf_syst_JP[4] = 0.116591; sf_eror_JP[4] = 0.117788; 
mupt_syst_JP[4] = 0.017; gluon_syst_JP[4] = 0.006; pu_syst_JP[4] = 0.001; total_syst_JP[4] = 0.126; 
cor_syst_JP[4] = 0.002; inc_syst_JP[4] = 0.075; bias_syst_JP[4] = 0.1; 
eff_JP[5] = 0.453029; eff_stat_JP[5] = 0.00722991; effMC_JP[5] = 0.489487; effMC_stat_JP[5] = 0.00135817; 
sf_JP[5] = 0.925516; sf_stat_JP[5] = 0.0147704; sf_syst_JP[5] = 0.116615; sf_eror_JP[5] = 0.117547; 
mupt_syst_JP[5] = 0.017; gluon_syst_JP[5] = 0.007; pu_syst_JP[5] = 0; total_syst_JP[5] = 0.126; 
cor_syst_JP[5] = 0.002; inc_syst_JP[5] = 0.075; bias_syst_JP[5] = 0.1; 
eff_JP[6] = 0.508532; eff_stat_JP[6] = 0.00885768; effMC_JP[6] = 0.527855; effMC_stat_JP[6] = 0.00130025; 
sf_JP[6] = 0.963387; sf_stat_JP[6] = 0.0167805; sf_syst_JP[6] = 0.0799611; sf_eror_JP[6] = 0.0817029; 
mupt_syst_JP[6] = 0.013; gluon_syst_JP[6] = 0.006; pu_syst_JP[6] = 0; total_syst_JP[6] = 0.083; 
cor_syst_JP[6] = 0.001; inc_syst_JP[6] = 0.056; bias_syst_JP[6] = 0.06; 
eff_JP[7] = 0.522915; eff_stat_JP[7] = 0.00792334; effMC_JP[7] = 0.543575; effMC_stat_JP[7] = 0.00130063; 
sf_JP[7] = 0.961996; sf_stat_JP[7] = 0.0145764; sf_syst_JP[7] = 0.0798457; sf_eror_JP[7] = 0.0811653; 
mupt_syst_JP[7] = 0.013; gluon_syst_JP[7] = 0.003; pu_syst_JP[7] = 0.001; total_syst_JP[7] = 0.083; 
cor_syst_JP[7] = 0.001; inc_syst_JP[7] = 0.056; bias_syst_JP[7] = 0.06; 
eff_JP[8] = 0.54764; eff_stat_JP[8] = 0.00808728; effMC_JP[8] = 0.569476; effMC_stat_JP[8] = 0.00124516; 
sf_JP[8] = 0.961658; sf_stat_JP[8] = 0.0142012; sf_syst_JP[8] = 0.0798176; sf_eror_JP[8] = 0.0810711; 
mupt_syst_JP[8] = 0.013; gluon_syst_JP[8] = 0.002; pu_syst_JP[8] = 0; total_syst_JP[8] = 0.083; 
cor_syst_JP[8] = 0.001; inc_syst_JP[8] = 0.056; bias_syst_JP[8] = 0.06; 
eff_JP[9] = 0.561938; eff_stat_JP[9] = 0.00871858; effMC_JP[9] = 0.58; effMC_stat_JP[9] = 0.00151095; 
sf_JP[9] = 0.968858; sf_stat_JP[9] = 0.015032; sf_syst_JP[9] = 0.0881661; sf_eror_JP[9] = 0.0894384; 
mupt_syst_JP[9] = 0.013; gluon_syst_JP[9] = 0.008; pu_syst_JP[9] = 0.001; total_syst_JP[9] = 0.091; 
cor_syst_JP[9] = 0.001; inc_syst_JP[9] = 0.066; bias_syst_JP[9] = 0.06; 
eff_JP[10] = 0.546425; eff_stat_JP[10] = 0.0100723; effMC_JP[10] = 0.58462; effMC_stat_JP[10] = 0.00268844; 
sf_JP[10] = 0.93467; sf_stat_JP[10] = 0.0172288; sf_syst_JP[10] = 0.113095; sf_eror_JP[10] = 0.1144; 
mupt_syst_JP[10] = 0.01; gluon_syst_JP[10] = 0.016; pu_syst_JP[10] = 0.001; total_syst_JP[10] = 0.121; 
cor_syst_JP[10] = 0; inc_syst_JP[10] = 0.066; bias_syst_JP[10] = 0.1; 
eff_JP[11] = 0.559882; eff_stat_JP[11] = 0.0123053; effMC_JP[11] = 0.582492; effMC_stat_JP[11] = 0.00452636; 
sf_JP[11] = 0.961187; sf_stat_JP[11] = 0.0211253; sf_syst_JP[11] = 0.117265; sf_eror_JP[11] = 0.119153; 
mupt_syst_JP[11] = 0.01; gluon_syst_JP[11] = 0.02; pu_syst_JP[11] = 0.005; total_syst_JP[11] = 0.122; 
cor_syst_JP[11] = 0; inc_syst_JP[11] = 0.066; bias_syst_JP[11] = 0.1; 
eff_JP[12] = 0.575243; eff_stat_JP[12] = 0.0183368; effMC_JP[12] = 0.600982; effMC_stat_JP[12] = 0.00748033; 
sf_JP[12] = 0.957174; sf_stat_JP[12] = 0.0305114; sf_syst_JP[12] = 0.128261; sf_eror_JP[12] = 0.131841; 
mupt_syst_JP[12] = 0.01; gluon_syst_JP[12] = 0.036; pu_syst_JP[12] = 0.007; total_syst_JP[12] = 0.134; 
cor_syst_JP[12] = 0; inc_syst_JP[12] = 0.08; bias_syst_JP[12] = 0.1; 
eff_JP[13] = 0.583626; eff_stat_JP[13] = 0.0219979; effMC_JP[13] = 0.592463; effMC_stat_JP[13] = 0.00647857; 
sf_JP[13] = 0.985084; sf_stat_JP[13] = 0.0371295; sf_syst_JP[13] = 0.130031; sf_eror_JP[13] = 0.135228; 
mupt_syst_JP[13] = 0.01; gluon_syst_JP[13] = 0.024; pu_syst_JP[13] = 0.017; total_syst_JP[13] = 0.132; 
cor_syst_JP[13] = 0; inc_syst_JP[13] = 0.08; bias_syst_JP[13] = 0.1; 
eff_JP[14] = 0.554872; eff_stat_JP[14] = 0.0286099; effMC_JP[14] = 0.588701; effMC_stat_JP[14] = 0.00615905; 
sf_JP[14] = 0.942533; sf_stat_JP[14] = 0.0485984; sf_syst_JP[14] = 0.123472; sf_eror_JP[14] = 0.132692; 
mupt_syst_JP[14] = 0.01; gluon_syst_JP[14] = 0.021; pu_syst_JP[14] = 0.01; total_syst_JP[14] = 0.131; 
cor_syst_JP[14] = 0.007; inc_syst_JP[14] = 0.08; bias_syst_JP[14] = 0.1; 
}

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  int kk;
  for (int k = 0; k <= 7; k++) {
    if ( k <=  2 )      kk = 0; // 20-50
    else if ( k <=  4 ) kk = 1; // 50-80
    else if ( k <=  6 ) kk = 2; // 80-120
    else                kk = 3; // >120
    mupt_syst_SY[k]  =      mu[kk]/100.;
    gluon_syst_SY[k] =  gsplit[kk]/100.;
    pu_syst_SY[k]    =      pu[kk]/100.;
    away_syst_SY[k]  =    away[kk]/100.;
    ptrel_syst_SY[k] =   ptrel[kk]/100.;
    clos_syst_SY[k]  = closure[kk]/100.;
    total_syst_SY[k] = combError(mupt_syst_SY[k],gluon_syst_SY[k],pu_syst_SY[k],away_syst_SY[k],ptrel_syst_SY[k],clos_syst_SY[k]);
    sf_syst_SY[k] = total_syst_SY[k] * sf_SY[k];
    sf_eror_SY[k] = combError(sf_stat_SY[k],sf_syst_SY[k]);
// cout << k << " " << sf_SY[k] << sf_stat_SY[k] << sf_syst_SY[k] 
//                  << sf_eror_SY[k] << endl;
  }

  for (int k = 0; k <= 8; k++) {
    sf_syst_PT[k] = total_syst_PT[k];
    sf_eror_PT[k] = combError(sf_stat_PT[k],sf_syst_PT[k]);
  }

  for (int k = 0; k <= 5; k++) { 
       eff_IP[k] = 0.; 
  eff_stat_IP[k] = 0.;
     effMC_IP[k] = 0.; 
effMC_stat_IP[k] = 0.;
        sf_IP[k] = 0.; 
   sf_stat_IP[k] = 0.;
   sf_syst_IP[k] = 0.;
   sf_eror_IP[k] = 0.;
 away_syst_IP[k] = 0.; 
 mupt_syst_IP[k] = 0.; 
gluon_syst_IP[k] = 0.;
   pu_syst_IP[k] = 0.; 
total_syst_IP[k] = 0.; 
  }
  for (int k = 6; k <= 10; k++) {
   sf_syst_IP[k] = total_syst_IP[k];
   sf_eror_IP[k] = combError(sf_stat_IP[k],sf_syst_IP[k]);
  }

  for (int k = 0; k <= 15; k++) exstat[k] = 0.0001;

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

     pad1->cd();

  TH2F* histo = new TH2F("histo","",65,20.,670.,100,0.7,1.3);

  histo->Draw(""); 
  histo->SetLabelSize(0.065, "XYZ");
  histo->SetTitleSize(0.07, "XYZ"); 
  histo->SetLabelFont(42, "XYZ"); 
  histo->SetTitleFont(42, "XYZ");
  histo->GetXaxis()->SetTitle("p_{T} (GeV)");
  histo->GetYaxis()->SetTitle("Data/Sim. b-tag SF");
  histo->SetTitleOffset(0.8,"X");
  histo->SetTitleOffset(0.8,"Y");
  histo->SetNdivisions(509, "XYZ");

  TGraphErrors* sf1stat = new TGraphErrors(8,x1,sf_PT,exstat,sf_stat_PT);
  sf1stat->Draw("P"); 
  sf1stat->SetFillStyle(3005);
  sf1stat->SetFillColor(kBlue);
  sf1stat->SetMarkerColor(kBlue);
  sf1stat->SetLineColor(kBlue);
  sf1stat->SetLineStyle(1);
  sf1stat->SetMarkerStyle(21);
  sf1stat->SetMarkerSize(1.3);
  sf1stat->SetLineWidth(4.);
  TGraphErrors* sf1 = new TGraphErrors(8,x1,sf_PT,ex1,sf_eror_PT);
  sf1->Draw("P"); 
  sf1->SetFillStyle(3005);
  sf1->SetFillColor(kBlue);
  sf1->SetMarkerColor(kBlue);
  sf1->SetLineColor(kBlue);
  sf1->SetLineStyle(1);
  sf1->SetMarkerStyle(21);
  sf1->SetMarkerSize(1.3);
  sf1->SetLineWidth(2.);

  TGraphErrors* sf2stat = new TGraphErrors(7,x2,sf_SY,exstat,sf_stat_SY);
  sf2stat->Draw("P"); 
  sf2stat->SetMarkerColor(kGreen);
  sf2stat->SetLineColor(kGreen);
  sf2stat->SetLineStyle(1);
  sf2stat->SetMarkerStyle(22);
  sf2stat->SetMarkerSize(1.3);
  sf2stat->SetLineWidth(4.);
  TGraphErrors* sf2 = new TGraphErrors(7,x2,sf_SY,ex2,sf_eror_SY);
  sf2->Draw("P"); 
  sf2->SetMarkerColor(kGreen);
  sf2->SetLineColor(kGreen);
  sf2->SetLineStyle(1);
  sf2->SetMarkerStyle(22);
  sf2->SetMarkerSize(1.3);
  sf2->SetLineWidth(2.);

  TGraphErrors* sf3stat = new TGraphErrors(11,x3,sf_IP,exstat,sf_stat_IP);
  sf3stat->Draw("P"); 
  sf3stat->SetMarkerColor(kViolet);
  sf3stat->SetLineColor(kViolet);
  sf3stat->SetLineStyle(1);
  sf3stat->SetMarkerStyle(23);
  sf3stat->SetMarkerSize(1.3);
  sf3stat->SetLineWidth(4.);
  TGraphErrors* sf3 = new TGraphErrors(11,x3,sf_IP,ex1,sf_eror_IP);
  sf3->Draw("P"); 
  sf3->SetMarkerColor(kViolet);
  sf3->SetLineColor(kViolet);
  sf3->SetLineStyle(1);
  sf3->SetMarkerStyle(23);
  sf3->SetMarkerSize(1.3);
  sf3->SetLineWidth(2.);

  TGraphErrors* sf4stat = new TGraphErrors(15,x,sf_JP,exstat,sf_stat_JP);
  sf4stat->Draw("P"); 
  sf4stat->SetMarkerColor(kRed);
  sf4stat->SetLineColor(kRed);
  sf4stat->SetLineStyle(1);
  sf4stat->SetMarkerStyle(20);
  sf4stat->SetMarkerSize(1.3);
  sf4stat->SetLineWidth(4.);
  TGraphErrors* sf4 = new TGraphErrors(15,x,sf_JP,ex,sf_eror_JP);
  sf4->Draw("P"); 
  sf4->SetMarkerColor(kRed);
  sf4->SetLineColor(kRed);
  sf4->SetLineStyle(1);
  sf4->SetMarkerStyle(20);
  sf4->SetMarkerSize(1.3);
  sf4->SetLineWidth(2.);

//   TLegend* leg = new TLegend(0.53,0.20,0.75,0.40);
  TLegend* leg = new TLegend(0.23,0.65,0.45,0.93);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->SetTextSize(0.06);   
  leg->SetHeader(title);
  leg->AddEntry(sf1,"PtRel","PL");
  leg->AddEntry(sf2,"System8","PL");
  leg->AddEntry(sf3,"IP3D","PL");
  leg->AddEntry(sf4,"JP","PL");
  leg->Draw();

  //include the official CMS label
//   TPaveText* pt = new TPaveText(0.55,0.21,0.90,0.31,"brNDC");   
  TPaveText* pt = new TPaveText(0.55,0.85,0.90,0.95,"brNDC");   
  pt->SetBorderSize(0);   
  pt->SetFillStyle(0);  
  pt->SetTextAlign(13);   
  pt->SetTextFont(22);   
  pt->SetTextSize(0.06);   
//   pt->AddText("CMS at #sqrt{s} = 7 TeV, 4.7 fb^{-1}");   
  pt->AddText("CMS prelim. at 7 TeV, 4.7 fb^{-1}");   
  pt->Draw(""); 

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  pair<double,double> combResult;

  sf_stat_SY[4] *= TMath::Sqrt(2);
  for (int k=7; k>=4; k--) {
    sf_SY[k+1]         = sf_SY[k];
    sf_stat_SY[k+1]    = sf_stat_SY[k];
    sf_syst_SY[k+1]    = sf_syst_SY[k];
    away_syst_SY[k+1]  = away_syst_SY[k];
    mupt_syst_SY[k+1]  = mupt_syst_SY[k];
    gluon_syst_SY[k+1] = gluon_syst_SY[k];
    pu_syst_SY[k+1]    = pu_syst_SY[k];
    ptrel_syst_SY[k+1] = ptrel_syst_SY[k];
    clos_syst_SY[k+1]  = clos_syst_SY[k];
  }

  for (int k=12; k<=14; k++) {
    sf_IP[k]      = sf_IP[10];
    sf_stat_IP[k] = sf_stat_IP[10]*TMath::Sqrt(3);
    sf_syst_IP[k] = sf_syst_IP[10];
    away_syst_IP[k]  = away_syst_IP[10];
    mupt_syst_IP[k]  = mupt_syst_IP[10];
    gluon_syst_IP[k] = gluon_syst_IP[10];
    pu_syst_IP[k]    = pu_syst_IP[10];
  }
  sf_stat_IP[9] *= TMath::Sqrt(3);
  for (int k=10; k<=11; k++) {
    sf_IP[k]      = sf_IP[9];
    sf_stat_IP[k] = sf_stat_IP[9];
    sf_syst_IP[k] = sf_syst_IP[9];
    away_syst_IP[k]  = away_syst_IP[9];
    mupt_syst_IP[k]  = mupt_syst_IP[9];
    gluon_syst_IP[k] = gluon_syst_IP[9];
    pu_syst_IP[k]    = pu_syst_IP[9];
  }

//$$  for (int k = 0; k <= 14; k++) {
  for (int k = 0; k <= 6; k++) {
    if ( k == 0 ) {
      if ( title == "TCHPT" || title == "SSVHPT" || title == "CSVT"
                            || title == "JPT" || title == "JBPT" ) { 
	combResult = combine(k,title.Data(),true,true,false,false);
      }
      else { 
	combResult = combine(k,title.Data(),true,true,false,true);
      }
    }
    else if ( k <= 5 ) {
	combResult = combine(k,title.Data(),true,true,false,true);
    }
    else if ( k <= 7 ) {
      if ( title == "TCHPT" || title == "SSVHPT" || title == "CSVT"
                            || title == "JPT" || title == "JBPT" ) { 
	combResult = combine(k,title.Data(),true,true,false,true);
      } 
      else  {
	combResult = combine(k,title.Data(),true,true,true,true);
cout << " PT " << sf_PT[k] << " +_ " << sf_stat_PT[k] << " +_ " << sf_syst_PT[k] << endl;
cout << " SY " << sf_SY[k] << " +_ " << sf_stat_SY[k] << " +_ " << sf_syst_SY[k] << endl;
cout << " IP " << sf_IP[k] << " +_ " << sf_stat_IP[k] << " +_ " << sf_syst_IP[k] << endl;
cout << " JP " << sf_JP[k] << " +_ " << sf_stat_JP[k] << " +_ " << sf_syst_JP[k] << endl;
      } 
    }
    else if ( k == 8 ) {
	combResult = combine(k,title.Data(),false,false,true,true);
    }
    else if ( k <= 11 ) {
	combResult = combine(k,title.Data(),false,false,true,true);
    }
    else {
      if ( title == "TCHPT" || title == "SSVHPT" || title == "CSVT"
                            || title == "JPT" || title == "JBPT" ) { 
	combResult = combine(k,title.Data(),false,false,false,true);
      }
      else { 
	combResult = combine(k,title.Data(),false,false,true,true);
      }
    }
    sf[k] = combResult.first;
    sf_eror[k] = combResult.second;
    sf_stat[k] = 0.;
cout << k << " " << sf[k] <<  " +_ " << sf_eror[k] << endl;
  }

  TGraphErrors* sf0 = new TGraphErrors(15,x,sf,ex,sf_eror);
  sf0->Draw("e2"); 
  sf0->SetFillStyle(3005);
  sf0->SetFillColor(kGray+3);
  sf3->Draw("P"); 
  sf4->Draw("P"); 
  sf2->Draw("P"); 
  sf1->Draw("P"); 

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

//   pad2->cd();
// 
// //$$  TF1* fun1 = new TF1("fun1","[0]*(1.+[1]*x)/(1.+[2]*x)",20.,670.);
// //$$  TF1* fun1 = new TF1("fun1","[0]+[1]*x",20.,670.);
//   TF1* fun1 = new TF1("fun1","[0]+(1.+[1]*x)/(1.+[2]*x)",20.,670.);
// 
//   histo = new TH2F("histo","",65,20.,670.,100,0.7,1.3);
// 
//   histo->Draw(""); 
//   histo->SetLabelSize(0.065, "XYZ");
//   histo->SetTitleSize(0.07, "XYZ"); 
//   histo->SetLabelFont(42, "XYZ"); 
//   histo->SetTitleFont(42, "XYZ");
//   histo->GetXaxis()->SetTitle("p_{T} (GeV)");
//   histo->GetYaxis()->SetTitle("Data/Sim. b-tag SF");
//   histo->SetTitleOffset(0.8,"X");
//   histo->SetTitleOffset(0.8,"Y");
//   histo->SetNdivisions(509, "XYZ");
// 
// //   TGraphErrors* sf0stat = new TGraphErrors(15,x,sf,ex,sf_stat);
// //   sf0stat->Draw("P"); 
// //   sf0stat->SetMarkerColor(kBlack);
// //   sf0stat->SetLineColor(kBlack);
// //   sf0stat->SetLineStyle(1);
// //   sf0stat->SetMarkerStyle(20);
// //   sf0stat->SetMarkerSize(1.3);
// //   sf0stat->SetLineWidth(2.);
// //   TGraphErrors* sf0 = new TGraphErrors(15,x,sf,ex,sf_eror);
//   sf0->Draw("e2"); 
// //   sf0->SetMarkerColor(kBlack);
// //   sf0->SetLineColor(kBlack);
// //   sf0->SetLineStyle(1);
// //   sf0->SetMarkerStyle(20);
// //   sf0->SetMarkerSize(1.3);
// //   sf0->SetLineWidth(2.);
//   sf0->SetFillStyle(3005);
//   sf0->SetFillColor(kBlack);
// 
//   fun1->SetLineColor(kBlack);
//   fun1->SetLineWidth(2.);
//   fun1->SetLineStyle(1);
//   sf0->Fit("fun1","0"); 
//   sf0->Fit("fun1","rve0"); 
// //$$  sf0->Fit("fun1","rvee0"); 
//   fun1->Draw("same"); 
// 
//   float fun_val[NBINS], fun_err[NBINS], fun_sys[NBINS];
//   cout << endl;
//   for (int k=0; k<=14; k++) {
//     fun_val[k] = fun1->Integral(x[k]-1.,x[k]+1.) / 2.;
//     fun_err[k] = fun1->IntegralError(x[k]-ex[k],x[k]+ex[k]) / ex[k];
//     cout << k << " " << fun_val[k] << " " << fun_err[k] << endl;
//     fun_sys[k] = sf_eror[k] * fun_val[k] / sf[k];
//   }
//   cout << endl;
//   TGraphErrors* fun0 = new TGraphErrors(15,x,fun_val,ex,fun_sys);
// //   fun0->Draw("Pe2"); 
//   fun0->Draw("P"); 
//   fun0->SetMarkerColor(kRed);
//   fun0->SetLineColor(kRed);
//   fun0->SetLineStyle(1);
//   fun0->SetMarkerStyle(24);
//   fun0->SetMarkerSize(0.001);
//   fun0->SetLineWidth(2.);
// //   fun0->SetFillStyle(3004);
// //   fun0->SetFillColor(kRed);
// 
//   leg = new TLegend(0.23,0.65,0.45,0.93);
//   leg->SetBorderSize(0);
//   leg->SetFillColor(kWhite);
//   leg->SetTextFont(42);
//   leg->SetTextSize(0.06);   
//   leg->SetHeader(title);
//   leg->AddEntry(sf0,"weighted average","PF");
//   leg->AddEntry(fun1,"fit","L");
//   leg->AddEntry(fun0,"fit #pm (stat #oplus syst)","LE");
//   leg->Draw();
// 
//   //include the official CMS label
// //   pt = new TPaveText(0.55,0.21,0.90,0.31,"brNDC");   
//   pt = new TPaveText(0.55,0.85,0.90,0.95,"brNDC");   
//   pt->SetBorderSize(0);   
//   pt->SetFillStyle(0);  
//   pt->SetTextAlign(13);   
//   pt->SetTextFont(22);   
//   pt->SetTextSize(0.06);   
// //   pt->AddText("CMS at #sqrt{s} = 7 TeV, 4.7 fb^{-1}");   
//   pt->AddText("CMS prelim. at 7 TeV, 4.7 fb^{-1}");   
//   pt->Draw(""); 
// 
//   float sf_ttbar = 0., sf_ttbar_eror = 0.; 
//   for (int k = 1; k <= 13; k++) {
//     sf_ttbar += fun_val[k]*ttbar_pt[k]; 
//     sf_ttbar_eror += sf_eror[k]*ttbar_pt[k] * fun_val[k]/sf[k]; 
// //   cout << " k " << k << " " << sf[k] << " " << ttbar_pt[k] << endl;
// //   cout << " k " << k << " " << sf_ttbar << " " << sf_ttbar_stat << " " << sf_ttbar_syst << endl;
// //   cout << endl;
//   }
// //   sf_ttbar_eror = TMath::Sqrt( sf_ttbar_eror );
// 
//   cout << endl;
// //   cout << " ttbar-like SF = " << sf_ttbar << " +_ " << sf_ttbar_stat << " +_ " << sf_ttbar_syst
// //        << " = " << sf_ttbar << " +_ " << sf_ttbar_eror << endl;
//   cout << " ttbar-like SF = " << sf_ttbar << " +_ " << sf_ttbar_eror << endl;
//   cout << endl;

// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  c1->Update();

  TString oname(title);
  oname.Append(".pdf");
  c1->SaveAs(oname);
}

