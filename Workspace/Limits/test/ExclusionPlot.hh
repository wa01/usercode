#ifndef MAKE_EXCLUSIONPLOT_HH
#define MAKE_EXCLUSIONPLOT_HH

#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TFile.h"
#include "TSpline.h"
#include "TGraphErrors.h"

#include "GridDump.C" //tan beta 40 squark, gluino lines

#include <vector>
#include <cassert>

void ExclusionPlot(TFile*);  

void CommandMSUGRA(TString plotName,Int_t tanBeta, Bool_t plotLO, Bool_t tb40_plotExpected);

void setPlottingStyle(TH1F& hsig);

 
//a little plotting routine to calculate the NLO cross-section
TH2F* sysPlot(TString mSuGraFile){
  //read In mSuGra Histo
  TFile* f = new TFile(mSuGraFile);
  TDirectory* dir = (TDirectory*)f->Get("mSuGraScan_beforeAll");
  TDirectory* dir2 = (TDirectory*)f->Get("mSuGraScan_350");
  
  TH2F* gg = (TH2F*)dir2->Get("m0_m12_gg_0");
  TH2F* gg_noweight = (TH2F*)dir->Get("m0_m12_gg_5");
  TH2F* sb = (TH2F*)dir2->Get("m0_m12_sb_0");
  TH2F* sb_noweight = (TH2F*)dir->Get("m0_m12_sb_5");
  TH2F* ss = (TH2F*)dir2->Get("m0_m12_ss_0");
  TH2F* ss_noweight = (TH2F*)dir->Get("m0_m12_ss_5");
  TH2F* sg = (TH2F*)dir2->Get("m0_m12_sg_0");
  TH2F* sg_noweight = (TH2F*)dir->Get("m0_m12_sg_5");
  TH2F* ll = (TH2F*)dir2->Get("m0_m12_ll_0");
  TH2F* ll_noweight = (TH2F*)dir->Get("m0_m12_ll_5");
  TH2F* nn = (TH2F*)dir2->Get("m0_m12_nn_0");
  TH2F* nn_noweight = (TH2F*)dir->Get("m0_m12_nn_5");
  TH2F* ns = (TH2F*)dir2->Get("m0_m12_ns_0");
  TH2F* ns_noweight = (TH2F*)dir->Get("m0_m12_ns_5");
  TH2F* ng = (TH2F*)dir2->Get("m0_mg12_ng_0");                                                               
  TH2F* ng_noweight = (TH2F*)dir->Get("m0_m12_ng_5");
  TH2F* bb = (TH2F*)dir->Get("m0_m12_bb_0");
  TH2F* bb_noweight = (TH2F*)dir->Get("m0_m12_bb_5");
  TH2F* tb = (TH2F*)dir->Get("m0_m12_tb_0");
  TH2F* tb_noweight = (TH2F*)dir->Get("m0_m12_tb_5");
  
  gg->Divide(gg_noweight);
  sb->Divide(sg_noweight);
  ss->Divide(ss_noweight);
  sg->Divide(sg_noweight);
  ll->Divide(ll_noweight);
  nn->Divide(nn_noweight);
  //ng->Divide(ng_noweight);
  bb->Divide(bb_noweight);
  tb->Divide(tb_noweight);
  ns->Divide(ns_noweight);

  TH2F* all = (TH2F*)gg->Clone();
  all->Add(sb);
  all->Add(ss);
  all->Add(sg);
  all->Add(ll);
  all->Add(nn);
  //all->Add(ng);
  all->Add(bb);
  all->Add(tb);
  all->Add(ns);

  all->Scale(100);
  
  return all;
}



TGraph* set_lep_ch(Int_t tanBeta);
TGraph* set_lep_ch_tanBeta10();
TGraph* set_lep_ch_tanBeta40();
TGraph* set_lep_sl(Int_t tanBeta);//slepton curve
TGraph* set_tev_sg_cdf(Int_t tanBeta);//squark gluino cdf
TGraph* set_tev_sg_d0(Int_t tanBeta);//squark gluino d0
//TGraph* set_tev_tlp_cdf(Int_t tanBeta);//trilepton cdf
//TGraph* set_tev_tlp_d0(Int_t tanBeta);//trilepton d0
TGraph* set_tev_stau(Int_t tanBeta);//stau 
TGraph* set_NoEWSB(Int_t tanBeta);//stau 

TGraph* set_sneutrino_d0_1(Int_t tanBeta);
TGraph* set_sneutrino_d0_2(Int_t tanBeta);


TF1* constant_squark(int tanBeta,int i);
TF1* constant_gluino(int tanBeta,int i);

TLatex* constant_squark_text(Int_t it,TF1& lnsq,Int_t tanBeta);
TLatex* constant_gluino_text(Int_t it,TF1& lngl);

TGraph* constant_mass(double mass,TGraph2D* massGrid);
TLatex* constant_squark_text_tanBeta40(int mass,TGraph* massLine);
TLatex* constant_gluino_text_tanBeta40(int mass,TGraph* massLine);


TLegend* makeStauLegend(Double_t txtsz,Int_t tanBeta_);
TLegend* makeNoEWSBLegend(Double_t txtsz,Int_t tanBeta_);
TLegend* makeExpLegend(TGraph& sg_gr, TGraph& sgd_gr,TGraph& ch_gr,TGraph& sl_gr,TGraph&,Double_t txtsz,Int_t tanbeta);



  
  vector<TH1F*> exclusionPlots;

  TH2F* yieldPlot(TString mSuGraFile,TString mSuGraDir, TString mSuGraHist){
  //read In mSuGra Histo
  TFile* f = new TFile(mSuGraFile);
  TDirectory* dir = (TDirectory*)f->Get(mSuGraDir);
  
  TH2F* hnev = (TH2F*)dir->Get(mSuGraHist);

  return hnev;
}

TH2F* yieldPlot_nosubdir(TString mSuGraFile,TString mSuGraHist){
  //read In mSuGra Histo
  TFile* f = new TFile(mSuGraFile);
  //  TDirectory* dir = (TDirectory*)f->Get(mSuGraDir);
  
  TH2F* hnev = (TH2F*)f->Get(mSuGraHist);

 

  return hnev;
}

TH1F* yieldPlot(TString mSuGraFile, TString mSuGraHist){
  //read In mSuGra Histo
  TFile* f = new TFile(mSuGraFile);
  // TDirectory* dir = (TDirectory*)f->Get(mSuGraDir);
  
  TH1F* hnev = (TH1F*)f->Get(mSuGraHist);

 

  return hnev;
}



TH2F* getHisto( TString path,
	       TString nameHist,
	       TString nameFile ) {
  TString name = path + nameFile;
  TFile* file =  new TFile(name);
  //  TDirectory* dir = (TDirectory*)file->Get(Dirname);
  TH2F* hist = (TH2F*)file->Get(nameHist);
  if (!hist) {
    std::cout << " name: " << nameHist
	      << " file: " << nameFile
	      << std::endl;
//     abort();

  }
  assert(hist);
  hist->SetLineWidth(1);

  hist->GetXaxis()->SetTitleSize(0.055);
  hist->GetYaxis()->SetTitleSize(0.055);
  hist->GetXaxis()->SetLabelSize(0.05);
  hist->GetYaxis()->SetLabelSize(0.05);
  hist->SetStats(kFALSE);
  return hist;
}







TGraphErrors* getLO_tanBeta3(){



  Int_t nl = 9;
  Double_t xl[9];
  Double_t yl[9];
  Double_t exl[9];
  Double_t eyl[9];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  
 
   xl[0] = 0;
  yl[0] = 265;
  xl[1] = 100;
  yl[1] = 258;
  xl[2] = 200;
  yl[2] = 250;
  xl[3] = 250;
  yl[3] = 240;
  xl[4] = 300;
  yl[4] = 210;
  xl[5] = 340;
  yl[5] = 177;
  xl[6] = 400;
  yl[6] = 140;
  xl[7] = 450;
  yl[7] = 120;
  xl[8] = 520;
  yl[8] =100;

  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kGreen+2);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kGreen+2);
  s->SetLineStyle(4);
  s->SetLineWidth(3);

  return gr1;
}

TGraphErrors* getLO_tanBeta10(){



  Int_t nl = 10;
  Double_t xl[10];
  Double_t yl[10];
  Double_t exl[10];
  Double_t eyl[10];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  
  xl[0] = 0;
  yl[0] = 270;
  xl[1] = 100;
  yl[1] = 260;
  xl[2] = 200;
  yl[2] = 250;
  xl[3] = 250;
  yl[3] = 240;
  xl[4] = 300;
  yl[4] = 210;
  xl[5] = 350;
  yl[5] = 174;
  xl[6] = 400;
  yl[6] = 147;
  xl[7] = 450;
  yl[7] = 127;
  xl[8] = 500;
  yl[8] =115;
  xl[9] = 520;
  yl[9] = 112;
  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kGreen+2);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kGreen+2);
  s->SetLineStyle(4);
  s->SetLineWidth(3);
  

  return gr1;



}

TGraphErrors* getLO_tanBeta50(){



  Int_t nl = 10;
  Double_t xl[10];
  Double_t yl[10];
  Double_t exl[10];
  Double_t eyl[10];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  

  xl[0] = 200;
  yl[0] = 239;
  xl[1] = 210;
  yl[1] = 249;
  xl[2] = 229;
  yl[2] = 260;
  xl[3] = 250;
  yl[3] = 245;
  xl[4] = 300;
  yl[4] = 210;
  xl[5] = 350;
  yl[5] = 180;
  xl[6] = 400;
  yl[6] = 160;
  xl[7] = 450;
  yl[7] = 150;
  xl[8] = 500;
  yl[8] =140;
  xl[9] = 520;
  yl[9] = 137;
  
  
  
  
  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kGreen+2);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kGreen+2);
  s->SetLineStyle(4);
  s->SetLineWidth(3);
  

  return gr1;



}



TGraphErrors* getExpected_NLO_tanBeta3(){

 Int_t nl = 11;
  Double_t xl[11];
  Double_t yl[11];
  Double_t exl[11];
  Double_t eyl[11];
  
    xl[0] = 0;
  yl[0] = 283;
  xl[1] = 100;
  yl[1] = 280;
  xl[2] = 150;
  yl[2] = 279;
  xl[3] = 200;
  yl[3] = 275;
  xl[4] = 250;
  yl[4] = 270;
  xl[5] = 300;
  yl[5] = 255;
  xl[6] = 350;
  yl[6] = 225;
  xl[7] = 400;
  yl[7] = 195;
  xl[8] = 450;
  yl[8] = 175;
  xl[9] = 500;
  yl[9] = 155;
  xl[10] = 550;
  yl[10] = 150;
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kBlue);
  s->SetLineStyle(2);
  s->SetLineWidth(3);
  

  return gr1;





}

TGraphErrors* getExpected_NLO_tanBeta10(){

 Int_t nl = 11;
  Double_t xl[11];
  Double_t yl[11];
  Double_t exl[11];
  Double_t eyl[11];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  
   xl[0] = 0;
  yl[0] = 283;
  xl[1] = 100;
  yl[1] = 280;
  xl[2] = 150;
  yl[2] = 279;
  xl[3] = 200;
  yl[3] = 275;
  xl[4] = 250;
  yl[4] = 270;
  xl[5] = 300;
  yl[5] = 255;
  xl[6] = 350;
  yl[6] = 225;
  xl[7] = 400;
  yl[7] = 195;
  xl[8] = 450;
  yl[8] = 175;
  xl[9] = 500;
  yl[9] = 165;
  xl[10] = 550;
  yl[10] = 150;

  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kBlue);
  s->SetLineStyle(2);
  s->SetLineWidth(3);
  

  return gr1;





}


TGraphErrors* getExpected_NLO_tanBeta50(){

 Int_t nl = 10;
  Double_t xl[10];
  Double_t yl[10];
  Double_t exl[10];
  Double_t eyl[10];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  
   xl[0] = 200;
  yl[0] = 287;
  xl[1] = 220;
  yl[1] = 287;
  xl[2] = 245;
  yl[2] = 287;
  xl[3] = 270;
  yl[3] = 265;
  xl[4] = 300;
  yl[4] = 245;
  xl[5] = 350;
  yl[5] = 222;
  xl[6] = 400;
  yl[6] = 197;
  xl[7] = 450;
  yl[7] = 180;
  xl[8] = 500;
  yl[8] = 168;
  xl[9] = 550;
  yl[9] = 145;



  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);

  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kBlue);
  s->SetLineStyle(2);
  s->SetLineWidth(3);
  

  return gr1;

}


TGraphErrors* getObserved_NLO_tanBeta3(){

  Int_t nl = 11;
  Double_t xl[11];
  Double_t yl[11];
  Double_t exl[11];
  Double_t eyl[11];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  

  xl[0] = 0;
  yl[0] = 274;
  xl[1] = 100;
  yl[1] = 270;
  xl[2] = 150;
  yl[2] = 268;
  xl[3] = 200;
  yl[3] = 265;
  xl[4] = 250;
  yl[4] = 255;
  xl[5] = 300;
  yl[5] = 230;
  xl[6] = 350;
  yl[6] = 195;
  xl[7] = 400;
  yl[7] = 160;
  xl[8] = 450;
  yl[8] = 140;
  xl[9] = 480;
  yl[9] = 130;
  xl[10] = 530;
  yl[10] = 120;
 
  
  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kRed);
  // s->SetLineStyle(2);
  s->SetLineWidth(3);
  

  return gr1;





}



TGraphErrors* getObserved_NLO_tanBeta10(){

 Int_t nl = 11;
  Double_t xl[11];
  Double_t yl[11];
  Double_t exl[11];
  Double_t eyl[11];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  
   xl[0] = 0;
  yl[0] = 278;
  xl[1] = 100;
  yl[1] = 270;
  xl[2] = 150;
  yl[2] = 267;
  xl[3] = 200;
  yl[3] = 262;
  xl[4] = 250;
  yl[4] = 250;
  xl[5] = 300;
  yl[5] = 225;
  xl[6] = 350;
  yl[6] = 192;
  xl[7] = 400;
  yl[7] = 163;
  xl[8] = 450;
  yl[8] = 148;
  xl[9] = 500;
  yl[9] = 140;
  xl[10] = 520;
  yl[10] = 137;
  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kRed);
  //  s->SetLineStyle(2);
  s->SetLineWidth(3);
  

  return gr1;

}


TGraphErrors* getObserved_NLO_tanBeta50(){

 Int_t nl = 10;
  Double_t xl[10];
  Double_t yl[10];
  Double_t exl[10];
  Double_t eyl[10];
  
  // cout << " n " << hist->GetXaxis()->GetNbins() << endl;
  

  xl[0] = 200;
  yl[0] = 243;
  xl[1] = 220;
  yl[1] = 264;
  xl[2] = 235;
  yl[2] = 278;
  xl[3] = 250;
  yl[3] = 267;
  xl[4] = 300;
  yl[4] = 230;
  xl[5] = 350;
  yl[5] = 205;
  xl[6] = 400;
  yl[6] = 184;
  xl[7] = 450;
  yl[7] = 168;
  xl[8] = 500;
  yl[8] = 156;
  xl[9] = 520;
  yl[9] = 148;
  
 
 
  
  
  TGraphErrors* gr1 = new TGraphErrors(nl,xl,yl,exl,eyl);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(21);
  
  
  //gr1->Draw("LP");

  TSpline3 *s = new TSpline3("grs",gr1);
  s->SetLineColor(kRed);
  //  s->SetLineStyle(2);
  s->SetLineWidth(3);
  

  return gr1;

}




TH1F* getHisto_1d( TString path,
	       TString nameHist,
	       TString nameFile ) {
  TString name = path + nameFile;
  TFile* file =  new TFile(name);
  //  TDirectory* dir = (TDirectory*)file->Get(Dirname);
  TH2F* hist = (TH2F*)file->Get(nameHist);
  if (!hist) {
    std::cout << " name: " << nameHist
	      << " file: " << nameFile
	      << std::endl;
//     abort();

  }
  assert(hist);
  

  TH1F* Onedhist = new TH1F(nameHist,nameHist,int(hist->GetXaxis()->GetNbins()+0.5),hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());

  for(int x = 0;x < hist->GetXaxis()->GetNbins();x++){

    bool firsthit = false;
    for(int y = hist->GetYaxis()->GetNbins(); y>0; y--){

      double y_height = hist->GetYaxis()->GetXmin()+hist->GetYaxis()->GetBinWidth(0)*y;

      if(firsthit == false && hist->GetBinContent(x+1,y) > 0){
	Onedhist->SetBinContent(x+1,y_height);
	firsthit = true;
	}


    }

  }


  Onedhist->SetLineWidth(1);

  Onedhist->GetXaxis()->SetTitleSize(0.055);
  Onedhist->GetYaxis()->SetTitleSize(0.055);
  Onedhist->GetXaxis()->SetLabelSize(0.05);
  Onedhist->GetYaxis()->SetLabelSize(0.05);
  Onedhist->SetStats(kFALSE);

  
  

  return Onedhist;
}


//utility function to convert two curves into a shaded region, TSpline3 version
TGraph* getShadedRegion( TSpline3 * a, TSpline3 * b) {

  const int n=50;

  double xmin = a->GetXmin();
  double xmax = a->GetXmax();

  TGraph *grshade = new TGraph(2*n);
  for ( int i=0; i<n; i++) {
    grshade->SetPoint(i, xmin+i*(xmax-xmin)/double(n), a->Eval(xmin+i*(xmax-xmin)/double(n)));
    grshade->SetPoint(n+i, xmax-i*(xmax-xmin)/double(n), b->Eval(xmax-i*(xmax-xmin)/double(n)));
  }
  return grshade;
}

//utility function to convert two curves into a shaded region, TGraph version
TGraph* getShadedRegion( TGraph * a, TGraph * b) {

  if (a==0 || b==0) return 0;

  const int n=a->GetN() + b->GetN();
  TGraph *grshade = new TGraph(n);

  int offset=0;
  for (int i=0; i<a->GetN(); i++) {
    grshade->SetPoint(i, a->GetX()[i], a->GetY()[i]);
    offset=i;
  }

  for (int i=0; i<b->GetN(); i++) {
    grshade->SetPoint(i+offset+1, b->GetX()[b->GetN()-i-1], b->GetY()[b->GetN()-i-1]);
  }

  return grshade;
}


// NLO exclusion counter here
#include "results.hh"



#endif
