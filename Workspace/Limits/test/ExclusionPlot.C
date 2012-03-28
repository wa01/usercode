#include <sstream>
#include <string>
#include <iomanip>
#include "ExclusionPlot.hh"
 
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TMarker.h"
#include <vector>
#include "TMath.h"

//Produce the limit plot with the command: root -l ExclusionPlot.C+

TFile* ContourFile(0);
TFile* ContourFileM(0);
TFile* ContourFileP(0);

void ExclusionPlot(TFile* contours, TFile* contoursM, TFile* contoursP){
  ContourFile = contours;
  ContourFileM = contoursM;
  ContourFileP = contoursP;
  gStyle->SetPalette(1);

  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);

  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");



  Int_t tanBeta = 10; //can also be set to 40
  Bool_t plotLO = false;
  Bool_t tb40_plotExpected=true; //can be set to false or true for tan beta 40. no effect on tan beta 10
   
  CommandMSUGRA("Limits_2011.root",tanBeta, plotLO, tb40_plotExpected);
}


void CommandMSUGRA(TString plotName,Int_t tanBeta_, Bool_t plotLO_, Bool_t tb40_plotExpected) {
  
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1); 
  gStyle->SetTextFont(42);
  gStyle->SetFrameBorderMode(0);

  //convert tanb value to string
  std::stringstream tmp;
  tmp << tanBeta_;
  TString tanb( tmp.str() );
  
  
  // Output file
  std::cout << " create " << plotName << std::endl;
  TFile* output = new TFile( plotName, "RECREATE" );
  if ( !output || output->IsZombie() ) { std::cout << " zombie alarm output is a zombie " << std::endl; }
  

  //set old exclusion Limits
  TGraph* LEP_ch = set_lep_ch(tanBeta_);
  TGraph* LEP_sl = set_lep_sl(tanBeta_);//slepton curve
  TGraph* TEV_sg_cdf = set_tev_sg_cdf(tanBeta_);//squark gluino cdf
  TGraph* TEV_sg_d0 = set_tev_sg_d0(tanBeta_);//squark gluino d0
  //  TGraph* TEV_tlp_cdf = set_tev_tlp_cdf(tanBeta_);//trilepton cdf
  //  TGraph* TEV_tlp_d0 = set_tev_tlp_d0(tanBeta_);//trilepton d0
  TGraph* stau   = set_tev_stau(tanBeta_);//stau 
  TGraph* NoEWSB = set_NoEWSB(tanBeta_); 

  TGraph* TEV_sn_d0_1 = set_sneutrino_d0_1(tanBeta_);
  TGraph* TEV_sn_d0_2 = set_sneutrino_d0_2(tanBeta_);

  //some tan beta 40 stuff (load the squark and gluino mass lines)
  TGraph2D* squarkMasses=0;
  TGraph2D* gluinoMasses=0;
  if (tanBeta_==40) {
    const int nPoints = nSusyGridPoints();
    double m0[nPoints],m12[nPoints],squarkMass[nPoints],gluinoMass[nPoints];
    
    susyGrid(m0,m12,squarkMass,gluinoMass);
    
    squarkMasses = new TGraph2D("squarkMasses","",nPoints,m0,m12,squarkMass);
    gluinoMasses = new TGraph2D("gluinoMasses","",nPoints,m0,m12,gluinoMass);
    
    gluinoMasses->GetHistogram();
    squarkMasses->GetHistogram();
  }
  // end of tan beta 40 stuff

  //constant squark and gluino lines
  TF1* lnsq[10];
  TF1* lngl[10];
  TLatex* sq_text[10];
  TLatex* gl_text[10];

  //versions for tan beta 40
  TGraph* lnsq_40[15];
  TGraph* lngl_40[15];
  TLatex* sq_40_text[15];
  TLatex* gl_40_text[15];

  int loopmax = 6;
  if (tanBeta_==40) loopmax=15;
  for(int i = 0; i < loopmax; i++){
    if (tanBeta_==10) {
      lnsq[i] = constant_squark(tanBeta_,i);
      sq_text[i] = constant_squark_text(i,*lnsq[i],tanBeta_);
      lngl[i] = constant_gluino(tanBeta_,i);
      gl_text[i] = constant_gluino_text(i,*lngl[i]);
    }
    else if (tanBeta_==40) {
      lnsq_40[i] = constant_mass(i*250,squarkMasses);
      lngl_40[i] = constant_mass(i*250,gluinoMasses);
      sq_40_text[i] = constant_squark_text_tanBeta40(i*250,lnsq_40[i]);
      gl_40_text[i] = constant_gluino_text_tanBeta40(i*250,lngl_40[i]);;
    }
  }
  

  //Legends
  TLegend* legst  = makeStauLegend(0.05,tanBeta_);
//   TLegend* legNoEWSB  = makeNoEWSBLegend(0.05,tanBeta_);
  TLegend* legexp = makeExpLegend( *TEV_sg_cdf,*TEV_sg_d0,*LEP_ch,*LEP_sl,*TEV_sn_d0_1,0.035,tanBeta_);
  
 
  //make Canvas
  TCanvas* cvsSys = new TCanvas("cvsnm","cvsnm",0,0,800,600);
  gStyle->SetOptTitle(0);
  cvsSys->SetFillColor(0);
  cvsSys->GetPad(0)->SetRightMargin(0.07);
  cvsSys->Range(-120.5298,26.16437,736.0927,750);
  //  cvsSys->Range(-50.5298,26.16437,736.0927,500);
  cvsSys->SetFillColor(0);
  cvsSys->SetBorderMode(0);
  cvsSys->GetPad(0)->SetBorderMode(0);
  cvsSys->GetPad(0)->SetBorderSize(2);
  cvsSys->GetPad(0)->SetLeftMargin(0.1407035);
  cvsSys->GetPad(0)->SetTopMargin(0.08);
  cvsSys->GetPad(0)->SetBottomMargin(0.13);

  cvsSys->SetTitle("tan#beta="+tanb);
 
  output->cd();
  
  TDirectory* curDir = gDirectory;
//   TFile* f = new TFile("limits_binc_ht1000_met350-contours.root");
//   TFile* f = new TFile("limits_msugraNLO_multibtag_ht1000_met250_m12_0-550_HN_comb-contours.root");
  TFile* f = ContourFile;
  TGraph* ra4VieObs = (TGraph*)f->Get("gObs");
  TGraph* ra4VieExpM2 = (TGraph*)f->Get("gExpMinus2");
  TGraph* ra4VieExpM1 = (TGraph*)f->Get("gExpMinus1");
  TGraph* ra4VieExp = (TGraph*)f->Get("gExpMedian");
  TGraph* ra4VieExpP1 = (TGraph*)f->Get("gExpPlus1");
  TGraph* ra4VieExpP2 = (TGraph*)f->Get("gExpPlus2");

  TGraph* ra4VieObsThM(0);
  TGraph* ra4VieExpThM(0);
  if ( ContourFileM ) {
    ra4VieObsThM = (TGraph*)ContourFileM->Get("gObs");
    ra4VieExpThM = (TGraph*)ContourFileM->Get("gExpMedian");
  }
  TGraph* ra4VieObsThP(0);
  TGraph* ra4VieExpThP(0);
  if ( ContourFileP ) {
    ra4VieObsThP = (TGraph*)ContourFileP->Get("gObs");
    ra4VieExpThP = (TGraph*)ContourFileP->Get("gExpMedian");
  }

  curDir->cd();


  double m0min = 0;
  double m0max=1800;
  double xscale = m0max-m0min;
  if (tanBeta_ == 50) m0min=200;
  if (tanBeta_ == 40) {m0min=400;  m0max=2000;}
  xscale = (m0max-m0min)/xscale;
  TH2D* hist = new TH2D("h","h",100,m0min,m0max,100,120,700);
  hist->Draw();  
  hist->GetXaxis()->SetTitle("m_{0} (GeV/c^{2})");
  hist->GetYaxis()->SetTitle("m_{1/2} (GeV/c^{2})");
  hist->GetXaxis()->SetTitleOffset(.9);
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleOffset(1.0);
  hist->GetYaxis()->SetTitleSize(0.06);

  hist->GetXaxis()->SetNdivisions(506);
  //  if (tanBeta_ == 50)  hist->GetXaxis()->SetNdivisions(504);
  hist->GetYaxis()->SetNdivisions(506);

  int col[]={2,3,4};

  TSpline3 *sRA4_LP =0;
  TSpline3 *sRA1 = 0;
  TSpline3 *sRA2 = 0;
  TSpline3 *sRA4_old =0;
  TSpline3 *sRAZ =0;

  if (tanBeta_==10) {
    ra4VieObs->SetLineWidth(3);
    ra4VieObs->SetLineColor(2);
    ra4VieExp->SetLineWidth(3);
    ra4VieExp->SetLineStyle(2);
    ra4VieExp->SetLineColor(4);
//     ra4VieObs->RemovePoint(0);
//     ra4VieObs->RemovePoint(0);
//     double x,y;
//     ra4VieObs->GetPoint(36,x,y);
//     std::cout << x << " " << y << std::endl;
//     ra4VieObs->RemovePoint(37);
//     ra4VieObs->RemovePoint(37);
//     ra4VieObs->RemovePoint(37);
    if ( ra4VieExpM1 && ra4VieExpP1 ) {
      TGraph* ra4ExpArea = new TGraph();
      int np(0);
      double* xExp = ra4VieExpM1->GetX();
      double* yExp = ra4VieExpM1->GetY();
      for ( int i=0; i<ra4VieExpM1->GetN(); ++i ) 
	ra4ExpArea->SetPoint(np++,xExp[i],yExp[i]);
      xExp = ra4VieExpP1->GetX();
      yExp = ra4VieExpP1->GetY();
      for ( int i=ra4VieExpP1->GetN()-1; i>=0; --i ) 
	ra4ExpArea->SetPoint(np++,xExp[i],yExp[i]);
      ra4ExpArea->SetFillColor(7);
      ra4ExpArea->Draw("F");
      ra4VieExpM1->SetLineWidth(2);
      ra4VieExpM1->SetLineStyle(3);
      ra4VieExpM1->SetLineColor(4);
      ra4VieExpP1->SetLineWidth(2);
      ra4VieExpP1->SetLineStyle(3);
      ra4VieExpP1->SetLineColor(4);
      ra4VieExpM1->Draw();
      ra4VieExpP1->Draw();
    }
    if ( ra4VieObsThM && ra4VieObsThP ) {
      ra4VieObsThM->SetLineWidth(2);
      ra4VieObsThM->SetLineColor(2);
      ra4VieObsThM->SetLineStyle(2);
      ra4VieObsThM->Draw();
      ra4VieObsThP->SetLineWidth(2);
      ra4VieObsThP->SetLineColor(2);
      ra4VieObsThP->SetLineStyle(2);
      ra4VieObsThP->Draw();
      
      ra4VieExpThM->SetLineWidth(2);
      ra4VieExpThM->SetLineColor(4);
      ra4VieExpThM->SetLineStyle(2);
      ra4VieExpThM->Draw();
      ra4VieExpThP->SetLineWidth(2);
      ra4VieExpThP->SetLineColor(4);
      ra4VieExpThP->SetLineStyle(2);
      ra4VieExpThP->Draw();
      
    }
    ra4VieExp->Draw();
    ra4VieObs->Draw();
//     if ( ra4VieExpM2 && ra4VieExpP2 ) {
//       ra4VieExpM2->SetLineWidth(2);
//       ra4VieExpM2->SetLineStyle(3);
//       ra4VieExpM2->SetLineColor(4);
//       ra4VieExpP2->SetLineWidth(2);
//       ra4VieExpP2->SetLineStyle(3);
//       ra4VieExpP2->SetLineColor(4);
//       ra4VieExpM2->Draw();
//       ra4VieExpP2->Draw();
//     }
  }
  else if (tanBeta_==40 ) {

  }
  
  TLegend* myleg;

  float leg_x1=0.39+0.23;
  float leg_y1=0.65+0.05;
  float leg_x2= 0.55+0.25;
  float leg_y2= 0.84+0.05;

  if( plotLO_ ) myleg = new TLegend(0.3,0.65,0.65,0.8,NULL,"brNDC");
  else if (tb40_plotExpected) myleg = new TLegend(0.25,0.76,0.44,0.91,NULL,"brNDC"); // copied from else block below
  else if (tanBeta_==40) myleg = new TLegend(leg_x1,leg_y1,leg_x2,leg_y2,NULL,"brNDC");
  else          myleg = new TLegend(0.25,0.76,0.44,0.91,NULL,"brNDC");


  myleg->SetFillColor(0); 
  myleg->SetShadowColor(0);
  myleg->SetTextSize(0.04);
  myleg->SetBorderSize(0);

  TLegendEntry *entry=0;

  if (tanBeta_ == 10 ) {
//     myleg->SetHeader("RA4Tmpl, (NLO, exp. unc.)");
    myleg->AddEntry(ra4VieObs,"obs","l");
    myleg->AddEntry(ra4VieExp,"exp","l");
    if ( ra4VieExpM1 && ra4VieExpP1 ) {
      myleg->AddEntry(ra4VieExpP1,"exp #pm 1#sigma","l");
    }
    //     entry=myleg->AddEntry("RA1","2011 Limits","l");
    //     entry->SetLineColor(1);
    //     entry->SetLineStyle(1);
    //     entry->SetLineWidth(3);
    
    //     entry=myleg->AddEntry("sRA1","2010 Limits","l");
    //     entry->SetLineColor(1);
    //     entry->SetLineStyle(2);
    //     entry->SetLineWidth(3);
  }
  else if (tanBeta_==40) {

    //     entry=myleg->AddEntry("obs","Observed Limit","l");
    //     entry->SetLineColor(1);
    //     entry->SetLineStyle(1);
    //     entry->SetLineWidth(3);
    //     entry->SetTextColor(1);

    //     entry=myleg->AddEntry("exp","Expected Limit #pm 1#sigma","lf");
    //     entry->SetFillColor(kGray);
    //     entry->SetTextColor(1);
    //     entry->SetLineColor(1);
    //     entry->SetLineStyle(7);
    //     entry->SetLineWidth(3);
    //     entry->SetFillStyle(3002);


  }
  
  //constant squark and gluino mass contours
  if (tanBeta_==10) {
    for (int it=0;it<5;it++) {   
      lngl[it]->Draw("same");   
      lnsq[it]->Draw("same");
      sq_text[it]->Draw();
      gl_text[it]->Draw();
    }
  }
  else if (tanBeta_==40) {
    for (int it=2;it<9;it++) {  
      if(it<7){
	if(lngl_40[it]!=0)lngl_40[it]->Draw("samec");   
	if(gl_40_text[it]!=0)gl_40_text[it]->Draw();
      }
      if(lnsq_40[it]!=0)lnsq_40[it]->Draw("samec");
      if(it<6){
	if(sq_40_text[it]!=0)sq_40_text[it]->Draw();
      }
    }
  }

  if (tanBeta_==10) {
    //     SSdilep->Draw("samec");
    //     OSdilep->Draw("samec");
    //     Multilep->Draw("samec");
    //     RA1->Draw("samec");
    //     RA2->Draw("samec");
    //     MT2->Draw("samec");
    //     sRA4_LP->Draw("samec");
    // 	RAZ->Draw("samec");
    
    //     sRA1->Draw("same"); 
    //     sRA2->Draw("same"); 
    //     sRA4_old->Draw("samec");
    // 	sRAZ->Draw("samec");
    //     RA5_old->Draw("c same");
    //     RA6_old->Draw("c same");
  }
  else if (tanBeta_==40) {

    //     //expected curves and errors bands
    //     if ( tb40_plotExpected) {
    //       //not drawing the expected curves for better legibility
    // //       RA1_tb40_exp_p->Draw("samel");
    // //       RA1_tb40_exp_m->Draw("samel");
    // //       RA2b_1b_tight_exp_p->Draw("samel");
    // //       RA2b_1b_tight_exp_m->Draw("samel");
    // //keep the expected +/-1 sigma band
    //       RA1_tb40_exp_band->Draw("f") ;
    //       RA2b_1b_tight_exp_band->Draw("f") ;

    //       //important to draw lines after fill
    //       RA1_tb40_exp->Draw("samel");
    //       RA2b_1b_tight_exp->Draw("samel") ;

  }

  //     RA2b_1b_tight->Draw("samel");
  //     RA1_tb40->Draw("samel");

  // }
  
  TLegend* leg2=0;
  if (tanBeta_==10) {
    //     TLatex* RA1label = new TLatex(670,370.,"#alpha_{T}");
    //     RA1label->SetTextFont(42);
    //     RA1label->SetTextSize(0.05);
    //     RA1label->SetTextColor(kRed+2);
    //     RA1label->Draw("same");
    
    //     TLatex* RA2label = new TLatex(640,465.,"Jets+MHT");
    //     RA2label->SetTextFont(42);
    //     RA2label->SetTextSize(0.04);
    //     RA2label->SetTextColor(kBlue+2);
    //     RA2label->Draw("same");
	  
    //     TLatex* RAZlabel = new TLatex(740,415.,"Razor (0.8 fb^{-1})");
    //     RAZlabel->SetTextFont(42);
    //     RAZlabel->SetTextSize(0.04);
    //     RAZlabel->SetTextColor(kMagenta+1);
    //     RAZlabel->Draw("same");
    
    //     TLatex* RA5label = new TLatex(300,350.,"SS Dilepton");
    //     RA5label->SetTextFont(42);
    //     //RA5label->SetTextAngle(-10);
    //     RA5label->SetTextSize(0.04);
    //     RA5label->SetTextColor(kGreen+2);
    //     RA5label->Draw("same");
    
    //     //TLatex* RA6label = new TLatex(400,280.,"OS Dilepton");
    //     TLatex* RA6label = new TLatex(355,285.,"OS Dilepton");
    //     RA6label->SetTextFont(42);
    //     //RA6label->SetTextAngle(-15);
    //     RA6label->SetTextSize(0.04);
    //     RA6label->SetTextColor(kCyan+2);
    //     RA6label->Draw("same");

    //     TLatex* RA7label = new TLatex(338,215.,"Multi-Lepton");
    //     RA7label->SetTextFont(42);
    //     RA7label->SetTextSize(0.04);
    //     RA7label->SetTextColor(kYellow+2);
    //     RA7label->Draw("same");
    //     TLatex* RA7lumi = new TLatex(450,185.,"(2.1 fb^{-1})");
    //     RA7lumi->SetTextFont(42);
    //     RA7lumi->SetTextSize(0.04);
    //     RA7lumi->SetTextColor(kYellow+2);
    //     RA7lumi->Draw("same");

    
    //     TLatex* MT2label = new TLatex(400,440.,"MT2");
    //     MT2label->SetTextFont(42);
    //     MT2label->SetTextSize(0.04);
    //     MT2label->SetTextColor(kRed);
    //     MT2label->Draw("same");
    
    //     TLatex* RA4label = new TLatex(250,400.,"1 Lepton");
    //     RA4label->SetTextFont(42);
    //     RA4label->SetTextSize(0.04);
    //     RA4label->SetTextColor(kBlue);
    //     RA4label->Draw("same");
  }
  else if (tanBeta_==40  && !tb40_plotExpected) {
    //     int xposRA1 = 480;
    //     int yposRA1 = 510;

    //     int xposRA2b=450;
    //     int yposRA2b=275;

    //     TLatex* RA1label = new TLatex(xposRA1,yposRA1,"#alpha_{T}");
    //     RA1label->SetTextFont(42);
    //     RA1label->SetTextSize(0.05);
    //     RA1label->SetTextColor(kRed+2);
    //     RA1label->Draw("same");

    //     TLatex* RA2blabel = new TLatex(xposRA2b,yposRA2b,"Jets+MET+b");
    //     RA2blabel->SetTextFont(42);
    //     RA2blabel->SetTextSize(0.05);
    //     RA2blabel->SetTextColor(kBlue+1);
    //     RA2blabel->Draw("same");
    //   }
    //   else if (tanBeta_==40 && tb40_plotExpected) {

    //     leg2 = new TLegend(0.63,0.76,0.94,0.91,NULL,"brNDC");
    //     leg2->SetFillColor(0); 
    //     leg2->SetShadowColor(0);
    //     leg2->SetTextSize(0.04);
    //     leg2->SetBorderSize(0);
    
    //     entry=leg2->AddEntry("RA1_tb40","#alpha_{T}","l");
    //     entry->SetLineColor(kRed+2);
    //     entry->SetLineStyle(1);
    //     entry->SetLineWidth(3);
    //     entry->SetTextColor(kRed+2);

    //     entry=leg2->AddEntry("RA2b_1btight","Jets+MET+b","l");
    //     entry->SetLineColor(kBlue+1);
    //     entry->SetLineStyle(1);
    //     entry->SetLineWidth(3);
    //     entry->SetTextColor(kBlue+1);

  }

  //exclusion limits previous experiments
  if(tanBeta_ == 3){
    TEV_sn_d0_1->Draw("fsame");
    TEV_sn_d0_2->Draw("fsame");
  }
  if (tanBeta_==10)   LEP_ch->Draw("fsame");
  if (tanBeta_ != 50 && tanBeta_!=40) LEP_sl->Draw("fsame");

  //remove CDF/D0 excluded regions
  if (tanBeta_==10) {
    TEV_sg_cdf->Draw("fsame");
    TEV_sg_d0->Draw("same");  
    TEV_sg_d0->Draw("fsame");
  }

  //other labels
  Double_t xpos = 0;
  Double_t xposi = 0;
  Double_t ypos = 0;
  if(tanBeta_ == 50) xposi = 100;
  if(tanBeta_ == 50) xpos = 200;
  if(tanBeta_ == 50) ypos = -10;

  if(tanBeta_ == 40) xposi = 180+160;
  if(tanBeta_ == 40) xpos = 400;//240;
  if(tanBeta_ == 40) ypos = 100;

 
  //TLatex* lumilabel = new TLatex(750 +xposi + 100,767.-154,"#sqrt{s} = 7 TeV, #scale[0.65]{#int}Ldt = 0.98 fb^{-1}");
  TLatex* lumilabel = new TLatex(450*xscale+xpos,767.-154+100,"#sqrt{s} = 7 TeV,   Ldt 4.98 fb^{ -1}");
  TLatex* integral_symbol = new TLatex((577+100)*xscale+xpos,767.-145+100,"#int");

  lumilabel->SetTextSize(0.03);
  integral_symbol->SetTextSize(0.015);
  lumilabel->Draw("same");
  integral_symbol->Draw("same");

  TLatex* cmslabel = new TLatex(10.+xpos,767.-154+100,"CMS Preliminary");
  cmslabel->SetTextSize(0.03);
  cmslabel->Draw("same");

  TString text_tanBeta;
  TString a0str="0";
  if (tanBeta_==40) a0str = "-500 GeV";
  text_tanBeta =  "tan#beta = "+tanb+",  A_{0} = "+a0str+",  #mu > 0";
  int anotherOffset = (tb40_plotExpected && tanBeta_==40) ? -100 : 0;
//   TLatex* cmssmpars = new TLatex(/*530.+xpos,690.+ypos-130*/120+xpos,555+ypos+anotherOffset,text_tanBeta);
  TLatex* cmssmpars = new TLatex(0.61,0.60,text_tanBeta);
  cmssmpars->SetNDC(1);

  cmssmpars->SetTextSize(0.04);
  cmssmpars->Draw("same");

  TLatex* lep_chargino = new TLatex(250,135,"LEP2 #tilde{#chi}_{1}^{#pm}");
  lep_chargino->SetTextSize(0.03);
  lep_chargino->SetTextFont(42);
  //  lep_chargino->Draw("same");

  TLatex* lep_slepton = new TLatex(26,190,"LEP2 #tilde{#font[12]{l}}^{#pm}");
  lep_slepton->SetTextSize(0.03);
  lep_slepton->SetTextAngle(-83);
  lep_slepton->SetTextFont(42);
  //  lep_slepton->Draw("same");



  //LM points
  TMarker* LM0 = new TMarker(200.,160.,20);
  TMarker* LM1 = new TMarker(60.,250.,20);
  TMarker* LM3 = new TMarker(330.,240.,20);
  TMarker* LM6 = new TMarker(80.,400.,20);
    
  LM0->SetMarkerSize(1.2);
  LM1->SetMarkerSize(1.2);
    
  TLatex* tLM0 = new TLatex(205.,160.," LM0");
  tLM0->SetTextSize(0.035);
    
  TLatex* tLM1 = new TLatex(80.,245.,"LM1");
  tLM1->SetTextSize(0.035);
  
  //TLatex* tLM3 = new TLatex(350.,235.,"LM3 (tan#beta=20)");
  TLatex* tLM3 = new TLatex(350.,235.,"LM3");
  tLM3->SetTextSize(0.035);
  
  TLatex* tLM6 = new TLatex(100.,395.,"LM6");
  tLM6->SetTextSize(0.035);
  
  //  if (tanBeta_ != 50){
  //  LM0->Draw("same");   
  //  tLM0->Draw("same");
  //  LM1->Draw("same");   
  //  tLM1->Draw("same");
  // }

  /*
    if (tanBeta_ == 10){ 
    LM1->Draw("same");
    tLM1->Draw("same");
    LM3->Draw("same");
    tLM3->Draw("same");
    LM6->Draw("same");
    tLM6->Draw("same");
    }
  */



  //stau=LSP contour
  if (tanBeta_==10) {
    stau->Draw("fsame");
//     NoEWSB->Draw("fsame");
    
    //legends
    legexp->Draw();
    legst->Draw();
    //legNoEWSB->Draw();
  }
  if (tanBeta_!=40 || tb40_plotExpected)  myleg->Draw();
  if (tanBeta_==40 && tb40_plotExpected) leg2->Draw();

  hist->Draw("sameaxis");
  cvsSys->RedrawAxis();
  cvsSys->Update();
  cvsSys->Write();
  
  if( plotLO_ ){
    cvsSys->SaveAs("ExclusionLimit_tanb"+tanb+"_LO.pdf");
    cvsSys->SaveAs("ExclusionLimit_tanb"+tanb+"_LO.png");
  }else{
    cvsSys->SaveAs("ExclusionLimit_tanb"+tanb+".eps");
    cvsSys->SaveAs("ExclusionLimit_tanb"+tanb+".ps");
    cvsSys->SaveAs("ExclusionLimit_tanb"+tanb+".pdf");
    cvsSys->SaveAs("ExclusionLimit_tanb"+tanb+".png");
  }
  
  output->Write();
  //output->Close();
  //delete output; 

}


void setPlottingStyle(TH1F& hsig){
  
  hsig.SetStats(kFALSE);
  
  hsig.SetAxisRange(80,750,"Y");
  hsig.SetAxisRange(0,520,"X");
  hsig.SetAxisRange(200,520,"X");

  hsig.GetXaxis()->SetTitle("m_{0} (GeV)");
  hsig.GetYaxis()->SetTitle("m_{1/2} (GeV)");
  hsig.GetYaxis()->SetTitleOffset(0.8);
  hsig.GetYaxis()->SetTitleSize(0.06);
  hsig.GetYaxis()->SetLabelSize(0.06);
  hsig.GetXaxis()->SetTitleOffset(0.9);
  hsig.GetXaxis()->SetTitleSize(0.06);
  hsig.GetXaxis()->SetLabelSize(0.06);

  hsig.SetLineWidth(1);  
  hsig.SetLineColor(kBlue);  
  
}


TGraph* set_sneutrino_d0_1(Int_t tanBeta){
  double sn_m0[14]= {0,  0, 48, 55, 80, 90,100,105,109,105,100, 72, 55,0};
  double sn_m12[14]={0,140,210,220,237,241,242,241,230,220,210,170,150,0};

  TGraph* sn_d0_gr = new TGraph(14,sn_m0,sn_m12);

  sn_d0_gr->SetFillColor(kGreen+3);
  sn_d0_gr->SetFillStyle(1001);

  return sn_d0_gr;
}

TGraph* set_sneutrino_d0_2(Int_t tanBeta){
  double sn_m0[9]= {0, 45, 75,115,130,150,163,185,0};
  double sn_m12[9]={0,140,170,213,202,183,168,140,0};

  TGraph* sn_d0_gr_2 = new TGraph(9,sn_m0,sn_m12);

  sn_d0_gr_2->SetFillColor(kGreen+3);
  sn_d0_gr_2->SetFillStyle(1001);

  return sn_d0_gr_2;
}

TGraph* set_lep_ch(Int_t tanBeta) {
  if(tanBeta == 10) return set_lep_ch_tanBeta10();
  if(tanBeta == 40) return set_lep_ch_tanBeta40();

  return 0;
}

TGraph* set_lep_ch_tanBeta10(){

double ch_m0[] ={0  ,100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,2000,0.};  
double ch_m12[]={162,162, 161, 160, 160, 159, 158, 157, 156, 155, 154, 154, 153, 152, 150, 149, 148,  147 , 146 , 146 , 146 , 147 , 148 , 149 , 151 , 154 ,159, 0., 0.}; 

  TGraph* ch_gr = new TGraph(29,ch_m0,ch_m12);

  ch_gr->SetFillColor(3);
  ch_gr->SetLineColor(3);
  //  ch_gr->SetLineWidth(3);
  ch_gr->SetFillStyle(1001);

  return ch_gr;

}

TGraph* set_lep_ch_tanBeta40(){

  double ch_m0[] = {0,   240, 400, 500.,700.,800.,1000,1200., 1300., 1400., 1450., 1500., 1550., 1600., 1650., 1700., 1750., 1800., 2150., 2150.,0};
  double ch_m12[]= {140, 140, 139, 138, 137, 136, 136,  137,   138,  139,   140,   142,   143,   145,   147,   150,   154,   158,  134,   0,   0};

  TGraph* ch_gr = new TGraph(21,ch_m0,ch_m12);

  ch_gr->SetFillColor(3);
  ch_gr->SetLineColor(3);
  ch_gr->SetFillStyle(1001);

  return ch_gr;

}



TGraph* set_tev_stau(Int_t tanBeta){

// taken from Frederic Ronga's calculation

    double st_m0_tanBeta10[] =  {0,   10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110 , 130,  147, 0, 170, 190, 210, 230, 250, 270, 290, 310, 330, 350, 370, 390, 410, 430, 450, 470, 490, 510, 530, 550, 570, 590, 610, 630, 6500};
    double st_m12_tanBeta10[] = {213,220,240,275,312,351,393,435,476,518, 559, 600., 682., 750.,750, 842., 921., 999., 1076, 1152, 1228, 1304, 1378, 1453, 1527, 1600, 1673, 1746, 1818, 1890, 1962, 2034, 2105, 2175, 2246, 2316, 2386, 2456, 2526, 2595}; 


    double st_m0_tanBeta40[] = {240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 820, 840, 860, 880};
    double st_m12_tanBeta40[] = {186, 256, 329, 400, 470, 537, 603, 666, 727, 787, 845, 902, 958, 1013, 1067, 1121, 1174, 1226, 1278, 1330, 1381, 1431, 1481, 1531, 1581, 1630, 1679, 1728, 1779, 1825, 1874, 1920, 1971};

    TGraph* st_gr_tanBeta10 = new TGraph(15,st_m0_tanBeta10,st_m12_tanBeta10);
    TGraph* st_gr_tanBeta40 = new TGraph(10,st_m0_tanBeta40,st_m12_tanBeta40);

    
    st_gr_tanBeta40->SetFillColor(40);
    st_gr_tanBeta40->SetFillStyle(1001);
    
    st_gr_tanBeta10->SetFillColor(40);
    st_gr_tanBeta10->SetFillStyle(1001);


    if(tanBeta == 10)return st_gr_tanBeta10;
    if(tanBeta == 40)return st_gr_tanBeta40;
    return 0;
}

TGraph* set_NoEWSB(Int_t tanBeta){

    double st_m0_tanBeta10[]  = { 0, 70, 90,  80,  60, 1000, 1090, 1170, 1240, 1320, 1370, 1440, 1500, 1550, 1610, 1660, 1720, 1780, 1830, 1860, 1920, 1970 , 2000, 2000, 0};
    double st_m12_tanBeta10[] = {10.,10.,20., 30., 40., 50.,  60.,  70.,  80.,  90.,  100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220,   0,   0}; 

//// Needs to be modified
    double st_m0_tanBeta40[] = { 0, 70, 90,  80,  60, 1000, 1090, 1170, 1240, 1320, 1370, 1440, 1500, 1550, 1610, 1660, 1720, 1780, 1830, 1860, 1920, 1970 , 2000, 2000, 0};
    double st_m12_tanBeta40[]= {10.,10.,20., 30., 40., 50.,  60.,  70.,  80.,  90.,  100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220,   0,   0};

    TGraph* st_gr_tanBeta10 = new TGraph(25,st_m0_tanBeta10,st_m12_tanBeta10);
    TGraph* st_gr_tanBeta40 = new TGraph(25,st_m0_tanBeta40,st_m12_tanBeta40);

    
    st_gr_tanBeta40->SetFillColor(40);
    st_gr_tanBeta40->SetFillStyle(1001);
    
    st_gr_tanBeta10->SetFillColor(40);
    st_gr_tanBeta10->SetFillStyle(1001);

    if(tanBeta == 10)return st_gr_tanBeta10;
    if(tanBeta == 40)return st_gr_tanBeta40;
    return 0;
}





TGraph* set_lep_sl(Int_t tanBeta){


  //contour from D0 trilepton paper (PLB 680 (2009) 34-43)

  double *sl_m0 = 0;
  double *sl_m12 = 0;
  int n = 0;

  double sl_m0_3[] ={0,  0, 10, 20, 30, 40, 50, 60, 70, 77,88,95};
  double sl_m12_3[]={0,245,242,239,232,222,209,189,165,140,60,0};
  int n_3 = 12;

  double sl_m0_10[]={ 0,  0, 11, 20, 24, 49, 70, 82,88,90};
  double sl_m12_10[]={0,240,237,233,230,200,150,100,50,0};
  int n_10 = 10;

  if (tanBeta==3){
    sl_m0 = sl_m0_3;
    sl_m12 = sl_m12_3;
    n = n_3;
  }
  //CMS PTDR-II
  //* Selectron_R line mass=99, ISASUGRA7.69, A0=0, m_top=175, tan(beta]=10
  if (tanBeta==10 || tanBeta==50){
    sl_m0 = sl_m0_10;
    sl_m12 = sl_m12_10;
    n = n_10;
  }

  TGraph* lep_sl = new TGraph(n,sl_m0,sl_m12);

  lep_sl->SetFillColor(5);
  lep_sl->SetLineColor(5);
  lep_sl->SetFillStyle(1001);
  
  return lep_sl;
}


TGraph* set_tev_sg_cdf(Int_t tanBeta){

  //New CHF from CDF plot in ICHEP2010 talk (E. Halkiadakis)
  double sg_m0[]= {0,  0, 30, 75,150,185,225,310,360,400,430,500,600,600};
  double sg_m12[]={0,162,168,170,160,150,130,120,109,108,100, 96, 95,  0};
  int np=14;

  TGraph* sg_gr = new TGraph(np,sg_m0,sg_m12);

  //  gStyle->SetHatchesLineWidth(3);

  sg_gr->SetFillColor(2);
  sg_gr->SetLineColor(2);
  //  sg_gr->SetLineWidth(3);
  sg_gr->SetFillStyle(1001); 

  return sg_gr;

}

TGraph* set_tev_sg_d0(Int_t tanBeta){

  double sgd_m0[]= {0,  0, 30, 80,150,240,320,400,500,600,600,0};
  double sgd_m12[]={0,167,166,162,156,138,121,109,105,105,  0,0};
  int npd=12;

  TGraph* sgd_gr = new TGraph(npd,sgd_m0,sgd_m12);

  gStyle->SetHatchesLineWidth(3);

  sgd_gr->SetFillColor(kMagenta+3);
  sgd_gr->SetLineColor(kMagenta+3);
  sgd_gr->SetLineWidth(3);
  sgd_gr->SetFillStyle(3335);

  return sgd_gr;

}



//From Sanjay
TF1* constant_squark(int tanBeta,int i){
//---lines of constant gluino/squark.
// Min squark mass from 1st and 2nd generations using fit for tanbeta = 10.

  double coef1[] = {2.67058e+04, 6.39642e+04, 1.16565e+05, 1.95737e+05, 2.86190e+05};
  double coef2[] = {1.98772e-01, 2.11242e-01, 2.17734e-01, 2.39535e-01, 2.39768e-01};
  double coef3[] = {2.67058e+04, 6.39641e+04, 1.16565e+05, 1.95736e+05, 2.86189e+05};
 
  char hname[200];

  sprintf(hname,"lnsq_%i",i);
  TF1* lnsq = new TF1(hname,"sqrt([0]-x*x*[1]+[2])",0,1800);
  lnsq->SetParameter(0,coef1[i-1]);
  lnsq->SetParameter(1,coef2[i-1]);
  lnsq->SetParameter(2,coef3[i-1]);
  lnsq->SetLineWidth(1);
  lnsq->SetLineColor(kGray);

  return lnsq;
}


TF1* constant_gluino(int tanBeta,int i){
//---lines of constant gluino/squark
  char hname[200];
  sprintf(hname,"lngl_%i",i);

  double coef1[] = {201.77, 311.027, 431.582, 553.895, 676.137};
  double coef2[] = {-0.0146608, -0.01677, -0.022244, -0.0271851, -0.0292212};
   
  TF1* lngl = new TF1(hname,"[0]+x*[1]",0,1800);
  lngl->SetParameter(0,coef1[i-1]);
  lngl->SetParameter(1,coef2[i-1]);
  lngl->SetLineWidth(1);
  lngl->SetLineColor(kGray);

  return lngl;
}


TLatex* constant_squark_text(Int_t it,TF1& lnsq,Int_t tanBeta_){
  char legnm[200];
  it--; //For Sanjay's code

  sprintf(legnm,"#font[92]{#tilde{q}(%i)GeV}",500+250*it);
  Double_t place_x = 170;
  if(tanBeta_ == 50)place_x = 290;
  TLatex* t3 = new TLatex(place_x+10*it,lnsq.Eval(-50+place_x+10*it)+5,legnm);
  t3->SetTextSize(0.03);
  t3->SetTextAngle(-15+it*2);
  t3->SetTextColor(kGray+2);

  return t3;
}


TLatex* constant_gluino_text(Int_t it,TF1& lngl){ //, Int_t tanBeta_){
  char legnm[200];
  Int_t tanBeta_ = 10;
  it--; //For Sanjay's code

  sprintf(legnm,"#font[12]{#tilde{g}}#font[92]{(%i)GeV}",500+250*it);

  Double_t place_x = 423;
  Double_t place_y = 18;
  if (tanBeta_ == 10 ) {
    place_x = 825;
    place_y = -15;
  }
  if (tanBeta_ == 50 ) {
    place_x = 543;
    place_y = 13;
  }
  TLatex* t4 = new TLatex(place_x,place_y+lngl.Eval(480),legnm);
  t4->SetTextSize(0.03);
  t4->SetTextAlign(13);
  t4->SetTextColor(kGray+2);

  return t4;
}

//part of the tan beta 40 code from Luke
TGraph* constant_mass(double mass,TGraph2D* massGrid)
{
  if(!(mass<massGrid->GetZmax() && mass>massGrid->GetZmin())) return 0;
  TGraph* theGraph = (TGraph*)massGrid->GetContourList(mass)->At(0);
  theGraph->SetLineWidth(1);
  theGraph->SetLineColor(kGray);
  return theGraph;
}


TLatex* constant_squark_text_tanBeta40(int mass,TGraph* massLine){
  if(massLine==0) return 0;
  char legnm[200];
  sprintf(legnm,"#font[92]{#tilde{q}(%i)GeV}",mass);
  Double_t place_x = 180 + 100;
  if(place_x+0.16*mass < 400 || massLine->Eval(place_x+0.16*(mass))+10 <100) return 0;
  TLatex* t3 = new TLatex(place_x+0.16*mass,massLine->Eval(place_x+0.16*(mass))+10,legnm);
  t3->SetTextSize(0.03);
  t3->SetTextAngle(-25);
  t3->SetTextColor(kGray+2);
  return t3;
}

TLatex* constant_gluino_text_tanBeta40(int mass,TGraph* massLine){
  if(!massLine) return 0;
  char legnm[200];
  sprintf(legnm,"#font[12]{#tilde{g}}#font[92]{(%i)GeV}",mass);

  Double_t place_x = 1600;
  if(25+massLine->Eval(place_x) <100) return 0;
  TLatex* t4 = new TLatex(place_x,25+massLine->Eval(place_x),legnm);
  t4->SetTextSize(0.03);
  t4->SetTextAlign(13);
  t4->SetTextColor(kGray+2);

  return t4;
}


TLegend* makeStauLegend(Double_t txtsz,Int_t tanBeta_){
  Double_t ypos_1 = 0.78;
  Double_t ypos_2 = 0.80;
  Double_t xpos_1 = 0.16;
  Double_t xpos_2 = 0.17;
  if(tanBeta_ == 40){
    xpos_1 = 0.17;
    xpos_2 = 0.18;
    ypos_1 = 0.76;
    ypos_2 = 0.78;

  }
    
  TLegend* legst = new TLegend(xpos_1+0.025,ypos_1,xpos_2+0.025,ypos_2);
  legst->SetHeader("#tilde{#tau} = LSP");
  legst->SetFillStyle(0);
  legst->SetBorderSize(0);
  legst->SetTextSize(0.03);
  legst->SetTextAngle(80);

  return legst;
}

TLegend* makeNoEWSBLegend(Double_t txtsz,Int_t tanBeta_){
  Double_t ypos_1 = 0.10+0.02;
  Double_t ypos_2 = 0.20+0.02;
  Double_t xpos_1 = 0.82;
  Double_t xpos_2 = 0.92;
  if(tanBeta_ == 40){
    xpos_1 = 0.10;
    xpos_2 = 0.20;
    ypos_1 = 0.85;
    ypos_2 = 0.95;

  }

  TLegend* legst = new TLegend(xpos_1,ypos_1,xpos_2,ypos_2);
  legst->SetHeader("NoEWSB");
  legst->SetFillStyle(0);
  legst->SetBorderSize(0);
  legst->SetTextSize(0.03);
  legst->SetTextAngle(20);

  return legst;
}


TLegend* makeExpLegend(TGraph& sg_gr, TGraph& sgd_gr,TGraph& ch_gr,TGraph& sl_gr,TGraph& tev_sn,Double_t txtsz,Int_t tanbeta){

  //TLegend* legexp = new TLegend(0.61,0.65,0.91,0.9,NULL,"brNDC");
  //TLegend* legexp = new TLegend(0.70,0.70,0.91,0.9,NULL,"brNDC");
  TLegend* legexp = new TLegend(0.61,0.65,0.99,0.9,NULL,"brNDC");


  legexp->SetFillColor(0);
  legexp->SetShadowColor(0);
  legexp->SetTextSize(txtsz);
  legexp->SetBorderSize(0);

  sg_gr.SetLineColor(1);
  
  legexp->AddEntry(&sg_gr,"CDF  #tilde{#font[12]{g}}, #tilde{#font[12]{q}}, #scale[0.8]{tan#beta=5, #mu<0}","f"); 

  legexp->AddEntry(&sgd_gr,"D0   #tilde{#font[12]{g}}, #tilde{#font[12]{q}}, #scale[0.8]{tan#beta=3, #mu<0}","f");  

  ch_gr.SetLineColor(1);
  legexp->AddEntry(&ch_gr,"LEP2   #tilde{#chi}_{1}^{#pm}","f");  
  
  sl_gr.SetLineColor(1);
  if(tanbeta != 50) legexp->AddEntry(&sl_gr,"LEP2   #tilde{#font[12]{l}}^{#pm}","f"); 
  if(tanbeta == 3) legexp->AddEntry(&tev_sn,"D0  #chi^{#pm}_{1}, #chi^{0}_{2}","f");  
 

  return legexp;

}

