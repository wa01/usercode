void scaleGraph (TGraphErrors* graph, float scale)
{
  double x,ex;
  double y,ey;
  unsigned int np = graph->GetN();
  for ( unsigned int i=0; i<np; ++i ) {
    graph->GetPoint(i,x,y);
    ey = graph->GetErrorY(i);
    graph->SetPoint(i,x,scale*y);
    graph->SetPointError(i,ex,scale*ey);
  }
}

void drawBeamSpotGraph (TDirectory* directory, TH1* refHisto, const char* name,
			const char* fname, float* runSummary = 0)
{
  TGraphErrors* graph = directory->Get(name);
  if ( graph==0 )  return;

  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);

  char newName[64];
  newName[0] = 'c';
  strcpy(&newName[1],name);
  std::string fullName(name);
  if ( fname )  fullName += fname;
  else  fullName += directory->GetName();
  TCanvas* c = new TCanvas(fullName.c_str(),fullName.c_str());
  newName[0] = 'h';
  TH1* h = refHisto->Clone(newName);
  h->Reset();
  h->SetTitle(name);
  int nb = h->GetNbinsX();
  if ( nb>50 ) {
    int iscale = nb/50+1;
    TAxis* xaxis = h->GetXaxis();
    for ( int i=1; i<=h->GetNbinsX(); ++i ) {
      if ( (i-1)%iscale )  xaxis->SetBinLabel(i,"");
    }
  }
  TAxis* yaxis = h->GetYaxis();
  float scale(1);
  if ( strcmp(name,"x")==0 ) {
    yaxis->SetTitle("PV x position [mm]");
    scale = 10.;
  }
  else if ( strcmp(name,"y")==0 ) {
    yaxis->SetTitle("PV y position [mm]");
    scale = 10.;
  }
  else if ( strcmp(name,"z")==0 ) {
    yaxis->SetTitle("PV z position [cm]");
  }
  else if ( strcmp(name,"ex")==0 ) {
    yaxis->SetTitle("PV x width [#mum]");
    scale = 10000.;
  }
  else if ( strcmp(name,"ey")==0 ) {
    yaxis->SetTitle("PV y width [#mum]");
    scale = 10000.;
  }
  else if ( strcmp(name,"ez")==0 ) {
    yaxis->SetTitle("PV z width [cm]");
  }
  else if ( strcmp(name,"corrxy")==0 ) {
    yaxis->SetTitle("PV x-y correlation");
  }
  else if ( strcmp(name,"dxdz")==0 ) {
    yaxis->SetTitle("PV slope dx/dz [10^{-3}]");
    scale = 1000.;
  }
  else if ( strcmp(name,"dydz")==0 ) {
    yaxis->SetTitle("PV slope dy/dz [10^{-3}]");
    scale = 1000.;
  }
  scaleGraph(graph,scale);

  double xmin,xmax,ymin,ymax;
  graph->ComputeRange(xmin,ymin,xmax,ymax);
  h->SetMinimum((ymax+ymin)/2.-2.*(ymax-ymin)/2.);
  h->SetMaximum((ymax+ymin)/2.+2.*(ymax-ymin)/2.);
  h->Draw();
  graph->SetMarkerStyle(20);
//   graph->SetMarkerColor(2);
//   graph->SetLineColor(2);
  graph->Draw("P");
  graph->Fit("pol1","same");
  graph->GetFunction("pol1")->SetLineStyle(2);
  graph->GetFunction("pol1")->SetLineWidth(2);

  TF1* fit = graph->GetFunction("pol1");
  cout << "Pol1 fit chi2 = " << fit->GetChisquare() 
       << " " << fit->GetNDF() << endl;

  string epsName = fullName + ".eps";
  c->SaveAs(epsName.c_str());
  string pngName = fullName + ".png";
  c->SaveAs(pngName.c_str());

  if ( runSummary ) {
    TF1* fit = graph->GetFunction("pol1");
    runSummary[0] = fit->GetChisquare();
    runSummary[1] = fit->GetNDF();
    runSummary[2] = fit->GetParameter(0);
    runSummary[3] = fit->GetParError(0);
    runSummary[4] = fit->GetParameter(1);
    runSummary[5] = fit->GetParError(1);
  }
}

void drawBeamSpotGraphs (TDirectory* directory, const char* fname=0, float* runSummary = 0)
{
  TH1* refHisto = directory->Get("pvcounts");
  if ( refHisto==0 )  return;

  drawBeamSpotGraph(directory,refHisto,"x",fname,runSummary);
  drawBeamSpotGraph(directory,refHisto,"y",fname,runSummary);
  drawBeamSpotGraph(directory,refHisto,"z",fname,runSummary);
  drawBeamSpotGraph(directory,refHisto,"ex",fname,runSummary);
  drawBeamSpotGraph(directory,refHisto,"ey",fname,runSummary);
  drawBeamSpotGraph(directory,refHisto,"ez",fname,runSummary);
  drawBeamSpotGraph(directory,refHisto,"corrxy",fname,runSummary);
  drawBeamSpotGraph(directory,refHisto,"dxdz",fname,runSummary);
  drawBeamSpotGraph(directory,refHisto,"dydz",fname,runSummary);
}

void drawBeamSpotGraphs (TFile* file, const char* dirname,
			 const char* fname=0, float* runSummary = 0)
{
  TDirectory* dir = (TDirectory*)file->FindObjectAny(dirname);
  drawBeamSpotGraphs(dir,fname,runSummary);
}

void drawBeamSpotGraphsAll (TFile* file, const char* fname=0)
{
  TDirectory* module = (TDirectory*)_file0->FindObjectAny("beamSpotFitPV");
  if ( module==0 )  return;

//   unsigned int nrun = module->GetListOfKeys()->GetSize();
//   TDirectory* currDir = gDirectory;
//   gROOT->cd();
//   TH1* hChi2 = new TH1F("runChi2","runChi2",nrun,0.5,nrun+0.5);
//   TH1* hConst = new TH1F("runConst","runConst",nrun,0.5,nrun+0.5);
//   TH1* hSlope = new TH1F("runSlope","runSlope",nrun,0.5,nrun+0.5);
//   currDir->cd();

  TIter iter(module->GetListOfKeys());
  TKey* key;
//   unsigned int irun(0);
//   float runSummary[6];
  while ( (key=(TKey*)iter()) ) {
    std::cout << key->GetName() << std::endl;
    drawBeamSpotGraphs(_file0,key->GetName(),fname);
    TIter ic(gROOT->GetListOfCanvases());
    TCanvas* c;
    while ( (c=(TCanvas*)ic()) )  delete c;
//     ++irun;
//     if ( runSummary[1]>0.1 )  hChi2->SetBinContent(irun,runSummary[0]/runSummary[1]);
//     hConst->SetBinContent(irun,runSummary[2]);
//     hConst->SetBinError(irun,runSummary[3]);
//     hSlope->SetBinContent(irun,runSummary[4]);
//     hSlope->SetBinError(irun,runSummary[5]);
//     hChi2->GetXaxis()->SetBinLabel(irun,key->GetName());
//     hConst->GetXaxis()->SetBinLabel(irun,key->GetName());
//     hSlope->GetXaxis()->SetBinLabel(irun,key->GetName());
//     if ( irun==3 )  return;
  }

}
