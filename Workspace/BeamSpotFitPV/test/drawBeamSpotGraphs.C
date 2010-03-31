void drawBeamSpotGraph (TDirectory* directory, TH1* refHisto, const char* name,
			const char* fname)
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
  if ( strcmp(name,"x")==0 )  yaxis->SetTitle("PV x position [cm]");
  else if ( strcmp(name,"y")==0 )  yaxis->SetTitle("PV y position [cm]");
  else if ( strcmp(name,"z")==0 )  yaxis->SetTitle("PV z position [cm]");
  else if ( strcmp(name,"ex")==0 )  yaxis->SetTitle("PV x width [cm]");
  else if ( strcmp(name,"ey")==0 )  yaxis->SetTitle("PV y width [cm]");
  else if ( strcmp(name,"ez")==0 )  yaxis->SetTitle("PV z width [cm]");
  else if ( strcmp(name,"corrxy")==0 )  yaxis->SetTitle("PV x-y correlation");
  else if ( strcmp(name,"dxdz")==0 )  yaxis->SetTitle("PV slope dx/dz");
  else if ( strcmp(name,"dydz")==0 )  yaxis->SetTitle("PV slope dy/dz");
  double xmin,xmax,ymin,ymax;
  graph->ComputeRange(xmin,ymin,xmax,ymax);
//   h->SetMinimum(ymin);
//   h->SetMaximum(ymax);
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

  string epsName = fullName + ".eps";
  c->SaveAs(epsName.c_str());
  string pngName = fullName + ".png";
  c->SaveAs(pngName.c_str());

//   double range(0.);
//   if ( strcmp(name,"x")==0 || strcmp(name,"y")==0 )  range = 0.005;
//   else if ( strcmp(name,"z")==0 )  range = 0.75;
//   else if ( strcmp(name,"ex")==0 || strcmp(name,"ey")==0 )  range = 0.003;
//   else if ( strcmp(name,"ez")==0 )  range = 0.5;
//   else if ( strcmp(name,"corrxy")==0 )  range = 0.4;
//   else if ( strcmp(name,"dxdz")==0 || strcmp(name,"dydz")==0 )  range = 0.001;
//   double mean = graph->GetMean();
//   h->SetMinimum(mean-range);
//   h->SetMaximum(mean+range);
//   double x,y;
//   double ex,ey;
//   for ( unsigned int ip=0; ip<graph->GetN(); ++ip ) {
//     graph->GetPoint(ip,x,y);
//     int ib = int(x+0.5);
//     std::cout << h->GetXaxis()->GetBinLabel(ib) << " " << y 
// 	      << " " << graph->GetErrorY(ip) << std::endl;
//   }

}

void drawBeamSpotGraphs (TDirectory* directory, const char* fname=0)
{
  TH1* refHisto = directory->Get("pvcounts");
  if ( refHisto==0 )  return;

  drawBeamSpotGraph(directory,refHisto,"x",fname);
  drawBeamSpotGraph(directory,refHisto,"y",fname);
  drawBeamSpotGraph(directory,refHisto,"z",fname);
  drawBeamSpotGraph(directory,refHisto,"ex",fname);
  drawBeamSpotGraph(directory,refHisto,"ey",fname);
  drawBeamSpotGraph(directory,refHisto,"ez",fname);
  drawBeamSpotGraph(directory,refHisto,"corrxy",fname);
  drawBeamSpotGraph(directory,refHisto,"dxdz",fname);
  drawBeamSpotGraph(directory,refHisto,"dydz",fname);
}

void drawBeamSpotGraphs (TFile* file, const char* dirname,
			 const char* fname=0)
{
  TDirectory* dir = (TDirectory*)file->FindObjectAny(dirname);
  drawBeamSpotGraphs(dir,fname);
}
