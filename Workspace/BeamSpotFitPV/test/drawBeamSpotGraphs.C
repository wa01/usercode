void drawBeamSpotGraph (TDirectory* directory, TH1* refHisto, const char* name,
			const char* fname)
{
  TGraphErrors* graph = directory->Get(name);
  if ( graph==0 )  return;

  char newName[64];
  newName[0] = 'c';
  strcpy(&newName[1],name);
  std::string fullName(name);
  if ( fname )  fullName += fname;
  TCanvas* c = new TCanvas(fullName.c_str(),fullName.c_str());
  newName[0] = 'h';
  TH1* h = refHisto->Clone(newName);
  h->Reset();
  h->SetTitle(name);
  double xmin,xmax,ymin,ymax;
  graph->ComputeRange(xmin,ymin,xmax,ymax);
  h->SetMinimum(ymin);
  h->SetMaximum(ymax);
  h->Draw();
  graph->SetMarkerStyle(20);
  graph->Draw("P");
  graph->Fit("pol1","same");
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
