void drawBeamSpotGraph (TDirectory* directory, TH1* refHisto, const char* name)
{
  TGraphErrors* graph = directory->Get(name);
  if ( graph==0 )  return;

  char newName[64];
  newName[0] = 'c';
  strcpy(&newName[1],name);
  TCanvas* c = new TCanvas(newName,newName);
  newName[0] = 'h';
  TH1* h = refHisto->Clone(newName);
  h->Reset();
  h->SetTitle(name);
  double xmin,xmax,ymin,ymax;
  graph->ComputeRange(xmin,ymin,xmax,ymax);
  h->SetMinimum(ymin);
  h->SetMaximum(ymax);
  h->Draw();
  graph->Draw("P");
}

void drawBeamSpotGraphs (TDirectory* directory)
{
  TH1* refHisto = directory->Get("pvcounts");
  if ( refHisto==0 )  return;

  drawBeamSpotGraph(directory,refHisto,"x");
  drawBeamSpotGraph(directory,refHisto,"y");
  drawBeamSpotGraph(directory,refHisto,"z");
  drawBeamSpotGraph(directory,refHisto,"ex");
  drawBeamSpotGraph(directory,refHisto,"ey");
  drawBeamSpotGraph(directory,refHisto,"ez");
  drawBeamSpotGraph(directory,refHisto,"corrxy");
  drawBeamSpotGraph(directory,refHisto,"dxdz");
  drawBeamSpotGraph(directory,refHisto,"dydz");
}
