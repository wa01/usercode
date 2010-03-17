void drawBeamSpotPdfs (TDirectory* directory, const char* coord,
			const char* fname)
{
  TH3* obsHisto = directory->Get("PVobs");
  TH3* estHisto = directory->Get("PVest");
  if ( obsHisto==0 || estHisto==0 )  return;
  
  std::string fullName("PV");
  fullName += coord;
  if ( fname )  fullName += fname;
  TCanvas* c = new TCanvas(fullName.c_str(),fullName.c_str());
  TH1* obsHisto1D = obsHisto->Project3D(coord);
  TH1* estHisto1D = estHisto->Project3D(coord);
  obsHisto1D->SetMarkerStyle(21);
  obsHisto1D->SetMarkerColor(2);
  obsHisto1D->Draw("E2");
  estHisto1D->Draw("hist same c");
}


void drawBeamSpotPdfs (TDirectory* directory, const char* fname=0)
{
  TH1* chi2Histo = directory->Get("chi2");
  if ( chi2Histo==0 )  return;
  std::string cname("Chi2");
  if ( fname )  cname += fname;
  TCanvas* c = new TCanvas(cname.c_str(),cname.c_str());
  chi2Histo->Draw();

  drawBeamSpotPdfs(directory,"x",fname);
  drawBeamSpotPdfs(directory,"y",fname);
  drawBeamSpotPdfs(directory,"z",fname);
}

void drawBeamSpotPdfs (TFile* file, const char* dirname,
		       const char* fname=0)
{
  TDirectory* dir = (TDirectory*)file->FindObjectAny(dirname);
  drawBeamSpotPdfs(dir,fname);
}
