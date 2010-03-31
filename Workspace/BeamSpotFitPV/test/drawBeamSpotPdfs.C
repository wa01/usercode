void drawBeamSpotPdfs (TDirectory* directory, const char* coord,
			const char* fname)
{
  TH3* obsHisto = directory->Get("PVobs");
  TH3* estHisto = directory->Get("PVest");
  if ( obsHisto==0 || estHisto==0 )  return;
  
  std::string fullName("PV");
  fullName += coord;
  if ( fname )  fullName += fname;
  else  fullName += directory->GetName();
  TCanvas* c = new TCanvas(fullName.c_str(),fullName.c_str());
  c->SetLogy(true);
  TH1* obsHisto1D = obsHisto->Project3D(coord);
  TH1* estHisto1D = estHisto->Project3D(coord);
  TAxis* xaxis = obsHisto1D->GetXaxis();
  string atitle("Primary vertex ");
  atitle += coord;
  atitle += " [cm]";
  xaxis->SetTitle(atitle.c_str());
  if ( c->GetLogy() )  obsHisto1D->SetMinimum(0.5);
  obsHisto1D->SetMarkerStyle(21);
//   obsHisto1D->SetLineColor(2);
//   obsHisto1D->SetMarkerColor(2);
  obsHisto1D->Draw("E0");
  estHisto1D->SetLineColor(2);
  estHisto1D->SetLineWidth(2);
  estHisto1D->Draw("hist same c");

  if ( strcmp(coord,"x") || strcmp(coord,"y") )  xaxis->SetNdivisions(508);

  string epsName = fullName + ".eps";
  c->SaveAs(epsName.c_str());
  string pngName = fullName + ".png";
  c->SaveAs(pngName.c_str());
  
}


void drawBeamSpotPdfs (TDirectory* directory, const char* fname=0)
{
  TH1* chi2Histo = directory->Get("chi2");
  if ( chi2Histo==0 )  return;
  std::string cname("Chi2");
  if ( fname )  cname += fname;
  else  cname += directory->GetName();
  TCanvas* c = new TCanvas(cname.c_str(),cname.c_str());
  chi2Histo->SetMinimum(0.);
  chi2Histo->GetXaxis()->SetTitle("#chi^{2} probability vertex - beamspot");
  chi2Histo->Draw();
  string epsName = cname + ".eps";
  c->SaveAs(epsName.c_str());
  string pngName = cname + ".png";
  c->SaveAs(pngName.c_str());

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
