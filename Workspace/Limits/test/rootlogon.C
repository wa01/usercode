{
  gROOT->ProcessLine(".L useNiceColorPalette.C");
  useNiceColorPalette();
  gROOT->ProcessLine(".L PlotLimits.C+");
}
