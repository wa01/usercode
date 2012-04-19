{
  gROOT->ProcessLine(".L useNiceColorPalette.C");
  useNiceColorPalette();
//   gROOT->ProcessLine(".L PlotLimits.C+");
  gROOT->ProcessLine(".L PlotLimits2.C+");
  gROOT->ProcessLine(".L ExclusionPlot.C+");
}
