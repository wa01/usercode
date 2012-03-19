import cPickle
import ROOT
import math
import sys
#
# read JES systematics for signal from input dictionary, smooth and
#   write smoothed systematics to a new dictionary

from ROOT import gROOT
gROOT.ProcessLine(".L SmoothingUtils.C+")

from signalUtils import buildSignalString

hts = [ 750, 1000 ]
mets = [ 250, 350, 450, 550 ]

jesDict = {}
for ht in hts:
  if not ht in jesDict:  jesDict[ht] = {}
  for met in mets:
    if not met in jesDict[ht]:  jesDict[ht][met] = {}
    filename = "msugra_"+str(ht)+"_"+str(met)+".root"
    file = ROOT.TFile(filename)
#    hRaw = ROOT.findJesHisto(file)
#    print hRaw.GetName()
    hSmoothed = ROOT.doJES(file,0)
    nbx = hSmoothed.GetNbinsX()
    xl = hSmoothed.GetXaxis().GetXmin()
    xh = hSmoothed.GetXaxis().GetXmax()
    dx = (xh-xl) / nbx
    nby = hSmoothed.GetNbinsY()
    yl = hSmoothed.GetYaxis().GetXmin()
    yh = hSmoothed.GetYaxis().GetXmax()
    dy = (yh-yl) / nby
    for ix in range(nbx):
      for iy in range(nby):
        v = hSmoothed.GetBinContent(ix+1,iy+1)
        if v < 0.000001:  continue
        m0 = int(xl+ix*dx+0.5)
        m12 = int(yl+iy*dy+0.5)
        msugra = buildSignalString("msugra",m0,m12)
        jesDict[ht][met][msugra] = v
    c = ROOT.TCanvas()
    hSmoothed.Draw("zcol")
    try:
      input("abc")
    except:
      pass

ofile = open("msugraJES-smoothed.pkl","wb")
cPickle.dump(jesDict,ofile)
ofile.close()
