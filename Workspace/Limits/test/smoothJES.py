import cPickle
import ROOT
import math
import sys
import os
#
# read JES systematics for signal from input dictionary, smooth and
#   write smoothed systematics to a new dictionary

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-M", "--model", dest="model", default="msugra", type="string", action="store", help="signal model")
parser.add_option("--m3Ratio", dest="m3Ratio", default=-1., type="float", action="store", help="ratio for intermediate mass in SMS models")
(options, args) = parser.parse_args()

from ROOT import gROOT
gROOT.ProcessLine(".L SmoothingUtils.C+")
gROOT.ProcessLine(".L useNiceColorPalette.C")
gROOT.ProcessLine("useNiceColorPalette()")
ROOT.gStyle.SetOptStat(0)

from signalUtils import buildSignalString

hts = [ 750, 1000 ]
mets = [ 250, 350, 450, 550 ]

jesDict = {}
for ht in hts:
  if not ht in jesDict:  jesDict[ht] = {}
  for met in mets:
    if not met in jesDict[ht]:  jesDict[ht][met] = {}
    filename = "/afs/hephy.at/user/s/schoefbeck/www/pngJES/"+options.model
    if options.model.startswith("T3w") and options.m3Ratio > 0:
      filename += "_"+str(options.m3Ratio)
    filename += "_"+str(ht)+"_"+str(met)+".root"
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
        if options.model == "msugra":
          msugra = buildSignalString("msugra",m0,m12)
        elif options.model == "T1tttt":
          msugra = options.model+"_"+str(m0)+"_"+str(m12)+"_-1"
        elif options.model.startswith("T3w") and options.m3Ratio > 0:
          m3 = int(options.m3Ratio*(m0-m12)+m12+0.5)
          msugra = options.model+"_"+str(m0)+"_"+str(m12)+"_"+str(m3)
        else:
          print "Unknown model",options.model
          sys.exit(1)
        jesDict[ht][met][msugra] = v
    c = ROOT.TCanvas()
    hSmoothed.Draw("zcol")
    try:
      input("abc")
    except:
      pass

oname = options.model+"JES-smoothed.pkl"
if options.model.startswith("T3w") and options.m3Ratio > 0:
  oname += str(int(100*options.m3Ratio+0.5))
oname += "JES-smoothed.pkl"
if os.path.exists(oname):
  print "Output file",oname,"exists"
  sys.exit(1)
ofile = open(oname,"wb")
cPickle.dump(jesDict,ofile)
ofile.close()
