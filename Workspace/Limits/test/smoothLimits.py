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
parser.add_option("--suffix", dest="suffix", default="It1", type="string", action="store")
(options, args) = parser.parse_args()

from ROOT import gROOT
gROOT.ProcessLine(".L SmoothingUtils.C+")
gROOT.ProcessLine(".L useNiceColorPalette.C")
gROOT.ProcessLine("useNiceColorPalette()")
ROOT.gStyle.SetOptStat(0)

from signalUtils import buildSignalString

assert len(args) > 1
lfname = args[0]
islash = lfname.rfind("/")
if islash >= 0:  lfname = lfname[2:]
parts = lfname.split("_")
btag = None
ht = None
met = None
for i,p in enumerate(parts):
  if p.startswith("ht"):
    ht = int(p[2:])
    if i > 0:  btag = parts[i-1]
  elif p.startswith("met"):
    met = int(p[3:])
assert ht != None and met != None and btag != None
refDict = cPickle.load(file(args[1]))
assert ht in refDict and met in refDict[ht]
bt = btag if btag != 'binc' else 'inc'

file = ROOT.TFile(args[0])
hSmoothed = ROOT.doLimits(file,"hObs",0)
nbx = hSmoothed.GetNbinsX()
xl = hSmoothed.GetXaxis().GetXmin()
xh = hSmoothed.GetXaxis().GetXmax()
dx = (xh-xl) / nbx
nby = hSmoothed.GetNbinsY()
yl = hSmoothed.GetYaxis().GetXmin()
yh = hSmoothed.GetYaxis().GetXmax()
dy = (yh-yl) / nby
newDict = {}
newDict[ht] = {}
newDict[ht][met] = {}
for ix in range(nbx):
  for iy in range(nby):
    v = hSmoothed.GetBinContent(ix+1,iy+1)
    if v < 0.000001:  continue
    m0 = int(xl+(ix+0.5)*dx+0.5)
    m12 = int(yl+(iy+0.5)*dy+0.5)
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
    if not msugra in refDict[ht][met]:
      print msugra,"missing in refDict"
      continue
    if not bt in refDict[ht][met][msugra]:
      print bt,"missing in refDict"
      continue
    if not msugra in newDict[ht][met]: newDict[ht][met][msugra] = {}
    if m0 == 850 and m12 == 450:
      print ix,iy,v,refDict[ht][met][msugra][bt]
    newDict[ht][met][msugra][bt] = v*refDict[ht][met][msugra][bt]
c = ROOT.TCanvas()
hSmoothed.Draw("zcol")
try:
  input("abc")
except:
  pass

oname = args[1]
oname = oname.replace(".pkl","_"+options.suffix+".pkl")
if os.path.exists(oname):
  print "Output file",oname,"exists"
  sys.exit(1)
ofile = open(oname,"wb")
cPickle.dump(newDict,ofile)
ofile.close()
