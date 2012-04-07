import cPickle
import ROOT
import math
import sys
import os

#
# find range and spacing of a 1D grid of masses
#   returns tuple of ( delta, min, max )
#
def findRange (ms):
  dm0 = None
  m0 = None
  m0min = None
  m0max = None
  for im in ms:
    if m0 == None:
      m0min = im
      m0max = im
    else:
      dm = abs(im-m0)
      if dm > 0 and ( dm0 == None or dm < dm0 ):  dm0 = dm
      if im < m0min:  m0min = im
      if im > m0max:  m0max = im
    m0 = im;
  return ( dm0,m0min,m0max)

#
# read contamination from ROOT files (in WKs format) and interpolate
#   write result to a new dictionary
#
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--prefix", dest="prefix", default=None, type="string", action="store")
parser.add_option("-M", "--model", dest="model", default="msugra", type="string", action="store", help="signal model")
parser.add_option("--m3Ratio", dest="m3Ratio", default=-1., type="float", action="store", help="ratio for intermediate mass in SMS models")
parser.add_option("-f", dest="force", default=False, action="store_true", help="replace output file")
(options, args) = parser.parse_args()


from ROOT import gROOT
gROOT.ProcessLine(".L TriangularInterpolation.C+")

from signalUtils import buildSignalString
from signalUtils import buildSignalString

assert len(args) > 0

refFileName = args[0]
sigDict = cPickle.load(file(refFileName))

#btags = [ 'binc', 'b0', 'b1', 'b2' ]
#hts = [ 750, 1000 ]
#mets = [ 250, 350, 450, 550 ]

contDict = {}
print sigDict.keys()
for ht in sigDict:
  if not ht in contDict:  contDict[ht] = {}
  for met in sigDict[ht]:
    if not met in contDict[ht]:  contDict[ht][met] = {}
    m0s = []
    m12s = []
    btags = None
    for msugra in sigDict[ht][met]:
      fields = msugra.split('_')
      m0 = int(fields[1])
      m12 = int(fields[2])
      if not m0 in m0s:  m0s.append(m0)
      if not  m12 in m12s:  m12s.append(m12)
      if btags == None:  btags = sigDict[ht][met][msugra].keys()
    m0s.sort()
    m12s.sort()
    m0range = findRange(m0s)
    m12range = findRange(m12s)
    nb0 = int((m0range[2]-m0range[1])/float(m0range[0])+1.5)
    nb12 = int((m12range[2]-m12range[1])/float(m12range[0])+1.5)
    hInter = ROOT.TH2F("hCont","hCont", \
                       nb0,m0range[1]-m0range[0]/2.,m0range[2]+m0range[0]/2., \
                       nb12,m12range[1]-m12range[0]/2.,m12range[2]+m12range[0]/2.)

    print btags    
    for bt in btags:
      bt2 = bt if bt != 'binc' else 'inc'
      filename = "Contamination/"
      if not options.prefix == None:
        filename += options.prefix
      filename += "_"+options.model
      if options.m3Ratio > 0:  filename += str(options.m3Ratio).replace(".","")
      filename += "_"+bt+"_"+str(ht)+"_"+str(met)+".root"
      if not os.path.exists(filename):
        print "Missing file with contamination information",filename
        continue
      print filename
      file = ROOT.TFile(filename)
      hInter.Reset()
      hInter = ROOT.triangular(hInter,file)
      nbx = hInter.GetNbinsX()
      xl = hInter.GetXaxis().GetXmin()
      xh = hInter.GetXaxis().GetXmax()
      dx = (xh-xl) / nbx
      nby = hInter.GetNbinsY()
      yl = hInter.GetYaxis().GetXmin()
      yh = hInter.GetYaxis().GetXmax()
      dy = (yh-yl) / nby
      for ix in range(nbx):
        print ix
        for iy in range(nby):
          print ix,iy
          v = hInter.GetBinContent(ix+1,iy+1)
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
          if not msugra in contDict[ht][met]: contDict[ht][met][msugra] = {}
          contDict[ht][met][msugra][bt] = v
      c = ROOT.TCanvas("c","c")
      hInter.Draw("zcol")
      try:
        input("abc")
      except:
        pass

oname = "sigCont_"+options.model+"NLO"
if options.model == "msugra":
  oname += "0"
elif options.m3Ratio > 0:
  oname += str(options.m3Ratio).replace(".","")
oname += ".pkl"
if os.path.exists(oname):
  if options.force:
    print "Replacing",oname
    os.remove(oname)
  else:
    print oname,"exists"
    sys.exit(1)
    
ofile = open(oname,"wb")
cPickle.dump(contDict,ofile)
ofile.close()
