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
from ROOT import gROOT
gROOT.ProcessLine(".L TriangularInterpolation.C+")

from signalUtils import buildSignalString
from signalUtils import buildSignalString

assert len(sys.argv) > 1

refFileName = sys.argv[1]
sigDict = cPickle.load(file(refFileName))

#btags = [ 'binc', 'b0', 'b1', 'b2' ]
#hts = [ 750, 1000 ]
#mets = [ 250, 350, 450, 550 ]

contDict = {}
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
    
    for bt in btags:
      bt2 = bt if bt != 'binc' else 'inc'
      filename = "Contamination/h_cont_msugra_"+bt+"_"+str(ht)+"_"+str(met)+".root"
      if not os.path.exists(filename):  continue
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
        for iy in range(nby):
          v = hInter.GetBinContent(ix+1,iy+1)
          if v < 0.000001:  continue
          m0 = int(xl+ix*dx+0.5)
          m12 = int(yl+iy*dy+0.5)
          msugra = buildSignalString("msugra",m0,m12)
          if not msugra in contDict[ht][met]: contDict[ht][met][msugra] = {}
          contDict[ht][met][msugra][bt] = v
      c = ROOT.TCanvas("c","c")
      hInter.Draw("zcol")
      try:
        input("abc")
      except:
        pass

ofile = open("sigCont_msugraNLO0.pkl","wb")
cPickle.dump(contDict,ofile)
ofile.close()
