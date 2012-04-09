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

def sigString(m0,m12,model,m3Ratio=-1):
  msugra = None  
  if model == "T1tttt":
    msugra = model+"_"+str(m0)+"_"+str(m12)+"_-1"
  elif model.startswith("T3w") and m3Ratio > 0:
    m3 = int(m3Ratio*(m0-m12)+m12+0.5)
    msugra = model+"_"+str(m0)+"_"+str(m12)+"_"+str(m3)
  else:
    print "Unknown model",model
    sys.exit(1)
  return msugra
#
# extend JES systematics to low m12 (copy for equal dm)

from optparse import OptionParser
parser = OptionParser()
parser.add_option("-M", "--model", dest="model", default="msugra", type="string", action="store", help="signal model")
parser.add_option("--m3Ratio", dest="m3Ratio", default=-1., type="float", action="store", help="ratio for intermediate mass in SMS models")
(options, args) = parser.parse_args()


from signalUtils import buildSignalString

hts = [ 750, 1000 ]
mets = [ 250, 350, 450, 550 ]

assert len(args) > 1
inFileName = args[0]
outFileName = args[1]
if os.path.exists(outFileName):
  print "File exists",outFileName
  sys.exit(1)
  
jesDict = cPickle.load(file(inFileName))


for ht in jesDict:
  for met in jesDict[ht]:

    mDict = {}
    m0s = [ ]
    m12s = [ ]
    m0Max = None
    m12Min = None
    for msugra in jesDict[ht][met]:
      parts = msugra.split("_")
      m0 = int(parts[1])
      if not m0 in m0s:  m0s.append(m0)
      if m0Max == None or m0 > m0Max:  m0Max = m0
      m12 = int(parts[2])
      if not m12 in m12s:  m12s.append(m12)
      if m12Min == None or m12 < m12Min:  m12Min = m12
      if not m12 in mDict:  mDict[m12] = [ ]
      if not m0 in mDict[m12]:  mDict[m12].append(m0)
    for m12 in mDict:  mDict[m12].sort()
    m0s.sort()
    m12s.sort()
    m0range = findRange(m0s)
    m12range = findRange(m12s)

    print m0range
    print m12range
    
    m12a = m12Min
    while m12a > m12range[0]:
      m12b = m12a - m12range[0]
      mDict[m12b] = [ ]
      for m0a in mDict[m12a]:
        msugraA = sigString(m0a,m12a,options.model,options.m3Ratio)
        m0b = m0a - m0range[0]
        if m0b < 0: continue
        msugraB = sigString(m0b,m12b,options.model,options.m3Ratio)
        jesDict[ht][met][msugraB] = jesDict[ht][met][msugraA]
        mDict[m12b].append(m0b)
        if m0a == m0Max:
          m0b = m0a
          msugraB = sigString(m0b,m12b,options.model,options.m3Ratio)
          jesDict[ht][met][msugraB] = jesDict[ht][met][msugraA]
          mDict[m12b].append(m0b)
      mDict[m12b].sort()
      m12a = m12b
    
    outFile = open(outFileName,"wb")
    cPickle.dump(jesDict,outFile)
    outFile.close()
