import pickle
import ROOT
import math
import sys
#
# smoothed b-tag systematic efficiencies
#

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
# create dictionary with ratio of smoothed / raw efficiencies
#   argument: dictionary with raw efficiencies
#
filename = sys.argv[1]
mydic = pickle.load(file(filename))

from ROOT import gROOT
gROOT.ProcessLine(".L EfficiencyUtils.C+")

# output dictionary
newdic = {}

m0s = [ ]
m12s = [ ]
labels = [ ]
# loop over all ht and met values
# in input dictionary
for btag in mydic:
  for ht in mydic[btag]:
    for met in mydic[btag][ht]:
      # create lists of M0 and M12 values and
      #   set of process names, if applicable
      for msugra in mydic[btag][ht][met]:
        # extract M0 and M12 from msugra string
        parts = msugra.split('_')
        m0 = int(parts[1])
        if not m0 in m0s:  m0s.append(m0)
        m12 = int(parts[2])
        if not m12 in m12s:  m12s.append(m12)
#
# extraction of ranges and bin widths
#
m0s.sort()
m12s.sort()
if len(m0s) < 2 or len(m12s) < 2:
  print "too few points"
  sys.exit(1)
#  print m0s
#  print m12s        
m0Range = findRange(m0s)
m12Range = findRange(m12s)
if m0Range[0] == 0 or m12Range[0] == 0:
  print "dm0 or dm12 == 0 ??"
  sys.exit(1)

print m0Range
print m12Range

# loop over all btag bins, ht and met values in input dictionary
#   and create corresponding levels in output dictionary
for btag in mydic:
  if not btag in newdic:  newdic[btag] = {}
  for ht in mydic[btag]:
    if not ht in newdic[btag]:  newdic[btag][ht] = {}
    for met in mydic[btag][ht]:
      if not met in newdic[btag][ht]:  newdic[btag][ht][met] = {}
      #
      # create histograms (grid points at bin centre)
      # 
      print "ht / met = ",ht,met
      nb0 = (m0Range[2]-m0Range[1])/m0Range[0] + 1
      if (m0Range[2]-m0Range[1])%m0Range[0] != 0:  nb0 += 1
      fm0Min = m0Range[1] - m0Range[0]/2.
      fm0Max = fm0Min + nb0*m0Range[0]
      nb12 = (m12Range[2]-m12Range[1])/m12Range[0] + 1
      if (m12Range[2]-m12Range[1])%m12Range[0] != 0:  nb12 += 1
      fm12Min = m12Range[1] - m12Range[0]/2.
      fm12Max = fm12Min + nb12*m12Range[0]
      print "m0 ",m0Range[0]," ",m0Range[1]," ",m0Range[2]," ; ",nb0," ",fm0Min," ",fm0Max
      print "m12 ",m12Range[0]," ",m12Range[1]," ",m12Range[2]," ; ",nb12," ",fm12Min," ",fm12Max

      #  print "h","h",nb0,fm0Min,fm0Max,nb12,fm12Min,fm12Max
      #  fill histogram with valid points (add entries / process to output dictionary, if necessary)
      hRaws = {}
      hRatios = {}
      for msugra in mydic[btag][ht][met]:
        for label in mydic[btag][ht][met][msugra]:
          # create histogram
          if not label in hRaws:
            hname = "hraw" + label
            hRaws[label] = ROOT.TH2F(hname,hname,nb0,fm0Min,fm0Max,nb12,fm12Min,fm12Max)
          parts = msugra.split('_')
          v = mydic[btag][ht][met][msugra][label]
          if not v == None and not math.isnan(v):
            hRaws[label].Fill(float(parts[1]),float(parts[2]),v)
      #
      # create histogram with ratios and fill bin contents into dictionary
      #
      labels = {}
      for label in hRaws:
        hRatio = ROOT.doEff(hRaws[label],0)
        nbx = hRatio.GetNbinsX()
        nby = hRatio.GetNbinsY()
        for ix in range(1,nbx):
          for iy in range(1,nby):
            ratio = hRatio.GetBinContent(ix,iy)
            if abs(ratio) < 0.00000001:  continue
            if not label in labels:  labels[label] = 0
            labels[label] += 1
            m0 = int(hRatio.GetXaxis().GetBinCenter(ix)+0.5)
            m12 = int(hRatio.GetYaxis().GetBinCenter(iy)+0.5)
            msugra = 'msugra_'+str(m0)+'_'+str(m12)+'_10_0_1'
            if not msugra in newdic[btag][ht][met]: newdic[btag][ht][met][msugra] = {}
            newdic[btag][ht][met][msugra][label] = ratio  #*hRaws[label].GetBinContent(ix,iy)
      # clear root memory before creating the next histogram
#      print btag,ht,met
#      gROOT.ls()
      gROOT.Clear()
#      print btag,ht,met,labels
#        gROOT.ls()
#
# write output dictionary (derive name from input file name)
#
import os.path
fname = filename
ind = fname.rfind("/")
if ind != -1:  fname = fname[ind+1:]
ind = fname.find(".")
if ind != -1:  fname = fname[0:ind]
fname += "-smoothed.pkl"
if os.path.exists(fname):
  print fname," exists"
  sys.exit(1)
  
f = open(fname,"w")
pickle.dump(newdic,f)
f.close()

