import cPickle
import ROOT
import math
import sys
#
# read efficiencies from dictionary, smooth using
#   the SmoothingUtils.C macros and write ratios
#   to new dictionary. If 2 arguments: rescale
#   event counts in 2nd file with ratios
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
# dictionary with efficiencies
fnameEff = sys.argv[1]
effdic = cPickle.load(file(fnameEff))
# dictionary with yields
fnameEvt = None
if len(sys.argv) > 2:
  fnameEvt = sys.argv[2]
#
# load macros
#
from ROOT import gROOT
gROOT.ProcessLine(".L SmoothingUtils.C+")

# output dictionary
allRatios = {}

# loop over all btag bins, ht and met values
# in input dictionary
for btag in effdic:
  for ht in effdic[btag]:
    for met in effdic[btag][ht]:
      # create lists of M0 and M12 values and
      #   set of process names, if applicable
      m0s = [ ]
      m12s = [ ]          
      processes = set([ ])
      for msugra in effdic[btag][ht][met]:
        # extract M0 and M12 from msugra string
        parts = msugra.split('_')
        m0 = int(parts[1])
        if not m0 in m0s:  m0s.append(m0)
        m12 = int(parts[2])
        if not m12 in m12s:  m12s.append(m12)
        # if dictionary contains processes: add to list
        if type(effdic[btag][ht][met][msugra]) is dict:
          processes.update(effdic[btag][ht][met][msugra].keys())
          
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
print processes

# loop over all btag bins, ht and met values in input dictionary
#   and create corresponding levels in output dictionary
for btag in effdic:
  if not btag in allRatios:  allRatios[btag] = {}
  for ht in effdic[btag]:
    if not ht in allRatios[btag]:  allRatios[btag][ht] = {}
    for met in effdic[btag][ht]:
      if not met in allRatios[btag][ht]:  allRatios[btag][ht][met] = {}
      #
      # create histograms (grid points at bin centre)
      # 
      print "btag / ht / met = ",btag,ht,met
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

      if len(processes) == 0:  processes.add( None )
      rawHistos = []
      ratioHistos = []
      # loop over processes
      for p in processes:
        # create histogram
        hname = "hraw"
        if p != None:  hname += p
        hRaw = ROOT.TH2F(hname,hname,nb0,fm0Min,fm0Max,nb12,fm12Min,fm12Max)
        #  fill histogram with valid points (add entries / process to output dictionary, if necessary)
        msugraStrings = {}
        for msugra in effdic[btag][ht][met]:
          parts = msugra.split('_')
          m0m12 = ( int(parts[1]), int(parts[2]) )
          v = None
          if p == None:
            v = effdic[btag][ht][met][msugra]
          elif p in effdic[btag][ht][met][msugra]:
            v = effdic[btag][ht][met][msugra][p]
            allRatios[btag][ht][met][msugra][p] = {}
          if not v == None and not math.isnan(v):
            hRaw.Fill(float(parts[1]),float(parts[2]),v)
          msugraStrings[m0m12] = msugra
        #
        # create histogram with ratios and fill bin contents into dictionary
        #
        hRatio = ROOT.doEff(hRaw)
        nbx = hRatio.GetNbinsX()
        nby = hRatio.GetNbinsY()
        for ix in range(1,nbx):
          for iy in range(1,nby):
            ratio = hRatio.GetBinContent(ix,iy)
            if abs(ratio) < 0.000001:  continue
            m0 = int(hRatio.GetXaxis().GetBinCenter(ix)+0.5)
            m12 = int(hRatio.GetYaxis().GetBinCenter(iy)+0.5)
#            msugra = 'msugra_'+str(m0)+'_'+str(m12)+'_10_0_1'
            msugra = msugraStrings[ (m0, m12) ]
            if p == None:
              allRatios[btag][ht][met][msugra] = ratio
            else:
              allRatios[btag][ht][met][msugra][p] = ratio
        # clear root memory before creating the next histogram
        gROOT.Clear()
#        gROOT.ls()
#
# rescale event yields
#
if fnameEvt != None:
  evdic = cPickle.load(file(fnameEvt))
  allEvts = {}
  for btag in allRatios:
    if not btag in evdic:  continue
    allEvts[btag] = {}
    for ht in allRatios[btag]:
      if not ht in evdic[btag]:  continue
      allEvts[btag][ht] = {}
      for met in allRatios[btag][ht]:
        if not met in evdic[btag][ht]:  continue
        allEvts[btag][ht][met] = {}
        for msugra in allRatios[btag][ht][met]:
          if not msugra in evdic[btag][ht][met]:  continue
          allEvts[btag][ht][met][msugra] = {}
          for p in processes:
            if p == None:
              allEvts[btag][ht][met][msugra] = evdic[btag][ht][met][msugra]*allRatios[btag][ht][met][msugra]
            else:
              if not p in evdic[btag][ht][met][msugra]:  continue
              allEvts[btag][ht][met][msugra][p] = evdic[btag][ht][met][msugra][p]* \
                                                  allRatios[btag][ht][met][msugra][p]
            
#
# write output dictionary (derive name from input file name)
#
import os.path
if fnameEvt == None:
  fname = fnameEff
  ind = fname.rfind("/")
  if ind != -1:  fname = fname[ind+1:]
  ind = fname.find(".")
  if ind != -1:  fname = fname[0:ind]
  fname += "-Ratio.pkl"
else:
  fname = fnameEvt
  ind = fname.rfind("/")
  if ind != -1:  fname = fname[ind+1:]
  ind = fname.find(".")
  if ind != -1:  fname = fname[0:ind]
  fname += "-smoothed.pkl"
if os.path.exists(fname):
  print fname," exists"
  sys.exit(1)
  
f = open(fname,"w")
if fnameEvt == None:
  cPickle.dump(allRatios,f)
else:
  cPickle.dump(allEvts,f)
f.close()

