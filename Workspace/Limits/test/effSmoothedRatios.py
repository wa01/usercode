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
#
def findM0M12 (m0,m12,msugras):
  matches = [ ]
  for key in msugras:
    parts = key.split("_")
    if int(parts[1]) == m0 and int(parts[2]) == m12:
      matches.append(key)
  if len(matches) != 1:
    return None
  else:
    return matches[0]
  
#
# create dictionary with ratio of smoothed / raw efficiencies
#   argument: dictionary with raw efficiencies
#
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--m3Ratio", dest="m3Ratio", default=-1., type="float", action="store", help="ratio for intermediate mass in SMS models")
parser.add_option("-f", dest="force", default=False, action="store_true", help="replace output file")
(options, args) = parser.parse_args()

# dictionary with efficiencies
fnameEff = args[0]
effdic = cPickle.load(file(fnameEff))
# dictionary with yields
fnameEvt = args[1] if len(args) > 1 else None
#
# load macros
#
from ROOT import gROOT
gROOT.ProcessLine(".L SmoothingUtils.C+")
gROOT.ProcessLine(".L useNiceColorPalette.C")
ROOT.useNiceColorPalette()
ROOT.gStyle.SetOptStat(0)

# output dictionary
allEvts = {}

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

if fnameEvt != None:
  evdic = cPickle.load(file(fnameEvt))

# loop over all btag bins, ht and met values in input dictionary
#   and create corresponding levels in output dictionary
model = None
suffix = None
for btag in effdic:
  if not btag in allEvts:  allEvts[btag] = {}
  for ht in effdic[btag]:
    if not ht in allEvts[btag]:  allEvts[btag][ht] = {}
    for met in effdic[btag][ht]:
      if not met in allEvts[btag][ht]:  allEvts[btag][ht][met] = {}
      #
      # create histograms (grid points at bin centre)
      # 
      # print "btag / ht / met = ",btag,ht,met
      nb0 = (m0Range[2]-m0Range[1])/m0Range[0] + 1
      if (m0Range[2]-m0Range[1])%m0Range[0] != 0:  nb0 += 1
      fm0Min = m0Range[1] - m0Range[0]/2.
      fm0Max = fm0Min + nb0*m0Range[0]
      nb12 = (m12Range[2]-m12Range[1])/m12Range[0] + 1
      if (m12Range[2]-m12Range[1])%m12Range[0] != 0:  nb12 += 1
      fm12Min = m12Range[1] - m12Range[0]/2.
      fm12Max = fm12Min + nb12*m12Range[0]
      # print "m0 ",m0Range[0]," ",m0Range[1]," ",m0Range[2]," ; ",nb0," ",fm0Min," ",fm0Max
      # print "m12 ",m12Range[0]," ",m12Range[1]," ",m12Range[2]," ; ",nb12," ",fm12Min," ",fm12Max

      if len(processes) == 0:  processes.add( None )
      rawHistos = []
      ratioHistos = []
      evtHistos = []
      # loop over processes
      for p in processes:
        # create histogram
        hname = "hraw"
        if p != None:  hname += p
        hRaw = ROOT.TH2F(hname,hname,nb0,fm0Min,fm0Max,nb12,fm12Min,fm12Max)
        hEvt = None
        if fnameEvt != 0 and btag in evdic  and ht in evdic[btag] and met in evdic[btag][ht]:
          hname = "hevt"
          if p != None:  hname += p
          hEvt = ROOT.TH2F(hname,hname,nb0,fm0Min,fm0Max,nb12,fm12Min,fm12Max)
        if hEvt == None:  continue
        #  fill histogram with valid points (add entries / process to output dictionary, if necessary)
        for msugra in effdic[btag][ht][met]:
          parts = msugra.split('_')
          # check model name
          if model == None:
            model = parts[0]
          else:
            assert parts[0] == model
          # check remaining fields
          if suffix == None:  suffix = parts[3:]
          if model == 'msugra':
            assert len(parts) == 6 and parts[3:] == suffix
          elif model.startswith('T'):
            assert len(parts) == 4
            if options.m3Ratio < 0:
              if parts[3:] != suffix:
                print "Varying suffix ",parts[3:],' / ',suffix,': probably missing m3Ratio argument'
                sys.exit(1)
            else:
              m0 = int(parts[1])
              m12 = int(parts[2])
              m3 = options.m3Ratio*(m0-m12) + m12
              if abs(int(m3+0.5)-int(parts[3])) > 1:
                print "Intermediate mass for ",msugra,"is inconsistent with m3Ratio =",options.m3Ratio
                sys.exit(1)
          #
          v = None
          if p == None:
            v = effdic[btag][ht][met][msugra]
          elif p in effdic[btag][ht][met][msugra]:
            v = effdic[btag][ht][met][msugra][p]
            allEvts[btag][ht][met][msugra][p] = {}
          if not v == None and not math.isnan(v):
            if v<0.000001: v = 0.000001
            hRaw.Fill(float(parts[1]),float(parts[2]),v)
          v = None
          if hEvt != None:
            if p == None:
              v = evdic[btag][ht][met][msugra]
            elif p in evdic[btag][ht][met][msugra]:
              v = evdic[btag][ht][met][msugra][p]
              allEvts[btag][ht][met][msugra][p] = {}
            if not v == None and not math.isnan(v):
              hEvt.Fill(float(parts[1]),float(parts[2]),v)
        #
        # create histogram with ratios and fill bin contents into dictionary
        #
        hEvt = ROOT.doEffFit(hRaw,hEvt)
#        try:
#          input("abc")
#        except:
#          pass
        nbx = hEvt.GetNbinsX()
        nby = hEvt.GetNbinsY()
        for ix in range(1,nbx+1):
          for iy in range(1,nby+1):
            ratio = hEvt.GetBinContent(ix,iy)
            if abs(ratio) < 0.000001:  continue
            m0 = int(hEvt.GetXaxis().GetBinCenter(ix)+0.5)
            m12 = int(hEvt.GetYaxis().GetBinCenter(iy)+0.5)
#            msugra = 'msugra_'+str(m0)+'_'+str(m12)+'_10_0_1'
            msugra = model+"_"+str(m0)+"_"+str(m12)
            if model.startswith('T3w') and options.m3Ratio > 0:
              m3 = options.m3Ratio*(m0-m12) + m12
              msugra += "_"+str(int(m3+0.5))
#              msugra = findM0M12(m0,m12,effdic[btag][ht][met])
            else:
              for part in suffix:  msugra += "_"+part
            if p == None:
              allEvts[btag][ht][met][msugra] = ratio
            else:
              allEvts[btag][ht][met][msugra][p] = ratio
        # clear root memory before creating the next histogram
        gROOT.Clear()
#        gROOT.ls()
##
## rescale event yields
##
#if fnameEvt != None:
#  evdic = cPickle.load(file(fnameEvt))
#  allEvts = {}
#  for btag in allRatios:
#    if not btag in evdic:  continue
#    allEvts[btag] = {}
#    for ht in allRatios[btag]:
#      if not ht in evdic[btag]:  continue
#      allEvts[btag][ht] = {}
#      for met in allRatios[btag][ht]:
#        if not met in evdic[btag][ht]:  continue
#        allEvts[btag][ht][met] = {}
#        for msugra in allRatios[btag][ht][met]:
#          if model.startswith('T3w') and options.m3Ratio > 0:
#            parts = msugra.split("_")
#            m0 = int(parts[1])
#            m12 = int(parts[2])
#            msugraIn = findM0M12(m0,m12,evdic[btag][ht][met])
#          else:
#            msugraIn = msugra
#          if msugraIn == None or not msugraIn in evdic[btag][ht][met]:  continue
#          allEvts[btag][ht][met][msugra] = {}
#          for p in processes:
#            if p == None:
#              allEvts[btag][ht][met][msugra] = evdic[btag][ht][met][msugraIn]*allRatios[btag][ht][met][msugra]
#            else:
#              if not p in evdic[btag][ht][met][msugraIn]:  continue
#              allEvts[btag][ht][met][msugra][p] = evdic[btag][ht][met][msugraIn][p]* \
#                                                  allRatios[btag][ht][met][msugra][p]
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
  if options.force:
    print "Replacing",fname
    os.remove(fname)
  else:
    print fname," exists"
    sys.exit(1)
  
f = open(fname,"w")
if fnameEvt == None:
  cPickle.dump(allRatios,f)
else:
  cPickle.dump(allEvts,f)
f.close()

