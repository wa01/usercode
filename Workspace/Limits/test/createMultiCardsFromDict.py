#!/usr/bin/env python
import cPickle
import math
import sys
import re


#
# writes Higgs combination program compatible cards files
#   for multi-channel analysis
#
leptons = [ 'Mu', 'Ele' ]
def sqr (val):
    return val*val
def getCountsLep (countDict,lep,btag):
    result = 0
    bt = btag    
    if bt == 'binc' or bt == 'b0' or bt == 'b1' or bt == 'b2':
        result += countDict[options.ht][options.met][bt][lep]
    elif bt == 'b1p':
        for bt1 in [ 'b1', 'b2' ]:
            result += countDict[options.ht][options.met][bt1][lep]
    else:
        print "Unknown b-tag bin ",btag
        sys.exit(1)
    return result
def getCounts (countDict,btag):
    result = 0
    bt = btag
    for lep in leptons:
        result += getCountsLep(countDict,lep,btag)
    return result
def getErrors (errorDict,btag):
    result = 0.
    bt = btag    
    if bt == 'binc' or bt == 'b0' or bt == 'b1' or bt == 'b2':
        for lep in leptons:  result += sqr(errorDict[options.ht][options.met][bt][lep])
    elif bt == 'b1p':
        for bt1 in [ 'b1', 'b2' ]:
            for lep in leptons:  result += sqr(errorDict[options.ht][options.met][bt1][lep])
    else:
        print "Unknown b-tag bin ",btag
        sys.exit(1)
    return math.sqrt(result)
def sumLepErrors (key,bt,predLep):
    sumPred = 0
    sumPredErr = 0
    for lep in leptons:
        err = largestAbsDoubleRatioDeviation[key][lep][bt]
        err2 = largestAbsSingleRatioDeviation[key][lep][bt]
        if err2 < err:
            print "single < double ratio for ",key,lep,btag,err2,err
            err = err2
        sumPred += predLep[lep]
        sumPredErr += predLep[lep]*err
    return sumPredErr/sumPred
def sumBTErrors (keyUp,keyDown,bt,predLep):
    # (signed) variation / lepton channel
    sumPred = 0.
    sumPredErr = 0.
    for lep in leptons:
        # get variations w.r.t. 1. (revert sign for down)
        dUp = singleRatio[keyUp][lep][bt] - 1.
        dDown = -(singleRatio[keyDown][lep][bt]-1.)
        # take the average of up/down and the maximum of Mu/Ele
        errAve = (dUp+dDown)/2.
        sumPred += predLep[lep]
        sumPredErr += predLep[lep]*errAve
    return sumPredErr/sumPred

#from sys import argv
import os.path
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--ht", dest="ht", type="int", action="store", help="HT cut")
parser.add_option("--met", dest="met", type="int", action="store", help="MET cut")
#parser.add_option("--btag", dest="btag", default="binc", type="string", action="store")
(options, args) = parser.parse_args()

btags = [ "b0", "b1", "b2" ]

#
# output file name
#
oname = "multibtag-ht" + str(options.ht) + "-met" + str(options.met) + ".txt"
ofile = open(oname,"w")
#
# numbers of channels and nuisance parameters
#
ofile.write("imax"+"%3d" % len(btags) +"\n")
ofile.write("jmax *\n")
ofile.write("kmax *\n")
#ofile.write("jmax"+"%3d" % 1 +"\n")
#ofile.write("kmax"+"%3d" % 4 +"\n")
#
# read observed and predicted background numbers
#
execfile("eventCounts.py")
obs = {}
pred = {}
predLep = {}
errpred = {}
for btag in btags:
    obs[btag] = getCounts(countsObs,btag)
    pred[btag] = getCounts(countsPred,btag)
    predLep[btag] = {}
    for lep in leptons:  predLep[btag][lep] = getCountsLep(countsPred,lep,btag)
    errpred[btag] = getErrors(errorsPred,btag)
#
# observations
#
lbin = "bin".ljust(15)
lobs  = "observation".ljust(15)
for btag in btags:
    lbin = lbin + btag.rjust(10)
    lobs = lobs + "%10d" % obs[btag]
ofile.write(lbin+"\n")
ofile.write(lobs+"\n")
ofile.write("#\n")

# very simple estimate for btag correlations assuming 2 bjets and fixed eff
beff = 0.60
dbeff = 0.06
beffsyst = {}
beffsyst["b0"] = -2*dbeff/(1-beff)
beffsyst["b1p"] = 2*(1-beff)*dbeff/beff/(2-beff)
beffsyst["b1"] = -dbeff*(2*beff-1)/beff/(1-beff)
beffsyst["b2"] = 2*dbeff/beff
#
# predictions (signal rates just for template - they will be replaced
#   when creating the final files)
#
lbin = "bin".ljust(15)
lproc1 = "process".ljust(15)
lproc2 = "process".ljust(15)
lrate = "rate".ljust(15)
for btag in btags:
    lbin = lbin + btag.rjust(10) + btag.rjust(10)
    lproc1 = lproc1 + "susy".rjust(10) + "bkg".rjust(10)
    lproc2 = lproc2 + "%10d" % 0 + "%10d" % 1
    lrate = lrate + "%10.3f" % 1. + "%10.3f" % pred[btag]
ofile.write(lbin+"\n")
ofile.write(lproc1+"\n")
ofile.write(lproc2+"\n")
ofile.write(lrate+"\n")
ofile.write("#\n")
#
# uncertainties: luminosity, other signal systematics, background statistics
#
lumiSyst = 1.045
sigSyst  = 1.20
llumi = "lumi".ljust(10) + "lnN".ljust(5)
lsigsyst = "sigSyst".ljust(10) + "lnN".ljust(5)
lsigbtag = "sigBTag".ljust(10) + "lnN".ljust(5)
lbkgstats = {}
for btag in btags:
    # background statistics (uncorrelated)
    lbkgstats[btag] = (btag+"Stat").ljust(10) + "lnN".ljust(5)
for btag in btags:
    # luminosity (correlated in all signal channels)
    llumi = llumi + "%10.3f" % lumiSyst + "-".rjust(10)
    # global signal systematics (correlated in all signal channels)
    lsigsyst = lsigsyst + "%10.3f" % sigSyst + "-".rjust(10)
    # btag (signal)
    lsigbtag = lsigbtag + "%10.3f" % (1+beffsyst[btag]) + "-".rjust(10)
    for btag1 in btags:
        lbkgstats[btag1] = lbkgstats[btag1] + "-".rjust(10)
        if btag1 == btag:
            lbkgstats[btag] = lbkgstats[btag] + "%10.3f" % (1+errpred[btag]/pred[btag])
        else:
            lbkgstats[btag1] = lbkgstats[btag1] + "-".rjust(10)
ofile.write(llumi+"\n")
ofile.write(lsigsyst+"\n")
ofile.write(lsigbtag+"\n")
for btag in btags:  ofile.write(lbkgstats[btag]+"\n")
#
# systematics (non-b related)
#
inf = float('inf')
sname = "systematics_htSig-" + str(options.ht) + "_metSig-" + str(options.met) + ".py"
execfile(sname)

sumerr = {}
lbkgsyst = "bkgSyst".ljust(10) + "lnN".ljust(5)
for btag in btags:
    # use 'b1p' for 'b1' and 'b2'
    bt = btag
    if btag == 'b1' or btag == 'b2':  bt = 'b1p'
    sumerr[btag] = 0
    for key in largestAbsDoubleRatioDeviation:
        err = sumLepErrors(key,bt,predLep[btag])
        sumerr[btag] += err*err
        print key,btag,err

#
# systematics (b-tag related)
#
sname = "systematics_BT_htSig-" + str(options.ht) + "_metSig-" + str(options.met) + ".py"
execfile(sname)
#
# W/tt ratio
#
for btag in btags:
    key = 'ScaleFrac'
    err = sumLepErrors(key,btag,predLep[btag])
    sumerr[btag] += err*err
    lbkgsyst = lbkgsyst + "-".rjust(10) + "%10.3f" % (1+math.sqrt(sumerr[btag]))
ofile.write(lbkgsyst+"\n")
#
# b-tag systs
#
# keys for up and down variations / source
#
btKeys = { 'beff' : [ 'btagEff3_Up_b_sf0', 'btagEff3_Down_b_sf0' ], \
           'leff' : [ 'btagEff3_Up_l_sf0', 'btagEff3_Down_l_sf0' ] }
for vari in btKeys:
    lbsyst = vari.ljust(10) + "lnN".ljust(5)
    for btag in btags:
        # keys
        keyUp = btKeys[vari][0]
        keyDown = btKeys[vari][1]
        # (signed) variation / lepton channel
        errLep = sumBTErrors(keyUp,keyDown,btag,predLep[btag])
        lbsyst = lbsyst + "-".rjust(10) + "%10.3f" % (1+errLep)
    ofile.write(lbsyst+"\n")
    

ofile.close()

