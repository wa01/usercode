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
        err = doubleRatio[key][lep][bt] - 1
        err2 = singleRatio[key][lep][bt] - 1
        if abs(err2) < abs(err):
            print "single < double ratio for ",key,lep,btag,err2,err
            err = err2
        sumPred += predLep[lep]
        sumPredErr += predLep[lep]*err
        print "JES correction ",key,lep,btag,err
    print "JES correction ",key,btag,sumPredErr/sumPred
    return sumPredErr/sumPred
def sumAbsLepErrors (key,bt,predLep):
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
        dUp1 = singleRatio[keyUp][lep][bt] - 1.
        dUp2 = doubleRatio[keyUp][lep][bt] - 1.
        dUp = dUp1 if abs(dUp1) < abs(dUp2) else dUp2
        dDown1 = singleRatio[keyDown][lep][bt] - 1.
        dDown2 = doubleRatio[keyDown][lep][bt] - 1.
        dDown = dDown1 if abs(dDown1) < abs(dDown2) else dDown2
        # take the average of up/down and the maximum of Mu/Ele
        errAve = (dUp-dDown)/2.
        sumPred += predLep[lep]
        sumPredErr += predLep[lep]*errAve
    return sumPredErr/sumPred

#from sys import argv
import os.path
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--ht", dest="ht", type="int", action="store", help="HT cut")
parser.add_option("--met", dest="met", type="int", action="store", help="MET cut")
parser.add_option("--btag", dest="btag", default="", type="string", action="store")
(options, args) = parser.parse_args()

btagname = "multibtag"
btags = [ "b0", "b1", "b2" ]
if options.btag != "":
    btagname = options.btag
    btags = [ options.btag ]
#
# MC closure
#
errMCClosure = { 750 : { 250 : 0.072, 350 : 0.041, 450 : 0.073, 550 : 0.224 }, \
                 1000 : { 250 : 0.056, 350 : 0.074, 450  : 0.098, 550 : 0.198 } }
#
# output file name
#
oname = "Backgrounds/" + btagname + "-ht" + str(options.ht) + "-met" + str(options.met) + ".txt"
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
llumi = "lumi".ljust(15) + "lnN".ljust(5)
lsigsyst = "sigSyst".ljust(15) + "lnN".ljust(5)
if options.btag != 'binc':
    lsigbtag = "sigBTag".ljust(15) + "lnN".ljust(5)
lbkgstats = {}
for btag in btags:
    # background statistics (uncorrelated)
    lbkgstats[btag] = (btag+"Stat").ljust(15) + "lnN".ljust(5)
for btag in btags:
    # luminosity (correlated in all signal channels)
    llumi = llumi + "%10.3f" % lumiSyst + "-".rjust(10)
    # global signal systematics (correlated in all signal channels)
    lsigsyst = lsigsyst + "%10.3f" % sigSyst + "-".rjust(10)
    # btag (signal)
    if options.btag != 'binc':
        lsigbtag = lsigbtag + "%10.3f" % (1+beffsyst[btag]) + "-".rjust(10)
    for btag1 in btags:
        lbkgstats[btag1] = lbkgstats[btag1] + "-".rjust(10)
        if btag1 == btag:
            lbkgstats[btag] = lbkgstats[btag] + "%10.3f" % (1+errpred[btag]/pred[btag])
        else:
            lbkgstats[btag1] = lbkgstats[btag1] + "-".rjust(10)
ofile.write(llumi+"\n")
ofile.write(lsigsyst+"\n")
if options.btag != 'binc':
    ofile.write(lsigbtag+"\n")
for btag in btags:  ofile.write(lbkgstats[btag]+"\n")
#
# systematics (non-b related)
#
inf = float('inf')
sname = "Systematics/systematics_htSig-" + str(options.ht) + "_metSig-" + str(options.met) + ".py"
execfile(sname)

sumerr = {}
jeserr = {}
lbkgsyst = "bkgSyst".ljust(15) + "lnN".ljust(5)
lbkgjes = "bkgSystJES".ljust(15) + "lnN".ljust(5)
for btag in btags:
    # use 'b1p' for 'b1' and 'b2'
    bt = btag
    if bt == 'binc':
        bt = 'inc'
    elif btag == 'b1' or btag == 'b2':
        bt = 'b1p'
    sumerr[btag] = 0
    for key in largestAbsDoubleRatioDeviation:
        if key == 'JES':  continue
        err = sumAbsLepErrors(key,bt,predLep[btag])
        sumerr[btag] += err*err
        print key,btag,err
    jeserr[btag] = 0
    jeserrpm = []
    for key in [ 'pfRA4TupelizerJESMinus', 'pfRA4TupelizerJESPlus' ]:
        jeserrpm.append(sumLepErrors(key,bt,predLep[btag]))
    if abs(jeserrpm[0]) > abs(jeserrpm[1]):
        jeserr[btag] = jeserrpm[0]
    else:
        jeserr[btag] = jeserrpm[1]
    print "JES correction ",btag,jeserr[btag]
        
#
# add non-closure to other errors
#
errClos = errMCClosure[options.ht][options.met]
print "Adding MC closure syst ",errClos
for btag in btags:
    sumerr[btag] += errClos*errClos
#
# write JES errors
#
for btag in btags:
    bt = btag
    if bt == 'binc':  bt = 'inc'
    lbkgjes = lbkgjes + "-".rjust(10) + "%10.3f" % (1+jeserr[btag])
ofile.write(lbkgjes+"\n")
#
# systematics (b-tag related)
#
sname = "Systematics/systematics_BT_htSig-" + str(options.ht) + "_metSig-" + str(options.met) + ".py"
execfile(sname)
#
# W/tt ratio
#
for btag in btags:
    bt = btag
    if bt == 'binc':  bt = 'inc'
    key = 'ScaleFrac'
    err = sumAbsLepErrors(key,bt,predLep[btag])
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
    lbsyst = vari.ljust(15) + "lnN".ljust(5)
    for btag in btags:
        bt = btag
        if bt == 'binc':  bt = 'inc'
        # keys
        keyUp = btKeys[vari][0]
        keyDown = btKeys[vari][1]
        # (signed) variation / lepton channel
        errLep = sumBTErrors(keyUp,keyDown,bt,predLep[btag])
        lbsyst = lbsyst + "-".rjust(10) + "%10.3f" % (1+errLep)
    ofile.write(lbsyst+"\n")
    

ofile.close()

