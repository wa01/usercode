#!/usr/bin/env python
import cPickle
import math
import sys
#
# writes Higgs combination program compatible cards files
#
leptons = [ 'Mu', 'Ele' ]
def sqr (val):
    return val*val
def getCountsLep (countDict,lep):
    result = 0
    bt = options.btag    
    if bt == 'binc' or bt == 'b0' or bt == 'b1' or bt == 'b2':
        result += countDict[options.ht][options.met][bt][lep]
    elif bt == 'b1p':
        for bt1 in [ 'b1', 'b2' ]:
            result += countDict[options.ht][options.met][bt1][lep]
    else:
        print "Unknown b-tag bin ",options.btag
        sys.exit(1)
    return result
def getCounts (countDict):
    result = 0
    bt = options.btag
    for lep in leptons:
        result += getCountsLep(countDict,lep)
    return result
def getErrors (errorDict):
    result = 0.
    bt = options.btag    
    if bt == 'binc' or bt == 'b0' or bt == 'b1' or bt == 'b2':
        for lep in leptons:  result += sqr(errorDict[options.ht][options.met][bt][lep])
    elif bt == 'b1p':
        for bt1 in [ 'b1', 'b2' ]:
            for lep in leptons:  result += sqr(errorDict[options.ht][options.met][bt1][lep])
    else:
        print "Unknown b-tag bin ",options.btag
        sys.exit(1)
    return math.sqrt(result)
def sumLepErrors (key,bt,predLep):
    sumPred = 0
    sumPredErr = 0
    for lep in leptons:
        err = largestAbsDoubleRatioDeviation[key][lep][bt]
        err2 = largestAbsSingleRatioDeviation[key][lep][bt]
        if err2 < err:
            print "single < double ratio for ",key,lep,options.btag,err2,err
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

import re
#from sys import argv
import os.path
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--ht", dest="ht", type="int", action="store", help="HT cut")
parser.add_option("--met", dest="met", type="int", action="store", help="MET cut")
parser.add_option("--btag", dest="btag", default="binc", type="string", action="store", help="btag bin name")
(options, args) = parser.parse_args()
#
# MC closure
#
errMCClosure = { 750 : { 250 : 0.072, 350 : 0.041, 450 : 0.073, 550 : 0.224 }, \
                 1000 : { 250 : 0.056, 350 : 0.074, 450  : 0.098, 550 : 0.198 } }
#
#
# output file name
#
oname = options.btag + "-ht" + str(options.ht) + "-met" + str(options.met) + ".txt"
ofile = open(oname,"w")
#
# numbers of channels and nuisance parameters
#
ofile.write("imax"+"%3d" % 1 +"\n")
ofile.write("jmax"+"%3d" % 1 +"\n")
if options.btag == 'binc':
    ofile.write("kmax"+"%3d" % 4 +"\n")
else:
    ofile.write("kmax"+"%3d" % 7 +"\n")
#
# read observed and predicted background numbers
#
execfile("eventCounts.py")
obs = getCounts(countsObs)
pred = getCounts(countsPred)
predLep = {}
for lep in leptons:  predLep[lep] = getCountsLep(countsPred,lep)
errpred = getErrors(errorsPred)
    
# very simple estimate for btag correlations assuming 2 bjets and fixed eff
beff = 0.60
dbeff = 0.06
beffsyst = {}
beffsyst["b0"] = -2*dbeff/(1-beff)
beffsyst["b1p"] = 2*(1-beff)*dbeff/beff/(2-beff)
beffsyst["b1"] = -dbeff*(2*beff-1)/beff/(1-beff)
beffsyst["b2"] = 2*dbeff/beff
#
# one bin with signal and one background
#   signal rate just for template - it will be replaced when creating the final files
#
ofile.write("bin 1\n")
ofile.write("observation "+str(obs)+"\n")
ofile.write("bin      " + "%9d" % 1 + "%9d" % 1 + "\n")
ofile.write("process  " + "%9s" % "susy" + "%9s" % "bkg" + "\n")
ofile.write("process  " + "%9d" % 0 + "%9d" % 1 + "\n")
ofile.write("rate     " + "%9.3f" % 1. + "%9.3f" % pred + "\n")
#
# uncertainties: luminosity, other signal systematics, background statistics
#
ofile.write("lumi    lnN    1.045   -     \n")
ofile.write("sigSyst lnN    1.20    -   \n")
ofile.write("bkgStat lnN     -   " + "%7.3f" % (1+errpred/pred) + "\n")
if options.btag != 'binc':
    ofile.write("sigBTag lnN" + "%9.3f" % (1+abs(beffsyst[options.btag])) + "   -     \n")
#
# background systematics from file
#
inf = float('inf')
#
# systematics for inclusive, 0, >=1
#
sname = "systematics_htSig-" + str(options.ht) + "_metSig-" + str(options.met) + ".py"
execfile(sname)

sumerr = 0
bt = options.btag
if bt == "binc":
    bt = "inc"
elif bt != "b0":
    bt = "b1p"
for key in largestAbsDoubleRatioDeviation:
    err = sumLepErrors(key,bt,predLep)
    sumerr = sumerr + err*err
    print "Adding syst ",bt,key,err
#ofile.write("bkgSyst lnN     -   " + "%7.2f" % (1+math.sqrt(sumerr)) + "\n")

errClos = errMCClosure[options.ht][options.met]
print "Adding MC closure syst ",errClos
sumerr += errClos*errClos
#
# b-tag related
#
sname = "systematics_BT_htSig-" + str(options.ht) + "_metSig-" + str(options.met) + ".py"
execfile(sname)

bt = options.btag
if bt == "binc":  bt = "inc"
key = 'ScaleFrac'
err = sumLepErrors(key,bt,predLep)
sumerr = sumerr + err*err
print "Adding BT syst ",bt,key,err
ofile.write("bkgSyst lnN     -   " + "%7.3f" % (1+math.sqrt(sumerr)) + "\n")
print "total syst = ",math.sqrt(sumerr)
        
#
# b-tag systs
#
# keys for up and down variations / source
#
if options.btag != 'binc':
    btKeys = { 'beff' : [ 'btagEff3_Up_b_sf0', 'btagEff3_Down_b_sf0' ], \
               'leff' : [ 'btagEff3_Up_l_sf0', 'btagEff3_Down_l_sf0' ] }
    for vari in btKeys:
        lbsyst = vari.ljust(8) + "lnN".ljust(5)
        # keys
        keyUp = btKeys[vari][0]
        keyDown = btKeys[vari][1]
        # (signed) variation / lepton channel
        errLep = sumBTErrors(keyUp,keyDown,options.btag,predLep)
        print "Adding BT syst ",options.btag,vari,errLep
        ofile.write(vari.ljust(8) + "lnN     -   " + "%7.3f" % (1+abs(errLep)) + "\n")
    

ofile.close()

