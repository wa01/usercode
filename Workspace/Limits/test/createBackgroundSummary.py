#!/usr/bin/env python
import cPickle
import math
import sys
#import re


leptons = [ 'Mu', 'Ele' ]
def sqr (val):
    return val*val
def getCountsLep (countDict,lep,btag):
    result = 0
    bt = btag    
    if bt == 'binc' or bt == 'b0' or bt == 'b1' or bt == 'b2':
        result += countDict[ht][met][bt][lep]
    elif bt == 'b1p':
        for bt1 in [ 'b1', 'b2' ]:
            result += countDict[ht][met][bt1][lep]
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
        for lep in leptons:  result += sqr(errorDict[ht][met][bt][lep])
    elif bt == 'b1p':
        for bt1 in [ 'b1', 'b2' ]:
            for lep in leptons:  result += sqr(errorDict[ht][met][bt1][lep])
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

btags = [ "binc", "b0", "b1", "b1p", "b2" ]
hts = [ 750, 1000 ]
mets = [ 250, 350, 450, 550 ]

#
# MC closure
#
errMCClosure = { 750 : { 250 : 0.072, 350 : 0.041, 450 : 0.073, 550 : 0.224 }, \
                 1000 : { 250 : 0.056, 350 : 0.074, 450  : 0.098, 550 : 0.198 } }
#
# output file name
#
bkgDict = { }
oname = "backgrounds.pkl"
if os.path.exists(oname):
    print "output file ",oname," exists"
    sys.exit(1)
ofile = open(oname,"wb")
    
#
# read observed and predicted background numbers
#
execfile("eventCounts.py")


for btag in btags:
    bkgDict[btag] = {}
    for ht in hts:
        bkgDict[btag][ht] = {}
        for met in mets:
            bkgDict[btag][ht][met] = {}
            bkgDict[btag][ht][met]['obs'] = getCounts(countsObs,btag)
            pred = getCounts(countsPred,btag)
            bkgDict[btag][ht][met]['pred'] = pred
            predLep = {}
            for lep in leptons:  predLep[lep] = getCountsLep(countsPred,lep,btag)
            errpred = getErrors(errorsPred,btag)

            # background statistics (uncorrelated)
            bkgDict[btag][ht][met]['stats'] = errpred/pred
            #
            # systematics (non-b related)
            #
            inf = float('inf')
            sname = "Systematics/systematics_htSig-" + str(ht) + "_metSig-" + str(met) + ".py"
            execfile(sname)

            bt = btag
            if bt == 'binc':
                bt = 'inc'
            elif btag == 'b1' or btag == 'b2':
                bt = 'b1p'

            sumerr = 0
            for key in largestAbsDoubleRatioDeviation:
                if key == 'JES':  continue
                print key,bt,predLep
                err = sumAbsLepErrors(key,bt,predLep)
                sumerr += err*err
                print key,btag,err
            jeserr = 0
            jeserrpm = []
            for key in [ 'pfRA4TupelizerJESMinus', 'pfRA4TupelizerJESPlus' ]:
                jeserrpm.append(sumLepErrors(key,bt,predLep))
            if abs(jeserrpm[0]) > abs(jeserrpm[1]):
                jeserr = jeserrpm[0]
            else:
                jeserr = jeserrpm[1]
            print "JES correction ",btag,jeserr
            bkgDict[btag][ht][met]['systJES'] = jeserr
            #
            # add non-closure to other errors
            #
            errClos = errMCClosure[ht][met]
            print "Adding MC closure syst ",errClos
            sumerr += errClos*errClos
            #
            # systematics (b-tag related)
            #
            sname = "Systematics/systematics_BT_htSig-" + str(ht) + "_metSig-" + str(met) + ".py"
            execfile(sname)

            bt = btag
            if bt == 'binc':  bt = 'inc'
            #
            # W/tt ratio
            #
            key = 'ScaleFrac'
            err = sumAbsLepErrors(key,bt,predLep)
            sumerr += err*err
            bkgDict[btag][ht][met]['systOther'] = math.sqrt(sumerr)
            #
            # b-tag systs
            #
            # keys for up and down variations / source
            #
            btKeys = { 'beff' : [ 'btagEff3_Up_b_sf0', 'btagEff3_Down_b_sf0' ], \
                       'leff' : [ 'btagEff3_Up_l_sf0', 'btagEff3_Down_l_sf0' ] }
            for vari in btKeys:
                # keys
                keyUp = btKeys[vari][0]
                keyDown = btKeys[vari][1]
                # (signed) variation / lepton channel
                errLep = sumBTErrors(keyUp,keyDown,bt,predLep)
                bkgDict[btag][ht][met]['syst'+vari] = errLep
    

print bkgDict
cPickle.dump(bkgDict,ofile)
ofile.close()

