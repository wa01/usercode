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
#    if bt == 'binc' or bt == 'b0' or bt == 'b1' or bt == 'b2':
#        result += countDict[ht][met][bt][lep]
#    elif bt == 'b1p':
#        for bt1 in [ 'b1', 'b2' ]:
#            result += countDict[ht][met][bt1][lep]
    if bt == 'binc' or bt == 'b0' or bt == 'b1' or bt == 'b2' or bt == 'b1p':
        result += countDict[ht][met][bt][lep]
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
#    if bt == 'binc' or bt == 'b0' or bt == 'b1' or bt == 'b2':
#        for lep in leptons:  result += sqr(errorDict[ht][met][bt][lep])
#    elif bt == 'b1p':
#        for bt1 in [ 'b1', 'b2' ]:
#            for lep in leptons:  result += sqr(errorDict[ht][met][bt1][lep])
    if bt == 'binc' or bt == 'b0' or bt == 'b1' or bt == 'b2' or bt == 'b1p':
        for lep in leptons:  result += sqr(errorDict[ht][met][bt][lep])
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
#parser.add_option("--ht", dest="ht", type="int", action="store", help="HT cut")
#parser.add_option("--met", dest="met", type="int", action="store", help="MET cut")
#parser.add_option("--btag", dest="btag", default="", type="string", action="store")
(options, args) = parser.parse_args()

btags = [ "binc", "b0", "b1", "b1p", "b2" ]
hts = [ 750, 1000 ]
mets = [ 250, 350, 450, 550 ]

systGroups = { 'WPol' : [ 'WPol1', 'WPol2+', 'WPol2-', 'WPol3' ], \
               'Eff' : [ 'HighEtaLepEff', 'LowPtLepEff', 'PU' ], \
               'Xsec' : [ 'DiLep', 'otherBkgs', 'scaleT', 'scaleW' ], \
               'Beta' : [ 'ScaleExpTT', 'ScaleExpWm', 'ScaleExpWp' ], \
               'Alpha' : [ 'alphaSlopeTT', 'alphaSlopeWm', 'alphaSlopeWp' ], \
               'Erf' : [ 'ErfcVar' ], 'JES' : [ 'JES' ] }
systGroupsInv = {}
for k1 in systGroups:
    for k2 in systGroups[k1]:
        assert not k2 in systGroupsInv
        systGroupsInv[k2] = k1

#
# b-tag consistency
#
btagConsistency = { 'b2':  {1000: {250: 0.085, 550: 0.180, 450: 0.154, 350: 0.114}, \
                             750: {250: 0.081, 550: 0.175, 450: 0.150, 350: 0.115}}, \
                    'b1':  {1000: {250: 0.020, 550: 0.015, 450: 0.020, 350: 0.021}, \
                             750: {250: 0.022, 550: 0.012, 450: 0.019, 350: 0.023}}, \
                    'b1p': {1000: {250: 0.042, 350: 0.048, 450: 0.053, 550: 0.050}, \
                             750: {250: 0.042, 350: 0.050, 450: 0.051, 550: 0.046} } }
#
# MC closure
#
errMCClosure = { 750 : { 250 : 0.073, 350 : 0.049, 450 : 0.071, 550 : 0.222 }, \
                 1000 : { 250 : 0.059, 350 : 0.077, 450  : 0.148, 550 : 0.197 } }
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

if len(btags) > 1:
    bkgCorr = cPickle.load(file("bkgCorrelations.pkl"))

for btag in btags:
    btagD = btag
    if btagD == 'binc':  btagD = 'inc'
    bkgDict[btagD] = {}
    for ht in hts:
        bkgDict[btagD][ht] = {}
        for met in mets:
            bkgDict[btagD][ht][met] = {}
            bkgDict[btagD][ht][met]['obs'] = getCounts(countsObs,btag)
            pred = getCounts(countsPred,btag)
            bkgDict[btagD][ht][met]['pred'] = pred
            predLep = {}
            for lep in leptons:  predLep[lep] = getCountsLep(countsPred,lep,btag)
            errpred = getErrors(errorsPred,btag)

            #            # background statistics (uncorrelated)
            #            bkgDict[btagD][ht][met]['stats'] = errpred/pred
            #
            # correlated background error from fit via Cholesky decomposition
            #   (redundant calculations because of inverse order btags <=> sig. regions
            #    but here CPU is not concern)
            #
            bkgDict[btagD][ht][met]['stats'] = { }
            for btag2 in btags:
                btag2D = btag2
                if btag2D == 'binc':  btag2D = 'inc'
                if btag2 == btag:
                    sig = errpred/pred
                else:
                    if btag == 'binc' or btag2 == 'binc': continue
                    name = btag+"-"+btag2
                    if name in bkgCorr[ht][met]:
                        sig = bkgCorr[ht][met][name]
                    else:
                        name = btag2+"-"+btag
                        sig = bkgCorr[ht][met][name]
                bkgDict[btagD][ht][met]['stats'][btag2D] = sig
            #
            # background statistics (Poisson)
            #
            bkgDict[btagD][ht][met]['countsNorm'] = getCounts(countsNorm,btag)
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

            sumerr = {}
            for key in largestAbsDoubleRatioDeviation:
                if key == 'JES':  continue
                assert key in systGroupsInv
                sg = systGroupsInv[key]
                if not sg in sumerr:  sumerr[sg] = 0.
                print key,bt,predLep
                err = sumAbsLepErrors(key,bt,predLep)
                sumerr[sg] += err*err
#                bkgDict[btagD][ht][met]['syst'+key] = err
#                sumerr += err*err
#                print key,btag,err
            jeserr = 0
            jeserrpm = []
            for key in [ 'pfRA4TupelizerJESMinus', 'pfRA4TupelizerJESPlus' ]:
                jeserrpm.append(sumLepErrors(key,bt,predLep))
            if abs(jeserrpm[0]) > abs(jeserrpm[1]):
                jeserr = jeserrpm[0]
            else:
                jeserr = jeserrpm[1]
            print "JES correction ",btag,jeserr
            bkgDict[btagD][ht][met]['systJES'] = jeserr
            #
            # add non-closure to other errors
            #
            errClos = errMCClosure[ht][met]
            print "Adding MC closure syst ",errClos
            assert not 'Closure' in sumerr
            sumerr['Closure'] = errClos*errClos
#            bkgDict[btagD][ht][met]['systClosure'] = err
#            sumerr += errClos*errClos
            # if available: add Met vs. b-tag uncertainty
            if btagD in btagConsistency:
                errMetBTag = btagConsistency[btagD][ht][met]
                sumerr['Closure'] += errMetBTag*errMetBTag
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
            assert 'Xsec' in sumerr
            sumerr['Xsec'] += err*err
#            bkgDict[btagD][ht][met]['systScaleFrac'] = err
#            sumerr += err*err
#            bkgDict[btagD][ht][met]['systOther'] = math.sqrt(sumerr)

            for key in sumerr:
                bkgDict[btagD][ht][met]['syst'+key] = math.sqrt(sumerr[key])
            #
            # b-tag systs
            #
            # keys for up and down variations / source
            #
            if met <= 350:
                btKeys = { 'beff' : [ 'btagEff4_Up_b_sf0', 'btagEff4_Down_b_sf0' ], \
                           'leff' : [ 'btagEff4_Up_l_sf0', 'btagEff4_Down_l_sf0' ] }
            else:
                btKeys = { 'beff' : [ 'btagEff3_Up_b_sf0', 'btagEff3_Down_b_sf0' ], \
                           'leff' : [ 'btagEff3_Up_l_sf0', 'btagEff3_Down_l_sf0' ] }
            for vari in btKeys:
                # keys
                keyUp = btKeys[vari][0]
                keyDown = btKeys[vari][1]
                # (signed) variation / lepton channel
                errLep = sumBTErrors(keyUp,keyDown,bt,predLep)
                bkgDict[btagD][ht][met]['syst'+vari] = errLep
    

print bkgDict
cPickle.dump(bkgDict,ofile)
ofile.close()

