#!/usr/bin/env python
import cPickle
import re
import sys
import os
import math
#import subprocess

def sqrSum (conts):
    sum = 0.
    for cont in conts:  sum += cont*cont
    return math.sqrt(sum)

def getRange (rangeString):
    result = ( 0, -1 )
    parts = rangeString.split(",")
    if len(parts) > 2:
        print "Ignoring invalid range specification: ",rangeString
        return result
    for part in parts:
        if not part.isdigit():
            print "Ignoring invalid range specification: ",rangeString
            return result
    if len(parts) == 1:  parts.append(parts[0])
    result = ( int(parts[0]), int(parts[1]) )
    return result

def addSignalTuplesFromEff (effDict,model,btags,hts,mets,sigTuples):
    for btag in btags:
        if not btag in effDict:  continue
        for ht in hts:
            if not ht in effDict[btag]: continue
            for met in mets:
                if not met in effDict[btag][ht]: continue
                for msugra in effDict[btag][ht][met]:
                    m0 = getFromSignalString(msugra,"m0")
                    m12 = getFromSignalString(msugra,"m12")
                    tup = buildSignalTuple(model,m0,m12)
                    if not tup in sigTuples:  sigTuples.append(tup)
    return

from signalUtils import *

from optparse import OptionParser
parser = OptionParser()
#parser.add_option("--ht", dest="ht", type="int", action="store", help="HT cut")
#parser.add_option("--met", dest="met", type="int", action="store", help="MET cut")
parser.add_option("--m0s", dest="m0s", default="999,0", type="string", action="store", help="range for M0")
parser.add_option("--m12s", dest="m12s", default="999,0", type="string", action="store", help="range for M12")
parser.add_option("-M", "--model", dest="model", default="msugra", type="string", action="store", help="signal model")
parser.add_option("--nlo", dest="nlo", default=True, action="store_true", help="use NLO")
parser.add_option("--nloVariation", dest="nloVar", default="0", type="choice", action="store", choices=["", "0", "-", "+"], help="NLO variation")
(options, args) = parser.parse_args()
#options.bin = True # fake that is a binary output, so that we parse shape lines

btags = [ 'inc', 'b0', 'b1', 'b1p', 'b2' ]
hts = [ 750, 1000 ]
mets = [ 250, 350, 450, 550 ]
#hts = [ 1000 ]
#mets = [ 250 ]
#
# mass ranges
#
m0range = getRange(options.m0s)
m12range = getRange(options.m12s)
    
cwd = os.getcwd()
#
# output file name
#
#basename = "ht" + str(ht) + "_met" + str(met)
oname = "sig_" + options.model
if options.nlo:
    oname = oname + "NLO"
    if options.nloVar == '0':  oname = oname + "0"
    elif options.nloVar == '-':  oname = oname + "m"
    elif options.nloVar == '+':  oname = oname + "p"
#oname = oname + "_" + basename
if m0range[0] <= m0range[1]:
    oname = oname + "_m0_" + str(m0range[0])
    if m0range[0] != m0range[1]:
        oname = oname + "-" + str(m0range[1])
if m12range[0] <= m12range[1]:
    oname = oname + "_m12_" + str(m12range[0])
    if m12range[0] != m12range[1]:
        oname = oname + "-" + str(m0range[1])
oname = oname + ".pkl"
if os.path.exists(oname):
    print "output file ",oname," exists"
    sys.exit(1)
fout = open(oname,"wb")
#
loString = "LO"
loVar = ""
if options.nlo:
    loString = "NLO"
    loVar = options.nloVar
#
# Mu efficiencies
#
fEffMuName = effFileName("Mu",options.model,loString,loVar,True)
fEffMu = open(fEffMuName,"rb")
effsMu = cPickle.load(fEffMu)
fEffMu.close()
#fEffRatioMuName = ratioFileName("Mu",options.model,loString,loVar)
#fEffRatioMu = open(fEffRatioMuName,"rb")
#effRatiosMu = cPickle.load(fEffRatioMu)
#fEffRatioMu.close()
#
# Ele efficiencies
#
fEffEleName = effFileName("Ele",options.model,loString,loVar,True)
fEffEle = open(fEffEleName,"rb")
effsEle = cPickle.load(fEffEle)
fEffEle.close()
#fEffRatioEleName = ratioFileName("Ele",options.model,loString,loVar)
#fEffRatioEle = open(fEffRatioEleName,"rb")
#effRatiosEle = cPickle.load(fEffRatioEle)
#fEffRatioEle.close()
#
# JES
#
fJesSystsName = options.model+"JES-smoothed.pkl"
jesSysts = cPickle.load(file(fJesSystsName))
#
# B and udsg scale factors
#
if options.model == 'msugra':
    beffSysts = cPickle.load(file("msugraBeffSyst-smoothed.pkl"))
#
# cross sections
#
if options.nlo:
    fXsecName = xsecFileName(options.model,"NLO")
    fXsecNLO = open(fXsecName)
    xsecsNLO = cPickle.load(fXsecNLO)
    fXsecNLO.close()
    xsecsLO = []
#    for m12 in xsecsNLO:
#        for m0 in xsecsNLO[m12]:
#            xsecsLO.append( (  m0, m12, 10, 0, 1 ) )
    addSignalTuplesFromEff(effsMu,options.model,btags,hts,mets,xsecsLO)
    addSignalTuplesFromEff(effsEle,options.model,btags,hts,mets,xsecsLO)
else:
    fXsecName = xsecFileName(options.model,"LO")
    fXsecLO = open(fXsecName)
    xsecsLO = cPickle.load(fXsecLO)
    fXsecLO.close()

#
# list of all M0/M12 pairs
#
from getM0M12 import *
#m0m12s = getM0M12a(effsMu,effsEle,xsecsLO,"binc",hts[0],mets[0],1,1)
sigTuples = getM0M12b(xsecsLO,m0range,m12range,1,1)

sigDict = {}
for ht in hts:
    if not ht in sigDict:  sigDict[ht] = {}
    for met in mets:
        if not met in sigDict[ht]:  sigDict[ht][met] = {}
        #
        # loop on mass combinations
        #
        idx = 0
        for msugraTuple in sigTuples:
            m0 = msugraTuple[0]
            m12 = msugraTuple[1]

            nmass = 0
            lmass = ""
            

            # signal characterization (string / tuple)
            msugraString = buildSignalString(options.model,m0,m12)
            if options.nlo:
                # yield for each b-tag bin in list
                sigYields = getSigYieldsNLO(btags,ht,met,msugraString,effsMu,effsEle)
                #sigYields = getSmoothedSigYieldsNLO(btags,ht,met,msugraString,msugraTuple,1.,\
                #                                    effRatiosMu,effsMu,effRatiosEle,effsEle)
            else:
                # cross section
                xsLO = xsecsLO[msugraTuple]
                # yield for each b-tag bin in list
                sigYields = getSigYieldsLO(btags,ht,met,msugraString,msugraTuple,1.,xsecsLO,effsMu,effsEle)
            if msugraString in jesSysts[ht][met]:
                systJES = jesSysts[ht][met][msugraString]
            else:
                continue
            for btag in btags:
                if sigYields[btag] == None or sigYields[btag] < 0.001:  continue
                if btag != 'inc':
                    corrFactor = beffSysts[btag][ht][met][msugraString]['sfFactor']
                    sigYields[btag] *= corrFactor
                if btag != 'inc':
                    try:
                        systBeff = beffSysts[btag][ht][met][msugraString]['relVarB']
                        systLeff = beffSysts[btag][ht][met][msugraString]['relVarL']
                    except:
                        continue
                if not msugraString in sigDict[ht][met]:  sigDict[ht][met][msugraString] = {}
                if not btag in sigDict[ht][met][msugraString]:  sigDict[ht][met][msugraString][btag] = {}
                sigDict[ht][met][msugraString][btag] = {}
                sigDict[ht][met][msugraString][btag]['yield'] = sigYields[btag]
                sigDict[ht][met][msugraString][btag]['sigSystOther'] = 0.05
                sigDict[ht][met][msugraString][btag]['sigSystJES'] = systJES
#                sigDict[ht][met][msugraString][btag]['sigSyst'] = 0.20
                if btag != 'inc':
                    sigDict[ht][met][msugraString][btag]['beff'] = systBeff
                    sigDict[ht][met][msugraString][btag]['leff'] = systLeff

                
cPickle.dump(sigDict,fout)
fout.close()
