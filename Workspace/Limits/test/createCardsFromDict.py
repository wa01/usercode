#!/usr/bin/env python
import cPickle
import math

#
# writes Higgs combination program compatible cards files
#

import re
#from sys import argv
import sys
import os.path
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--ht", dest="ht", type="int", action="store", help="HT cut")
parser.add_option("--met", dest="met", type="int", action="store", help="MET cut")
parser.add_option("--btag", dest="btag", default="binc", type="string", action="store", help="btag bin name")
(options, args) = parser.parse_args()

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
if options.btag == 'binc' or options.btag == 'b0' or options.btag == 'b1' or options.btag == 'b2':
    obs = countsObs[options.ht][options.met][options.btag]
    pred = countsPred[options.ht][options.met][options.btag]
    errpred = errorsPred[options.ht][options.met][options.btag]
elif options.btag == 'b1p':
    obs1 = countsObs[options.ht][options.met]['b1']
    pred1 = countsPred[options.ht][options.met]['b1']
    errpred1 = errorsPred[options.ht][options.met]['b1']
    obs2 = countsObs[options.ht][options.met]['b2']
    pred2 = countsPred[options.ht][options.met]['b2']
    errpred2 = errorsPred[options.ht][options.met]['b2']
    obs = obs1 + obs2
    pred = pred1 + pred2
    errpred = math.sqrt(errpred1*errpred1+errpred2*errpred2)
else:
    print "Unknown b-tag bin ",options.btag
    sys.exit(1)
    
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
ofile.write("rate     " + "%9.2f" % 1. + "%9.2f" % pred + "\n")
#
# uncertainties: luminosity, other signal systematics, background statistics
#
ofile.write("lumi    lnN    1.045   -     \n")
ofile.write("sigSyst lnN    1.20    -   \n")
ofile.write("bkgStat lnN     -   " + "%7.2f" % (1+errpred/pred) + "\n")
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
    for lep in [ "Mu", "Ele" ]:
        err = largestAbsDoubleRatioDeviation[key][lep][bt]
        err2 = largestAbsSingleRatioDeviation[key][lep][bt]
        if err2 < err:
            print "single < double ratio for ",key,lep,options.btag,err2,err
            err = err2
#        print key,lep,err
        sumerr = sumerr + err*err
#ofile.write("bkgSyst lnN     -   " + "%7.2f" % (1+math.sqrt(sumerr)) + "\n")

#
# b-tag related
#
sname = "systematics_BT_htSig-" + str(options.ht) + "_metSig-" + str(options.met) + ".py"
execfile(sname)

#sumerr = 0
bt = options.btag
if bt == "binc":  bt = "inc"
key = 'ScaleFrac'
for lep in [ "Mu", "Ele" ]:
    err = largestAbsDoubleRatioDeviation[key][lep][bt]
    err2 = largestAbsSingleRatioDeviation[key][lep][bt]
    if err2 < err:
        print "single < double ratio for ",key,lep,options.btag,err2,err
        err = err2
#      print key,lep,err
    sumerr = sumerr + err*err
ofile.write("bkgSyst lnN     -   " + "%7.2f" % (1+math.sqrt(sumerr)) + "\n")

        
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
        errLep = 0.
        for lep in [ "Mu", "Ele" ]:
            # get variations w.r.t. 1. (revert sign for down)
            dUp = singleRatio[keyUp][lep][options.btag] - 1.
            dDown = -(singleRatio[keyDown][lep][options.btag]-1.)
            # take the average of up/down and the maximum of Mu/Ele
            errAve = (dUp+dDown)/2.
            if abs(errAve) > abs(errLep):  errLep = errAve
        ofile.write(vari.ljust(8) + "lnN     -   " + "%7.2f" % (1+abs(errLep)) + "\n")
    

ofile.close()

