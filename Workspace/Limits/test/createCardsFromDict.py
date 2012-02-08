#!/usr/bin/env python
import cPickle
import math

#
# writes Higgs combination program compatible cards files
#

import re
from sys import argv
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
ofile.write("kmax"+"%3d" % 4 +"\n")
#
# read observed and predicted background numbers
#
execfile("eventCounts.py")
obs = countsObs[options.ht][options.met][options.btag]
pred = countsPred[options.ht][options.met][options.btag]
errpred = errorsPred[options.ht][options.met][options.btag]

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

#
# background systematics from file
#
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
        print key,lep,err
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
for key in largestAbsDoubleRatioDeviation:
    for lep in [ "Mu", "Ele" ]:
        if key == 'ScaleFrac' or bt != 'inc':
            err = largestAbsDoubleRatioDeviation[key][lep][bt]
            err2 = largestAbsSingleRatioDeviation[key][lep][bt]
            if err2 < err:
                print "single < double ratio for ",key,lep,options.btag,err2,err
                err = err2
            print key,lep,err
            sumerr = sumerr + err*err
ofile.write("bkgSyst lnN     -   " + "%7.2f" % (1+math.sqrt(sumerr)) + "\n")

        


ofile.close()

