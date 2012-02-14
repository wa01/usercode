#!/usr/bin/env python
import cPickle
import re
import sys
import os

def which(program):
#    import os
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

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

def zeroSignalChannel (yields):
    noSig = False
    for btag in yields:
        if yields[btag] < 0.001: return True
    return False

from signalUtils import *

from optparse import OptionParser
parser = OptionParser()
parser.add_option("--btag", dest="btag", default="", type="string", action="store", help="btag bin (default = all)")
parser.add_option("--ht", dest="ht", type="int", action="store", help="HT cut")
parser.add_option("--met", dest="met", type="int", action="store", help="MET cut")
parser.add_option("--lumi", dest="lumi", default=4700., type="float", action="store", help="luminosity (pb-1)")
parser.add_option("--regroupM0", dest="regroupM0", default=1, type="int", action="store", help="step size in M0")
parser.add_option("--regroupM12", dest="regroupM12", default=1, type="int", action="store", help="step size in M12")
#parser.add_option("--m0", dest="m0", default=-1, type="int", action="store", help="")
#parser.add_option("--m12", dest="m12", default=-1, type="int", action="store", help="")
parser.add_option("--m0s", dest="m0s", default="999,0", type="string", action="store", help="range for M0")
parser.add_option("--m12s", dest="m12s", default="999,0", type="string", action="store", help="range for M12")
parser.add_option("-t", "--text", dest="text", default=False, action="store_true", help="keep output in text format")
parser.add_option("--njobs", dest="njobs", default=10, type="int", action="store", help="#points / job")
parser.add_option("-a", "--algo", dest="algo", default="HybridNew", type="string", action="store", help="limit algorithm (Asymptotic or HybridNew")
parser.add_option("--single", dest="single", default=False, action="store_true", help="single point evaluation (for HybridNew)")
parser.add_option("--exp", dest="exp", default=-1, type="float", action="store", help="quantile for expected limit")
parser.add_option("-M", "--model", dest="model", default="msugra", type="string", action="store", help="signal model")
parser.add_option("--nlo", dest="nlo", default=False, action="store_true", help="use NLO")
(options, args) = parser.parse_args()
#options.bin = True # fake that is a binary output, so that we parse shape lines

# default for b-tag bins (all)
btags = [ "b0", "b1", "b2" ]
if options.btag != "":
    btags = [ options.btag ]
    if options.btag == "binc":  btags = [ "inc"]

#if len(args) == 0:
#    sys.exit(1)
# template card file
if len(args) > 0:
    incards = args[0]
else:
    if options.btag != "":
        incards = options.btag
    else:
        incards = "multibtag"
    incards = incards+"-ht"+str(options.ht)+"-met"+str(options.met)+".txt"

# options for "combine"
combopt = "-M "+options.algo
if options.algo == "HybridNew":
    if not options.single:
        combopt = "-H Asymptotic " + combopt
    else:
        combopt = combopt + " --singlePoint 1.00"
    combopt = combopt + " --frequentist --testStat LHC"
    if options.exp > 0:
        combopt = combopt + " --expectedFromGrid " + str(options.exp)
elif options.algo == "Asymptotic":
    pass
else:
    print "Unknown algorithm",options.algo
    sys.exit(1)
#
# mass ranges
#
m0range = getRange(options.m0s)
m12range = getRange(options.m12s)
    
cwd = os.getcwd()
#
# output file name
#
if options.btag != "":
    basename = options.btag
else:
    basename = "multibtag"
basename = basename + "_ht" + str(options.ht) + "_met" + str(options.met)
dirname = "/tmp/adamwo/job_" + options.model
if options.nlo:  dirname = dirname + "NLO"
dirname = dirname + "_" + basename
if m0range[0] <= m0range[1]:
    dirname = dirname + "_m0_" + str(m0range[0])
    if m0range[0] != m0range[1]:
        dirname = dirname + "-" + str(m0range[1])
if m12range[0] <= m12range[1]:
    dirname = dirname + "_m12_" + str(m12range[0])
    if m12range[0] != m12range[1]:
        dirname = dirname + "-" + str(m12range[1])
if options.algo == "HybridNew":
    dirname = dirname + "_HN"
    if options.single:
        dirname = dirname + "SP"
    if options.exp > 0:
        dirname = dirname + "_exp" + "%3.3d" % int(100*options.exp)
elif options.algo == "Asymptotic":
    dirname = dirname + "_A"
os.mkdir(dirname)
os.system("cp "+incards+" "+dirname)
#
loString = "LO"
if options.nlo:  loString = "NLO"
#
# Mu efficiencies
#
fEffMuName = effFileName("Mu",options.model,loString)
fEffMu = open(fEffMuName,"rb")
effsMu = cPickle.load(fEffMu)
fEffMu.close()
#
# Ele efficiencies
#
fEffEleName = effFileName("Ele",options.model,loString)
fEffEle = open(fEffEleName,"rb")
effsEle = cPickle.load(fEffEle)
fEffEle.close()
#
# cross sections
#
if options.nlo:
    fXsecName = xsecFileName(options.model,"NLO")
    fXsecNLO = open(fXsecName)
    xsecsNLO = cPickle.load(fXsecNLO)
    fXsecNLO.close()
    xsecsLO = []
    for m12 in xsecsNLO:
        for m0 in xsecsNLO[m12]:
            xsecsLO.append( (  m0, m12, 10, 0, 1 ) )
else:
    fXsecName = xsecFileName(options.model,"LO")
    fXsecLO = open(fXsecName)
    xsecsLO = cPickle.load(fXsecLO)
    fXsecLO.close()

#
# list of all M0/M12 pairs
#
from getM0M12 import *
m0m12s = getM0M12a(effsMu,effsEle,xsecsLO,"binc",options.ht,options.met,options.regroupM0,options.regroupM12)

m0s = m0m12s.keys()
m0s.sort()

from createCards import createCards
from createCards import createMultiCards

#
# loop on mass combinations
#
idx = 0
massname = dirname + "/m0m12.lis"
fmass = open(massname,"w")
for m0 in m0s:

#    if options.m0 > 0 and m0 != options.m0:  continue
    # M0 within range?
    if m0range[1] >= m0range[0] and ( m0 < m0range[0] or m0 > m0range[1] ): continue

#    idx = idx + 1
    nmass = 0
    lmass = ""
    for m12 in m0m12s[m0]:

#        if options.m12 > 0 and m12 != options.m12:  continue
        # M12 within range?
        if m12range[1] >= m12range[0] and ( m12 < m12range[0] or m12 > m12range[1] ): continue

        # signal characterization (string)
        msugraString = signalString(options.model,m0,m12)
        # signal characterization (tuple)
        msugraTuple = signalTuple(options.model,m0,m12)
        if options.nlo:
            # yield for each b-tag bin in list
            sigYields = getSigYieldsNLO(btags,options.ht,options.met,msugraString,msugraTuple,options.lumi,xsecsNLO,effsMu,effsEle)
            if zeroSignalChannel(sigYields):
                # skip points with (at least one) channel without signal
                print "No signal for ",msugraString," ",sigYields
                continue
        else:
            # cross section
            xsLO = xsecsLO[msugraTuple]
            # yield for each b-tag bin in list
            sigYields = getSigYieldsLO(btags,options.ht,options.met,msugraString,msugraTuple,options.lumi,xsecsLO,effsMu,effsEle)
            if zeroSignalChannel(sigYields):
                # skip points with (at least one) channel without signal
                print "No signal for ",msugraString," ",sigYields
                continue
        print "sigYields = ",m0," ",m12," ",sigYields

        # output (text) file name
        modelname = dirname + "/model_" + str(m0) + "_" + str(m12)
        dcname = modelname + ".txt"
        dcfile = open(dcname,"w")
        save_stdout = sys.stdout
        sys.stdout = dcfile
        # creation of the final cards file (with signal yields)
        if len(btags) == 1:
            createCards(incards,btags[0],options.ht,options.met,sigYields[btags[0]])
        else:
            createMultiCards(incards,options.ht,options.met,sigYields)
        sys.stdout = save_stdout
        dcfile.close()

        # conversion to root workspace
        if not options.text:
            wsname = modelname + ".root"
            wscmd = "text2workspace.py " + dcname + " -m " + str(m0+m12/10000.) + " -o " + wsname
            os.system(wscmd)
            os.remove(dcname)

        # mass count and entry in mass list for CRAB jobs
        nmass = nmass + 1
        lmass = lmass + " " + str(m0) + " " + str(m12)
#        break

        if nmass >= options.njobs:
            fmass.write(lmass+"\n")
            nmass = 0
            lmass = ""
            idx = idx + 1

    if nmass != 0:
        fmass.write(lmass+"\n")
        nmass = 0
        lmass = ""
        idx = idx + 1

fmass.close()

#
# bash script for executing job
#
replacements = {}
replacements["${OPTIONS}"] = combopt

fshtmp = open("combine_template.sh","r")
shname = dirname + "/combine_"+options.algo+".sh"
fsh = open(shname,"w")
for line in fshtmp:
    for repl in replacements:
        line = line.replace(repl,replacements[repl])
    fsh.write(line)
fshtmp.close()
fsh.close()

#os.system("cp combine_"+options.algo+".sh "+dirname)

#
# CRAB cfg file
#
replacements = {}
replacements["${ALGORITHM}"] = options.algo
replacements["${NUMBER_OF_JOBS}"] = str(idx)

fcfgtmp = open("combine_template.cfg","r")
cfgname = dirname + "/combine_"+options.algo+".cfg"
fcrab = open(cfgname,"w")
for line in fcfgtmp:
    for repl in replacements:
        line = line.replace(repl,replacements[repl])
    fcrab.write(line)
fcfgtmp.close()
fcrab.close()
#
# link to executable
#
combine = which("combine")
if combine != "None":
    os.system("ln -s "+combine+" "+dirname+"/combine")


#if not ( options.m0 > 0 and options.m12 > 0 ):
#
# tar file of all models
#
if not ( m0range[0] == m0range[1] and m12range[0] == m12range[1] ):
    if options.text:
        os.system("cd "+dirname+"; tar -cf models.tar model_*.txt")
        os.system("rm "+dirname+"/model_*.txt")
    else:
        os.system("cd "+dirname+"; tar -cf models.tar model_*.root")
        os.system("rm "+dirname+"/model_*.root")
