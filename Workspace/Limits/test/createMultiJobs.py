#!/usr/bin/env python
import cPickle
import re
import sys
import os
import subprocess

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

def getSig (btags,ht,met,msugra,key,bdict):
    result = {}
    for btag in btags:
        if btag in bdict[ht][met][msugra]:
            result[btag] = bdict[ht][met][msugra][btag][key]
        else:
            result[btag] = None
    return result

def getBkg (btags,ht,met,key,bdict):
    result = {}
    for btag in btags:
        result[btag] = bdict[btag][ht][met][key]
    return result

def getBkgNumbers (btags,ht,met,bdict):
    result = {}
    for btag in btags:
        result[btag] = {}
        for key in bdict[btag][ht][met]:
            if key == 'obs' or key == 'pred' or key == 'countsNorm' or key == 'stats':
                result[btag][key] = bdict[btag][ht][met][key]
            elif key.startswith('syst'):
                if not 'syst' in result[btag]:  result[btag]['syst'] = {}
                k = key[4:]
                result[btag]['syst'][k] = bdict[btag][ht][met][key]
            else:
                print "unknown key in background numbers: ",key
                sys.exit(1)
    return result
                
def checkProcess (processes):
    result = {}
    for name in processes:
        p = processes[name]
        p.poll()
        if p.returncode == None:
            result[name] = p
        elif p.returncode != 0:
            print "job with pid ",p.pid," returned code ",p.returncode
        else:
            print "job with pid ",p.pid," terminated"
#            os.remove(name)
    return result

def getCorrStat (btags,bkgNumbers):
    result = { }
    # build list of sigmas / correlations
    covL = [ ]
    for i,bi in enumerate(btags):
        for j,bj in enumerate(btags):
            if j > i: continue
            covL.append(bkgNumbers[bi]['stats'][bj])
    # convert to TVector
    nb = len(btags)
    nb2 = len(covL)
    print nb,nb2
    covV = ROOT.TVectorD(nb2)
    for i,c in enumerate(covL):  covV[i] = c
    # get matrix
    matU = ROOT.Chol(covV)
    for i,bi in enumerate(btags):
        result[bi] = {}
        for j,bj in enumerate(btags):
            if j >= i:  result[bi][bj] = matU(i,j)
    return result

    
def waitProcess (processes,all=False):
    for name in processes:
        p = processes[name]
        p.wait()
        if not all:
            break
    return checkProcess(processes)

from signalUtils import *

from optparse import OptionParser
parser = OptionParser()
parser.add_option("--btag", dest="btag", default="", type="string", action="store", help="btag bin (default = all)")
parser.add_option("--ht", dest="ht", type="int", action="store", help="HT cut")
parser.add_option("--met", dest="met", type="int", action="store", help="MET cut")
#parser.add_option("--lumi", dest="lumi", default=4700., type="float", action="store", help="luminosity (pb-1)")
parser.add_option("--regroupM0", dest="regroupM0", default=1, type="int", action="store", help="step size in M0")
parser.add_option("--regroupM12", dest="regroupM12", default=1, type="int", action="store", help="step size in M12")
parser.add_option("--m0s", dest="m0s", default="999,0", type="string", action="store", help="range for M0")
parser.add_option("--m12s", dest="m12s", default="999,0", type="string", action="store", help="range for M12")
parser.add_option("-t", "--text", dest="text", default=False, action="store_true", help="keep output in text format")
parser.add_option("--njobs", dest="njobs", default=10, type="int", action="store", help="#points / job")
parser.add_option("-a", "--algo", dest="algo", default="HybridNew", type="string", action="store", help="limit algorithm (Asymptotic or HybridNew")
parser.add_option("--single", dest="single", default=False, action="store_true", help="single point evaluation (for HybridNew)")
parser.add_option("--exp", dest="exp", default=-1, type="float", action="store", help="quantile for expected limit")
parser.add_option("-M", "--model", dest="model", default="msugra", type="string", action="store", help="signal model")
parser.add_option("--lo", dest="lo", default=False, action="store_true", help="use LO")
parser.add_option("--nloVariation", dest="nloVar", default="0", type="choice", action="store", choices=["", "0", "-", "+"], help="NLO variation")
(options, args) = parser.parse_args()
#options.bin = True # fake that is a binary output, so that we parse shape lines

# default for b-tag bins (all)
btags = [ "b0", "b1", "b2" ]
if options.btag != "":
    btags = [ options.btag ]
    if options.btag == "binc":  btags = [ "inc"]
if len(btags) > 1:
    import ROOT
    from ROOT import gROOT
    gROOT.ProcessLine(".L Chol.C+")
#
# file with background information
#
bkgDict = cPickle.load(file("backgrounds.pkl"))
sigName = "sig_"+options.model
if not options.lo:
    sigName += "NLO"
if not options.lo and options.nloVar != "":
    sigName += options.nloVar
sigName += ".pkl"
sigDict = cPickle.load(file(sigName))


# options for "combine"
combopt = "-M "+options.algo
if options.algo == "HybridNew":
    if not options.single:
        combopt = "-H Asymptotic " + combopt
    else:
        combopt = combopt + " --singlePoint 1.00 --clsAcc 0 -T 500 -i 25 --saveToys --saveHybridResult -n Toys "
    combopt = combopt + " --frequentist --testStat LHC"
    if options.exp > 0:
        combopt = combopt + " --expectedFromGrid " + str(options.exp)
    combopt = combopt + " --fork 0"
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
if not options.lo:
    dirname = dirname + "NLO"
    if options.nloVar == '0':  dirname = dirname + "0"
    elif options.nloVar == '-':  dirname = dirname + "m"
    elif options.nloVar == '+':  dirname = dirname + "p"
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
#os.system("cp "+incards+" "+dirname)
fc = open(dirname+"/args.txt","w")
larg = ""
for a in sys.argv:  larg = larg + a + " "
fc.write(larg+"\n")
fc.close()
#
loString = "LO"
loVar = ""
if not options.lo:
    loString = "NLO"
    loVar = options.nloVar

#
# list of all M0/M12 pairs
#
from getM0M12 import *
m0m12s = getM0M12c(sigDict[options.ht][options.met].keys(),m0range,m12range,options.regroupM0,options.regroupM12)

m0s = m0m12s.keys()
m0s.sort()

#
# loop on mass combinations
#
idx = 0
massname = dirname + "/m0m12.lis"
fmass = open(massname,"w")
processes = {}
bkgCorrTags = [ ]
bkgCorrs = None
for m0 in m0s:

    nmass = 0
    lmass = ""
    for m12 in m0m12s[m0]:

        # signal characterization (string)
        msugraString = buildSignalString(options.model,m0,m12)
        # signal characterization (tuple)
        msugraTuple = buildSignalTuple(options.model,m0,m12)

        sigEvts = getSig(btags,options.ht,options.met,msugraString,'yield',sigDict)
        sigSyst = getSig(btags,options.ht,options.met,msugraString,'sigSystOther',sigDict)
        sigJES = getSig(btags,options.ht,options.met,msugraString,'sigSystJES',sigDict)
        if options.btag != 'binc':
            sigBeff = getSig(btags,options.ht,options.met,msugraString,'beff',sigDict)
            sigLeff = getSig(btags,options.ht,options.met,msugraString,'leff',sigDict)

        bkgNumbers = getBkgNumbers(btags,options.ht,options.met,bkgDict)
#        obsEvts = getBkg(btags,options.ht,options.met,'obs',bkgDict)
#        predEvts = getBkg(btags,options.ht,options.met,'pred',bkgDict)
#        bkgStats = getBkg(btags,options.ht,options.met,'stats',bkgDict)
#        bkgSyst = getBkg(btags,options.ht,options.met,'systOther',bkgDict)
#        if options.btag != 'binc':
#            bkgBeff = getBkg(btags,options.ht,options.met,'systbeff',bkgDict)
#            bkgLeff = getBkg(btags,options.ht,options.met,'systleff',bkgDict)
#        bkgJES = getBkg(btags,options.ht,options.met,'systJES',bkgDict)

        btagsFiltered = []
        for btag in btags:
            if sigEvts[btag] != None and sigEvts[btag] > 0.001:  btagsFiltered.append(btag)
        nc = len(btagsFiltered)
        if nc == 0:
            print "No signals for ",msugraString
            continue

        if len(btagsFiltered) > 1 and ( bkgCorrTags != btagsFiltered or bkgCorr == None ):
            bkgCorrTags = btagsFiltered
            bkgCorr = getCorrStat(btagsFiltered,bkgNumbers)

        # output (text) file name
        modelname = dirname + "/model_" + str(m0) + "_" + str(m12)
        dcname = modelname + ".txt"
        dcfile = open(dcname,"w")

        line = "imax".ljust(10) + str(nc)
        dcfile.write(line+"\n")
        dcfile.write("jmax".ljust(10)+"*\n")
        dcfile.write("kmax".ljust(10)+"*\n")

        dcfile.write("-"*(20+22*nc)+"\n")

        line = "bin".ljust(20)
        for btag in btagsFiltered:
            line += btag.rjust(22)
        dcfile.write(line+"\n")
        line = "observation".ljust(20)
        for btag in btagsFiltered:
            line += "%22d" % bkgNumbers[btag]['obs']
        dcfile.write(line+"\n")

        dcfile.write("-"*(20+22*nc)+"\n")
        
        line = "bin".ljust(20)
        for btag in btagsFiltered:
            line += btag.rjust(12) + btag.rjust(10)
        dcfile.write(line+"\n")

        line = "process".ljust(20)
        for btag in btagsFiltered:
            line += "sig".rjust(12) + "bkg".rjust(10)
        dcfile.write(line+"\n")

        line = "process".ljust(20)
        for btag in btagsFiltered:
            line += "0".rjust(12) + "1".rjust(10)
        dcfile.write(line+"\n")

        line = "rate".ljust(20)
        for btag in btagsFiltered:
            line += "%12.3f" % sigEvts[btag]
            line += "%10.3f" % bkgNumbers[btag]['pred']
        dcfile.write(line+"\n")
        
        dcfile.write("-"*(20+22*nc)+"\n")

        for btag in btagsFiltered:
            line = (btag+"StatF").ljust(15) + "lnN".ljust(5)
            if len(btagsFiltered) == 1:
                line += "-".rjust(12)
                line += "%10.3f" % (bkgNumbers[btag]['stats'][btag]+1)
            else:
                for btag1 in btagsFiltered:
                    line += "-".rjust(12)
                    if btag1 in bkgCorr[btag]:
                        line += "%10.3f" % (bkgCorr[btag][btag1]+1)
                    else:
                        line += "-".rjust(10)
            dcfile.write(line+"\n")

        for btag in btagsFiltered:
            line = (btag+"StatP").ljust(12) + "gmN".ljust(4)
            cnorm = bkgNumbers[btag]['countsNorm']
            line += "%4d" % cnorm
            for btag1 in btagsFiltered:
                line += "-".rjust(12)
                if btag1 == btag:
                    line += "%10.4f" % (bkgNumbers[btag]['pred']/cnorm)
                else:
                    line += "-".rjust(10)
            dcfile.write(line+"\n")

        for key in bkgNumbers[btagsFiltered[0]]['syst']:
            if key == 'beff' or key == 'leff' or key == 'JES': continue
            line = key.ljust(15) + "lnN".ljust(5)
            for btag in btagsFiltered:
                line += "-".rjust(12)
                line += "%10.3f" % (bkgNumbers[btag]['syst'][key]+1)
            dcfile.write(line+"\n")

        line = "JES".ljust(15) + "lnN".ljust(5)
        for btag in btagsFiltered:
            line += "%12.3f" % (sigJES[btag]+1)
            line += "%10.3f" % (bkgNumbers[btag]['syst']['JES']+1)
        dcfile.write(line+"\n")

        if options.btag != 'binc':

            line = "beff".ljust(15) + "lnN".ljust(5)
            for btag in btagsFiltered:
                line += "%12.3f" % (sigBeff[btag]+1)
                line += "%10.3f" % (bkgNumbers[btag]['syst']['beff']+1)          
            dcfile.write(line+"\n")

            line = "leff".ljust(15) + "lnN".ljust(5)
            for btag in btagsFiltered:
                line += "%12.3f" % (sigLeff[btag]+1)
                line += "%10.3f" % (bkgNumbers[btag]['syst']['leff']+1)
            dcfile.write(line+"\n")

        line = "sigSyst".ljust(15) + "lnN".ljust(5)
        for btag in btagsFiltered:
            line += "%12.3f" % (sigSyst[btag]+1)
            line += "-".rjust(10)
        dcfile.write(line+"\n")

        line = "lumi".ljust(15) + "lnN".ljust(5)
        for btag in btagsFiltered:
            line += "%12.3f" % 1.045
            line += "-".rjust(10)
        dcfile.write(line+"\n")

        dcfile.close()

        # asynchronous conversion to root workspace
        if not options.text:
            processes = checkProcess(processes)
            print "Currently ",len(processes.keys())," text2workspace commands running"
            if len(processes.keys()) > 5:
                print "Waiting for a text2workspace command to finish"
                processes = waitProcess(processes)
                print "Continuing"
            wsname = modelname + ".root"
            wscmd = [ "text2workspace.py", dcname , " -m ", str(m0+m12/10000.), " -o ", wsname ]
            processes[dcname] = subprocess.Popen(wscmd)
            
#        # conversion to root workspace
#        if not options.text:
#            wsname = modelname + ".root"
#            wscmd = "text2workspace.py " + dcname + " -m " + str(m0+m12/10000.) + " -o " + wsname
#            os.system(wscmd)
#            os.remove(dcname)

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
# waiting for all text2workspace commands to finish
if not options.text:
    waitProcess(processes,True)
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
    if line.find("return_data") != -1:
        return_data = 0
        fcrab.write("copy_data = 1\n")
        fcrab.write("publish_data = 0\n")
        fcrab.write("storage_element=srm-cms.cern.ch\n")
        castor = dirname
        ijob = castor.find("job_")
        if ijob != -1:  castor = castor[(ijob+4):] 
        fcrab.write("user_remote_dir=user/a/adamwo/CrabOutput/"+castor+"\n")
        fcrab.write("storage_path=/srm/managerv2?SFN=/castor/cern.ch/\n")

fcfgtmp.close()
fcrab.close()
#
# link to executable
#
combine = which("combine")
if combine != None:
    os.system("ln -s "+combine+" "+dirname+"/combine")

#if not ( options.m0 > 0 and options.m12 > 0 ):
#
# tar file of all models
#
if not ( m0range[0] == m0range[1] and m12range[0] == m12range[1] ):
    os.system("cd "+dirname+"; tar -zcf modelsText.tgz model_*.txt")
    os.system("rm "+dirname+"/model_*.txt")
    if not options.text:
        os.system("cd "+dirname+"; tar -cf models.tar model_*.root")
        os.system("rm "+dirname+"/model_*.root")
