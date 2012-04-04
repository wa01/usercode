import cPickle
import sys
import os
#
# convert the SMS efficiencies from files in /data/schoef/efficiencies/*_Efficiency.py
#   to "msugra-like" dictionaries
#
from optparse import OptionParser
parser = OptionParser()
parser.add_option("--xsName", dest="xsname", default=None, type="string", action="store")
parser.add_option("--xsHisto", dest="xshisto", default="gluino", type="string", action="store")
parser.add_option("--lumi", dest="lumi", default=4700., type="float", action="store", help="luminosity (pb-1)")
(options, args) = parser.parse_args()
# input file
fnameEff = args[0]
execfile(fnameEff)
#
# loop over input dictionary
#
effDicts = {}
# btag, ht and met bins
for btag in efficiency:
    for ht in [ 750, 1000 ]:
        for met in [ 250, 350, 450, 550 ]:
            # intermediate key(s) (for mass ratios)
            for mass in efficiency[btag][ht][met]:
                # create missing levels in output dictionary
                if not mass in effDicts:  effDicts[mass] = {}
                if not btag in effDicts[mass]:  effDicts[mass][btag] = {}
                if not ht in effDicts[mass][btag]:  effDicts[mass][btag][ht] = {}
                if not met in effDicts[mass][btag][ht]:  effDicts[mass][btag][ht][met] = {}
                # copy efficiency values for each signal point
                for key in efficiency[btag][ht][met][mass]:
                    effDicts[mass][btag][ht][met][key] = efficiency[btag][ht][met][mass][key]
#
# create base name for output file (strip directories and extension)
#
fname0 = fnameEff
ind = fname0.rfind("/")
if ind != -1:  fname0 = fname0[ind+1:]
ind = fname0.find(".")
if ind != -1:  fname0 = fname0[0:ind]
ind = fname0.rfind("_Efficiencies")
if ind != -1:  fname0 = fname0[0:ind]
#
# loop over intermediate keys
#
for mass in effDicts:
    suffix = ""
    # if more than one intermediate key or a float: create suffix for model name
    if len(effDicts.keys()) > 1 or type(mass) == float or type(mass) == int:
        suffix = str(mass).replace(".","")
    fname = fname0 + suffix + "_NLO_efficiency.pkl"
    # don't overwrite output file
    if os.path.exists(fname):
        print "output file ",fname," exists"
        continue
    # dump the dictionary for the intermediate key
    fout = open(fname,"wb")
    cPickle.dump(effDicts[mass],fout)
    fout.close()

#
#
#
if options.xsname != None:
    import ROOT
    xsfile = ROOT.TFile(options.xsname)
    xsTH1 = xsfile.Get(options.xshisto)
    xsecs = {}
    evtDicts = {}
    nb = xsTH1.GetNbinsX()
    for i in range(1,nb+1):
        x = int(xsTH1.GetXaxis().GetBinLowEdge(i)+0.5)
        xsecs[x] = options.lumi*xsTH1.GetBinContent(i)

    import copy
    for mass in effDicts:
        evtDict = copy.deepcopy(effDicts[mass])
        for btag in evtDict:
            for ht in evtDict[btag]:
                for met in evtDict[btag][ht]:
                    for key in evtDict[btag][ht][met]:
                        parts = key.split("_")
                        m0 = int(parts[1])
                        evtDict[btag][ht][met][key] *= xsecs[m0]
        
        suffix = ""
        # if more than one intermediate key or a float: create suffix for model name
        if len(effDicts.keys()) > 1 or type(mass) == float or type(mass) == int:
            suffix = str(mass).replace(".","")
        fname = fname0 + suffix + "_NLO_events.pkl"
        # don't overwrite output file
        if os.path.exists(fname):
            print "output file ",fname," exists"
            continue
        # dump the dictionary for the intermediate key
        fout = open(fname,"wb")
        cPickle.dump(evtDict,fout)
        fout.close()
