import cPickle
import sys
import os
#
# convert the SMS efficiencies from files in /data/schoef/efficiencies/*_Efficiency.py
#   to "msugra-like" dictionaries
#

# input file
fnameEff = sys.argv[1]
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

