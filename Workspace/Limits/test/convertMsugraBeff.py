import os
import sys
import cPickle
from string import Template
#
# read Tina's b efficiency variations for msugra on heplx* and assemble
# them in one pickle file
#

#
# get effect of scale factor (defined as SF1 / SF0)
#   if necessary : sum over b-tag bins
#
def getFactor (dic,signal,btags):
    norm_sf0 = 0.
    norm_sf1 = 0.
    for lep in dic:
        for btag in btags:
            norm_sf0 += dic[lep][signal]['Norm_sf0'][btag][0]
            norm_sf1 += dic[lep][signal]['Norm_sf1'][btag][0]
    if norm_sf0 < 0.000001:  return None
    return norm_sf1/norm_sf0
#
# get variation on SF1 (up-down)/(2*nominal)
#   if necessary : sum over b-tag bins
#
def getVar (dic,signal,btags,mode):
    norm_sf1 = 0.
    up = 0.
    down = 0.
    strup = 'Up_'+mode+'_sf1'
    strdown = 'Down_'+mode+'_sf1'
    s1 = 0.
    for lep in dic:
        for btag in btags:
            norm_sf1 += dic[lep][signal]['Norm_sf1'][btag][0]
            up += dic[lep][signal][strup][btag][0]
            down += dic[lep][signal][strdown][btag][0]
            if mode == 'b':
                s1 += dic[lep][signal]['Norm_sf1'][btag][0]* \
                      dic[lep][signal]['BDelta/(2Norm_sf1)'][btags[0]][0]
            else:
                s1 += dic[lep][signal]['Norm_sf1'][btag][0]* \
                      dic[lep][signal]['LDelta/(2Norm_sf1)'][btags[0]][0]
    var = (up-down)/2./norm_sf1
    if len(btags) == 1:
        if abs(var-s1/norm_sf1) > 0.0001:
            print "Mismatch in variations: ",btags[0],mode,var,s1/norm_sf1
            return None
    return var



#
# Template for directory / file name
#
template = Template('/data/trauner/Counts/Counts_MSUGRA_${lepton}_BGEff_newMCEff_Eff24Feb/${lepton}_MSUGRA_${ht}_ht_${met}_barepfmet_allJets_withEffgt500_eff24Feb_absoluteErr_BGEff.py')
#
# definition of flavours, HT and MET regions
#
leptons = [ 'Muon', 'Electron' ]
hts = [ 750, 1000 ]
mets = [ 250, 350, 450, 550 ]
#
# translate btag bin labels
#
btagLabels = { 'b0' : [ '0' ], 'b1' : [ '1' ], 'b1p' : [ '1' , '>=2' ],  'b2' : [ '>=2' ] }
#
# create dictionary and loop over flavours, HT and MET cuts
#
effDict = {}
for ht in hts:
#    if not ht in effDict:  effDict[ht] = {}
    for met in mets:
#        if not met in effDict[ht]:  effDict[ht][met] = {}
        inDict = { 'Muon' : {}, 'Electron' : {} }
        for lepton in inDict:
            ifname = template.substitute(lepton=lepton,ht=str(ht),met=str(met))
            if not os.path.exists(ifname):
                print "No such file ",ifname
            print ifname
            # execute input file
            execfile(ifname)
            labels = nbtags.keys()
            assert(len(labels)==1)
            inDict[lepton] = nbtags[labels[0]]
        # unused label
        for label in nbtags:
            # msugra points
            for signal in nbtags[label]:
                # translate to standard string
                msugraString = signal.replace("signal_","msugra_")
                msugraString += "_10_0_1"
#                if not msugraString in effDict[ht][met]: effDict[ht][met][msugraString] = {}
#                    btags = nbtags[label][signal]['originalMC'].keys()
#                    print lep,ht,met,msugraString
                # btag bins (output / input notation)
                for btagOut, btagsIn in btagLabels.iteritems():
                    # scaling after application of SF
                    sfFactor = getFactor(inDict,signal,btagsIn)
                    # skip empty bins
                    if sfFactor == None:  continue
                    if not btagOut in effDict:  effDict[btagOut] = {}
                    if not ht in effDict[btagOut]: effDict[btagOut][ht] = {}
                    if not met in effDict[btagOut][ht]: effDict[btagOut][ht][met] = {}
                    if not msugraString in effDict[btagOut][ht][met]:
                        effDict[btagOut][ht][met][msugraString] = {}
                    # add correction factor and variations to dictionary
#                    if not btagOut in effDict[ht][met][msugraString]:
#                        effDict[btagOut][ht][met][msugraString] = {}
                    effDict[btagOut][ht][met][msugraString]['sfFactor'] = sfFactor
                    effDict[btagOut][ht][met][msugraString]['relVarB'] = getVar(inDict,signal,btagsIn,'b')
                    effDict[btagOut][ht][met][msugraString]['relVarL'] = getVar(inDict,signal,btagsIn,'l')
#                print effDict
#                sys.exit(0)
#
# write dictionary
#
fout = open("msugraBeffSyst.pkl","wb")
cPickle.dump(effDict,fout)
fout.close()
