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
def getFactor (dic,btags):
    norm_sf0 = 0.
    norm_sf1 = 0.
    for btag in btags:
        norm_sf0 += dic['Norm_sf0'][btag][0]
        norm_sf1 += dic['Norm_sf1'][btag][0]
    if norm_sf0 < 0.000001:  return None
    return norm_sf1/norm_sf0
#
# get variation on SF1 (up-down)/(2*nominal)
#   if necessary : sum over b-tag bins
#
def getVar (dic,btags,mode):
    norm_sf1 = 0.
    up = 0.
    down = 0.
    strup = 'Up_'+mode+'_sf1'
    strdown = 'Down_'+mode+'_sf1'
    for btag in btags:
        norm_sf1 += dic['Norm_sf1'][btag][0]
        up += dic[strup][btag][0]
        down += dic[strdown][btag][0]
    var = (up-down)/2./norm_sf1
    if len(btags) == 1:
        if ( mode == 'b' and abs(var-dic['BDelta/(2Norm_sf1)'][btags[0]][0]) > 0.0001 ) or \
           ( mode == 'l' and abs(var-dic['LDelta/(2Norm_sf1)'][btags[0]][0]) > 0.0001 ):
            print "Mismatch in variations: ",btags[0],mode,var,dic['BDelta/(2Norm_sf1)'][btags[0]][0],dic['LDelta/(2Norm_sf1)'][btags[0]][0]
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
for lepton in leptons:
    if lepton == 'Muon':
        lep = 'Mu'
    else:
        lep = 'Ele'
    if not lep in effDict:  effDict[lep] = {}
    for ht in hts:
        if not ht in effDict[lep]:  effDict[lep][ht] = {}
        for met in mets:
            if not met in effDict[lep][ht]:  effDict[lep][ht][met] = {}
            ifname = template.substitute(lepton=lepton,ht=str(ht),met=str(met))
            if not os.path.exists(ifname):
                print "No such file ",ifname
            print ifname
            # execute input file
            execfile(ifname)
            # unused label
            for label in nbtags:
                # msugra points
                for signal in nbtags[label]:
                    # translate to standard string
                    msugraString = signal.replace("signal_","msugra_")
                    msugraString += "_10_0_1"
                    if not msugraString in effDict[lep][ht][met]: effDict[lep][ht][met][msugraString] = {}
#                    btags = nbtags[label][signal]['originalMC'].keys()
#                    print lep,ht,met,msugraString
                    # btag bins (output / input notation)
                    for btagOut, btagsIn in btagLabels.iteritems():
                        # scaling after application of SF
                        sfFactor = getFactor(nbtags[label][signal],btagsIn)
                        # skip empty bins
                        if sfFactor == None:  continue
                        # add correction factor and variations to dictionary
                        if not btagOut in effDict[lep][ht][met][msugraString]:
                            effDict[lep][ht][met][msugraString][btagOut] = {}
                        effDict[lep][ht][met][msugraString][btagOut]['sfFactor'] = sfFactor
                        effDict[lep][ht][met][msugraString][btagOut]['relVarB'] = getVar(nbtags[label][signal],btagsIn,'b')
                        effDict[lep][ht][met][msugraString][btagOut]['relVarL'] = getVar(nbtags[label][signal],btagsIn,'l')
#                print effDict
#                sys.exit(0)
#
# write dictionary
#
fout = open("msugraBeffSyst.pkl","wb")
cPickle.dump(effDict,fout)
fout.close()
