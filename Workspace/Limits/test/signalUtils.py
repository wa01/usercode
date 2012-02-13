from string import Template
#
# small utility functions for signal decoding
#
channelNames = [ "Ele", "Mu" ]
modelTemplates = { "msugra": [ Template("msugra_${m0}_${m12}_10_0_1"), ( "{m0}", "{m12}", 10, 0, 1 ) ],
                   "T1tttt": [ Template("T1tttt_${m0}_${m12}_-1"), ( "{m0}", "{m12}", -1 ) ] }
def effFileName (channel,model,order="LO"):
    assert channel in channelNames
    assert model in modelTemplates
    if model == "msugra":
#        return channel+"_"+model+"Efficiencies.pkl"
        if order == "LO":
            return channel+"_"+model+"_LO_efficiency.pkl"
        else:
            return channel+"_"+model+"_NLO_eventsPP.pkl"
    else:
        return channel+"_"+model+"_Efficiencies.pkl"
    
def xsecFileName (model,order="LO"):
    assert model in modelTemplates
    if model == "msugra":
        if order == "LO":
            return "goodModelNames_10_0_1.pkl"
        else:
            return "tanb10.msugra_xsecs.pc"
    else:
        return "xsec"+model+".pkl"
    
    
def signalString (model,m0_,m12_):
    assert model in modelTemplates
    return modelTemplates[model][0].substitute(m0=str(m0_),m12=str(m12_))

def signalTuple (model,m0_,m12_):
    assert model in modelTemplates
    return (m0_,m12_)+modelTemplates[model][1][2:]

def getSigYieldsLO (btags,ht,met,msugraString,msugraTuple,lumi,xsecs,effsMu,effsEle):
    sigYields = {}
    xsLO = xsecs[msugraTuple]
    for btag in btags:
        if btag != 'b1p':
            effMu = effsMu[btag][ht][met][msugraString]
            effEle = effsEle[btag][ht][met][msugraString]
        else:
            effMu1 = effsMu['b1'][ht][met][msugraString]
            effEle1 = effsEle['b1'][ht][met][msugraString]
            effMu2 = effsMu['b2'][ht][met][msugraString]
            effEle2 = effsEle['b2'][ht][met][msugraString]
            effMu = effMu1 + effMu2
            effEle = effEle1 + effEle2
        sigYields[btag] = lumi*(effEle+effMu)*xsLO
    return sigYields

def getSigYieldsNLO (btags,ht,met,msugraString,msugraTuple,lumi,xsecs,effsMu,effsEle):
    sigYields = {}
#    meanEffs = {}
    for btag in btags:
        sigYields[btag] = 0.
#        meanEffs[btag] = [ 0, 0, 0 ]
    for pp in xsecs[msugraTuple[1]][msugraTuple[0]]:
        if pp == 'total':  continue
        xsNLO = xsecs[msugraTuple[1]][msugraTuple[0]][pp]
#        print "xs for ",pp," = ",xsNLO
        if not ( pp in effsMu[btag][ht][met][msugraString] and \
                 pp in effsEle[btag][ht][met][msugraString] ):
            print "Missing process in efficiency files: ",pp
            continue
        for btag in btags:
            effMu = effsMu[btag][ht][met][msugraString][pp]
            effEle = effsEle[btag][ht][met][msugraString][pp]
#            meanEffs[btag][0] += xsNLO
#            meanEffs[btag][1] += xsNLO*effMu
#            meanEffs[btag][2] += xsNLO*effEle
            # temporary fix: using yields instead of efficiencies
            sigYields[btag] += effEle+effMu
#            print "contribution for ",btag,' is ',(effEle+effMu)
#            print "contribution for ",btag,' is ',lumi*(effEle+effMu)*xsNLO
#            sigYields[btag] += lumi*(effEle+effMu)*xsNLO
#    for btag in btags:
#        if meanEffs[btag][0] > 0:
#            meanEffs[btag][1] /= meanEffs[btag][0]
#            meanEffs[btag][2] /= meanEffs[btag][0]
#    print meanEffs
#    print sigYields
    return sigYields
