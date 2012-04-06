from string import Template
#
# small utility functions for signal decoding
#
channelNames = [ "Ele", "Mu" ]
modelTemplates = { "msugra": [ Template("msugra_${m0}_${m12}_10_0_1"), ( "{m0}", "{m12}", 10, 0, 1 ) ],
                   "T1tttt": [ Template("T1tttt_${m0}_${m12}_-1"), ( "{m0}", "{m12}", -1 ) ],
                   "T3w": [ Template("T3w_${m0}_${m12}_${m3}"), ( "{m0}", "{m12}", "{m3}" ) ],
                   "*": [ Template("${model}_${m0}_${m12}_${suffix}"), ( "{model}", "{m0}", "{m12}", "{suffix}" ) ] }
def effFileName (channel,model,order="LO",var="",smoothed=False,m3Ratio=-1):
    assert channel in channelNames
    assert model in modelTemplates
    name = channel + "_"
    if model == "msugra":
#        return channel+"_"+model+"Efficiencies"
        if order == "LO":
            name = name + model + "_LO_efficiency"
        else:
            if var != "":  name += var + "_"
            print name+model+"_NLO_events"
            name = name + model + "_NLO_events"
    else:
        name += model
        if m3Ratio > 0:
            name += str(m3Ratio).replace(".","")
#        name = name + model + "_Efficiencies"
        if order == "LO":
            name += "_LO_efficiency"
        else:
#            if var != "":  name += var + "_"
            name += "_NLO_events"
    if smoothed:  name += "-smoothed"
    return name+".pkl"
    
def ratioFileName (channel,model,order="LO",var=""):
    assert channel in channelNames
    assert model in modelTemplates
    name = channel + "_"
    if model == "msugra":
#        return channel+"_"+model+"Efficiencies.pkl"
        if order == "LO":
            return name+model+"_LO_efficiency.pkl"
        else:
            if var != "":  name += var + "_"
            print name+model+"_NLO_efficiency-EffRatio.pkl"
            return name+model+"_NLO_efficiency-EffRatio.pkl"
    else:
        return name+model+"_Efficiencies.pkl"
    
def xsecFileName (model,order="LO"):
    assert model in modelTemplates
    if model == "msugra":
        if order == "LO":
            return "goodModelNames_10_0_1.pkl"
        else:
            return "tanb10.msugra_xsecs.pc"
    else:
        return "xsec"+model+".pkl"
    
    
def buildSignalString (model_,m0_,m12_,m3Ratio=-1,suffices=None):
    assert model_ in modelTemplates
    if suffices != None:
        return modelTemplates['*'][0].substitute(model=model_,m0=str(m0_),m12=str(m12_),suffix="_".join(suffices))
    if m3Ratio > 0:
        m3_ = int(m3Ratio*(m0_-m12_)+m12_+0.5)
        return modelTemplates[model_][0].substitute(m0=str(m0_),m12=str(m12_),m3=str(m3_))
    else:
        return modelTemplates[model_][0].substitute(m0=str(m0_),m12=str(m12_))

def buildSignalTuple (model,m0_,m12_,m3Ratio=-1):
    assert model in modelTemplates
    if m3Ratio > 0:
        m3_ = int(m3Ratio*(m0_-m12_)+m12_+0.5)
        return (m0_,m12_,m3_)+modelTemplates[model][1][3:]
    else:
        return (m0_,m12_)+modelTemplates[model][1][2:]

def getFromSignalString (sigString,field):
    parts = sigString.split('_')
    assert parts[0] in modelTemplates
    partsT = modelTemplates[parts[0]][0].template.split('_')
    assert len(parts) == len(partsT)
    for i, p in enumerate(partsT):
        if p == "${"+field+"}":  return int(parts[i])
    return None

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

def getSigYieldsNLO (btags,ht,met,msugraString,evtsMu,evtsEle):
    sigYields = {}
    for btag in btags:
        if btag == 'b1p':  continue
        sigYields[btag] = None
        evtMu = None
        if msugraString in evtsMu[btag][ht][met]:
            evtMu = evtsMu[btag][ht][met][msugraString]
        evtEle = None
        if msugraString in evtsEle[btag][ht][met]:
            evtEle = evtsEle[btag][ht][met][msugraString]
        # temporary fix: using yields instead of efficiencies
        if evtMu != None and evtEle != None:
            sigYields[btag] = evtMu + evtEle
    if 'b1p' in btags:
        if sigYields['b1'] != None and sigYields['b2'] != None:
            sigYields['b1p'] = sigYields['b1'] + sigYields['b2']
        else:
            sigYields['b1p'] = None
            #    for btag in sigYields:
            #       if sigYields[btag] == None:  sigYields[btag] = 0.
    return sigYields

def getSmoothedSigYieldsNLO (btags,ht,met,msugraString,msugraTuple,lumi,ratiosMu,yieldsMu,ratiosEle,yieldsEle):
    sigYields = {}
    for btag in btags:
        sigYields[btag] = 0.
        # temporary fix: using yields instead of efficiencies
        yMu = None
        if msugraString in ratiosMu[btag][ht][met] and msugraString in yieldsMu[btag][ht][met]:
            yMu = ratiosMu[btag][ht][met][msugraString]*yieldsMu[btag][ht][met][msugraString]
            sigYields[btag] += yMu
        yEle = None
        if msugraString in ratiosEle[btag][ht][met] and msugraString in yieldsMu[btag][ht][met]:
            yEle = ratiosEle[btag][ht][met][msugraString]*yieldsEle[btag][ht][met][msugraString]
            sigYields[btag] += yEle
        if yMu == None or yEle == None:
            print btag,ht,met,msugraString,yMu,yEle
    return sigYields

def getSmoothedSigYieldsNLOtmp (btags,ht,met,msugraString,msugraTuple,lumi,effMu,yieldsMu,effEle,yieldsEle):
    sigYields = {}
    btagEffs = { 'b0' : 0.001, 'b1' : 0.0015, 'b2' : 0.002 }
    for btag in btags:
        sigYields[btag] = 0.
        # temporary fix: using yields instead of efficiencies
        if msugraString in effMu[btag][ht][met] and effMu[btag][ht][met][msugraString] > 0:
            print btagEffs[btag],effMu[btag][ht][met][msugraString],yieldsMu[btag][ht][met][msugraString]
            yMu = btagEffs[btag]/effMu[btag][ht][met][msugraString]*yieldsMu[btag][ht][met][msugraString]
            sigYields[btag] += yMu
        if msugraString in effEle[btag][ht][met] and effEle[btag][ht][met][msugraString] > 0:
            yEle = btagEffs[btag]/effEle[btag][ht][met][msugraString]*yieldsEle[btag][ht][met][msugraString]
            sigYields[btag] += yEle
    return sigYields
