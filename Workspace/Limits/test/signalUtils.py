from string import Template
#
# small utility functions for signal decoding
#
channelNames = [ "Ele", "Mu" ]
modelTemplates = { "msugra": [ Template("msugra_${m0}_${m12}_10_0_1"), ( "{m0}", "{m12}", 10, 0, 1 ) ],
                   "T1tttt": [ Template("T1tttt_${m0}_${m12}_-1"), ( "{m0}", "{m12}", -1 ) ] }
def effFileName (channel,model):
    assert channel in channelNames
    assert model in modelTemplates
    if model == "msugra":
#        return channel+"_"+model+"Efficiencies.pkl"
        return channel+"_"+model+"_LO_efficiency.pkl"
    else:
        return channel+"_"+model+"_Efficiencies.pkl"
    
def xsecFileName (model):
    assert model in modelTemplates
    if model == "msugra":
        return "goodModelNames_10_0_1.pkl"
    else:
        return "xsec"+model+".pkl"
    
    
def signalString (model,m0_,m12_):
    assert model in modelTemplates
    return modelTemplates[model][0].substitute(m0=str(m0_),m12=str(m12_))

def signalTuple (model,m0_,m12_):
    assert model in modelTemplates
    return (m0_,m12_)+modelTemplates[model][1][2:]
