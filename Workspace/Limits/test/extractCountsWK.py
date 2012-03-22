from string import Template
from math import sqrt
#
# extract observed counts, predicted counts and stat. errors from WKs summary files ("res.py")
#
#
# template for directory / file names
#
template = Template('/scratch/kwolf/fitTuple_120319/copy3j_highM0_dataset_${lep}_with_msugra/Data_${lep}_newModel_V120212_htSig-${ht}_metSig-${met}_RooMinuit_sampleOff/res_extended.py')
#
# lepton flavour, HT, MET and btag bins
#
leptons = [ 'Mu', 'Ele' ]
hts = [ 750, 1000 ]
mets = [ 250, 350, 450, 550 ]
btags = [ 'binc', 'b0', 'b1', 'b1p', 'b2' ]
#
# output dictionaries:
#
countsObs = {}          # observed counts in signal region
countsNorm = {}         # observed counts in normalization region
countsPred = {}         # predicted background in signal region
errorsPred = {}         # stat. error on pred. background in signal region
errorsPredUp = {}         # stat. error on pred. background in signal region
errorsPredDown = {}         # stat. error on pred. background in signal region
for lep in leptons:
    for ht in hts:
        if not ht in countsObs:
            countsObs[ht] = {}
            countsNorm[ht] = {}
            countsPred[ht] = {}
            errorsPred[ht] = {}
            errorsPredUp[ht] = {}
            errorsPredDown[ht] = {}
        for met in mets:
            if not met in countsObs[ht]:
                countsObs[ht][met] = {}
                countsNorm[ht][met] = {}
                countsPred[ht][met] = {}
                errorsPred[ht][met] = {}
                errorsPredUp[ht][met] = {}
                errorsPredDown[ht][met] = {}
            # read input file (one file / lepton flavour / HT / MET)
            ldict = {}
            fname = template.substitute(lep=lep,ht=str(ht),met=str(met))
            execfile(fname,ldict)
            for btag in btags:
                if not btag in countsObs[ht][met]:
                    countsObs[ht][met][btag] = { }
                    countsNorm[ht][met][btag] = { }
                    countsPred[ht][met][btag] = { }
                    errorsPred[ht][met][btag] = { }
                    errorsPredUp[ht][met][btag] = { }
                    errorsPredDown[ht][met][btag] = { }
                # translate name for inclusive btag bin
                #   all btag bins (except 'inc') are split by lepton charge
                if btag == 'binc':
                    bt = 'inc'
                    signs = [ "" ]
                elif btag == 'b1p':
                    bt = 'b12'
                    signs = [ "" ]
                else:
                    bt = btag
#                    signs = [ "Minus", "Plus" ]
                    signs = [ "PM" ]
                # initialise sums
                countsObs[ht][met][btag][lep] = 0
                countsNorm[ht][met][btag][lep]  = 0
                countsPred[ht][met][btag][lep] = 0. 
                errorsPred[ht][met][btag][lep] = 0.
                errorsPredUp[ht][met][btag][lep] = 0.
                errorsPredDown[ht][met][btag][lep] = 0.
                # upper / lower errors
                ePred = { "Hi" : 0., "Low" : 0. }
                # loop over lepton charges
                for sign in signs:
                    bts = bt+sign
                    # counts: add lepton charges
                    countsObs[ht][met][btag][lep] += ldict['count'][bts]
                    countsNorm[ht][met][btag][lep] += ldict['countNorm'][bts]
                    countsPred[ht][met][btag][lep] += ldict['prediction'][bts]
                    # errors: quadratic sum for upper / lower error
                    for var in ePred:
                        e = ldict['predictionError'+var][bts]
                        ePred[var] += e*e
                # final error: average of combined upper / lower errors
                for var in ePred:  ePred[var] = sqrt(ePred[var])
                errorsPredUp[ht][met][btag][lep] = ePred["Hi"]
                errorsPredDown[ht][met][btag][lep] = ePred["Low"]
                errorsPred[ht][met][btag][lep] = (ePred["Hi"]+ePred["Low"])/2.
                
print "countsObs =",countsObs
print "countsNorm =",countsNorm
print "countsPred =",countsPred
print "errorsPred =",errorsPred
print "errorsPredUp =",errorsPredUp
print "errorsPredDown =",errorsPredDown

