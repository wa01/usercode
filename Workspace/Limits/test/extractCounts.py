import re
import math
#
# extract observed and predicted numbers from latex tables
#
# column numbers and names for the different b-tag bins
#
btagPred = {}
btagPred[1] = "b0"
btagPred[3] = "b1"
btagPred[5] = "b2"
btagPred[7] = "binc"
btagObs = {}
btagObs[2] = "b0"
btagObs[4] = "b1"
btagObs[6] = "b2"
btagObs[8] = "binc"

inTable = False
label = ""
#
# regexps for begin / end of tables and the label
#
#reBegTable = re.compile(r"^\\begin.sidewaystable")
#reEndTable = re.compile(r"^\\end.sidewaystable")
reBegTable = re.compile(r"^\\begin.table")
reEndTable = re.compile(r"^\\end.table")
reLabel = re.compile(r"\\label")
#
# result dictionaries
#
countsPred = {}
countsObs = {}
errorsPredSign = {}
#
# Open tex file (TemplateFit from AN-11-379)
#
f = open("TemplateFit.tex","r")
for line in f:

    if not inTable:
        # Begin of table?
        mBegTable = re.search(reBegTable,line)
        if mBegTable:
            # set flag and clear working arrays
            inTable = True
            systkeys = []
            systs = {}
            continue
    
    else:
        # Found label? 
        mLabel = re.search(reLabel,line)
        if mLabel:
            # extract LaTex label name
            line = re.sub(r".*\\label.","",line)
            line = re.sub(r"}.*","",line)
            line = re.sub(r"tab:","",line)
            label = line
            idx = 0
            continue

        # Found end of table?
        mEndTable = re.search(reEndTable,line)
        if mEndTable:
            print label
            # convert variance into error
            if label.startswith("DataResultsDetails-"):
                print "converting for ",ht,met
                for sign in errPred:
                    for key in errPred[sign]:
                        if not key in errorsPredSign[ht][met]:  errorsPredSign[ht][met][key] = {}
                        errorsPredSign[ht][met][key][sign] = math.sqrt(errPred[sign][key])
#            sums = []
#            for key in systkeys:
#                if len(sums) == 0:
#                    for v in systs[key]: sums.append(v*v)
#                else:
#                    for i,v in enumerate(systs[key]):
#                        sums[i] = sums[i] + v*v
#                line = (key + ":").ljust(50)
#                for v in systs[key]:
#                    line= line + "%6.1f" % v
#                print line
#            line = "sum:".ljust(50)
#            for v in sums:
#                line= line + "%6.1f" % math.sqrt(v)
#            print line
            # reset variables
            label = ""
            inTable = False
            errPred = { 1 : {}, -1 : {} }
            ht = -1
            met = -1
            continue

        # In data table?
        if label.startswith("DataResultsDetails-"):
            # extract cuts from header line
            if line.find("multicolumn") != -1:
                # convert variance into error for previous ht/met pair
                for sign in errPred:
                    for key in errPred[sign]:
                        if not key in errorsPredSign[ht][met]:  errorsPredSign[ht][met][key] = {}
                        errorsPredSign[ht][met][key][sign] = math.sqrt(errPred[sign][key])
                ht = -1
                met = -1
                for p in line.split("$"):
                    if p.find("HT") != -1:
                        p = re.sub(r"[^0-9]+","",p)
                        ht = int(p)
#                        print "HT = ",ht
                    if p.find("ETmiss") != -1:
                        p = re.sub(r"[^0-9]+","",p)
                        met = int(p)
#                        print "MET = ",met
                if ht < 0 or met < 0: continue
                if not ht in countsPred:
                    print "adding ",ht," to countsPred"
                    countsPred[ht] = {}
                    errorsPredSign[ht] = {}
                    countsObs[ht] = {}
                if not met in countsPred[ht]:
                    countsPred[ht][met] = {}
                    errorsPredSign[ht][met] = {}
                    countsObs[ht][met] = {}
                errPred = { 1 : {} , -1 : {} }
                continue
            # split line in cells acc. to LaTex notation
            # assume 9 columns (label, 4x (pred+-err,obs) for b0, b1, b>=2 and inc)
            parts = line.split("&")
            if len(parts) != 9: continue
            # extract text in mathmode and remove (some) LaTex formatting
            for i in range(1,len(parts)):
                parts[i] = re.sub(r"^[^\$]*\$","",parts[i])
                parts[i] = re.sub(r"\$[^\$]+$","",parts[i])
#                parts[i] = re.subn(r"[ \^\_\+\-}]","",parts[i])[0]
                parts[i] = re.subn(r"[ \^\_}]","",parts[i])[0]
            # remove multirow environment for sum of +-
            parts[8] = re.sub(r".*{","",parts[8])
            parts[8] = re.sub(r"[^0-9]+$","",parts[8])
            # extract predictions
            for i in btagPred:
                # skip empty cells in inclusive column (only sum of +-)
                if i > 6 and parts[0].find("+") < 0: continue
                # split into value and errors
                vals = parts[i].split("{")
                if not btagPred[i] in countsPred[ht][met]:
                    countsPred[ht][met][btagPred[i]] = 0
                    errPred[1][btagPred[i]] = 0
                    errPred[-1][btagPred[i]] = 0
                # add prediction to count in HT/MET/btag bin
                countsPred[ht][met][btagPred[i]] = countsPred[ht][met][btagPred[i]] + float(vals[0])
                # add average error to variance in HT/MET/btag bin
                errP = float(vals[1])
                errM = float(vals[2])
                if errP < 0 and errM > 0:
                    errP = float(vals[2])
                    errM = float(vals[1])
                errPred[1][btagPred[i]] += errP*errP
                errPred[-1][btagPred[i]] += errM*errM
            # extract observed count
            for i in btagObs:
                # skip empty cells in inclusive column (only sum of +-)
                if i > 6 and parts[0].find("+") < 0: continue
                if not btagObs[i] in countsObs[ht][met]:
                    countsObs[ht][met][btagObs[i]] = 0
                # add observation to count in HT/MET/btag bin
                countsObs[ht][met][btagObs[i]] = countsObs[ht][met][btagObs[i]] + int(parts[i])
            
# close input file
f.close()

# combine positive and negative errors
errorsPred = {}
for ht in errorsPredSign:
    errorsPred[ht] = {}
    for met in errorsPredSign[ht]:
        errorsPred[ht][met] = {}
        for btag in errorsPredSign[ht][met]:
            errorsPred[ht][met][btag] = (errorsPredSign[ht][met][btag][1]+errorsPredSign[ht][met][btag][-1])/2.
# print result
print "countsPred = "+str(countsPred)
print "errorsPred = "+str(errorsPred)
print "countsObs = "+str(countsObs)

#
# cross checks
#
if len(countsObs.keys()) != len(countsPred.keys()) or  len(countsObs.keys()) != len(errorsPred.keys()):
    print "Inconsistency in #HT bins"
for ht in countsObs:
    if ( not ht in countsPred ) or ( not ht in errorsPred ):
        print "No HT=",ht," in predictions"
    if len(countsObs[ht].keys()) != len(countsPred[ht].keys()) or  len(countsObs[ht].keys()) != len(errorsPred[ht].keys()):
        print "Inconsistency in #MET bins for HT=",ht
    for met in countsObs[ht]:
        if ( not met in countsPred[ht] ) or ( not met in errorsPred[ht] ):
            print "No MET=",met," in predictions"
        for btag in btagPred:
            if ( not btagPred[btag] in countsObs[ht][met] ) or \
               ( not btagPred[btag] in countsPred[ht][met] ) or \
               ( not btagPred[btag] in errorsPred[ht][met] ):
                print "No btag=",btagPred[btag]," in a dictionary for HT=",ht,", MET=",met

for ht in countsObs:
    for met in countsObs[ht]:
        sumObs = 0
        sumPred = 0
        sumErr = 0
        for key in btagPred:
            btag = btagPred[key]
            if btag != 'binc':
                sumObs = sumObs + countsObs[ht][met][btag]
                sumPred = sumPred + countsPred[ht][met][btag]
                sumErr = sumErr + errorsPred[ht][met][btag]*errorsPred[ht][met][btag]
        if sumObs != countsObs[ht][met]['binc']:
            print "Error in observed sums for HT=",ht," MET=",met," : ",sumObs,countsObs[ht][met]['binc']
        if abs(sumPred-countsPred[ht][met]['binc'])>0.01*len(btagPred.keys()):
            print "Error in predicted sums for HT=",ht," MET=",met," : ",sumPred,countsPred[ht][met]['binc']
        sumErr = math.sqrt(sumErr)
        if abs(sumErr/errorsPred[ht][met]['binc']-1)>0.15:
            print "Error in predicted sums for HT=",ht," MET=",met," : ",sumErr,errorsPred[ht][met]['binc']
        
                
