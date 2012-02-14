import re
import math
#
# extract observed and predicted numbers from latex tables
#
#
# result dictionaries
#
errPred = { 'Mu' : { 1 : {}, -1 : {} }, 'Ele' : { 1 : {}, -1 : {} } }
countsPred = {}
countsObs = {}
errorsPredSign = {}
#
# HT and MET values
#
ht = -1
met = -1
#
def storeHtMet ():
    global ht,met,errPred,errorsPredSign
    if ht < 0 or met < 0:  return
    print "converting for ",ht,met
    for lep in errPred:
        for sign in errPred[lep]:
            for btag in errPred[lep][sign]:
                if not btag in errorsPredSign[ht][met]:  errorsPredSign[ht][met][btag] = {}
                if not lep in errorsPredSign[ht][met][btag]:  errorsPredSign[ht][met][btag][lep] = {}
                errorsPredSign[ht][met][btag][lep][sign] = math.sqrt(errPred[lep][sign][btag])
    ht = -1
    met = -1
    errPred = { 'Mu' : { 1 : {}, -1 : {} }, 'Ele' : { 1 : {}, -1 : {} } }
#
def findHtMet (line):
    global ht,met,countsObs,countsPred,errorsPredSign
    ht = -1
    met = -1
    for p in line.split("$"):
        if p.find("HT") != -1:
            p = re.sub(r"[^0-9]+","",p)
            ht = int(p)
#            print "HT = ",ht
        if p.find("ETmiss") != -1:
            p = re.sub(r"[^0-9]+","",p)
            met = int(p)
#                        print "MET = ",met
    if ht < 0 or met < 0: return
    if not ht in countsPred:
        print "adding ",ht," to countsPred"
        countsPred[ht] = {}
        errorsPredSign[ht] = {}
        countsObs[ht] = {}
    if not met in countsPred[ht]:
        print "adding ",met," to countsPred[",ht,"]"
        countsPred[ht][met] = {}
        errorsPredSign[ht][met] = {}
        countsObs[ht][met] = {}
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
        continue
    
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
        if label.startswith("DataResultsDetails-"): storeHtMet()
        # reset variables
        label = ""
        inTable = False
        continue

    # In data table?
    if label.startswith("DataResultsDetails-"):
        # extract cuts from header line
        if line.find("multicolumn") != -1:
            storeHtMet()
            findHtMet(line)
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
        # lepton type
        if parts[0].find("\mu") != -1:
            lep = "Mu"
        elif parts[0].find("\Pe") != -1:
            lep = "Ele"
        else:
            print "Unknown lepton channel ",parts[0]
            lep = "Unknown"
        # extract predictions
        for i in btagPred:
            # skip empty cells in inclusive column (only sum of +-)
            if i > 6 and parts[0].find("+") < 0: continue
            # split into value and errors
            vals = parts[i].split("{")
            if not btagPred[i] in countsPred[ht][met]:
                countsPred[ht][met][btagPred[i]] = {}
            if not lep in countsPred[ht][met][btagPred[i]]:
                countsPred[ht][met][btagPred[i]][lep] = 0.
            if not btagPred[i] in errPred[lep][1]:
                errPred[lep][1][btagPred[i]] = 0
                errPred[lep][-1][btagPred[i]] = 0
            # add prediction to count in HT/MET/btag bin
            countsPred[ht][met][btagPred[i]][lep] += float(vals[0])
            # add average error to variance in HT/MET/btag bin
            errP = float(vals[1])
            errM = float(vals[2])
            if errP < 0 and errM > 0:
                errP = float(vals[2])
                errM = float(vals[1])
            errPred[lep][1][btagPred[i]] += errP*errP
            errPred[lep][-1][btagPred[i]] += errM*errM
        # extract observed count
        for i in btagObs:
            # skip empty cells in inclusive column (only sum of +-)
            if i > 6 and parts[0].find("+") < 0: continue
            if not btagObs[i] in countsObs[ht][met]:
                countsObs[ht][met][btagObs[i]] = {}
            if not lep in countsObs[ht][met][btagObs[i]]:
                countsObs[ht][met][btagObs[i]][lep] = 0
            # add observation to count in HT/MET/btag bin
            countsObs[ht][met][btagObs[i]][lep] += int(parts[i])
            
# close input file
f.close()

# combine positive and negative errors
errorsPred = {}
for ht in errorsPredSign:
    errorsPred[ht] = {}
    for met in errorsPredSign[ht]:
        errorsPred[ht][met] = {}
        for btag in errorsPredSign[ht][met]:
            errorsPred[ht][met][btag] = {}
            for lep in errorsPredSign[ht][met][btag]:
                errorsPred[ht][met][btag][lep] = (errorsPredSign[ht][met][btag][lep][1]+errorsPredSign[ht][met][btag][lep][-1])/2.
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
            for lep in errPred:
                if ( not lep in countsObs[ht][met][btagPred[btag]] ) or \
                   ( not lep in countsPred[ht][met][btagPred[btag]] ) or \
                   ( not lep in errorsPred[ht][met][btagPred[btag]] ):
                    print "No lep=",lep," in a dictionary for HT=",ht,", MET=",met," btag=",btagPred[btag]
                
                
for lep in errPred:
    for ht in countsObs:
        for met in countsObs[ht]:
            sumObs = 0
            sumPred = 0
            sumErr = 0
            for key in btagPred:
                btag = btagPred[key]
                if btag != 'binc':
                    sumObs = sumObs + countsObs[ht][met][btag][lep]
                    sumPred = sumPred + countsPred[ht][met][btag][lep]
                    sumErr = sumErr + errorsPred[ht][met][btag][lep]*errorsPred[ht][met][btag][lep]
            if sumObs != countsObs[ht][met]['binc'][lep]:
                print "Error in observed sums for HT=",ht," MET=",met," : ",sumObs,countsObs[ht][met]['binc'][lep]
            if abs(sumPred-countsPred[ht][met]['binc'][lep])>0.01*len(btagPred.keys()):
                print "Error in predicted sums for HT=",ht," MET=",met," lep=",lep," : ",sumPred,countsPred[ht][met]['binc'][lep]
            sumErr = math.sqrt(sumErr)
            if abs(sumErr/errorsPred[ht][met]['binc'][lep]-1)>0.15:
                print "Error in predicted errs for HT=",ht," MET=",met," lep=",lep," : ",sumErr,errorsPred[ht][met]['binc'][lep]
        

for ht in countsObs:
    for met in countsObs[ht]:
        print "            ht = ",ht,"         met = ",met
        for lep in errPred:
            out = lep.ljust(5)
            for btag in [ 'b0', 'b1', 'b2', 'binc' ]:
                out = out + "%5.2f" % countsPred[ht][met][btag][lep] + "%5.2f" % errorsPred[ht][met][btag][lep] + "%5d" % countsObs[ht][met][btag][lep] + " | "
            print out
    
