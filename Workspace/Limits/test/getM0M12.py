import cPickle
import os
import sys

def numM0M12 (m0m12Dict):
    res = 0
    for m0 in m0m12Dict:
        res = res + len(m0m12Dict[m0])
    return res

def getM0M12c (sigTuples, rangeM0, rangeM12, regroupM0=1, regroupM12=1 ):

    from signalUtils import getFromSignalString

    if len(sigTuples) == 0:  return []
    isString = ( type(sigTuples[0]) == str )

    noRangeM0 = rangeM0[1] < rangeM0[0] 
    noRangeM12 = rangeM12[1] < rangeM12[0] 
#    if regroupM0 <= 1 and regroupM12 <= 1 and noRangeM0 and noRangeM12:  return sigTuples

    m0s = []
    m12s = []
    for tup in sigTuples:
        if isString:
            m0 = getFromSignalString(tup,"m0")
            m12 = getFromSignalString(tup,"m12")
        else:
            m0 = tup[0]
            m12 = tup[1]
        if noRangeM0 or ( m0 >= rangeM0[0] and m0 <= rangeM0[1] ):
            if not m0 in m0s:  m0s.append(m0)
        if noRangeM12 or ( m12 >= rangeM12[0] and m12 <= rangeM12[1] ):
            if not m12 in m12s:  m12s.append(m12)

    m0s.sort()
    if regroupM0 > 1:
        m0sToKeep = []
        for i in range(0,len(m0s),regroupM0):
            m0sToKeep.append(m0s[i])
    else:
        m0sToKeep = m0s

    m12s.sort()
    if regroupM12 > 1:
        m12sToKeep = []
        for i in range(0,len(m12s),regroupM12):
            m12sToKeep.append(m12s[i])
    else:
        m12sToKeep = m12s

    print "Total number of combinations = ",len(sigTuples)

    result = {}
    for tup in sigTuples:
        if isString:
            m0 = getFromSignalString(tup,"m0")
            m12 = getFromSignalString(tup,"m12")
        else:
            m0 = tup[0]
            m12 = tup[1]
        if m0 in m0sToKeep:
            if m12 in m12sToKeep:
                if not m0 in result:  result[m0] = []
                result[m0].append(m12)

    nacc = 0
    for m0 in result:
        nacc += len(result[m0])
        result[m0].sort()
    print "Kept number of combinations = ",nacc

    return result

def getM0M12b (sigTuples, rangeM0, rangeM12, regroupM0=1, regroupM12=1 ):

    if len(sigTuples) == 0:  return []
    isString = ( type(sigTuples[0]) == str )

    noRangeM0 = rangeM0[1] < rangeM0[0] 
    noRangeM12 = rangeM12[1] < rangeM12[0] 
    if regroupM0 <= 1 and regroupM12 <= 1 and noRangeM0 and noRangeM12:  return sigTuples

    m0s = []
    m12s = []
    for tup in sigTuples:
        if isString:
            m0 = getFromSignalString(tup,"m0")
            m12 = getFromSignalString(tup,"m12")
        else:
            m0 = tup[0]
            m12 = tup[1]
        if noRangeM0 or ( m0 >= rangeM0[0] and m0 <= rangeM0[1] ):
            if m0 in m0s:  m0s.append(m0)
        if noRangeM12 or ( m12 >= rangeM12[0] and m12 <= rangeM12[1] ):
            if m12 in m12s:  m12s.append(m12)
    
    m0s.sort()
    if regroupM0 > 1:
        m0sToKeep = []
        for i in range(0,len(m0s),regroupM0):
            m0sToKeep.append(m0s[i])
    else:
        m0sToKeep = m0s

    m12s.sort()
    if regroupM12 > 1:
        m12sToKeep = []
        for i in range(0,len(m12s),regroupM12):
            m12sToKeep.append(m12s[i])
    else:
        m12sToKeep = m12s

    print "Total number of combinations = ",len(sigTuples)

    result = []
    for tup in sigTuples:
        if isString:
            m0 = getFromSignalString(tup,"m0")
            m12 = getFromSignalString(tup,"m12")
        else:
            m0 = tup[0]
            m12 = tup[1]
        if m0 in m0sToKeep:
            if m12 in m12sToKeep:
                result.append(tup)

    print "Kept number of combinations = ",len(result)

    return result

def getM0M12a (effMu,effEle,xsecs,btag,ht,met,regroupM0=1, regroupM12=1):

#    from optparse import OptionParser
#    parser = OptionParser()
#    parser.add_option("--ht", dest="ht", type="int", action="store")
#    parser.add_option("--met", dest="met", type="int", action="store")
#    parser.add_option("--btag", dest="btag", default="binc", type="string", action="store")
#    (options, args) = parser.parse_args()

    m0s = []
    m0m12s = {}

    bt = btag
    if bt == 'binc':  bt = 'inc'
    for key in effMu[bt][ht][met]:
    #    print key
        parts = key.split('_')
    #    print len(parts)
        m0 = int(parts[1])
        m12 = int(parts[2])
        if not m0 in m0s:
            m0s.append(m0)
            m0m12s[m0] = []
        try:
            dmy = effEle[bt][ht][met][key]
        except:
    #        print "**** key not found in Ele ",key
            continue
#        if not  ( int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4]), int(parts[5]) ) in xsecs:
        tup = ( m0, m12 )
        for p in parts[3:]:
            tup = tup + ( int(p), )
#        print tup
        if not tup in xsecs:
#            print "**** key not found in xsLO ",key
            continue
        m0m12s[m0].append(m12)
        
    for m0 in m0m12s:
        m0m12s[m0].sort()

#
# reduce number of points by n
#
#    regroupM0 = 1
#    regroupM12 = 1
    m0s = m0m12s.keys()
    m0s.sort()
    if regroupM0 > 1:
        m0sToKeep = []
        for i in range(0,len(m0s),regroupM0):
            m0sToKeep.append(m0s[i])
    else:
        m0sToKeep = m0s
    m12s = []
    for key in m0m12s:
        for value in m0m12s[key]:
            if not value in m12s:  m12s.append(value)
    m12s.sort()
    if regroupM12 > 1:
        m12sToKeep = []
        for i in range(0,len(m12s),regroupM12):
            m12sToKeep.append(m12s[i])
    else:
        m12sToKeep = m12s

#    print m0sToKeep
#    print m12sToKeep
    print "Total number of combinations = ",numM0M12(m0m12s)

    m0m12sToKeep = {}
    for m0 in m0m12s:
        if not m0 in m0sToKeep: continue
        m0m12sToKeep[m0] = []
        for m12 in m0m12s[m0]:
            if m12 in m12sToKeep:  m0m12sToKeep[m0].append(m12)
        m0m12sToKeep[m0].sort()

    print "Kept number of combinations = ",numM0M12(m0m12sToKeep)

    return m0m12sToKeep

def getM0M12 (btag,ht,met,regroupM0=1, regroupM12=1):

#    from optparse import OptionParser
#    parser = OptionParser()
#    parser.add_option("--ht", dest="ht", type="int", action="store")
#    parser.add_option("--met", dest="met", type="int", action="store")
#    parser.add_option("--btag", dest="btag", default="binc", type="string", action="store")
#    (options, args) = parser.parse_args()

    fMu = open("Mu_msugraEfficiencies.pkl","rb")
    effMu = cPickle.load(fMu)
    fMu.close()

    fEle = open("Ele_msugraEfficiencies.pkl","rb")
    effEle = cPickle.load(fEle)
    fEle.close()

    fXS = open("goodModelNames_10_0_1.pkl","rb")
    xsecs = cPickle.load(fXS)
    fXS.close()

    return getM0M12a(effMu,effEle,xsecs,btag,ht,met,regroupM0,regroupM12)
