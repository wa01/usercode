import cPickle
import os
import sys

def numM0M12 (m0m12Dict):
    res = 0
    for m0 in m0m12Dict:
        res = res + len(m0m12Dict[m0])
    return res

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
            print "**** key not found in xsLO ",key
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
