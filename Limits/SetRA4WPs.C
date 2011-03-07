#include "RA4WorkingPoint.h"

RA4WorkingPoint* eleTightWP () {
  RA4WorkingPoint* result = new RA4WorkingPoint(76.8,6.5,29.5,2.9,0.214,
						"/afs/hephy.at/user/a/adamwo/www/pngMSUGRA/kinMetSig_2.5_5.5_5.5_ht_300_650_650_",
						"Ele-tanb-3_A0_0_signMu_1_signalYield.root",
						"Ele-tanb-3_A0_0_signMu_1_kfactors.root",0.20);
  result->setObserved(80,4,30,2);
  return result;
}

RA4WorkingPoint* muTightWP () {
  RA4WorkingPoint* result = new RA4WorkingPoint(93.1,8.7,37.6,3.4,0.158,
						"/afs/hephy.at/user/a/adamwo/www/pngMSUGRA/kinMetSig_2.5_5.5_5.5_ht_300_650_650_",
						"Mu-tanb-3_A0_0_signMu_1_signalYield.root",
						"Mu-tanb-3_A0_0_signMu_1_kfactors.root",0.20);
  result->setObserved(98,4,41,5);
  return result;
}
