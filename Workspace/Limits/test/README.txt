Sequence for limit calculation (all executed in /test):
#
# get observed counts and prediction (from tex -> to python)
#
python extractCounts.py
  Input: TemplateFit.tex
  Output: eventCounts.py

# generate template cards file for limit calculation 
#  (including all background info - signal is added when creating jobs)
#  only needs to be redone if background or systematics change

# for a single b-tag bin
python createCardsFromDict.py
 Mandatory options: ht and met cuts
 Optional: btag bin (default: binc)
 Input: eventCounts.py systematics_htSig*_metSig*.py systematics_BT_htSig*_metSig*.py
 Output: <btag>-ht*-met*.txt
# OR: for multiple b-tag bins (0, 1, >=2)
python createMultiCardsFromDict.py
 Mandatory options: ht and met cuts
 Input: eventCounts.py systematics_htSig*_metSig*.py systematics_BT_htSig*_metSig*.py
 Output: multibtag-ht*-met*.txt

#
# create directory for crab job submission
#
python createMultiJobs.py
  Mandatory options: ht and met cuts
  Optional: see createMultiJobs.py -h
  Input: efficiencies: <Ele/Mu>_msugra_LO_Efficiency.pkl or <Ele/Mu>_<SMS>_Efficiencies.pkl
         cross sections: goodModelNames_10_0_1.pkl or xsec<SMS>.pkl
  Other code: signalUtils.py getM0M12.py createCards.py HiggsAnalysis/CombinedLimit/python/DatacardParser.py
  Output: directory /tmp/adamwo/job_<parameters> 
          (to be moved to /test and used for submission; contains a cfg file)

# merge outputs from all jobs
./mergeLimits.csh <job_working_directory>
  Input: the tgz files in the <job_working_directory>/crab_*/res directory 
  Output: limits_<job_working_directory_without_prefix>.root

# draw results
root -l limits_<job>.root
p = new PlotLimits(_file0,1.); // choose 0.05 instead of 1 if HN "single point" was used
p->Loop();
p->drawHistograms(); // saves pdf version of canvas 
p->saveContours();   // saves contours for the final limit plot
