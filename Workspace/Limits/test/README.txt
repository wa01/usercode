Sequence for limit calculation (all executed in /test):

# get observed counts and prediction (from tex -> to python)
python extractCounts.py
  Input: TemplateFit.tex
  Output: eventCounts.py

# generate template cards file for limit calculation 
#  (including all background info - signal is added when creating jobs)
#  only needs to be redone if background or systematics change
python createCardsFromDict.py
 Mandatory options: ht and met cuts
 Optional: btag bin (default: binc)
 Input: eventCounts.py systematics_htSig*_metSig*.py systematics_BT_htSig*_metSig*.py
 Output: <btag>-ht*-met*.txt


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
