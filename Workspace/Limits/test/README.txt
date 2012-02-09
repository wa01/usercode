Sequence for limit calculation:

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

