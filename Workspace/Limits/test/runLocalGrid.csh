#! /bin/tcsh -f
#
# to be sourced in a job_ directory
# runs a grid based on the single signal point in model_760_400.txt
# (used to obtain a single limit on the observable signal for use with SMS)
#
set echo
set srcDir = ~/scratch0/CMSSW_4_2_8_patch7/src/Workspace/Limits/test
cd $srcDir
eval `scram runtime -csh`
python createMultiJobs.py --ht $1 --met $2 --btag $3 --m0s 760 --m12s 440 -t
if ( $status != 0 ) then
  exit 1
endif
#
# temporary directory
#
if ( $?WORKDIR ) then
  cd $WORKDIR
else
  set tmpdir = createGrid_$$
  if ( -d $tmpdir ) then
    echo "Removing $tmpdir"
    rm -r $tmpdir
  endif
  mkdir $tmpdir
  cd $tmpdir
endif

set jobDir = /tmp/adamwo/job_msugra_${3}_ht${1}_met${2}_m0_760_m12_440_HN
awk '{if ($1=="rate") {print $1,1.,$3} else{print $0}}' $jobDir/model_760_440.txt > model.txt
cat model.txt
rm -r $jobDir
#
# first step: asymptotic limit
#
combine model.txt -M Asymptotic >& asymptotic.log
if ( $status != 0 ) then
  echo "Asymptotic failed"
  cat asymptotic.log
  cd ..
  exit 1
endif
#
# read observed / expected limit from file and create
#  20 grid points
#
cat <<EOF > tmp.awk
BEGIN{
 obs=-1;
 ex=-1;
}
/Observed Limit:/{
obs=int(\$NF+0.5);
}
/Expected 50.*:/{
ex=int(\$NF+0.5)
}
END{
if ( obs<=0 || ex<=0 )  exit 1;
if ( ex>obs )  obs = ex;
rmax = 5 * obs;
dr = rmax / 20.;
line = "";
for ( i=0; i<20; ++i )  line = line " " (i+1)*dr;
print line;
}
EOF
set res = ( `awk -f tmp.awk asymptotic.log` )
if ( ( $status != 0 ) || ( $#res < 2 ) ) then
  echo "Failed to find result"
  cd ..
  exit 1
endif
echo "res = $res"
#
# generate toys and hybrid results for each point
#
foreach r ( $res )
  combine model.txt -H Asymptotic -M HybridNew --frequentist --testStat LHC --clsAcc 0 -s -1 --singlePoint $r --saveToys --saveHybridResult -T 500 -i 25 -n SP$r
end
# combine the grid outputs
hadd higgsCombineGrid.HybridNew.root higgsCombineSP*.root
rm higgsCombineSP*.root
#
# calculate full observed / expected limits from grid
# 
combine model.txt -H Asymptotic -M HybridNew --frequentist --testStat LHC -s -1 --grid higgsCombineGrid.HybridNew.root --readHybridResult -T 500 -i 10 >& obs.log
combine model.txt -H Asymptotic -M HybridNew --frequentist --testStat LHC -s -1 --grid higgsCombineGrid.HybridNew.root --readHybridResult -T 500 -i 10 --expectedFromGrid 0.16 >& exp16.log
combine model.txt -H Asymptotic -M HybridNew --frequentist --testStat LHC -s -1 --grid higgsCombineGrid.HybridNew.root --readHybridResult -T 500 -i 10 --expectedFromGrid 0.50 >& exp50.log
combine model.txt -H Asymptotic -M HybridNew --frequentist --testStat LHC -s -1 --grid higgsCombineGrid.HybridNew.root --readHybridResult -T 500 -i 10 --expectedFromGrid 0.84 >& exp84.log
combine model.txt -H Asymptotic -M HybridNew --frequentist --testStat LHC -s -1 --grid higgsCombineGrid.HybridNew.root --readHybridResult -T 500 -i 10 --expectedFromGrid 0.975 >& exp975.log
#mv -i higgsCombineGrid.HybridNew.root higgsCombineTest.HybridNew*.root *.log ..

hadd higgsCombineTest.HybridNew.root higgsCombineTest.HybridNew.mH120*.root
rm higgsCombineTest.HybridNew.mH120*.root

tar -zcvf $srcDir/hnGrid_${3}_ht${1}_met${2}.tgz *
#
# cleanup
#
if ( ! ( $?WORKDIR ) ) then
 cd ..
 rm -r $tmpdir
endif

