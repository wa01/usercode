#
# to be sourced in a job_ directory
# runs a grid based on the single signal point in model_760_400.txt
# (used to obtain a single limit on the observable signal for use with SMS)
#
#
# temporary directory
#
set tmpdir = createGrid_$$
if ( -d $tmpdir ) then
  echo "Removing $tmpdir"
  rm -r $tmpdir
endif
mkdir $tmpdir
cd $tmpdir

#
# first step: asymptotic limit
#
combine ../model_760_400.txt -M Asymptotic >& asymptotic.log
if ( $status != 0 ) then
  echo "Asymptotic failed"
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
for ( i=0; i<20; ++i )  printf("%6.2f",(i+1)*dr);
printf("\n");
}
EOF
set res = ( `awk -f tmp.awk asymptotic.log` )
if ( ( $status != 0 ) || ( $#res < 2 ) ) then
  echo "Failed to find result"
  cd ..
  exit 1
endif
#
# generate toys and hybrid results for each point
#
foreach r ( $res )
  combine ../model_760_400.txt -H Asymptotic -M HybridNew --frequentist --testStat LHC --clsAcc 0 -s -1 --singlePoint $r --saveToys --saveHybridResult -T 500 -i 10 -n SP$r
end
# combine the grid outputs
hadd higgsCombineGrid.HybridNew.root higgsCombineSP*.root
#
# calculate full observed / expected limits from grid
# 
combine ../model_760_400.txt -H Asymptotic -M HybridNew --frequentist --testStat LHC -s -1 --grid higgsCombineGrid.HybridNew.root --readHybridResult -T 500 -i 10 >& obs.log
combine ../model_760_400.txt -H Asymptotic -M HybridNew --frequentist --testStat LHC -s -1 --grid higgsCombineGrid.HybridNew.root --readHybridResult -T 500 -i 10 --expectedFromGrid 0.16 >& exp16.log
combine ../model_760_400.txt -H Asymptotic -M HybridNew --frequentist --testStat LHC -s -1 --grid higgsCombineGrid.HybridNew.root --readHybridResult -T 500 -i 10 --expectedFromGrid 0.50 >& exp50.log
combine ../model_760_400.txt -H Asymptotic -M HybridNew --frequentist --testStat LHC -s -1 --grid higgsCombineGrid.HybridNew.root --readHybridResult -T 500 -i 10 --expectedFromGrid 0.84 >& exp84.log
combine ../model_760_400.txt -H Asymptotic -M HybridNew --frequentist --testStat LHC -s -1 --grid higgsCombineGrid.HybridNew.root --readHybridResult -T 500 -i 10 --expectedFromGrid 0.975 >& exp975.log
mv -i higgsCombineGrid.HybridNew.root higgsCombineTest.HybridNew*.root *.log ..
#
# cleanup
#
cd ..
rm -r $tmpdir
