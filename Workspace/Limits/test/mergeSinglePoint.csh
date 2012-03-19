#! /bin/tcsh -f
#
if ( $#argv < 1 ) then
  echo "Missing argument"
  exit 1
endif
if ( $?CMSSW_RELEASE_BASE == 0 ) then
  cmsenv
endif

set jobdir = $1
if ( $jobdir =~ */ ) set jobdir = $jobdir:h
if ( ! ( -d $jobdir ) ) then
  echo "No such directory $jobdir"
  exit 1
endif

set limargs = ( ' ' '--expectedFromGrid=0.16' '--expectedFromGrid=0.50' '--expectedFromGrid=0.84'  '--expectedFromGrid=0.025' '--expectedFromGrid=0.975' )
set ofile = `echo $jobdir | sed 's/job_/limits_/'`
set ofile = /tmp/adamwo/${ofile}_comb.root
if ( -e $ofile ) then
  echo "Output file $ofile exists"
  exit 1
endif

set files = ( )

if ( $?localDir )  unset localDir
set castorDir = `echo $jobdir | sed 's/job_//'`
set castorDir = "${CASTOR_HOME}/CrabOutput/${castorDir}"
nsls $castorDir >& /dev/null
if ( ( $status == 0 ) && ( $#argv == 1 ) ) then
#  echo "No such CASTOR directory: $castorDir"
#  exit 1
#endif
  set files = ( `nsls $castorDir | grep -E 'outputToy.*\.tgz' | sort -t'_' -r -n -k2,2 -k3,3` )
else
  echo "No such CASTOR directory: $castorDir - trying job area"
  ls $jobdir/crab_0_* >& /dev/null
  if ( $status != 0 ) then
    echo "No crab directory in $jobdir"
    exit 1
  endif
  if ( $#argv < 2 ) then
    set crabs = `ls -d $jobdir/crab_0_*`
  else
    set crabs = $2
  endif
  if ( $#crabs != 1 ) then
    echo "More than one crab directory in $jobdir"
    exit 1
  endif
  if ( !( -d $crabs ) ) then
    echo "No such directory $crabs"
    exit 1
  endif
  pushd $crabs/res
  set localDir = $PWD
  ls outputToy*.tgz >& /dev/null
  if ( $status == 0 ) then
    set files = ( `ls outputToy*.tgz` )
  endif
  popd
endif

if ( $#files == 0 ) then
  echo "No output files"
  exit 1
endif
#echo $files
#exit 0

set curdir = $PWD
set tmpdir = /tmp/adamwo/mergeSinglePoint_$$
mkdir $tmpdir
pushd $tmpdir

@ ntot = 0
@ imax = 0
@ icur = 0
foreach file ( $files ) 
    set parts = ( `echo $file | sed 's/[_\.]/ /g'` )
    if ( $#parts != 5 ) then
      echo "Wrong format of filename : $file"
      exit 1
    endif
    if ( $imax == 0 )  @ imax = $parts[2]
    if ( $parts[2] == $icur )  continue
    @ icur = $parts[2]
    @ ntot = $ntot + 1
    echo $file

    if ( $?localDir ) then
      cp $localDir/$file .
    else
     rfcp $castorDir/$file .
    endif
    tar -zxf $file
    if ( $status != 0 ) then
      echo "failed to unpack output file $file"
      break
    endif
    pushd outputToy
    pwd
    ls
    ls higgsCombineToys*.root >& /dev/null
    if ( $status != 0 ) then
      echo "No output root files for $file \!\!\!\!\!\!\!"
      popd
      rm -r outputToy
      continue
    endif
    foreach tfile ( higgsCombineToys*.root )
      set fields = ( `echo $tfile | sed 's/\./ /g'` )
      @ m0m12 = $fields[4]
      @ m0 = $m0m12 / 10000
      @ m12 = $m0m12 % 10000
      set model = "model_${m0}_${m12}.root"
      echo $m0,$m12
      tar -xvf $curdir/$jobdir/models.tar $model
      if ( !( -e $model ) ) then
        echo "Could not unpack $model"
        exit 1
       endif
       @ ilim = 0
       while ( $ilim < $#limargs )
         @ ilim = $ilim + 1
         combine $model -M HybridNew --frequentist --testStat LHC --singlePoint 1 -s $fields[4] --clsAcc 0 -T 100 -i 50 \
           --toysFile $tfile --readHybridResult $limargs[$ilim]
         if ( $status != 0 ) then
           echo "FAILED :  combine $model ... --singlePoint 1 -s $fields[4] ... --toysFile $tfile --readHybridResult $limargs[$ilim]"
           exit 1
         endif
       end
    end
    set cfiles = ( higgsCombineTest.*.root )
    if ( $#cfiles == 0 ) then  
      echo "No combine output files?"
      exit 1
    endif
    if ( -e $ofile ) then
       cp $ofile tmp.root
       set cfiles = ( tmp.root $cfiles )
       rm $ofile
    endif
    hadd $ofile $cfiles
    popd
    rm -r outputToy
    rm $file
#    if ( $ntot > 2 ) break

end
echo "Nr. of files: ${ntot} ; highest index = ${imax}"

popd
rm -r $tmpdir

exit
