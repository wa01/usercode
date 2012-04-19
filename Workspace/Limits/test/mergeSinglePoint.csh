#! /bin/tcsh -f
#
# create files with observed and expected limits from root files with toys (HybridNew)
# input are outputToy*.tgz tar files with the results (log & root files) of the jobs
# arguments: 
#  #1 : name of the job
#  #2 : crab working directory (optional)
# If only arg1 is given the job name will first be looked up in CASTOR. If no CASTOR directory
#   with name == job name - the "job_" prefix is found or if arg2 is given the files will be
#   taken from the "res" subdirectory
#
#
# check arguments and CMSSW environment
#
#set echo
if ( $#argv < 1 ) then
  echo "Missing argument"
  exit 1
endif
if ( $?CMSSW_RELEASE_BASE == 0 ) then
  cmsenv
endif
#
# check directory (and remove trailing "/", if necessary)
#
set jobdir = $1
if ( $jobdir =~ */ ) set jobdir = $jobdir:h
if ( ! ( -d $jobdir ) ) then
  echo "No such directory $jobdir"
  exit 1
endif
#
# list of limit types
#
set limargs = ( ' ' '--expectedFromGrid=0.16' '--expectedFromGrid=0.50' '--expectedFromGrid=0.84'  '--expectedFromGrid=0.025' '--expectedFromGrid=0.975' )
#
# output file name
#
set ofile = `echo $jobdir | sed 's/job_/limits_/'`
set ofile = /tmp/adamwo/${ofile}_comb.root
if ( -e $ofile ) then
  echo "Output file $ofile exists"
  exit 1
endif
#
# look for input files in CASTOR or in a local directory
#
set files = ( )
if ( $?localDir )  unset localDir
set castorDir = `echo $jobdir | sed 's/job_//'`
set castorDir = "${CASTOR_HOME}/CrabOutput/${castorDir}"
nsls $castorDir >& /dev/null
if ( ( $status == 0 ) && ( $#argv == 1 ) ) then
#  echo "No such CASTOR directory: $castorDir"
#  exit 1
#endif
# get list of files from CASTOR
  set files = ( `nsls $castorDir | grep -E 'outputToy.*\.tgz' | sort -t'_' -r -n -k2,2 -k3,3` )
###  set files = ( `nsls -l $castorDir | grep -F 'Apr 16' | awk '/outputToy/{print $NF}' | sort -t'_' -r -n -k2,2 -k3,3` )
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
  # get list of files from "res" directory
  pushd $crabs/res
  set localDir = $PWD
  ls outputToy*.tgz >& /dev/null
  if ( $status == 0 ) then
    set files = ( `ls outputToy*.tgz | sort -t'_' -r -n -k2,2 -k3,3` )
  endif
  popd
endif
#
# any file found?
#
if ( $#files == 0 ) then
  echo "No output files"
  exit 1
endif
#echo $files
#exit 0
#
# move to temporary working directory
#
set curdir = $PWD
set tmpdir = /tmp/adamwo/mergeSinglePoint_$$
mkdir $tmpdir
pushd $tmpdir
#
# create file for running the limits
#
cat <<EOF >runCombine.csh
#\! /bin/tcsh -f
set echo
set curDir = \$PWD
set tmpDir = /tmp/adamwo/runCombine_\$\$
mkdir \$tmpDir
pushd \$tmpDir
combine \$curDir/\$1 -M HybridNew --frequentist --testStat LHC --singlePoint 1 -s \$2 --clsAcc 0 \
  -T 100 -i 50 --toysFile \$curDir/\$3 --readHybridResult \$4
if ( \$status != 0 ) then
  echo "runCombine FAILED :  combine \$1 ... --singlePoint 1 -s \$2 ... --toysFile \$3 --readHybridResult \$4"
   exit 1
else
  cp higgsCombineTest.*.root \$curDir
endif
popd
rm -r \$tmpDir
exit 0
EOF
chmod u+x runCombine.csh
#
#
# loop of tar files
#
@ ntot = 0
@ imax = 0
@ icur = 0
foreach file ( $files ) 
    date
    #
    # decode file name: expect 5 fields with job# in the 2nd field
    #
    echo ";$file;"
    echo $file | sed 's/[_\.]/ /g'
    if ( $?parts )  unset parts
    set parts = ( `echo $file | sed 's/[_\.]/ /g'` )
    echo $parts
    echo $#parts
    if ( $#parts != 5 ) then
      echo "Wrong format of filename : $file"
      exit 1
    endif
    # skip multiple output files for the same job
    if ( $parts[2] == $icur ) then
      echo "Found multiple outputs for job $parts[2], skipping $file"
      continue
#      exit 1
    endif
    @ icur = $parts[2]
    @ ntot = $ntot + 1
    echo $file
    date
    # copy file to working directory
    if ( $?localDir ) then
      cp $localDir/$file .
    else
     rfcp $castorDir/$file .
    endif
    date
    # untar file
    tar -zxf $file
    if ( $status != 0 ) then
      echo "failed to unpack output file $file"
      if ( -d outputToy ) then
        rm -r outputToy
      endif
      rm $file
      @ icur = 0
      continue
    endif
    date
    # move to directory with output files and get list of root files
    pushd outputToy
#    pwd
#    ls
    ls higgsCombineToys*.root >& /dev/null
    if ( $status != 0 ) then
      echo "No output root files for $file \!\!\!\!\!\!\!"
      popd
      rm -r outputToy
      rm $file
      @ icur = 0
      continue
    endif
    #
    # loop over root files (first try to get list of root files with limits 
    # from all points of the job in case they were done by the job)
    #
#    rm higgsCombineTest.HybridNew*.root
    ls higgsCombineTest.HybridNew*.root >& /dev/null
    if ( $status != 0 ) then
      foreach tfile ( higgsCombineToys*.root )
        # decode file name: expect 10000*m0+m12 in 4th field
        set fields = ( `echo $tfile | sed 's/\./ /g'` )
        @ m0m12 = $fields[4]
        @ m0 = $m0m12 / 10000
        @ m12 = $m0m12 % 10000
        # build name of input workspace and retrieve the input file
        set model = "model_${m0}_${m12}.root"
        echo $m0,$m12
        date
        tar -xvf $curdir/$jobdir/models.tar $model
        if ( !( -e $model ) ) then
          echo "Could not unpack $model"
          exit 1
         endif
         # check in log file, if toys were written (i.e., no crash in combine)
         set errflg
         set log = "combine_${m0}_${m12}.log"
         if ( -e $log ) then
           grep -qF 'Hybrid result saved as ' $log
           if ( $status == 0 )  unset errflg
         endif
         if ( $?errflg ) then
           echo "No hybrid result for ${tfile}"
           continue
         endif
         #
         # loop over all types of limits
         #
         date
         @ ilim = 0
         while ( $ilim < $#limargs )
           #  calculate limit using the toys
           @ ilim = $ilim + 1
           time $tmpdir/runCombine.csh $model $fields[4] $tfile $limargs[$ilim] &
#           combine $model -M HybridNew --frequentist --testStat LHC --singlePoint 1 -s $fields[4] --clsAcc 0 -T 100 -i 50 \
#             --toysFile $tfile --readHybridResult $limargs[$ilim]
#           if ( $status != 0 ) then
#             echo "FAILED :  combine $model ... --singlePoint 1 -s $fields[4] ... --toysFile $tfile --readHybridResult $limargs[$ilim]"
#             exit 1
#           endif
         end
         wait
      end
      date
    endif
    # get list of root files with limits from all points of the job
    ls higgsCombineTest.HybridNew*.root >& /dev/null
    if ( $status != 0 ) then
      echo "No combine output files?"
      popd
      rm -r outputToy
      @ icur = 0
      continue
#      exit 1
    endif
    set cfiles = ( higgsCombineTest.HybridNew*.root )
    if ( $imax == 0 )  @ imax = $parts[2]
#    if ( $#cfiles == 0 ) then  
#    endif
    # prepare for hadd: if output file exists rename it and
    #   include it in the list of files to be added
    if ( -e $ofile ) then
       cp $ofile tmp.root
       set cfiles = ( tmp.root $cfiles )
       rm $ofile
    endif
    # combine previous results and all limits from the current job in the output file
    if ( $#cfiles == 1 ) then
      mv $cfiles $ofile
    else
      echo hadd $ofile $cfiles
      hadd $ofile $cfiles
    endif
    date
    # clean temporary files from this job
    popd
    rm -r outputToy
    rm $file
#    if ( $ntot > 2 ) break

end
echo "Nr. of files: ${ntot} ; highest index = ${imax}"

popd
rm -r $tmpdir

exit
