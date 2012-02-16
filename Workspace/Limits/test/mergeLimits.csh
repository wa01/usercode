#! /bin/tcsh
#
# merge output files from higgs limits jobs
#   unpacks each job output archive and hadd's the root files
#   argument: working directory of job (without trailing / !)
#   output: limits_<jobname>.root
#   
#
#cd ~/scratch0/CMSSW_4_2_8_patch7/src/HiggsAnalysis/CombinedLimit/test 
eval `scram runtime -csh`
#which root
#cd /tmp/adamwo

#set echo
if ( $#argv < 1 )  then
  exit 1
endif

set d = `echo $1 | sed 's:/::g'`
if ( !(-d $d) ) then
  echo "Directory does not exist"
  exit 1
endif


set crabs = ( `ls -d $d/crab_0_*` )
if ( $#crabs != 1 ) then
  echo "No or several crab directories"
  echo $crabs
  exit 1
endif

if ( !(-d $crabs/res) ) then
  echo "No res directory in $crabs"
  exit 1
endif

set oname = `echo $d | sed 's/job_//'`
set oname = "limits_${oname}"

pushd $crabs/res

ls -d outputToy_[0-9]*_*_*.tgz >& /dev/null
if ( $status != 0 ) then
  echo "No job output files"
  exit 1
endif

@ njobs = `ls -d outputToy_[0-9]*_*_*.tgz | cut -d'_' -f2 | sort -n | tail -1`
@ ijob = 1
while ( $ijob <= $njobs ) 

  ls outputToy_${ijob}_*_*.tgz >& /dev/null
  if ( $status == 0 ) then
    set alltars = ( outputToy_${ijob}_*_*.tgz )
    if ( $#alltars != 1 ) then
      echo "Missing or multiple job files : $alltars"
      exit 1
    endif
    set tgz = $alltars
    tar -zxf $tgz
    set odir = $tgz:r
    pushd outputToy
    set files1 = ( )
    echo "outputToy_${ijob}_*_*.tgz"
    ls
    ls higgs*.root >& /dev/null
    if ( $status == 0 ) then
      set jobfiles = ( higgs*.root )
      if ( $#jobfiles > 0 ) then
        echo "hadd ../${oname}_${odir}.root $jobfiles"
        hadd ../${oname}_${odir}.root $jobfiles
        if ( $status != 0 ) then
          exit 1
        endif
        set files1 = ( $files1 ${oname}_${odir}.root )
      endif
    endif
    popd
    rm -r outputToy/

    if ( $#files1 > 0 ) then
      set files2 = ( $files1 )
      if ( -e ${oname}.root ) then
        set files2 = ( ${oname}.root $files2 )
      endif
      hadd ${oname}_tmp.root $files2
      if ( $status != 0 ) then
        exit 1
      endif
      mv -f ${oname}_tmp.root ${oname}.root
      rm $files1
    endif
  endif

  @ ijob = $ijob + 1
end

popd
mv -i $crabs/res/${oname}.root .

exit
