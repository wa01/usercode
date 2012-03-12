#! /bin/tcsh -f
#
#
# extract result lines from log files in hnGrid*.tgz
#
set nlog = "exp50.log"

foreach f ( $* )

 if ( !( -e $1 ) ) then
   exit 1
 endif

 set srcDir = $PWD

 set fn = $f:t:r
 set parts = ( `echo $fn | sed 's/_/ /g'` )
 @ ht = `echo $parts[3] | sed 's/ht//'`
 @ met = `echo $parts[4] | sed 's/met//'`

 set tmpDir = /tmp/adamwo/$fn
 if ( -d $tmpDir )  rm -r $tmpDir
 mkdir $tmpDir
 pushd $tmpDir >& /dev/null
 tar -zxf $srcDir/$f
 set lim = ( `awk 'BEGIN{l=-1;el=-1;}/^Done in/{print l,el}/^Limit:/{l=$4;el=$6}' ${nlog}` )
 if ( $#lim == 2 ) then
   printf "%-6s %5d %4d %s %s %s\n" $parts[2] $ht $met $lim[1] "+-" $lim[2]
 else
   printf "%-6s %5d %4d %s\n" $parts[2] $ht $met "**** missing ****"
 endif

 popd >& /dev/null
 rm -r $tmpDir
end

