#!/bin/sh
#########################
#
# Driver script for Toy Monte Carlo submission with CRAB 
#
# author: Luca Lista, INFN
#                      
#########################

if [ -e outputToy ]; then 
  rm -rf outputToy 
fi
mkdir outputToy

i="$1"
if [ "$i" == "help" ]; then
  echo "usage: combine_crab.sh <job index> [<max events>]"
  exit 0;
fi
if [ "$i" = "" ]; then
  echo "Error: missing job index"
  exit 1;
fi
echo "max events from CRAB: $MaxEvents"
n="$MaxEvents" 
if [ "$n" = "" ]; then
  n="$2"
fi
if [ "$n" = "" ]; then
  echo "Error: missing number of experiments"
  exit 2;
fi

# first, link locally the executable:
# ln -s ../../../../bin/slc5_amd64_gcc434/combine .

set -x

n=1
while read -a line
do
if [ "$n" = "$i" ]; then break; fi
let "n += 1"
done < m0m12.lis
echo $line
echo ${#line}

n=0
until [ "${line[$n]}" = "" ]
do
  m0=${line[$n]}
  let "n += 1"
  m12=${line[$n]}

  tar -xf models.tar model_${m0}_${m12}.root
  let "seed = 10000 * ${m0} + ${m12}"
  ./combine model_${m0}_${m12}.root ${OPTIONS} -s ${seed} >& combine_${m0}_${m12}.log
  rm -f model_${m0}_${m12}.root

  let "n += 1"
done

ls
echo "job number: seed #$i with $n toys"
mv higgs*.root outputToy/
cp *.txt outputToy/
mv *.log outputToy/
echo "pack the results"
tar cvfz outputToy.tgz outputToy/
