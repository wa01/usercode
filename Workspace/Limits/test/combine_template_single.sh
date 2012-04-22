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
  if [ -e higgsCombineToys.*.root ]; then
    quantiles=( '--expectedFromGrid=0.16' '--expectedFromGrid=0.50' '--expectedFromGrid=0.84'  '--expectedFromGrid=0.025' '--expectedFromGrid=0.975' )
    infile=( higgsCombineToys.*.root )
    ./combine model.txt -H Asymptotic -M HybridNew --frequentist --testStat LHC -s ${seed} --grid ${infile} --readHybridResult -T 500 -i 10
    for quant in ${quantiles[@]}; do
      ./combine model.txt -H Asymptotic -M HybridNew --frequentist --testStat LHC -s ${seed} --grid ${infile} --readHybridResult -T 500 -i 10 ${quant}
    done

    ls higgsCombineTest.*.${seed}*root > /dev/null 2>&1
    if [ $? -eq 0 ]; then
      hadd tmp.root higgsCombineTest.*.${seed}*root
      rm higgsCombineTest.*.${seed}*root
      mv tmp.root higgsCombineTest.HybridNew.mH120.${seed}.root
    fi
  fi
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
