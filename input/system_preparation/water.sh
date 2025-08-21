#!/bin/bash

pdb=$1
numres=`grep OH2 $pdb | wc -l`

x=`python -c "print($numres/9999+1)"`

echo $numres $x

echo -n > temp.pdb
for n in $(seq 1 10); do
  m=`python -c "print($n - 1)"`
  t1=`python -c "print($m*10000)"`
  t2=`python -c "print($n*10000-1)"`
  seg=`printf "%02d" $m`
  echo $n $m $t1 $t2 "$t1"-"$t2" WT$seg
  convpdb.pl -nsel 0-99999 $pdb | convpdb.pl -nsel "$t1"-"$t2" $pdb | convpdb.pl -setseg "WT$seg" | grep -v END | grep -v TER >> temp.pdb  
done
#exit 0
for n in $(seq 11 19); do
  n2=`python -c "print($n - 10)"`
  m=`python -c "print($n - 11)"`
  t1=`python -c "print($m*10000)"`
  t2=`python -c "print($n2*10000-1)"`
  m=`python -c "print($n - 1)"`
  seg=`printf "%02d" $m`
  echo $n $m $t1 $t2 "$t1"-"$t2" WT$seg

  convpdb.pl -nsel 100000-"$numres" $pdb | convpdb.pl -renumber 1 | convpdb.pl -nsel "$t1"-"$t2" | convpdb.pl -setseg "WT$seg" | grep -v END | grep -v TER >> temp.pdb 
done
echo TER >> temp.pdb
echo END >> temp.pdb

