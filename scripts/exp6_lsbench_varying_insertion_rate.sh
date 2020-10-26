#!/bin/bash

basedir=../results/exp6_lsbench_varying_insertion_rate

rm -rf $basedir
mkdir $basedir

e=../symbi

datadir=../datasets/lsbench_x1
querydir=../querysets/lsbench_x1

size=10
for rate in 2 4 6 8 10
do
	if [ "${rate}" -eq 2 ]
	then
		numq=466291
	elif [ "${rate}" -eq 4 ]
	then
		numq=932582
	elif [ "${rate}" -eq 6 ]
	then
		numq=1398874
	elif [ "${rate}" -eq 8 ]
	then
		numq=1865165
	else
		numq=2331457
	fi
	resultdir=$basedir/insertion_$rate
	mkdir $resultdir
   	for q in {1..100}
   	do
       	resultfile=$resultdir/Q_$q
       	timeout 2h $e $datadir/lsbench.initial $datadir/lsbench.stream.insertion $querydir/G_$size/Q_$q $numq > $resultfile	
   	done
done
