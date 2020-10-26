#!/bin/bash

basedir=../results/exp2_lsbench_varying_query_size

rm -rf $basedir
mkdir $basedir

e=../symbi

datadir=../datasets/lsbench_x1
querydir=../querysets/lsbench_x1

# Tree queries
for size in 3 6 9 12
do
	resultdir=$basedir/T_$size
	mkdir $resultdir
	for q in {1..100}
	do
		resultfile=$resultdir/Q_$q
		timeout 2h $e $datadir/lsbench.initial $datadir/lsbench.stream.insertion $querydir/T_$size/Q_$q > $resultfile
	done
done

# Graph queries
for size in 6 9 12 10 15 20 25
do
	resultdir=$basedir/G_$size
	mkdir $resultdir
	for q in {1..100}
	do
		resultfile=$resultdir/Q_$q
		timeout 2h $e $datadir/lsbench.initial $datadir/lsbench.stream.insertion $querydir/G_$size/Q_$q > $resultfile
	done
done
		
