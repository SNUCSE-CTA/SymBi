#!/bin/bash

basedir=../results/exp4_lsbench_varying_deletion_rate

rm -rf $basedir
mkdir $basedir

e=../symbi

datadir=../datasets/lsbench_x1
querydir=../querysets/lsbench_x1

size=10
for rate in 2 4 6 8 10
do
	resultdir=$basedir/deletion_$rate
	mkdir $resultdir
	for q in {1..100}
	do
		resultfile=$resultdir/Q_$q
		timeout 2h $e $datadir/lsbench.initial $datadir/lsbench.stream.deletion.$rate $querydir/G_$size/Q_$q > $resultfile
	done
done
