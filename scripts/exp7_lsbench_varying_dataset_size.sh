#!/bin/bash

basedir=../results/exp7_lsbench_varying_dataset_size

rm -rf $basedir
mkdir $basedir

e=../symbi

size=10
for datasize in 1 5 25
do
	datadir=../datasets/lsbench_x$datasize
	querydir=../querysets/lsbench_x$datasize
	resultdir=$basedir/x$datasize
	mkdir $resultdir
	for q in {1..100}
	do
		resultfile=$resultdir/Q_$q
		timeout 2h $e $datadir/lsbench.initial $datadir/lsbench.stream.insertion $querydir/G_$size/Q_$q > $resultfile
	done
done
