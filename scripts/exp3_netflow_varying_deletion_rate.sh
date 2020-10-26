#!/bin/bash

basedir=../results/exp3_netflow_varying_deletion_rate

rm -rf $basedir
mkdir $basedir

e=../symbi

datadir=../datasets/netflow
querydir=../querysets/netflow

size=10
for rate in 2 4 6 8 10
do
	resultdir=$basedir/deletion_$rate
	mkdir $resultdir
	for q in {1..100}
	do
		resultfile=$resultdir/Q_$q
		timeout 2h $e $datadir/netflow.initial $datadir/netflow.stream.deletion.$rate $querydir/G_$size/Q_$q > $resultfile
	done
done
