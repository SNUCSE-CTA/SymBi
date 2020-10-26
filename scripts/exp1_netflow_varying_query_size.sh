#!/bin/bash

basedir=../results/exp1_netflow_varying_query_size

rm -rf $basedir
mkdir $basedir

e=../symbi

datadir=../datasets/netflow
querydir=../querysets/netflow

# Tree queries
for size in 3 6 9 12
do
	resultdir=$basedir/T_$size
	mkdir $resultdir
	for q in {1..100}
	do
		resultfile=$resultdir/Q_$q
		timeout 2h $e $datadir/netflow.initial $datadir/netflow.stream.insertion $querydir/T_$size/Q_$q > $resultfile
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
		timeout 2h $e $datadir/netflow.initial $datadir/netflow.stream.insertion $querydir/G_$size/Q_$q > $resultfile
	done
done
		
