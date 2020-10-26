#!/bin/bash

basedir=../results/exp5_netflow_varying_insertion_rate

rm -rf $basedir
mkdir $basedir

e=../symbi

datadir=../datasets/netflow
querydir=../querysets/netflow

size=10
for rate in 2 4 6 8 10
do
	if [ "${rate}" -eq 2 ]
	then
		numq=370415
	elif [ "${rate}" -eq 4 ]
	then
		numq=740830
	elif [ "${rate}" -eq 6 ]
	then
		numq=1111245
	elif [ "${rate}" -eq 8 ]
	then
		numq=1481660
	else
		numq=1852076
	fi
	resultdir=$basedir/insertion_$rate
	mkdir $resultdir
   	for q in {1..100}
   	do
       	resultfile=$resultdir/Q_$q
       	timeout 2h $e $datadir/netflow.initial $datadir/netflow.stream.insertion $querydir/G_$size/Q_$q $numq > $resultfile	
   	done
done
