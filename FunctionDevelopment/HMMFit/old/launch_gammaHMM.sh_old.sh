#!/bin/bash
var=$1
for i in `seq 1 $var`;
do
	sbatch ./pybatch_gammaHMM.sh $i
done
