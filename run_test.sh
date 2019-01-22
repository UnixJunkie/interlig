#!/bin/bash

T0=$1
ch=$2
seed=$3

date
for target in test/shells/*; do
	for interface in test/interfaces/*; do

	res=`./InterfaceComparison -PDB $target $interface -ch $ch -T0 $T0 -seed $seed`

	target_base=`basename $target .pdb.ca`
	interf_base=`basename $interface .pdb`
	echo "$target_base $interf_base $res"

	done
done
date
