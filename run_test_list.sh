#!/bin/bash

echo "Testing file list"

date
./InterfaceComparison2 -list $1 ./list_of_interfaces -ch 10 -T0 0.2 -seed 1 > test_list

date
echo "Testing file by file now"

date
for file in test/interfaces/*; do

	./InterfaceComparison2 -PDB $1 $file -ch 10 -T0 0.2 -seed 1 #-anneal 3

done > test_one_by_one

date
