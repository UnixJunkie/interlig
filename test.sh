#!/bin/bash

echo "Testing PDB alignment..."

./InterfaceComparison -PDB examples/target1.pdb examples/target2.pdb > output_example_target1_target2_pdb
./InterfaceComparison -list examples/target1.pdb examples/pdb_list > output_example_target1_list_pdb

dif1=`diff output_example_target1_target2_pdb examples/output_example_target1_target2_pdb`
dif2=`diff output_example_target1_list_pdb examples/output_example_target1_list_pdb`

if [ -n "$dif1" ] || [ -n "$dif2" ]; then

	echo " x Error! PDB alignment test failed."

else

	echo " * PDB alignment test successful."

fi

echo "Testing PDB superposition..."

./InterfaceComparison -PDB examples/target1.pdb examples/target2.pdb -super -apply examples/1mrw.pdb > /dev/null

dif1=`diff target1.pdb_target2.pdb examples/target1.pdb_target2.pdb`
dif2=`diff target2.pdb_target1.pdb examples/target2.pdb_target1.pdb`
dif3=`diff 1mrw.pdb_target1.pdb examples/1mrw.pdb_target1.pdb`

if [ -n "$dif1" ] || [ -n "$dif2" ] || [ -n "$dif3" ]; then

        echo " x Error! PDB superposition test failed."

else

        echo " * PDB superposition test successful."

fi

echo "Testing mol2 alignments..."

./InterfaceComparison -mol2 examples/target.mol2 examples/database.mol2 > output_example_target_database_mol2

dif=`diff output_example_target_database_mol2 examples/output_example_target_database_mol2`

if [ -n "$dif" ]; then

        echo " x Error! mol2 alignment test failed."

else

        echo " * mol2 alignment test successful."

fi

rm output_example_target1_target2_pdb output_example_target1_list_pdb output_example_target_database_mol2 target1.pdb_target2.pdb target2.pdb_target1.pdb 1mrw.pdb_target1.pdb

echo "Done."

