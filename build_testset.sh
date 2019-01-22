BASE="/proj/wallner/users/x_clami/interface_comparison/interface_comparison_nullmatching/interface_comparison/git_version/ialign_test/data/results_5_eps_0.5_dW_0.5_ch_10_d0_0.5_null_0/"
INTERFACES="/proj/wallner/users/x_clami/interface_comparison/interface_comparison_nullmatching/interface_comparison/interfaces_bio_5/"
SHELLS="/proj/wallner/users/x_clami/interface_comparison/interface_comparison_nullmatching/interface_comparison/git_version/ialign_test/set/shells/"

for file in $BASE/1*.scores_sorted_qvalues_ppv; do

	medium=`grep medium $file`
	target=`basename $file .scores_sorted_qvalues_ppv`

	cp $SHELLS/$target.shell.pdb.ca test/shells/
	
	while read line; do
		interface=`awk '{ print $1 }' <<< $line`
		cp $INTERFACES/$interface.pdb test/interfaces/
	done<<<"$medium"

done
