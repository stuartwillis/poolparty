#!/bin/bash

#This script plots results from PPanalyze, runs the FLK test (midpoint root; fails with <3 samples), and optionally runs the CMH test

if [ -z "$5" ]
then
echo "Not enough arguments."
echo "Suggested usage: bash SCRIPT.sh analysis_prefix CMH?[yes/no] FLK?[yes/no] PPmanhat_path Popoolation2_path {SCAFF_PREFIX}"
echo "command line arguments in {curly brackets} are optional"
echo "example: bash SCRIPT.sh phenoAB_analyze no yes ~/poolparty-master/ ~/bin/popoolation2_1201 JAAXML"
exit 1
#else
fi

MHPATH=$4
POPOOL=$5
SCAFF=$6

echo "Regarding CMH test you have selected: "$2
echo "Regarding FLK test you have selected: "$3
echo "The path to the PPmanhat script is: "$MHPATH
echo "The path to the Popoolation2 folder is: "$POPOOL
echo "Scaffold designation is: "$SCAFF

if  ( [[ $2 =~(yes)$ ]] )  ; then

	echo -e "\nPlease provide the population comparisons you wish to make and press [ENTER]"
	echo -e "Make sure this is consistent with the order of populations as indicated by PPanalyze"
	echo -e "PPanalyze reorders populations in the sync file, ordering by phenotype (e.g. all A then all B)"
	echo -e "So, 1:3,2:4 where 1 and 2 are from Pop 1 and 1 and 3 have phenotype A,"
	echo -e "Results in PPanalyze indicating \"Pops 1 3 2 4 are now in the order of 1 2 3 4 in the subset sync file\""
	echo -e "So, to look for differences among phenotypes consistent across populations"
	echo -e "the specification for the CMH test should be:"
	echo -e "e.g. 1-3,2-4 where 1&3 IN THE ANALYZE SUBSET SYNC WERE 1&2 BEFORE and are from the same locale "
	echo -e "OR E.g. \"1-3\" or \"1-2,3-4\" or \"1-4,2-5,3-6\" etc.\n"

	read Comparison

	if [ "$Comparison" == "" ]; then
        echo "Ok, I guess not."
        exit 1
	fi

fi

echo -e "\nNow running SNP density plot"
bash $MHPATH/PPmanhat.sh -i $1_density.txt -o Density_plot -a SNP_Density_10Kb -s $SCAFF -t line

echo -e "\nNow running Fst plot"
bash $MHPATH/PPmanhat.sh -i $1.fst -o FST_plot -a FST -s $SCAFF -L 0.5

echo -e "\nNow running sliding window Fst plot"
bash $MHPATH/PPmanhat.sh -i $1.Sfst -o SFst_plot -a SFst -s $SCAFF

echo -e "\nNow running Fisher's exact test plot"
bash $MHPATH/PPmanhat.sh -i $1.fet -o FET_plot -a -log10p -s $SCAFF

echo -e "\nNow running local score of FET and plot"
bash $MHPATH/PPrunls.sh -i $1.fet -o LS_Fet_run -s $SCAFF
sig=$(cut -d" " -f3 LS_Fet_run_mean_sig.txt | tail -n1)
bash $MHPATH/PPmanhat.sh -i LS_Fet_run.ls -o LS_Fet_plot -a -Local_score -s $SCAFF -L $sig

if  ( [[ $3 =~(yes)$ ]] )  ; then

	echo -e "\nNow running FLK test and plot"
	bash $MHPATH/PPrunflk.sh -i $1.fz -o FLK_run
	bash $MHPATH/PPmanhat.sh -i FLK_run.flk -o FLK_plot -a -log10p -s $SCAFF -l TRUE

fi

if  ( [[ $2 =~(yes)$ ]] )  ; then

	echo -e "\nNow running CMH test and plot"
	bash $MHPATH/PPruncmh.sh -i $1.sync -o CMH_out -p $Comparison -P $POPOOL
	bash $MHPATH/PPmanhat.sh -i CMH_out.cmh -o CMH_plot -a -log10p -s $SCAFF -l TRUE

	echo -e "\nNow running local score of CMH and plot"
	bash $MHPATH/PPrunls.sh -i CMH_out.cmh -o LS_CMH_run -s $SCAFF
	sig=$(cut -d" " -f3 LS_CMH_run_mean_sig.txt | tail -n1)
	bash $MHPATH/PPmanhat.sh -i LS_CMH_run.ls -o LS_CMH_plot -a -Local_score -s $SCAFF -L $sig

fi


