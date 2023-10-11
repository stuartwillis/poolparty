#!/bin/bash

#PoolParty v04.01.2022
#PPanalyze

source $1
BASEDIR=$(dirname "$0")

if [[ ${BASEDIR} != *poolparty* ]];then
	SLOC=$(readlink -f "$0")
	BASEDIR=${SLOC%/*}

set -o pipefail

fi

if  ( [[ $(echo $1)  = "" ]] )  ; then
	echo "ERROR: You must provide a config file after the PPanalyze command"
	echo "for example, PPanalyze example_analyze.config"
	exit
fi


#Declare integers for numerical parameters
declare -i MINCOV
declare -i MAXCOV

#Get date of run and save an additional config file with run parameters
RUNDATE=$(date +'%m_%d_%Y')

echo "ALERT: Beginning PoolParty Analyze at $(date)"

#make output dir
	if [ ! -d "$OUTDIR/" ]; then
		mkdir $OUTDIR/
	fi

#make temporary file location
	if [ ! -d "${OUTDIR}/temp" ]; then
		mkdir ${OUTDIR}/temp
	fi

echo "ALERT: Performing error checking"
	#Check for required files, quit if it one does not exist
		if [[ ! -f $COVFILE ]] ; then
			 echo 'ERROR: COVFILE file is missing, aborting.'
			 exit
		fi
		if [[ ! -f $SYNC ]] ; then
			 echo 'ERROR: SYNC file is missing, aborting. (Are you in the correct working directory??)'
			 exit
		fi
		if [[ ! -f $FZFILE ]] ; then
			 echo 'ERROR: FZFILE file is missing, aborting.'
			 exit
		fi
		if [ ! -z $BLACKLIST ]; then
			if [[ ! -f $BLACKLIST ]] ; then
			 	echo 'ERROR: BLACKLIST file is missing, aborting.'
				exit
			else
				echo "ALERT: Blacklist has been specified: " $BLACKLIST
				#echo "$BLACKLIST"
			fi
		fi

	#Ensure poolparty directory is intact
		if [ ! -d "$BASEDIR/rscripts" ] ; then
			echo "ERROR: Poolparty directory has been tampered with"
			exit
		fi
	#Ensure base directories all exist
		if  [ ! -d "$POOL2" ] ; then
			echo "ERROR: Wrong Directories: POOL2 does not exist"
			exit
		fi
	#Ensure required parameters are within a valid range
		if  [[ ! "$MINCOV" =~ ^[0-9]+$ ]] || [[ ! "$MAXCOV" =~ ^[0-9]+$ ]];  then
        	echo "ERROR: MINCOV and MAXCOV must be positive integers. Check parameters"
			echo "MINCOV currently read as $MINCOV, and MAXCOV currently read as $MAXCOV"
			exit
		fi
		if  [[ ! "$MAF" =~ ^[+-]?[0-9]+\.?[0-9]*$  ]] ; then 
			echo "ERROR: MAF must be a value between 0 and 1"
			echo "MAF currently read as $MAF"
			exit
		fi
	#Check for analyses turned on or off
		if [[ "$FST" != "off" ]] && [[ "$FST" != "on" ]] ; then
			echo "ERROR: Wrong Parameters: FST must be either on or off"
			exit
		fi
		if [[ "$SLIDINGFST" != "off" ]] && [[ "$SLIDINGFST" != "on" ]] ; then
			echo "ERROR: Wrong Parameters: SLIDINGFST must be either on or off"
			exit
		fi
		if [[ "$FET" != "off" ]] && [[ "$FET" != "on" ]] ; then
			echo "ERROR: Wrong Parameters: FET must be either on or off"
			exit
		fi
		if [[ "$NJTREE" != "off" ]] && [[ "$NJTREE" != "on" ]] ; then
			echo "ERROR: Wrong Parameters: NJTREE must be either on or off"
			exit
		fi
	#Check for analysis-specific paramters
		if [[ "$FST" == "on" ]] ; then
			declare -i NIND
			if [[ "$FSTTYPE" != "traditional" ]] && [[ "$FSTTYPE" != "karlsson" ]] ; then
				echo "ERROR: Wrong Parameters: FST type must be either traditional or karlsson"
				echo "FSTTYPE currently read as $FSTTYPE"
				exit
			fi
			if  [[ ! "$NIND" =~ ^[0-9]+$ ]] ; then
				echo "ERROR: Wrong Parameters: NIND must be a positive integer"
				echo "NIND currently read as $NIND"
				exit
			fi
		fi
		if [[ "$SLIDINGFST" == "on" ]]  ; then
			declare -i WINDOW
			declare -i STEP
			declare -i NIND
			if [[ "$FSTTYPE" != "traditional" ]] && [[ "$FSTTYPE" != "karlsson" ]] ; then
				echo "ERROR: Wrong Parameters: FST type must be either traditional or karlsson"
				echo "FSTTYPE currently read as $FSTTYPE"
				exit
			fi
			if  [[ ! "$NIND" =~ ^[0-9]+$ ]]  ; then
				echo "ERROR: Wrong Parameters: NIND must be a positive integer"
				echo "NIND currently read as $NIND"
				exit
			fi
			if [[ "$SLIDINGFST" == "on" ]]  ; then
				if [[ ! "$WINDOW" =~ ^[0-9]+$ ]] || [[ ! "$STEP" =~ ^[0-9]+$ ]] ; then
					echo "ERROR: Wrong Parameters: WINDOW AND STEP must be positive integers"
					echo "WINDOW currently read as $WINDOW, STEP currently read as $STEP"
					exit
				fi
			fi
		fi
		if [[ "$NJTREE" == "on" ]]  ; then
			declare -i STRWINDOW
			declare -i BSTRAP
			declare -i AFFILT
			if [[ "$METHOD" != "mean" ]] &&  [[ "$METHOD" != "random" ]] \
				&& [[ "$METHOD" != "rangemax" ]] && [[ "$METHOD" != "rangemin" ]] \
				&& [[ "$METHOD" != "first" ]] && [[ "$METHOD" != "last" ]]	; then
				echo "ERROR: Wrong Parameters: NJ Combine method type must be mean, random, rangemax, rangemin, first, or last"
				echo "METHODS currently read as $METHOD"
				exit
			fi
			if  [[ ! "$AFFILT" =~ ^[+-]?[0-9]+\.?[0-9]*$  ]]  ; then
				echo "ERROR: Wrong Parameters: AFFILT must be between 0 and 1"
				echo "AFFILT currently read as $AFFILT "
				exit
			fi
			if  [[ ! "$STRWINDOW" =~ ^[0-9]+$ ]]  ; then
				echo "ERROR: Wrong Parameters: STRWINDOW must be a positive integer"
				echo "STRWINDOW currently read as $STRWINDOW"
				exit
			fi
			if  [[ ! "$BSTRAP" =~ ^[0-9]+$ ]]  ; then
				echo "ERROR: Wrong Parameters: BSTRAP must be a positive integer"
				echo "BSTRAP currently read as $BSTRAP"
				exit
			fi
		fi
		#Check for R and PERL
		if  [[ $(command -v Rscript)  = "" ]] ; then
			echo "ERROR: R not detected on system. R is required for analysis "
			echo "R should initiate when 'Rscript' is typed in terminal"
			exit
		fi
		if  [[ $(command -v perl)  = "" ]] ; then
			echo "ERROR: perl not detected on system. perl is required for Popoolation2"
			echo "perl should initiate when 'perl' is typed in terminal"
			exit
		fi
		#Load R to check for dependencies
		if  [[ -f $OUTDIR/R_ERROR.txt ]] ; then
			rm $OUTDIR/R_ERROR.txt 
		fi	
		
		Rscript $BASEDIR/rscripts/r_analyze_check.R $OUTDIR
		
		if  [[ -f $OUTDIR/R_ERROR.txt ]] ; then
			echo "ERROR: R dependency check failed, install dependencies manually"
			exit
		fi

#Copy configuration file so you know what you did later

	cp $1 ${OUTDIR}/${PREFIX}_${RUNDATE}.config | awk '{ sub("\r$", ""); print $0 }'
	
#Get number of populations specified by POPS 
	POPS2=$(echo $POPS | sed 's/[,:]/ /g')
		declare -a arr2=($POPS2)
		declare -i NPOPS=${#arr2[@]}

	echo "ALERT: There are $NPOPS populations specified for this analysis"

	if [ -z "$(echo $POPS | sed -n 's/\([,:]\)/\1/p')" ];
		then echo "ERROR: POPS specification contains invalid characters"
		exit
	fi
#Create matching sequence for column name adjustment
	pseq=$(printf '%0.s3 ' $(seq 1 $NPOPS))

#Prepare expression for keeping columns
		symb=$(printf '%s\n' "${POPS//[[:digit:]]/}")
		#echo $symb
		#echo "ZERO"
		seq 1 $NPOPS > ${OUTDIR}/temp/${PREFIX}_seq.txt
		newORD=$(awk '{print}' ORS=' ' ${OUTDIR}/temp/${PREFIX}_seq.txt)
		#echo "ONE"
		echo -n $symb | sed 's/./&\n/g' > ${OUTDIR}/temp/${PREFIX}_char.txt
		#echo "TWO"
		paste -d \\n  ${OUTDIR}/temp/${PREFIX}_seq.txt ${OUTDIR}/temp/${PREFIX}_char.txt > ${OUTDIR}/temp/${PREFIX}_charseq.txt
		#echo "THREE"
		stored=$(cat ${OUTDIR}/temp/${PREFIX}_charseq.txt | tr -d '\r\n')
		#echo "FOUR"
		#sed 's/\S*\(:\)\S*//g' <<< $stored
		#echo "FIVE"
		echo "$stored" | sed -r 's/[,]+/ /g' > ${OUTDIR}/temp/${PREFIX}_charseq2.txt
		#echo "SIX"
		COMBLIST=$(cat ${OUTDIR}/temp/${PREFIX}_charseq2.txt | tr ' ' '\n' | grep ":")
		#split each term in comblist into its own sequence
		#echo "SEVEN"
		echo $COMBLIST | tr ' ' '\n' | grep ":" | awk  '{gsub(":"," ",$0); print;}'  > ${OUTDIR}/temp/${PREFIX}_expressionlist.txt
		#echo "EIGHT"
		echo $POPS2 > ${OUTDIR}/temp/${PREFIX}_popnames.txt
		#echo "NINE"
		
	echo "ALERT: Pops $POPS2 are now in the order of $newORD in the subset sync file"
	
		if  [[ -f ${OUTDIR}/${PREFIX}_avoidcols.txt ]] ; then
			rm ${OUTDIR}/${PREFIX}_avoidcols.txt
		fi	
	
		#Get names of columns to ignore by finding all possible combinations of similar populations
		while read p; do
			set $p
			for a; do
				shift
			for b; do
				printf "%s:%s\n" "$a" "$b" >> ${OUTDIR}/${PREFIX}_avoidcols.txt
			done
			done
		done <${OUTDIR}/temp/${PREFIX}_expressionlist.txt

if [[ -f $OUTDIR/${PREFIX}.sync ]] ; then
	echo "ALERT: ${PREFIX}.sync exists; will not subset file"
	echo "ALERT: Grabbing information from .sync file"
	# Freq file, run through R to get final MAF for populations being compared(black list)
	TSO3=$(date +%s | sha256sum | base64 | head -c 8 ; echo)
	declare -a fzarray=()
		for i in "${arr2[@]}" ; do
			fzarray+=( "$(( 4 + $i ))" )
		done
			bar3=$(printf "\n%s\n" "${fzarray[@]}")
			ONE=$(echo $bar3| sed -r 's/([^ ]+)/$\1/g')
					TWO=$(echo $ONE | awk '$NF=$NF "}" ' OFS=",")
					THREE='{print $1,$2,$3,'
					FOUR=$(echo ${THREE}${TWO})	
			awk "$FOUR" $FZFILE > $OUTDIR/temp/${PREFIX}_$TSO3
		
	#Call Rscript
		echo "ALERT: Calculating comparison-specific MAF at $(date)"
				rin=$OUTDIR/temp/${PREFIX}_$TSO3
				rout=$OUTDIR/temp/
				Rscript $BASEDIR/rscripts/r_anal_maf.R  $rin $rout $MAF
else
	#Create array adding column offset to each to filter sync file by specified populations
	TSO1=$(date +%s | sha256sum | base64 | head -c 10 ; echo)
	declare -a colarray=()
		for i in "${arr2[@]}" ; do
			colarray+=( "$(( 3 + $i ))" )
		done
			bar=$(printf "\n%s\n" "${colarray[@]}")
					ONE=$(echo $bar| sed -r 's/([^ ]+)/$\1/g')
					TWO=$(echo $ONE | awk '$NF=$NF "}" ' OFS=",")
					THREE='{print $1,$2,$3,'
					FOUR=$(echo ${THREE}${TWO})	
			awk "$FOUR" ${SYNC} > $OUTDIR/temp/${PREFIX}_$TSO1
			#subset sync created

	# Get coverage whitelist using coverage range for all files specified 
	TSO2=$(date +%s | sha256sum | base64 | head -c 9 ; echo)
	declare -a covarray=()
		for i in "${arr2[@]}" ; do
			covarray+=( "$(( 2 + $i ))" )  
		done
			bar2=$(printf "\n%s\n" "${covarray[@]}")
				ONE=$(echo $bar2| sed -r 's/([^ ]+)/$\1/g')
				TWO=$(echo $ONE | awk '$NF=$NF " >= '$MINCOV' && "' OFS=" >= $MINCOV && ")
				THREE=$(echo $ONE | awk '$NF=$NF " <= '$MAXCOV'"' OFS=" <= $MAXCOV && ")
				FOUR=$(echo $TWO$THREE)
			tail -n +2 $FZFILE | awk  "$FOUR" $COVFILE | awk '{print $1,$2}' > $OUTDIR/temp/${PREFIX}_$TSO2

			#subset coverage sync file created
					FLEN=$(wc -l $FZFILE )
 					echo "ALERT: There were $FLEN SNPs before filters"
					FLEN=$(wc -l < $OUTDIR/temp/${PREFIX}_$TSO2 )
 					echo "ALERT: There are $FLEN SNPs above coverage filters"

	# Get indel/N blacklist from coverage file for all files specified 
	TSO2B=$(date +%s | sha256sum | base64 | head -c 17 ; echo)
	declare -a covarray=()
                declare -i  size=$(awk '{print NF}' $SYNC| sort -nu | head -n 1)
	
		for i in "${arr2[@]}" ; do
	
			covarray+=( "$(( $size + $i ))" )
		done
			bar2=$(printf "\n%s\n" "${covarray[@]}")
				ONE=$(echo $bar2| sed -r 's/([^ ]+)/$\1/g')
				TWO=$(echo $ONE | awk '$NF=$NF " >= 1 || "' OFS=" >= 1 || ")
				THREE=$(echo $ONE | awk '$NF=$NF " > 0"' OFS=" > 0 || ")
				FOUR=$(echo $TWO$THREE)
			tail -n +2 $FZFILE | awk  "$FOUR" $COVFILE | awk '{print $1,$2}' > $OUTDIR/temp/${PREFIX}_$TSO2B
			#subset coverage sync file created
					FLEN=$(wc -l < $OUTDIR/temp/${PREFIX}_$TSO2B )
 					echo "ALERT: There are $FLEN SNPs marked by N/indel filters"

			
	# Freq file, run through R to get final MAF for populations being compared(black list)
	TSO3=$(date +%s | sha256sum | base64 | head -c 8 ; echo)
	declare -a fzarray=()
		for i in "${arr2[@]}" ; do
			fzarray+=( "$(( 4 + $i ))" )
		done
			bar3=$(printf "\n%s\n" "${fzarray[@]}")
			ONE=$(echo $bar3| sed -r 's/([^ ]+)/$\1/g')
					TWO=$(echo $ONE | awk '$NF=$NF "}" ' OFS=",")
					THREE='{print $1,$2,$3,'
					FOUR=$(echo ${THREE}${TWO})	
			awk "$FOUR" $FZFILE > $OUTDIR/temp/${PREFIX}_$TSO3
		
	#Call Rscript; creates a MAF black list or "mafBL" text file
		echo "ALERT: Calculating comparison-specific MAF at $(date)"
				rin=$OUTDIR/temp/${PREFIX}_$TSO3
				rout=$OUTDIR/temp/
				Rscript $BASEDIR/rscripts/r_anal_maf.R  $rin $rout $MAF

				MLEN=$(wc -l < $OUTDIR/temp/${PREFIX}_${TSO3}_mafBL )
				echo "ALERT: There are $MLEN SNPs marked for removal from MAF filters"

	#First Filter by coverage whitelist and remove by blacklist
	TSO4=$(date +%s | sha256sum | base64 | head -c 5 ; echo)
		echo "ALERT: Filtering blacklisted markers at $(date)"
		gawk 'NR==FNR{a[$1,$2]=$3;next} ($1,$2) in a{print $0, a[$1,$2]}' $OUTDIR/temp/${PREFIX}_$TSO2 $OUTDIR/temp/${PREFIX}_$TSO1 > $OUTDIR/temp/${PREFIX}_$TSO4
		awk 'NR==FNR{a[$1,$2]=$3;next} ($1,$2) in a{next}{print $0, a[$1,$2]}'  $OUTDIR/temp/${PREFIX}_$TSO2B $OUTDIR/temp/${PREFIX}_$TSO4 > $OUTDIR/temp/${PREFIX}_${TSO4}_B
					FLEN=$(wc -l < $OUTDIR/temp/${PREFIX}_${TSO4}_B )
 					echo "ALERT: There are $FLEN SNPs after coverage and N/indel filters"
		rm $OUTDIR/temp/${PREFIX}_$TSO2 &> /dev/null
		rm $OUTDIR/temp/${PREFIX}_$TSO1 &> /dev/null
		rm $OUTDIR/temp/${PREFIX}_$TSO2B &> /dev/null
			
	#Combine blacklist info into one temp file
	echo "The following has been specified to identify blacklisted positions: $BLACKLIST"
	TSO5=$(date +%s | sha256sum | base64 | head -c 12 ; echo)
					BLEN=$(wc -l < "$BLACKLIST" )
					echo "ALERT: There are $BLEN SNPs marked for removal from user-supplied blacklist"
					cat $OUTDIR/temp/${PREFIX}_${TSO3}_mafBL $BLACKLIST | sort | uniq > $OUTDIR/temp/${PREFIX}_$TSO5
					BLEN=$(wc -l < $OUTDIR/temp/${PREFIX}_$TSO5 )
					echo "ALERT: There are $BLEN SNPs marked for removal from blacklist and MAF filters"
					gawk 'NR==FNR{a[$1,$2]=$3;next} ($1,$2) in a{next}{print $0, a[$1,$2]}' $OUTDIR/temp/${PREFIX}_$TSO5  $OUTDIR/temp/${PREFIX}_${TSO4}_B |  awk  '{gsub(" ","\t",$0); print;}' > $OUTDIR/${PREFIX}.sync
				rm $OUTDIR/temp/${PREFIX}_$TSO4  &> /dev/null
				rm $OUTDIR/temp/${PREFIX}_$TSO4_B &> /dev/null
				rm $OUTDIR/temp/${PREFIX}_$TSO5 &> /dev/null
fi


FLEN=$(wc -l < $OUTDIR/${PREFIX}.sync )
echo "ALERT: There are $FLEN SNPs being analyzed after coverage, N/indel, MAF, and user-blacklist filters"

if [ "$FLEN" -eq 0 ] ; then
	echo "There are no SNPs left to analyze; something failed!"
	exit 1
fi

#Run analyses 

echo "ALERT: Running pairwise analyses at $(date)"

if [[ "$NJTREE" =~(on)$ ]] ; then 
	if [[ -f $OUTDIR/${PREFIX}.fz ]] ; then
		echo "ALERT: $OUTDIR/${PREFIX}.fz exists; will not redo FZ subsetting or NJ trees"
	else
		echo "ALERT: Running R for structure analyses $(date)"
		gawk 'NR==FNR{a[$1,$2]=$3;next} ($1,$2) in a{print $0, a[$1,$2]}' <(awk '{print $1,$2}' $OUTDIR/${PREFIX}.sync)  <( awk '{print $0}' $OUTDIR/temp/${PREFIX}_$TSO3)  >  $OUTDIR/${PREFIX}.fz
		#Call Rscript
		rin=$OUTDIR/${PREFIX}.fz
		rout=$OUTDIR/
		cfile=${OUTDIR}/temp/${PREFIX}_popnames.txt
		Rscript $BASEDIR/rscripts/r_structure.R $rin $rout $AFFILT $STRWINDOW $METHOD $BSTRAP $cfile
	fi
fi

	#FST, either traditional or karlsson
if [[ "$FST" =~(on)$ ]] && [[ "$FSTTYPE" =~(traditional)$ ]]; then 
	if [[ -f $OUTDIR/${PREFIX}_raw.fst ]] ; then
		echo "ALERT: $OUTDIR/${PREFIX}_raw.fst exists; will not redo FST calculation"
	else
		echo "ALERT: FST is running"
		perl /${POOL2}/fst-sliding.pl --input $OUTDIR/${PREFIX}.sync --output  $OUTDIR/${PREFIX}_raw.fst  --min-count 1 \
			--min-coverage 2 --max-coverage ${MAXCOV} --min-covered-fraction 1 --window-size 1 --step-size 1 --pool-size ${NIND} &
	fi
fi

if [[ "$FST" =~(on)$ ]] && [[ "$FSTTYPE" =~(karlsson)$ ]]; then 
	if [[ -f $OUTDIR/${PREFIX}_raw.fst ]] ; then
		echo "ALERT: $OUTDIR/${PREFIX}_raw.fst exists; will not redo FST calculation"
	else
		echo "ALERT: FST is running"
		perl /${POOL2}/fst-sliding.pl --input $OUTDIR/${PREFIX}.sync --output  $OUTDIR/${PREFIX}_raw.fst  --min-count 1 \
			--min-coverage 2 --max-coverage ${MAXCOV} --min-covered-fraction 1 --window-size 1 --step-size 1 --pool-size ${NIND} --karlsson-fst &
	fi
fi	

	#SFST, either traditional or karlsson
if [[ "$SLIDINGFST" =~(on)$ ]] && [[ "$FSTTYPE" =~(traditional)$ ]] ; then 
	if [[ -f $OUTDIR/${PREFIX}_raw.Sfst ]] ; then
		echo "ALERT: $OUTDIR/${PREFIX}_raw.Sfst exists; will not redo sliding FST calculation"
	else
		echo "ALERT: Sliding FST is running"
		perl /${POOL2}/fst-sliding.pl --input $OUTDIR/${PREFIX}.sync --output $OUTDIR/${PREFIX}_raw.Sfst  --min-count 1 \
			--min-coverage 2 --max-coverage ${MAXCOV} --min-covered-fraction 0 --window-size ${WINDOW} --step-size ${STEP} --pool-size ${NIND}  &
	fi
fi
	
if [[ "$SLIDINGFST" =~(on)$ ]] && [[ "$FSTTYPE" =~(karlsson)$ ]] ; then 
	if [[ -f $OUTDIR/${PREFIX}_raw.Sfst ]] ; then
		echo "ALERT: $OUTDIR/${PREFIX}_raw.Sfst exists; will not redo sliding FST calculation"
	else
		echo "ALERT: Sliding FST is running"
		perl /${POOL2}/fst-sliding.pl --input $OUTDIR/${PREFIX}.sync --output $OUTDIR/${PREFIX}_raw.Sfst  --min-count 1 \
			--min-coverage 2 --max-coverage ${MAXCOV} --min-covered-fraction 0 --window-size ${WINDOW} --step-size ${STEP} --karlsson-fst --pool-size ${NIND}  &
	fi
fi

if [[ "$FET" =~(on)$ ]] ; then 
	if [[ -f $OUTDIR/${PREFIX}_raw.fet ]] ; then
		echo "ALERT: $OUTDIR/${PREFIX}_raw.fet exists; will not redo FET calculation"
	else
		echo "ALERT: FET is running"
		perl /${POOL2}/fisher-test.pl --input $OUTDIR/${PREFIX}.sync --output $OUTDIR/${PREFIX}_raw.fet  \
			--min-count 1 --min-coverage 2 --max-coverage ${MAXCOV}
	fi
fi

wait
echo "ALERT: All pairwise analyses complete at $(date)"


#Chop up results
echo "ALERT: Formatting results"


if [[ "$FST" =~(on)$ ]]; then
#Cleanup FST OUTDIR format, average rows, and produce intuitive format
	if [[ -f $OUTDIR/${PREFIX}.fst ]] ; then
		echo "ALERT: $OUTDIR/${PREFIX}.fst exists; will not redo FST reformatting"
	else
		echo "ALERT: now doing FST reformatting"
		awk '$5 >= '$MINCOV' &&  $5  <= '$MAXCOV' ' $OUTDIR/${PREFIX}_raw.fst >  $OUTDIR/temp/${PREFIX}_reclassing.fst

		#Cut up Sfst file
		cut -d$'\t' -f 1-2 $OUTDIR/temp/${PREFIX}_reclassing.fst | gawk '$3=(FNR FS $3)' > $OUTDIR/temp/${PREFIX}_head.fst
		cut -d$'\t' -f 6- $OUTDIR/temp/${PREFIX}_reclassing.fst > $OUTDIR/temp/${PREFIX}_body.fst
		awk 'NR==1{print $0}' $OUTDIR/temp/${PREFIX}_body.fst | awk  '{gsub("=","\t",$0); print;}' |  awk '{ for (i=2;i<=NF;i+=2) $i="" } 1' | awk  '{gsub("  ","\t",$0); print;}' | awk  '{gsub(" ","",$0); print;}'   > $OUTDIR/temp/${PREFIX}_body_heading.fst
		awk '{gsub("=","\t",$0); print;}' $OUTDIR/temp/${PREFIX}_body.fst | awk '{ for (i=1;i<=NF;i+=2) $i="" } 1' | awk  '{gsub("  ","\t",$0); print;}' | awk  '{gsub(" ","",$0); print;}'  > $OUTDIR/temp/${PREFIX}_body2.fst
		cat $OUTDIR/temp/${PREFIX}_body_heading.fst $OUTDIR/temp/${PREFIX}_body2.fst > $OUTDIR/temp/${PREFIX}_prep.fst

		if  [[ -f ${OUTDIR}/${PREFIX}_avoidcols.txt ]] ; then
			COMBLIST2=$(awk '{print}' ORS=' ' ${OUTDIR}/${PREFIX}_avoidcols.txt) 
			COMBLIST3=$(echo $COMBLIST2 | sed 's/[^ ][^ ]*/"&"/g')
				ZERO=$(echo $COMBLIST3 | sed -r 's/([^ ]+)/$i!=\1/g') 
				TWO=$(echo $ZERO | sed -e 's/ /\  \&\&\ /g')
				ONE=$(echo " NR==1{for(i=1; i<=NF; i++) if (")
				THREE=$(echo ' ) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"} ')
				FOUR=$(echo ${ONE}${TWO}${THREE})
				#Use above expression to remove Sfst columns in similar populations
				awk "$FOUR" $OUTDIR/temp/${PREFIX}_prep.fst > $OUTDIR/temp/${PREFIX}_keep.fst
				NCOLZ=$(awk '{print NF}' $OUTDIR/temp/${PREFIX}_keep.fst | sort -nu | tail -n 1)
				if [ "$NCOLZ" -gt 1 ]; then
					awk '{if (NR!=1) {print}}' $OUTDIR/temp/${PREFIX}_keep.fst | awk '{ s = 0; for (i = 1; i <= NF; i++) s += $i; print $1, (NF > 1) ? s / (NF - 0) : 0; }' | awk '{print $2}' > $OUTDIR/temp/${PREFIX}_keep2.fst
				else
					awk '{if (NR!=1) {print}}' $OUTDIR/temp/${PREFIX}_keep.fst > $OUTDIR/temp/${PREFIX}_keep2.fst
				fi

				#Create final Sfst file for this analysis.
				paste $OUTDIR/temp/${PREFIX}_head.fst $OUTDIR/temp/${PREFIX}_keep2.fst | awk '{gsub("\t","",$0); print;}'  > $OUTDIR/temp/${PREFIX}_SfstNA.fst
		else 
			if [ "$NPOPS" -gt 2 ]; then
				awk '{if (NR!=1) {print}}' $OUTDIR/temp/${PREFIX}_prep.fst | awk '{ s = 0; for (i = 1; i <= NF; i++) s += $i; print $1, (NF > 1) ? s / (NF - 0) : 0; }' | awk '{print $2}' > $OUTDIR/temp/${PREFIX}_keep2.fst
			else
				cp <(awk '{if (NR!=1) {print}}' $OUTDIR/temp/${PREFIX}_prep.fst) $OUTDIR/temp/${PREFIX}_keep2.fst
			fi
			
				#Create final Sfst file for this analysis.
				paste $OUTDIR/temp/${PREFIX}_head.fst $OUTDIR/temp/${PREFIX}_keep2.fst | awk '{gsub("\t","",$0); print;}'  > $OUTDIR/temp/${PREFIX}_SfstNA.fst
		fi

					awk '!/na/' $OUTDIR/temp/${PREFIX}_SfstNA.fst > $OUTDIR/${PREFIX}.fst
					declare -i before=$(wc -l $OUTDIR/${PREFIX}_raw.fst |  cut -f1 -d' ')
					declare -i after=$(wc -l $OUTDIR/${PREFIX}.fst | cut -f1 -d' ')
					declare -i loss=$(($before-$after))
					echo "ALERT: $loss SNPs were removed from fst analysis due to Ns, indels, or uninformative comparisons"
					echo "ALERT: $after fst SNPs  remain"
					rm $OUTDIR/temp/${PREFIX}**.fst*
	fi

fi

if [[ "$SLIDINGFST" =~(on)$ ]]; then
#Cleanup FST sliding window OUTDIR format, average rows, and produce intuitive format
	if [[ -f $OUTDIR/${PREFIX}.Sfst ]] ; then
		echo "ALERT: $OUTDIR/${PREFIX}.Sfst exists; will not redo sliding FST reformatting"
	else
		echo "ALERT: now doing SFST reformatting"
		awk '$5 >= '$MINCOV' &&  $5  <= '$MAXCOV' ' $OUTDIR/${PREFIX}_raw.Sfst >  $OUTDIR/temp/${PREFIX}_reclassing.Sfst

		#Cut up Sfst file
		cut -d$'\t' -f 1-2 $OUTDIR/temp/${PREFIX}_reclassing.Sfst | gawk '$3=(FNR FS $3)' > $OUTDIR/temp/${PREFIX}_head.Sfst
		cut -d$'\t' -f 6- $OUTDIR/temp/${PREFIX}_reclassing.Sfst > $OUTDIR/temp/${PREFIX}_body.Sfst
		awk 'NR==1{print $0}' $OUTDIR/temp/${PREFIX}_body.Sfst | awk  '{gsub("=","\t",$0); print;}' |  awk '{ for (i=2;i<=NF;i+=2) $i="" } 1' | awk  '{gsub("  ","\t",$0); print;}' | awk  '{gsub(" ","",$0); print;}'   > $OUTDIR/temp/${PREFIX}_body_heading.Sfst
		awk '{gsub("=","\t",$0); print;}' $OUTDIR/temp/${PREFIX}_body.Sfst | awk '{ for (i=1;i<=NF;i+=2) $i="" } 1' | awk  '{gsub("  ","\t",$0); print;}' | awk  '{gsub(" ","",$0); print;}'  > $OUTDIR/temp/${PREFIX}_body2.Sfst
		cat $OUTDIR/temp/${PREFIX}_body_heading.Sfst $OUTDIR/temp/${PREFIX}_body2.Sfst > $OUTDIR/temp/${PREFIX}_prep.Sfst

		if  [[ -f ${OUTDIR}/${PREFIX}_avoidcols.txt ]] ; then
			COMBLIST2=$(awk '{print}' ORS=' ' ${OUTDIR}/${PREFIX}_avoidcols.txt) 
			COMBLIST3=$(echo $COMBLIST2 | sed 's/[^ ][^ ]*/"&"/g')
				ZERO=$(echo $COMBLIST3 | sed -r 's/([^ ]+)/$i!=\1/g') 
				TWO=$(echo $ZERO | sed -e 's/ /\  \&\&\ /g')
				ONE=$(echo " NR==1{for(i=1; i<=NF; i++) if (")
				THREE=$(echo ' ) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"} ')
				FOUR=$(echo ${ONE}${TWO}${THREE})
				#Use above expression to remove Sfst columns in similar populations
				awk "$FOUR" $OUTDIR/temp/${PREFIX}_prep.Sfst > $OUTDIR/temp/${PREFIX}_keep.Sfst
				NCOLZ=$(awk '{print NF}' $OUTDIR/temp/${PREFIX}_keep.Sfst | sort -nu | tail -n 1)
				if [ "$NCOLZ" -gt 1 ]; then
					awk '{if (NR!=1) {print}}' $OUTDIR/temp/${PREFIX}_keep.Sfst | awk '{ s = 0; for (i = 1; i <= NF; i++) s += $i; print $1, (NF > 1) ? s / (NF - 0) : 0; }' | awk '{print $2}' > $OUTDIR/temp/${PREFIX}_keep2.Sfst
				else
					awk '{if (NR!=1) {print}}' $OUTDIR/temp/${PREFIX}_keep.Sfst > $OUTDIR/temp/${PREFIX}_keep2.Sfst
				fi
				#Create final Sfst file for this analysis.
				paste $OUTDIR/temp/${PREFIX}_head.Sfst $OUTDIR/temp/${PREFIX}_keep2.Sfst | awk '{gsub("\t","",$0); print;}'  > $OUTDIR/temp/${PREFIX}_SfstNA.Sfst
		else 
			if [ "$NPOPS" -gt 2 ]; then
				awk '{if (NR!=1) {print}}' $OUTDIR/temp/${PREFIX}_prep.Sfst | awk '{ s = 0; for (i = 1; i <= NF; i++) s += $i; print $1, (NF > 1) ? s / (NF - 0) : 0; }' | awk '{print $2}' > $OUTDIR/temp/${PREFIX}_keep2.Sfst
			else
				cp <(awk '{if (NR!=1) {print}}' $OUTDIR/temp/${PREFIX}_prep.Sfst) $OUTDIR/temp/${PREFIX}_keep2.Sfst
			fi
			
				#Create final Sfst file for this analysis.
				paste $OUTDIR/temp/${PREFIX}_head.Sfst $OUTDIR/temp/${PREFIX}_keep2.Sfst | awk '{gsub("\t","",$0); print;}'  > $OUTDIR/temp/${PREFIX}_SfstNA.Sfst
		fi

					awk '!/na/' $OUTDIR/temp/${PREFIX}_SfstNA.Sfst > $OUTDIR/${PREFIX}.Sfst
					declare -i before=$(wc -l $OUTDIR/${PREFIX}_raw.Sfst |  cut -f1 -d' ')
					declare -i after=$(wc -l $OUTDIR/${PREFIX}.Sfst | cut -f1 -d' ')
					declare -i loss=$(($before-$after))
					echo "ALERT: $loss SNPs were removed from Sfst analysis due to Ns, indels, or uninformative comparisons"
					echo "ALERT: $after Sfst SNP positions remain"
					rm $OUTDIR/temp/${PREFIX}**.Sfst*
	fi
fi


if [[ "$FET" =~(on)$ ]]; then
	if [[ -f $OUTDIR/${PREFIX}.fet ]] ; then
		echo "ALERT: $OUTDIR/${PREFIX}.fet exists; will not redo FET reformatting"
	else
		echo "ALERT: now doing FET reformatting Step 1"
		awk '$5 >= '$MINCOV' &&  $5  <= '$MAXCOV' ' $OUTDIR/${PREFIX}_raw.fet >  $OUTDIR/temp/${PREFIX}_reclassing.fet

		#Cut up fet file
		echo "ALERT: now doing FET reformatting Step 2"
		cut -d$'\t' -f 1-2 $OUTDIR/temp/${PREFIX}_reclassing.fet | gawk '$3=(FNR FS $3)' > $OUTDIR/temp/${PREFIX}_head.fet
		echo "ALERT: now doing FET reformatting Step 3"
		cut -d$'\t' -f 6- $OUTDIR/temp/${PREFIX}_reclassing.fet > $OUTDIR/temp/${PREFIX}_body.fet
		echo "ALERT: now doing FET reformatting Step 4"
		awk 'NR==1{print $0}' $OUTDIR/temp/${PREFIX}_body.fet | awk  '{gsub("=","\t",$0); print;}' |  awk '{ for (i=2;i<=NF;i+=2) $i="" } 1' | awk  '{gsub("  ","\t",$0); print;}' | awk  '{gsub(" ","",$0); print;}'   > $OUTDIR/temp/${PREFIX}_body_heading.fet
		echo "ALERT: now doing FET reformatting Step 5"
		awk '{gsub("=","\t",$0); print;}' $OUTDIR/temp/${PREFIX}_body.fet | awk '{ for (i=1;i<=NF;i+=2) $i="" } 1' | awk  '{gsub("  ","\t",$0); print;}' | awk  '{gsub(" ","",$0); print;}'  > $OUTDIR/temp/${PREFIX}_body2.fet
		echo "ALERT: now doing FET reformatting Step 6"
		cat $OUTDIR/temp/${PREFIX}_body_heading.fet $OUTDIR/temp/${PREFIX}_body2.fet > $OUTDIR/temp/${PREFIX}_prep.fet
		

		if  [[ -f ${OUTDIR}/${PREFIX}_avoidcols.txt ]] ; then
			echo "ALERT: now doing FET reformatting Step 7 - the long one"
			COMBLIST2=$(awk '{print}' ORS=' ' ${OUTDIR}/${PREFIX}_avoidcols.txt) 
			COMBLIST3=$(echo $COMBLIST2 | sed 's/[^ ][^ ]*/"&"/g')
			ZERO=$(echo $COMBLIST3 | sed -r 's/([^ ]+)/$i!=\1/g') 
			TWO=$(echo $ZERO | sed -e 's/ /\  \&\&\ /g')
			ONE=$(echo " NR==1{for(i=1; i<=NF; i++) if (")
			THREE=$(echo ' ) {a[i]++;} } { for (i in a) printf "%s\t", $i; printf "\n"} ')
			FOUR=$(echo ${ONE}${TWO}${THREE})
			#Use above expression to remove fet columns in similar populations
			awk "$FOUR" $OUTDIR/temp/${PREFIX}_prep.fet > $OUTDIR/temp/${PREFIX}_keep.fet
			awk '{if (NR!=1) {print}}' $OUTDIR/temp/${PREFIX}_keep.fet > $OUTDIR/temp/${PREFIX}_keep2.fet
			declare -i COL=$( head -n 1 $OUTDIR/temp/${PREFIX}_keep2.fet | wc -w )

			#Create final fet file for this analysis.
		
# 				rin=$OUTDIR/temp/${PREFIX}_keep2.fet
# 				rout=$OUTDIR/temp/
# 				Rscript $BASEDIR/rscripts/r_combine_p.R $rin $rout

	##chunking code
#			declare -i LEN=( `wc -l $OUTDIR/temp/${PREFIX}_keep2.fet `)
			declare -i LEN=$( cat $OUTDIR/temp/${PREFIX}_keep2.fet | wc -l )

			if [ $LEN -le 1000000 ]; then 

				echo "Combining FET p-values without division"

				rin=$OUTDIR/temp/${PREFIX}_keep2.fet
				rout=$OUTDIR/temp/
				Rscript $BASEDIR/rscripts/r_combine_p.R $rin $rout $COL

			##reduce and tally subsets of individuals, then tally population together
			elif [ $LEN -gt 1000000 ]; then 

				echo "Combining FET p-values by dividing into subsets of 1,000,000 lines"

				#split the FET pairwise p-value file into multiple files
				split --numeric-suffixes=1 -l 1000000 --additional-suffix='_temp.fet' $OUTDIR/temp/${PREFIX}_keep2.fet $OUTDIR/temp/keep2_
			
				FILES=( `ls $OUTDIR/temp/keep2_*_temp.fet `)
				LENb=( `ls $OUTDIR/temp/keep2_*_temp.fet | wc -l `)
				echo "There appear to be "$LENb" set FET files to process and combine."
				LENb=$(($LENb - 1))

				for ((j = 0; j <= $LENb; j++));
				do

					rin=${FILES[j]}
					rout=$OUTDIR/temp/
					Rscript $BASEDIR/rscripts/r_combine_p.R $rin $rout $COL

				done
			 
				cat $OUTDIR/temp/keep2_*_temp2.fet > $OUTDIR/temp/${PREFIX}_keep22.fet
			fi

	##end chunking code

		else 
				#Create final fet file for this analysis.
			echo "ALERT: now doing FET reformatting Step 7 (the long one)"

			awk '{if (NR!=1) {print}}' $OUTDIR/temp/${PREFIX}_prep.fet > $OUTDIR/temp/${PREFIX}_keep2.fet


	##chunking code
			declare -i LEN=$( cat $OUTDIR/temp/${PREFIX}_keep2.fet | wc -l )
			declare -i COL=$( head -n 1 $OUTDIR/temp/${PREFIX}_keep2.fet | wc -w )

			if [ $LEN -le 1000000 ]; then 

				echo "Combining FET p-values without division"

				rin=$OUTDIR/temp/${PREFIX}_keep2.fet
				rout=$OUTDIR/temp/
				Rscript $BASEDIR/rscripts/r_combine_p.R $rin $rout $COL

			##reduce and tally subsets of individuals, then tally population together
			elif [ $LEN -gt 1000000 ]; then 

				echo "Combining FET p-values by dividing into subsets of 1,000,000 lines"

				#split the FET pairwise p-value file into multiple files
				split --numeric-suffixes=1 -l 1000000 --additional-suffix='_temp.fet' $OUTDIR/temp/${PREFIX}_keep2.fet $OUTDIR/temp/keep2_
			
				FILES=( `ls $OUTDIR/temp/keep2_*_temp.fet `)
				LENb=( `ls $OUTDIR/temp/keep2_*_temp.fet | wc -l `)
				echo "There appear to be "$LENb" set FET files to process and combine."
				LENb=$(($LENb - 1))

				for ((j = 0; j <= $LENb; j++));
				do

					rin=${FILES[j]}
					rout=$OUTDIR/temp/
					Rscript $BASEDIR/rscripts/r_combine_p.R $rin $rout $COL

				done
			 
				cat $OUTDIR/temp/keep2_*_temp2.fet > $OUTDIR/temp/${PREFIX}_keep22.fet
			fi

	##end chunking code

		fi
	
		echo "ALERT: now doing FET reformatting Step 8"
		paste $OUTDIR/temp/${PREFIX}_head.fet $OUTDIR/temp/${PREFIX}_keep22.fet | awk '{gsub("\t","",$0); print;}'  > $OUTDIR/temp/${PREFIX}_fetNA.fet
		echo "ALERT: now doing FET reformatting Step 9"
		awk '!/na/' $OUTDIR/temp/${PREFIX}_fetNA.fet > $OUTDIR/${PREFIX}.fet
		declare -i before=$(wc -l $OUTDIR/${PREFIX}_raw.fet |  cut -f1 -d' ')
		declare -i after=$(wc -l $OUTDIR/${PREFIX}.fet | cut -f1 -d' ')
		declare -i loss=$(($before-$after))
		echo "ALERT: $loss SNPs were removed from FET analysis due to Ns, indels, or uninformative comparisons"
		echo "ALERT: $after FET SNPs remain"

	fi
fi


#########################

#	ls $OUTDIR/temp/${PREFIX}*
#	rm $OUTDIR/temp/${PREFIX}*

echo "ALERT: PPanalyze done at $(date)"