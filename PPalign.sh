#!/bin/bash

#PoolParty v05.09.2022
#PPalign

set -o pipefail

BASEDIR=$(dirname "$0")

if [[ ${BASEDIR} != *poolparty* ]];then
	SLOC=$(readlink -f "$0")
	BASEDIR=${SLOC%/*}
fi

echo $BASEDIR

if  ( [[ $(echo $1)  = "" ]] )  ; then
	echo "ERROR: You must provide a config file after the PPalign command. Example usage"
	echo "bash PPalign.sh pp_align.config [Optional: sample list file]"
	exit
fi
	echo "ALERT: $1 has been specified as the configuration file"

source $1

#Declare integers for numerical parameters
declare -i THREADZ
declare -i BQUAL
declare -i MINLENGTH
declare -i INWIN
declare -i MAPQ
declare -i SNPQ
declare -i MINDP

		rm $OUTDIR/${OUTPOP}_sample_pops.txt  &> /dev/null 
		rm $OUTDIR/${OUTPOP}_sample_files.txt &> /dev/null 
		rm $OUTDIR/${OUTPOP}_prefixes.txt  &> /dev/null 



#Get date of run and save an additional config file with run parameters
RUNDATE=$(date +'%m_%d_%Y')

echo "ALERT: Beginning PoolParty run at $(date)"

########################################################################################################
####-----------------------------------------RUNNING-----------------------------------------------#####
########################################################################################################

	echo "ALERT: Performing error checking"

#Check for samplelist file, quit if it does not exist
if  ( [[ $(echo $2)  = "" ]] )  ; then
	if  ( [[ $(echo "${SAMPLELIST}")  = "" ]] )  ; then
		if [[ ! -f $INDIR/"samplelist.txt" ]] ; then
 			echo 'ERROR: No sample list is specified or available. Aborting.'
   			exit 1
   		else
   			SAMPLE=("samplelist.txt")
		fi
	else
		SAMPLE=${SAMPLELIST}
	fi
else
	SAMPLE=$2
fi

#Confirm with user that correct sample list has been identified [optional]
if [ ! -z $CONFIRM ] &&  [[ "$CONFIRM" =~(on)$ ]] ; then
	if [[ -f $INDIR/"$SAMPLE" ]] ; then
		echo -e "\nSample list file "$SAMPLE" was found. Is this correct? Type 'yes' or 'no'"
		read Sampcorrect

		if [ "$Sampcorrect" == "no" ]; then
        	echo "Please double check that your sample list has been specified correctly."
        	exit 1
		elif [ "$Sampcorrect" == "yes" ]; then
            echo "Proceeding with "$SAMPLE
		else
        	echo "Incorrect Input"
        	exit 1
		fi
	else
   		echo "ERROR: File "$SAMPLE" could not be found. Aborting."
   		exit 1
	fi
fi

	#Check samplelist for invalid characters
		if grep -q / $INDIR/"$SAMPLE" ; then
  			echo 'ERROR: "$SAMPLE" contains invalid characters. Do not include directory names or backslashes in samplelist'
			exit
		fi

	#Check for at least two entries in sample list (since pipeline input is paired-end data)
		leN=$(wc -l $INDIR/"$SAMPLE" | cut -f1 -d' ') 
		if  [[ $leN -lt 1 ]]  ; then
 			 echo 'ERROR: File "$SAMPLE" has less than two rows. Something wrong with file or not paired-end reads.'
   			 exit
		fi
	#Check for genome fasta file
		if [[ ! -f ${GENOME} ]] ; then
 	 		echo "ERROR: ${GENOME} genome file is missing"
   		 	exit
		fi
	#Check for index file
		if [[ ! -f ${GENOME}.fai ]] ; then
 			 echo "ERROR:: ${GENOME}.fai index file is missing, place this in the same dir as the genome .fasta"
   			 exit
		fi
	#Check skip break is either on or off
		if [[ "$SPLITDISC" != "off" ]] && [[ "$SPLITDISC" != "on" ]] ; then
			echo "ERROR: Wrong Parameters: SPLITDISC must be either on or off"
			exit
		fi
	#Check individual analysis is either on or off
		if [[ "$INDCONT" != "off" ]] && [[ "$INDCONT" != "on" ]] ; then
			echo "ERROR: Wrong Parameters: INDCONT must be either on or off"
			exit
		fi
	#Check quality report is either on or off
		if [[ "$QUALREPORT" != "off" ]] && [[ "$QUALREPORT" != "on" ]] ; then
			echo "ERROR: Wrong Parameters: QUALREPORT must be either on or off"
			exit
		fi
	#Ensure base directories all exist
		if [ ! -d "$INDIR" ] || [ ! -d "$BBMAPDIR" ]  || [ ! -d "$POOL2" ] ; then
			echo "ERROR: Wrong Directories: Either INDIR, POOL2, or BBMAPDIR don't exist"
			exit
		fi
	#Ensure poolparty directory is intact
		if [ ! -d "$BASEDIR/rscripts" ] ; then
			echo "ERROR: Poolparty directory has been tampered with"
			exit
		fi
	#Ensure Popoolation2 directory structure has not been tampered with
		if [[ ! -f ${POOL2}/indel_filtering/"identify-indel-regions.pl" ]] ||  [[ ! -f ${POOL2}/indel_filtering/"filter-sync-by-gtf.pl" ]] \
			||  [[ ! -f ${POOL2}/"mpileup2sync.jar" ]]  ; then
			echo "ERROR: A component of Popoolation2 is missing. Ensure the directory structure in Popoolation2 is unaltered."
			echo "Ensure mpileup2sync.jar, filter-sync-by-gtf.pl, and identify-indel-regions.pl are all present"
			exit
		fi
	#Ensure BBMAP directory structure has not been tampered with
		if [[ ! -f ${BBMAPDIR}/"bbduk.sh" ]] ||  [[ ! -f ${BBMAPDIR}/resources/"adapters.fa" ]]  ; then
			echo "ERROR: A component of bbmap is missing. Ensure the directory structure in bbmap is unaltered."
			echo "Ensure that bbduk.sh and /resources/adapters.fa exist in the bbmap directory"
			exit
		fi
	#Ensure that all dependencies exist on the system (both the path and the command have to fail)
		if  ( [[ $(command -v "$FASTQC")  = "" ]]  &&  [[ ! -f $FASTQC ]] ) || ( [[ $(command -v "$BWA")  = "" ]]  &&  [[ ! -f $BWA ]] ) \
				|| ( [[ $(command -v "$SAMTOOLS")  = "" ]]  &&  [[ ! -f $SAMTOOLS ]] ) ||  ( [[ $(command -v "$PICARDTOOLS")  = "" ]]  &&  [[ ! -f $PICARDTOOLS ]] ) \
				|| ( [[ $(command -v "$SAMBLASTER")  = "" ]]  &&  [[ ! -f $SAMBLASTER ]] ) || ( [[ $(command -v "$BCFTOOLS")  = "" ]]  &&  [[ ! -f $BCFTOOLS ]] ) ; then
			echo "ERROR: One or more dependencies are incorrect, double check dependency locations and names"
			echo "Type each dependency into the terminal, as it is listed in the config file, and ensure that it initiates"		
			exit
		fi
	#Make sure java memory is set properly
		if ! printf '%s\n' "$KMEM" | grep -Fqe "Xmx"; then
			echo "ERROR: Java memory set incorrectly, check value and ensure it followes Xmx#g format"
			exit
		fi
	#Ensure other parameters are within a valid range
		if  [[ ! "$THREADZ" =~ ^[0-9]+$ ]] || [[ ! "$BQUAL" =~ ^[0-9]+$ ]] || [[ ! "$MINLENGTH" =~ ^[0-9]+$ ]] || [[ ! "$INWIN" =~ ^[0-9]+$ ]] \
				|| [[ ! "$MAPQ" =~ ^[0-9]+$ ]] || [[ ! "$MINDP" =~ ^[0-9]+$ ]] ;  then
        		echo "ERROR: THREADZ, BQUAL, MINLENGTH, INWIN, and MAPQ must all be positive integers. Check parameters"
			exit
		fi
		if  [[ ! "$MAF" =~ ^[+-]?[0-9]+\.?[0-9]*$  ]] ; then 
				echo "ERROR: MAF must be a value between 0 and 1"
				exit
		fi
	#Check for R, perl, java on the system
		if  [[ $(command -v java)  = "" ]] ; then
			echo "ERROR: java not detected on system. java is required for multiple packages"
			echo "java should initiate when 'java' is typed in terminal"
			exit
		fi
		if  [[ $(command -v perl)  = "" ]] ; then
			echo "ERROR: perl not detected on system. perl is required for Popoolation2"
			echo "perl should initiate when 'perl' is typed in terminal"
			exit
		fi
		if  [[ $(command -v Rscript)  = "" ]] ; then
			echo "ERROR: R not detected on system. R is required for analysis "
			echo "R should initiate when 'Rscript' is typed in terminal"
			exit
		fi
	#Check for piping ability 
		if  [[ $(command -v mkfifo)  = "" ]] ; then
			echo "WARNING: piping 'mkfifo' not detected on system or available on drive. This may cause issues in downstream analyses"
			echo "Edit PPalign.sh and redirect mkfifo to another drive"
		fi
	#check for gawk
		if  [[ $(command -v gawk)  = "" ]] ; then
			echo "ERROR: gawk not detected on system. gawk is usually standard on Linux systems. Install then retry"
			exit
		fi
	#Check for process substitution 
		if  [[ $(command -v cat <(date); echo $? )  = "" ]] ; then
			echo "WARNING: process substitution does not appear to be working on this system. Errors may arise"
		fi
	#Load R to check for dependencies 
		if  [[ -f $OUTDIR/R_ERROR.txt ]] ; then
			rm $OUTDIR/R_ERROR.txt 
		fi		

		Rscript $BASEDIR/rscripts/r_align_check.R $OUTDIR

		if  [[ -f $OUTDIR/R_ERROR.txt ]] ; then
			echo "ERROR: R dependency check failed, install dependencies manually"
			exit
		fi

	echo "ALERT: Parameter check passed. moving on..."

#########################################################################################################
	#Make the output directory if it does not exist 
	if [ ! -d "$OUTDIR/" ]; then
		mkdir $OUTDIR/
	fi

	if [ ! -d "$OUTDIR/tmp" ]; then
		mkdir $OUTDIR/tmp
	fi

	#Copy configuration file so you know what you did later
	cp $1 ${OUTDIR}/${OUTPOP}_${RUNDATE}.config

	#Remove any potential extra characters (/r) and sort samplelist
	awk '{ sub("\r$", ""); print $0 }'  ${INDIR}/"$SAMPLE" > ${INDIR}/samplelist.tmp && mv ${INDIR}/samplelist.tmp ${INDIR}/$SAMPLE
	sort -f ${INDIR}/"$SAMPLE" | awk '{print $1}' | grep -P -v '^\s*$' >  $OUTDIR/${OUTPOP}_sample_files.txt
	sort -f ${INDIR}/"$SAMPLE" | awk '{print $2}' | grep -P -v '^\s*$' > $OUTDIR/${OUTPOP}_sample_pops.txt
	FILES=$(cat $OUTDIR/${OUTPOP}_sample_files.txt)
	POPS=$(cat $OUTDIR/${OUTPOP}_sample_pops.txt)
		aCHCK1=$(wc -l $OUTDIR/${OUTPOP}_sample_files.txt | awk '{print $1}' )
		aCHCK2=$(wc -l $OUTDIR/${OUTPOP}_sample_pops.txt | awk '{print $1}' )
			#Ensure each sample has been given a population designation
			if [[ $aCHCK1 != $aCHCK2  ]] ; then
				echo "ERROR: Something wrong with $SAMPLE. $aCHCK1 samples but only $aCHCK2 population designations.."
				echo "Check for empty lines in $SAMPLE"
				exit
			fi

	#Check mate pairs; fail if order is incorrect
	cut -d_ -f1  $OUTDIR/${OUTPOP}_sample_files.txt >  $OUTDIR/${OUTPOP}_prefixes.txt

	CHECK=$(cat $OUTDIR/${OUTPOP}_prefixes.txt)
	CHECK2=($CHECK)

	for (( i=0; i<${#CHECK2[@]} ; i+=2 )) ; do
		echo  Checking pairs "${CHECK2[i+1]}" and "${CHECK2[i]}"
		if [ "${CHECK2[i+1]}" == "${CHECK2[i]}" ]; then
		echo "ALERT: Correct Mates, proceeding"
		single=off
	else
		echo "ALERT: Naming convention suggests single-end reads. Proceeding using single-end fqs."
		echo "If reads are paired-end, check naming convention and try again."
		single=on
		break
	fi
	done

		#Check that populations are designated in the samplelist
		PnUm=$(awk '{ a[$2]++ } END { for (b in a) { print b } }' $INDIR/"$SAMPLE")
		if [[ $PnUm  = "" ]] ; then
			echo 'ERROR: Something wrong with "$SAMPLE" . Are population numbers specified in column 2?'
			exit
		fi
		parray=($PnUm)
		parraynum=$(echo ${#parray[@]})
		echo -e "\nALERT: You have specified $parraynum unique populations. If this sounds wrong then you done messed up..."


#Checks that, for paired read files, populations specified for each read file are the same; fail if not
if  ( [[ $single = "off" ]] )  ; then
	POPS2=( `cut -f1 $OUTDIR/${OUTPOP}_sample_pops.txt `)

	for (( i=0; i<${#POPS2[@]} ; i+=2 )) ; do
		echo  Checking pairs "${POPS2[i+1]}" and "${POPS2[i]}"
		if [ "${POPS2[i+1]}" == "${POPS2[i]}" ]; then
			echo "ALERT: Consistent population designations, proceeding"
		else
			echo "ALERT: Population designations are not consistent."
			echo "If reads are paired-end, check population designations and try again."
			exit 1
		fi
	done


fi



#Check PPalign has detected the proper number of read files to individuals and exits if incorrect
if  ( [[ $single = "off" ]] )  ; then
	INDnum=$((${#CHECK2[@]}/2))
else
	INDnum=${#CHECK2[@]}
fi

if [ ! -z $CONFIRM ] &&  [[ "$CONFIRM" =~(on)$ ]] ; then
echo -e "\n"${#CHECK2[@]} "read files are detected for "$INDnum" individuals in "$parraynum" populations. Is this correct? Enter 'yes' or 'no' and press [ENTER]"

read Indcorrect

if [ "$Indcorrect" == "no" ]; then
        echo "Please double check that all fastq files are named PopA001_R1.fq.gz and PopA001_R2.fq.gz or equivalent"
        exit 1
elif [ "$Indcorrect" == "yes" ]; then
            echo "Proceeding with "$INDnum" individuals"
else
        echo "Incorrect Input"
        exit 1
fi
fi



if [[ "$single" =~(off)$ ]] ; then
#Create prefix file; this indicates prefix before the first "_" which should be library identifiers
	cut -d_ -f1 $OUTDIR/${OUTPOP}_sample_files.txt | awk '!seen[$0]++' | sed '/^$/d'  >  $OUTDIR/${OUTPOP}_prefixes.txt
	printf "ALERT: Using libraries:\n$(cat $OUTDIR/${OUTPOP}_prefixes.txt)\n"
		declare -i CHCK1=$(wc -l $OUTDIR/${OUTPOP}_prefixes.txt | cut -f1 -d' ' )
		declare -i CHCK2=$(wc -l $OUTDIR/${OUTPOP}_sample_files.txt | cut -f1 -d' '  )
		declare -i CHCK3=($CHCK2)/2
		if [[ $CHCK1 != $CHCK3  ]] ; then
        	echo "ERROR: number of samples and prefixes don't match up ($CHCK1 vs $CHCK3) ; there is likely something wrong with the filenames or samplelist"
			echo "Note that file name prefixes (text before the first underscore) must be unique for the paired files!"
			echo "Check _sample_files.txt, paired-end libraries should be stacked on top of one another; if this isn't the case change the naming convention"
			exit
		fi
fi

if [[ "$single" =~(on)$ ]] ; then
#Create prefix file; this indicates prefix before the first "_" which should be library identifiers
	cut -d_ -f1 $OUTDIR/${OUTPOP}_sample_files.txt | awk '!seen[$0]++' | sed '/^$/d'  >  $OUTDIR/${OUTPOP}_prefixes.txt
	printf "ALERT: Using single-ended libraries:\n$(cat $OUTDIR/${OUTPOP}_prefixes.txt)\n"
fi

#Get number of populations that are represented by the libraries in your samplelist
	if  [[ -f $OUTDIR/${OUTPOP}_poplist.txt ]] ; then
		skipmerge=on
		rm  $OUTDIR/${OUTPOP}_poplist.txt 
	fi

	array=($FILES)
	array2=($POPS)
	echo "ALERT: Checking ${#array[@]} file names in samplelist" 
	if [[ "$single" =~(off)$ ]] ; then
		for (( i=0; i<${#array[@]} ; i+=2 )) ; do
			e=${array[i]%%.*}
			f=${array2[i]}
			oZ=$(echo $e $f)
			echo ${oZ} >> $OUTDIR/${OUTPOP}_poplist.txt
			#remove any introduced characters
			awk '{ sub("\r$", ""); print $0 }'  $OUTDIR/${OUTPOP}_poplist.txt > $OUTDIR/${OUTPOP}_poplist.tmp && mv $OUTDIR/${OUTPOP}_poplist.tmp $OUTDIR/${OUTPOP}_poplist.txt
		done
	fi
	if [[ "$single" =~(on)$ ]] ; then
		for (( i=0; i<${#array[@]} ; i+=1 )) ; do
			e=${array[i]%%.*}
			f=${array2[i]}
			oZ=$(echo $e $f)
			echo ${oZ} >> $OUTDIR/${OUTPOP}_poplist.txt
			#remove any introduced characters
			awk '{ sub("\r$", ""); print $0 }'  $OUTDIR/${OUTPOP}_poplist.txt > $OUTDIR/${OUTPOP}_poplist.tmp && mv $OUTDIR/${OUTPOP}_poplist.tmp $OUTDIR/${OUTPOP}_poplist.txt
		done
	fi

#Get anchored chromosome lengths in bp. Prints a copy to the rundir for later analyses
	echo "ALERT: Getting genome anchored index..."
		HNAM=$(cut -f1 ${GENOME}.fai |  sed -re "s/[^a-zA-Z]*//g"  |  awk '!x[$0]++')
		HNAM2=$(echo $HNAM)
		echo "ALERT: $HNAM2 are your chromosome and/or scaffold headings"
		if [ -z "$SCAHEAD" ]; then
			cut -f1,2 ${GENOME}.fai >  $OUTDIR/${OUTPOP}_CHRbp1.txt
		else
			cut -f1,2 ${GENOME}.fai  | grep -v "^${SCAHEAD}"  >  $OUTDIR/${OUTPOP}_CHRbp1.txt
		fi
	awk '{print $1}' $OUTDIR/${OUTPOP}_CHRbp1.txt  | awk  '$2="1"' | awk '{gsub(" ","\t",$0); print;}' > $OUTDIR/${OUTPOP}_CHRbp2.txt
	cat $OUTDIR/${OUTPOP}_CHRbp1.txt $OUTDIR/${OUTPOP}_CHRbp2.txt > $OUTDIR/${OUTPOP}_CHRbp.txt
	rm $OUTDIR/${OUTPOP}_CHRbp1.txt ; rm $OUTDIR/${OUTPOP}_CHRbp2.txt

# Make separate files for pops, indicating which library belongs to which population
	if [ ! -d "$OUTDIR/pops" ]; then
		mkdir $OUTDIR/pops
	fi

#Split into one file per population to combine
	awk '{print $1 > "'${OUTDIR}/pops/pop_'"$2".txt"}' $OUTDIR/${OUTPOP}_poplist.txt
	popfiles=$(awk '{ a[$2]++ } END { for (b in a) { print b } }' $INDIR/"$SAMPLE")
	printf '%s\n' "${popfiles[@]}" | awk '{print "pop_" $0;}' | sort > ${OUTDIR}/pops/${OUTPOP}_files_for_pops.txt
	awk '{ sub("\r$", ""); print $0 }'  ${OUTDIR}/pops/${OUTPOP}_files_for_pops.txt > ${OUTDIR}/pops/${OUTPOP}_files_for_pops.tmp && mv ${OUTDIR}/pops/${OUTPOP}_files_for_pops.tmp ${OUTDIR}/pops/${OUTPOP}_files_for_pops.txt


#Trim by quality score, uses BBMAP (java-based)
##! ADDITIONAL TRIM PARAMETERS CAN BE ADDED BELOW AFTER "${BBMAPDIR}/bbduk.sh" !##
##CHECK BBMAP DOCUMENTATION FOR ADDITIONAL OPTIONS AND MODIFY LINE BELOW ##
	if [ ! -d "$OUTDIR/trimmed" ]; then
		mkdir $OUTDIR/trimmed
	fi

	if  [[ -f $OUTDIR/${OUTPOP}_names.txt ]] ; then
		rm  $OUTDIR/${OUTPOP}_names.txt 
	fi

	array=($FILES)

	if [[ "$single" =~(off)$ ]] ; then
		echo "ALERT: Proceeding to trimming of paired-end reads."
		for (( i=0; i<${#array[@]} ; i+=2 )) ; do
			b=${array[i]%%.*}
			c=${array[i+1]%%.*}
			o1=${OUTDIR}/trimmed/${b}.trim_1
			o2=${OUTDIR}/trimmed/${b}.trim_2 

			if  [[ -f $o1 ]] || [[ -f $o1.gz ]] ; then
				echo "ALERT: $o1 pair exists; skipping"
				echo ${b} >> $OUTDIR/${OUTPOP}_names.txt
			fi

			if [[ ! -f $o1 ]] && [[ ! -f $o2 ]] && [[ ! -f $o1.gz ]] && [[ ! -f $o2.gz ]]; then
			#Add modifications below if needed#
				echo ${b} >> $OUTDIR/${OUTPOP}_names.txt
				if [[ -f ${OUTDIR}/BAM/${b}_filtered.bam ]] ; then
					echo "ALERT: ${b}_filtered.bam exists and will not be re-trimmed or re-aligned."
				else
					echo "ALERT: ${b}_filtered.bam does not exist. Trimming read files."
					${BBMAPDIR}/bbduk.sh -${KMEM} in=$INDIR/${array[i]} in2=$INDIR/${array[i+1]} out=$o1 out2=$o2 \
					ref=${BBMAPDIR}/resources/adapters.fa ktrim=r k=23 tpe tbo qtrim=r trimq=${BQUAL} minlength=${MINLENGTH} minavgquality=${BQUAL} threads=${THREADZ} stats=$OUTDIR/trimmed/${b}_trimstats.txt
#					gzip $o1; gzip $o2
				fi
			fi
		done
	fi
	
	if [[ "$single" =~(on)$ ]] ; then
		echo "ALERT: Proceeding to trimming of single-end reads."
		for (( i=0; i<${#array[@]} ; i+=1 )) ; do
			b=${array[i]%%.*}
			c=${array[i+1]%%.*}
			o1=${OUTDIR}/trimmed/${b}.trim_1

			if  [[ -f $o1 ]] || [[ -f $o1.gz ]] ; then
				echo "ALERT: $o1  exists; skipping"
				echo ${b} >> $OUTDIR/${OUTPOP}_names.txt
			fi

			if [[ ! -f $o1 ]] && [[ ! -f $o1.gz ]] ; then
				echo ${b} >> $OUTDIR/${OUTPOP}_names.txt
				if [[ -f ${OUTDIR}/BAM/${b}_filtered.bam ]] ; then
					echo "ALERT: ${b}_filtered.bam exists and will not be re-trimmed or re-aligned."
				else
					#Add modifications below if needed#
					${BBMAPDIR}/bbduk.sh -${KMEM} in=$INDIR/${array[i]} out=$o1  \
					ref=${BBMAPDIR}/resources/adapters.fa ktrim=r k=23 tpe tbo qtrim=r trimq=${BQUAL} minlength=${MINLENGTH} minavgquality=${BQUAL} threads=${THREADZ} stats=$OUTDIR/trimmed/${b}_trimstats.txt
#					gzip $o1
				fi
			fi
		done
	fi

		#remove weird characters from names
		awk '{ sub("\r$", ""); print $0 }'  $OUTDIR/${OUTPOP}_names.txt > $OUTDIR/${OUTPOP}_names.tmp && mv $OUTDIR/${OUTPOP}_names.tmp $OUTDIR/${OUTPOP}_names.txt
	
	#Check that trimmed folder was written to
	if [ -z "$(ls -A ${OUTDIR}/trimmed)" ]; then
#	if [ -z "$(ls -A ${OUTDIR}/trimmed)" ] && [ -z "$(ls -A ${OUTDIR}/BAM)" ]; then
		echo "WARNING: No trimmed files were produced"
		if [ -z "$(ls -A ${OUTDIR}/trimmed)" ] && [ -z "$(ls -A ${OUTDIR}/BAM)" ]; then
			echo "WARNING: No trimmed files were produced AND BAM folder is empty; exiting"
			exit 1
		fi
	else
		if [ ! -z "$(ls -A ${OUTDIR}/trimmed)" ]; then
			ls ${OUTDIR}/trimmed/*.trim_[12] | parallel -j ${THREADZ} "gzip {}" 
		fi
	fi

	if [[ "$QUALREPORT" =~(on)$ ]] && [ ! -z "$(ls -A ${OUTDIR}/trimmed)" ] ; then
	# Quality report of trimmed files (FASTQC). This runs in the background as alignments are produced
		if [ ! -d "$OUTDIR/quality" ]; then
			mkdir $OUTDIR/quality
		fi
			echo "ALERT: FASTQC has started running in background"
			for i in ${OUTDIR}/trimmed/*.trim_*; do
				ni=$(echo $i | sed 's/\.gz//' | awk -F/ '{print $NF}') 
				#echo $ni
				if [[ ! -f ${OUTDIR}/quality/${ni}_fastqc.zip ]] ; then
					${FASTQC} -q -o ${OUTDIR}/quality $i 
				fi
			done &
	fi


#############################Read Mapping section#############################
#Alignment and duplicate removal (BWA, SAMBLASTER, SAMTOOLS)
##!ADDITIONAL ALIGNMENT PARAMETERS CAN BE ADDED BELOW AFTER "${BWA} mem" !##
#CHECK bwa mem DOCUMENTATION FOR ADDITIONAL OPTIONS AND MODIFY LINE BELOW#
	
	#Make folder for BAMS
	echo "ALERT: Beginning BWA mem at $(date)"
	
	ITER=1
	if [ ! -d "$OUTDIR/BAM" ]; then
		mkdir $OUTDIR/BAM
	fi
	
	for i in $(cat $OUTDIR/${OUTPOP}_names.txt); do
		if [[ -f ${OUTDIR}/BAM/${i}_aligned.bam ]] || [[ -f ${OUTDIR}/BAM/${i}_sorted.bam ]] || [[ -f ${OUTDIR}/BAM/${i}_filtered.bam ]]; then
			echo "ALERT: ${i}_aligned.bam or ${i}_sorted.bam exists or ${i}_filtered.bam; skipping"
		else
			if [[ "$single" =~(off)$ ]] ; then
				TRIM1=(`ls "${OUTDIR}/trimmed/$i.trim_1"* `)
				TRIM2=(`ls "${OUTDIR}/trimmed/$i.trim_2"* `)
				
				if [[ "$SPLITDISC" =~(on)$ ]] ; then 
					echo "ALERT: Aligning with discordant and split read production"
					#Makes unique lane and pop name for each library
					a="@RG\tID:SMP"
					b="\tPL:ILLUMINA\tLB:Pooled\tSM:"
					c=$a$ITER$b$i
					nice -n 19 ${BWA} mem -M -t $THREADZ -R $c \
					"${GENOME}" \
					"${TRIM1}" "${TRIM2}" | ${SAMBLASTER} -M -r -d "${OUTDIR}/BAM/$i.disc.sam" -s "${OUTDIR}/BAM/$i.split.sam"  | ${SAMTOOLS} view -Sb -q ${MAPQ} - > "${OUTDIR}/BAM/${i}_aligned.bam"
				fi	
		
				if [[ "$SPLITDISC" =~(off)$ ]] ; then 
					echo "ALERT: Aligning without discordant and split read production"
					#Makes unique lane and pop name for each library
					a="@RG\tID:SMP"
					b="\tPL:ILLUMINA\tLB:Pooled\tSM:"
					c=$a$ITER$b$i
					nice -n 19 ${BWA} mem -M -t $THREADZ -R $c \
					"${GENOME}" \
					"${TRIM1}" "${TRIM2}" | ${SAMBLASTER} -M -r | ${SAMTOOLS} view -Sb -q ${MAPQ} - > "${OUTDIR}/BAM/${i}_aligned.bam"
				fi								
			fi

			if [[ "$single" =~(on)$ ]] ; then
				TRIM1=(`ls "${OUTDIR}/trimmed/$i.trim_1"*`)
				if [[ "$SPLITDISC" =~(on)$ ]] ; then 
					echo "ALERT: Cannot perform discordant/split-end analyses on single-end reads! ignoring"
				fi	
				echo "ALERT: Aligning without discordant and split read production"
				#Makes unique lane and pop name for each library
				a="@RG\tID:SMP"
				b="\tPL:ILLUMINA\tLB:Pooled\tSM:"
				c=$a$ITER$b$i
				nice -n 19 ${BWA} mem -M -t $THREADZ -R $c \
				"${GENOME}" \
				"${TRIM1}" | ${SAMBLASTER} -M -r | ${SAMTOOLS} view -Sb -q ${MAPQ} - > "${OUTDIR}/BAM/${i}_aligned.bam"
			fi
		fi

		ITER="$(($ITER + 1))"
	done

	echo "ALERT: Finished BWA mem at $(date)"
		#reports will have read alignment results for each bam file
		if [ ! -d "$OUTDIR/reports" ]; then
			mkdir $OUTDIR/reports
		fi

# Sorting BAM files and filtering unpaired reads (PICARDTOOLS, SAMTOOLS)
	for i in $(cat $OUTDIR/${OUTPOP}_names.txt); do
		if  [[ -f ${OUTDIR}/BAM/${i}_sorted.bam ]] || [[ -f ${OUTDIR}/BAM/${i}_filtered.bam ]] ; then
			echo "ALERT: ${i}_sorted.bam or filtered.bam exists; skipping"
		else
			echo "ALERT: Picardtools started sorting ${OUTDIR}/BAM/${i}_aligned.bam at $(date) "
			nice -n 19 java -XX:ParallelGCThreads=$THREADZ -${KMEM} -Djava.io.tmpdir=$OUTDIR/tmp -jar ${PICARDTOOLS} SortSam I= ${OUTDIR}/BAM/${i}_aligned.bam O= ${OUTDIR}/BAM/${i}_sorted.bam VALIDATION_STRINGENCY=SILENT QUIET=true SO=coordinate TMP_DIR=${OUTDIR}/tmp 
		fi
	done
	echo "ALERT: Picardtools finished sorting all BAMS at $(date) " 

	#Making reports for sorted BAM files
	echo "ALERT: Alignment reports started at $(date) for new samples"
	for i in $(cat $OUTDIR/${OUTPOP}_names.txt); do
		if [[ ! -f ${OUTDIR}/reports/${i}_aln_report.txt ]] ; then
				if [[ -f ${OUTDIR}/BAM/${i}_filtered.bam ]] ; then
					echo "ALERT: ${i}_filtered.bam exists; alignment report skipped."
				else
					nice -n 19 ${SAMTOOLS} flagstat ${OUTDIR}/BAM/${i}_sorted.bam > ${OUTDIR}/reports/${i}_aln_report.txt
				fi
		fi
	done 

#Create semaphore to parallel run at specified 'THREADZ'
# 	open_sem(){
# 			mkfifo pipe-$$
# 			exec 3<>pipe-$$
# 			rm pipe-$$
# 			local i=$1
# 			for((;i>0;i--)); do
# 				printf %s 000 >&3
# 			done
# 	}
# 	run_with_lock(){
# 	local x
# 	read -u 3 -n 3 x && ((0==x)) || exit $x
# 	(
# 	"$@" 
# 	printf '%.3d' $? >&3
# 	)&
# 	}
# 
# #Filtering bams; removing junk. Runs in parallel.
# ##! ADDITIONAL BAM FILTER PARAMETERS CAN BE ADDED BELOW AFTER "${SAMTOOLS} view " !##
# 		task(){
# 			if [[ -f ${OUTDIR}/BAM/${i}_filtered.bam ]] ; then
# 				echo "ALERT: ${i}_filtered.bam exists; skipping"
# 			else
# 				if [[ "$single" =~(off)$ ]] ; then
# 					echo "ALERT: samtools is filtering ${OUTDIR}/BAM/${i}_sorted.bam at $(date) "
# 					nice -n 19 ${SAMTOOLS} view -b -F 0x04 -f 0x02 ${OUTDIR}/BAM/${i}_sorted.bam > ${OUTDIR}/BAM/${i}_filtered.bam 
# 				else
# 					cp ${OUTDIR}/BAM/${i}_sorted.bam ${OUTDIR}/BAM/${i}_filtered.bam 
# 				fi
# 			fi
# 			if [[ ! -f ${OUTDIR}/BAM/${i}_filtered.bam ]] && [[ "$single" =~(on)$ ]] ; then
# 
# 				echo "ALERT: samtools is filtering ${OUTDIR}/BAM/${i}_sorted.bam at $(date) "
# 				nice -n 19 ${SAMTOOLS} view -b ${OUTDIR}/BAM/${i}_sorted.bam > ${OUTDIR}/BAM/${i}_filtered.bam 
# 			fi
# 			}
# 
# 		N=$THREADZ
# 		open_sem $N
# 			for i in $(cat $OUTDIR/${OUTPOP}_names.txt); do
# 				run_with_lock task $i
# 			done ; wait

		rm $OUTDIR/filterlist &> /dev/null 
		for i in $(cat $OUTDIR/${OUTPOP}_names.txt); do
			if [[ -f ${OUTDIR}/BAM/${i}_filtered.bam ]] ; then
				echo "ALERT: ${i}_filtered.bam exists; skipping BAM filtering"
			else
				if [[ "$single" =~(off)$ ]] ; then
					echo "ALERT: samtools is filtering ${OUTDIR}/BAM/${i}_sorted.bam at $(date) "
					#nice -n 19 ${SAMTOOLS} view -b -F 0x04 -f 0x02 ${OUTDIR}/BAM/${i}_sorted.bam > ${OUTDIR}/BAM/${i}_filtered.bam 
					echo "${i}" >> $OUTDIR/filterlist
				else
					cp ${OUTDIR}/BAM/${i}_sorted.bam ${OUTDIR}/BAM/${i}_filtered.bam 
				fi
			fi
			if [[ ! -f ${OUTDIR}/BAM/${i}_filtered.bam ]] && [[ "$single" =~(on)$ ]] ; then
				echo "ALERT: samtools is filtering ${OUTDIR}/BAM/${i}_sorted.bam at $(date) "
				#nice -n 19 ${SAMTOOLS} view -b ${OUTDIR}/BAM/${i}_sorted.bam > ${OUTDIR}/BAM/${i}_filtered.bam 
				echo "${i}" >> filterlist
			fi
		done

		if [[ "$single" =~(off)$ ]] && [[ -f $OUTDIR/filterlist ]] ; then
			cat $OUTDIR/filterlist | parallel -j ${THREADZ} "${SAMTOOLS} view -b -F 0x04 -f 0x02 ${OUTDIR}/BAM/{}_sorted.bam > ${OUTDIR}/BAM/{}_filtered.bam"
		fi
		if [[ "$single" =~(on)$ ]] && [[ -f $OUTDIR/filterlist ]] ; then
			cat $OUTDIR/filterlist | parallel -j ${THREADZ} "${SAMTOOLS} view -b ${OUTDIR}/BAM/{}_sorted.bam > ${OUTDIR}/BAM/{}_filtered.bam"
		fi	
		
#		rm $OUTDIR/filterlist

		echo "ALERT: Samtools finished filtering all BAMs at $(date) "

	#Remove initial bams to reduce storage (pretty useless at this point)
	for i in $(cat $OUTDIR/${OUTPOP}_names.txt); do
		if [[ -f ${OUTDIR}/BAM/${i}_aligned.bam ]] ; then
			rm ${OUTDIR}/BAM/${i}_aligned.bam
		fi
	done	



if [ ! -z $ALIGNONLY ] &&  [[ "$ALIGNONLY" =~(on)$ ]] ; then
		echo "ALERT: PPalign Alignment Only completed at $(date) "
		exit
fi



#############################BAM merging by Class/Population
#Combine BAMS into specified populations. Uses parallel
			POPZ=$(cat ${OUTDIR}/pops/${OUTPOP}_files_for_pops.txt)	
			declare -a farray=($POPZ)
			echo "ALERT: Samtools is combining ${#farray[@]} populations into BAMS at $(date)"



	if  [[ "$skipmerge" =~(on)$ ]] ; then
		echo "ALERT: Skipping population merge; identical poplist exists. Delete poplist.txt and population bams to make new populations "
	else

# 			open_sem(){
# 			mkfifo pipe-$$
# 			exec 3<>pipe-$$
# 			rm pipe-$$
# 			local i=$1
# 			for((;i>0;i--)); do
# 				printf %s 000 >&3
# 			done
# 			}
# 			run_with_lock(){
# 			local x
# 			read -u 3 -n 3 x && ((0==x)) || exit $x
# 			(
# 			"$@" 
# 			printf '%.3d' $? >&3
# 			)&
# 			}
# 
# 		task() {
# 			declare -i pnuM=$(wc -l ${OUTDIR}/pops/${i}.txt |  cut -f1 -d' ')
# 			declare -a barray=$(awk '{print "'${OUTDIR}/BAM/'"$0"_filtered.bam"}' ${OUTDIR}/pops/${i}.txt)
# 			#If population consists of one library, simply duplicate the bam file. 
# 			if [ $pnuM -lt 2 ] ; then 
# 				cp $barray ${OUTDIR}/BAM/${i}.bam
# 			else
# 				echo "ALERT: Samtools is merging $pnuM bams into populations ${i}.bam at $(date) "
# 				samtools merge -f ${OUTDIR}/BAM/${i}.bam ${barray[@]}
# 			fi
# 		}
# 
# 		N=$THREADZ
# 		open_sem $N
# 			for i in "${farray[@]}" ; do
# 				run_with_lock task $i
# 			done ; wait

		rm $OUTDIR/mergelist &> /dev/null
		for i in "${farray[@]}" ; do
			declare -i pnuM=$(wc -l ${OUTDIR}/pops/${i}.txt |  cut -f1 -d' ')
			declare -a barray=$(awk '{print "'${OUTDIR}/BAM/'"$0"_filtered.bam"}' ${OUTDIR}/pops/${i}.txt)
			#If population consists of one library, simply duplicate the bam file. 
			if [ $pnuM -lt 2 ] ; then 
				cp $barray ${OUTDIR}/BAM/${i}.bam
			else
				echo "ALERT: Samtools is merging $pnuM bams into population ${i}.bam at $(date) "
				paste -d'\0' <(for ((j = 1; j <= $pnuM; j++)); do echo "${OUTDIR}/BAM/"; done) <(paste -d'\0' ${OUTDIR}/pops/${i}.txt <(for ((j = 1; j <= $pnuM; j++)); do echo "_filtered.bam"; done) ) > $OUTDIR/BAM/"$i".mergelist
				echo "$i" >> $OUTDIR/BAM/mergelist.pops
			fi

		done

		cat $OUTDIR/BAM/mergelist.pops | parallel -j ${THREADZ} "samtools merge -f ${OUTDIR}/BAM/{}.bam -b $OUTDIR/BAM/{}.mergelist"

#		rm $OUTDIR/BAM/mergelist.pops &> /dev/null
#		rm $OUTDIR/BAM/"$i".mergelist &> /dev/null
		echo "ALERT: Samtools finished merging all bams at $(date) "

	fi


#############################Creating mpileup including every site covered, retains only a few columns, runs in background
	echo "ALERT: beginning stats mpileup creation"
	COMBINE=$(sed 's,^,'$OUTDIR/BAM/', ; s/$/.bam/' $OUTDIR/pops/${OUTPOP}_files_for_pops.txt)
	echo "ALERT: See ${OUTPOP}_files_for_pops.txt for population merge order"
	

	if [[ -f  $OUTDIR/${OUTPOP}_stats.mpileup ]] ; then
		echo "ALERT: ${OUTPOP}_stats.mpileup already exists, skipping mpileup creation "
	else

	echo "ALERT: stat mpileup is being created at  $(date)"
	
		#Get Columns of mpileup and determine population number
		declare -i NCOL=$(($parraynum * 3 + 3))
		# Start at 3, increment by NPOPS until NCOL, i.e, get column of each pop
		POPSEQ=$(seq 4 3 $NCOL)
		#Select columns to run analyses on 
		declare -a arr2=($POPSEQ)
		#Subset
		ONE=$(echo $POPSEQ | sed -r 's/([^ ]+)/$\1/g')
		TWO=$(echo $ONE | awk '$NF=$NF "}" ' OFS=",")			#echo $TWO
		THREE='{print $1,$2,$3,'
		FOUR=$(echo ${THREE}${TWO})	
		${SAMTOOLS} mpileup -f ${GENOME} -B $COMBINE | awk "$FOUR"  > $OUTDIR/${OUTPOP}_stats.mpileup &  PIDMIX=$!
	fi


	if [ ! -z $STATSONLY ] &&  [[ "$STATSONLY" =~(on)$ ]] ; then
			wait $PIDMIX
			echo "ALERT: PPalign Alignment and Stats Only completed at $(date) "
			exit
	fi




#############################SNP Calling section#############################
#Call SNPs and print variant sites

	if [[ -f $OUTDIR/${OUTPOP}_variants.txt ]]; then
		echo "ALERT: ${OUTPOP}_variants.txt already exists, skipping mpileup creation and variant calling"
	else

		echo "ALERT: beginning mpileup creation and variant calling"

	#Make the filter directory if it does not exist 
		if [ ! -d "$OUTDIR/filters" ]; then
			mkdir $OUTDIR/filters
		fi

#-z means 'has zero length'		
		if [ ! -z $USEVCF ]; then
#############################Use substitute "vcf" for further processing, produce variants file from this
			echo "ALERT: Using $USEVCF for SNPs, skipping variant calling"
			echo "ALERT: $USEVCF will NOT be filtered before use"
			# ${BCFTOOLS} view  -i  'MIN(DP)>'$MINDP' & MIN(QUAL)>'$SNPQ' ' $USEVCF  > $OUTDIR/${OUTPOP}_Qualtemp.VCF
				# declare -i after=$(grep -v "^#" $OUTDIR/${OUTPOP}_Qualtemp.VCF | wc -l) 
				# declare -i lost=$(($before-$after))
				# echo "ALERT: $lost SNPs removed due to QUAL < $SNPQ and total DP < $MINDP "
			# ${BCFTOOLS} view  -i  'MAF[0]> '$MAF' ' $OUTDIR/${OUTPOP}_Qualtemp.VCF > $OUTDIR/${OUTPOP}.VCF
				mv $OUTDIR/${USEVCF} $OUTDIR/${OUTPOP}.VCF
				declare -i afterm=$( grep -v "^#" $OUTDIR/${OUTPOP}.VCF | wc -l) 
#				declare -i lostm=$(($after-$afterm))
#				echo "ALERT: Additional $lostm SNPs removed due to global MAF < $MAF "
				echo "ALERT: $afterm total SNPs in $USEVCF, now called ${OUTPOP}.VCF"
				grep -v "^#"  $OUTDIR/${OUTPOP}.VCF | awk '{gsub("=|;", "\t", $0); print;}' \
						| awk '{print $1,$2,$6,$9}'  >  $OUTDIR/${OUTPOP}_variants.txt
		#end 'use VCF'
		else

#############################SNP Calling
			if [ -z $PARALLEL ] ||  [[ "$PARALLEL" =~(off)$ ]] ; then
			###non-parallel SNP calling (bcftools only uses additional processors for compressing output)
					echo -e "\nALERT: bcftools is calling variants at $(date) WITHOUT parallel"
					echo -e "This might take awhile...longer for no parallel :-/\n"
					#${SAMTOOLS} mpileup -uf ${GENOME} -B $COMBINE |  ${BCFTOOLS} call --threads ${THREADZ} -mv -Ov >  $OUTDIR/${OUTPOP}_full.VCF 
					${BCFTOOLS} mpileup -Ou -f ${GENOME} -B $COMBINE |  ${BCFTOOLS} call -mv -Ov >  $OUTDIR/${OUTPOP}_full.VCF 

					declare -i before=$(grep -v "^#" $OUTDIR/${OUTPOP}_full.VCF | wc -l)
					echo "ALERT: $before total SNPS called without filters"
				#############################SNP filtering, non-parallel
				${BCFTOOLS} view  -i  'MIN(DP)>'$MINDP' & MIN(QUAL)>'$SNPQ' ' $OUTDIR/${OUTPOP}_full.VCF  > $OUTDIR/${OUTPOP}_Qualtemp.VCF
					declare -i after=$(grep -v "^#" $OUTDIR/${OUTPOP}_Qualtemp.VCF | wc -l) 
					declare -i lost=$(($before-$after))
					echo "ALERT: $lost SNPs removed due to QUAL < $SNPQ and total DP < $MINDP "
				${BCFTOOLS} view  -i  'MAF[0]> '$MAF' ' $OUTDIR/${OUTPOP}_Qualtemp.VCF > $OUTDIR/${OUTPOP}.VCF
					declare -i afterm=$( grep -v "^#" $OUTDIR/${OUTPOP}.VCF | wc -l) 
					declare -i lostm=$(($after-$afterm))
					echo "ALERT: Additional $lostm SNPs removed due to global MAF < $MAF "
					echo "ALERT: $afterm total SNPs retained after SNP calling"

				if [ "$afterm" -gt 0 ] ; then
				#if [[ -f $OUTDIR/${OUTPOP}.VCF ]]; then
					rm $OUTDIR/${OUTPOP}_Qualtemp.VCF
					rm $OUTDIR/${OUTPOP}_full.VCF
				else
					echo "ALERT: variant calling failed: VCF file not found (non-parallel)"
					exit 1
				fi
			
			#end of 'cal SNPs without parallel'
			elif [ ! -z $PARALLEL ] &&  [[ "$PARALLEL" =~(on)$ ]]; then
			###parallel SNP call by 'chunking' genome
				
				if [[ -f $OUTDIR/refint0.txt.vcf.gz_Qualtemp.vcf.gz ]] || [[ -f $OUTDIR/refint0.txt.vcf ]]; then
					echo -e "\nALERT: temporary VCFs found; will not repeat variant calling "
				else
					echo -e "\nALERT: bcftools is calling variants at $(date) USING parallel "
					echo -e "This might take awhile...but shorter bc parallel ;-)\n"
					
					#index bam files
					ls $OUTDIR/BAM/pop*.bam > $OUTDIR/bamlist
					BAMS=(` cat $OUTDIR/bamlist ` )
					LENb=(` wc -l $OUTDIR/bamlist `)
					LENb=$(($LENb - 1 ))
					rm $OUTDIR/bailist &> /dev/null
					for ((i = 0; i <= $LENb; i++));
					do
						if [[ ! -f ${BAMS[i]}.bai ]]; then
							echo ${BAMS[i]} >> $OUTDIR/bailist
						fi
					done

					if [[ -f $OUTDIR/bailist ]]; then
						echo -e "\nALERT: Indexing population BAM files"
						cat $OUTDIR/bailist | parallel -j ${THREADZ} "samtools index {}"
					fi

					#create intervals for each chromosome and scaffold
					cut -f1-2 ${GENOME}.fai | awk '$1 = $1 FS "1"' | tr ' ' '\t' > $OUTDIR/reference_intervals.txt				
					
					#if scaffold prefix provided; evenly distribute chromosomes as well as scaffolds
					if [ -z "$SCAHEAD" ]; then
						#shuffle intervals to distribute across threads
						shuf $OUTDIR/reference_intervals.txt > $OUTDIR/reference_intervals_shuf.txt
						#split the intervals into files, one per thread
						rm $OUTDIR/refint*txt &> /dev/null
						split -e --numeric-suffixes=1 -n l/${THREADZ} --additional-suffix='.txt' $OUTDIR/reference_intervals_shuf.txt refint
						mv refint*txt $OUTDIR
					else
						#get chromosome interval names
						grep -v "^${SCAHEAD}" $OUTDIR/reference_intervals.txt >  $OUTDIR/chromosome_reference_intervals.txt
						#shuffle chromosome intervals to distribute across threads
						shuf $OUTDIR/chromosome_reference_intervals.txt > $OUTDIR/chromosome_reference_intervals_shuf.txt
						rm chr_refint*.txt &> /dev/null
						split -e --numeric-suffixes=1 -n l/${THREADZ} --additional-suffix='.txt' $OUTDIR/chromosome_reference_intervals_shuf.txt chr_refint
						#get scaffold interval names
						cat $OUTDIR/reference_intervals.txt | grep "^${SCAHEAD}"  >  $OUTDIR/scaffold_reference_intervals.txt
						#shuffle scaffold intervals to distribute across threads
						shuf $OUTDIR/scaffold_reference_intervals.txt > $OUTDIR/scaffold_reference_intervals_shuf.txt
						rm scaf_refint*.txt &> /dev/null
						split -e --numeric-suffixes=1 -n l/${THREADZ} --additional-suffix='.txt' $OUTDIR/scaffold_reference_intervals_shuf.txt scaf_refint

						#concatenate chromosome and scaffold interval files using loop
						#determine loop length by larger number of files
						#this is necessary for when #chromosome or #scaffold files is less than #threads specified
						CHR=( ` ls chr_refint*txt ` )
						LENchr=( ` ls chr_refint*txt | wc -l ` )
						SCAF=( ` ls scaf_refint*txt ` )
						LENscaf=( ` ls scaf_refint*txt | wc -l ` )

						if [ "$LENchr" -ge "$LENscaf" ] ; then
							LEN=$LENchr
							echo "ALERT: number of intervals will be "$LEN
						elif [ "$LENchr" -lt "$LENscaf" ] ; then
							LEN=$LENscaf
							echo "ALERT: number of intervals will be "$LEN
						else
							echo "ALERT: could not determine the number of interval files"
							exit 1
						fi


						LEN=$(($LEN - 1))

						for ((j = 0; j <= $LEN; j++));
						do
							k=$(($j+1))
							cat ${CHR[j]} ${SCAF[j]} > $OUTDIR/refint0$k.txt
						done

						rm chr_refint*
						rm scaf_refint*
					#end 'make reference intervals with or without scaffold prefix'
					fi



					#collect interval files, feed to parallel which implements bcftools
					ls $OUTDIR/refint*.txt > $OUTDIR/namelist
					LEN=( ` ls $OUTDIR/refint*.txt | wc -l ` )
					
					#set number of processors for SNP calling to lesser of specified threads or number of reference intervals
					if [ "$LEN" -ge "${THREADZ}" ] ; then
						THREADZb=${THREADZ}
						echo "ALERT: number of threads for SNP calling will be "$THREADZb
					elif [ "$LEN" -lt "${THREADZ}" ] ; then
						THREADZb=$LEN
						echo "ALERT: number of threads for SNP calling will be "$THREADZb
					else
						echo "ALERT: could not determine the number of threads for SNP calling"
						exit 1
					fi

					echo -e "\nALERT: calling variants in parallel "
					cat $OUTDIR/namelist | parallel -j ${THREADZb} "${BCFTOOLS} mpileup -Ou -f ${GENOME} -R{} -B -b $OUTDIR/bamlist | ${BCFTOOLS} call -mv -Ov > {}.vcf"
#					wc -l $OUTDIR/refint*.txt.vcf

					##on occasion, more processors are specified/requested than the system has to devote
					##this results in failed SNP calling for those intervals, and empty VCF files
					##here we check that the VCF files are not empty (0 lines), and first:
					##collect empty files for reference intervals and attempt to call ONCE more in parallel; then:
					##loop over empty files for reference intervals and attempt to call in SERIAL (sequentially); then:
					##if any are still empty, EXIT

					##loop over all reference files to check that output is not empty; if it is:
					##make new list; try to call in PARALLEL
					rm $OUTDIR/namelist2 &> /dev/null
					for ((j = 1; j <= $LEN; j++));
					do
						#k=$(($j+1))
						k=$j
						NUM=( ` wc -l $OUTDIR/refint0$k.txt.vcf `)
						if [ "$NUM" -eq 0 ] ; then
							echo "$OUTDIR/refint0$k.txt.vcf" | sed 's/\.vcf//' >> $OUTDIR/namelist2
						fi
					done

					##if any files were found to be empty, try to call them again in PARALLEL
					if [[ -f $OUTDIR/namelist2 ]]; then
 						echo -e "\n\nALERT: Some reference intervals were not scored in parallel"
 						echo -e "\n\nALERT: Attempting to score again in parallel"
						cat $OUTDIR/namelist2 | parallel -j ${THREADZb} "${BCFTOOLS} mpileup -Ou -f ${GENOME} -R{} -B -b $OUTDIR/bamlist | ${BCFTOOLS} call -mv -Ov > {}.vcf"
					fi
#					wc -l $OUTDIR/refint*.txt.vcf


					##loop over all reference files to check that output is not empty; if it is:
					##make new list again; try to call in SERIAL
					rm $OUTDIR/namelist3 &> /dev/null
					for ((j = 1; j <= $LEN; j++));
					do
						#k=$(($j+1))
						k=$j
						NUM=( ` wc -l $OUTDIR/refint0$k.txt.vcf `)
						if [ "$NUM" -eq 0 ] ; then
							echo "$OUTDIR/refint0$k.txt.vcf" | sed 's/\.vcf//' >> $OUTDIR/namelist3
						fi
					done

					##if any files were found to be empty, try to call them again in SERIAL
					if [[ -f $OUTDIR/namelist3 ]]; then
 						echo -e "\n\nALERT: Some reference intervals were not scored in parallel a second time"
 						echo -e "\n\nALERT: Attempting to score again in SERIAL (NOT in parallel)"
						for i in $(cat $OUTDIR/namelist3);
						do
							echo "Scoring variants in ${i}.\n"
							${BCFTOOLS} mpileup -Ou -f ${GENOME} -R${i} -B -b $OUTDIR/bamlist | ${BCFTOOLS} call -mv -Ov > ${i}.vcf 
						done
					fi
#					wc -l $OUTDIR/refint*.txt.vcf


					##loop over all reference files to check that output is not empty; if it is:
					for ((j = 0; j <= $LEN; j++));
					do
						k=$(($j+1))
						NUM=( ` wc -l $OUTDIR/refint0$k.txt.vcf `)
						if [ "$NUM" -eq 0 ] ; then
							echo -e "\n\nOne or more temporary VCF files is empty; SNP calling failed!"
							exit 1
						fi
					done

				#end 'make intervals, call SNPs unless raw refint VCFs found'	
				fi
				
				##concatenate VCFs, sort alphanumerically
				echo -e "\nALERT: SNP calling completed; proceeding to filtering, concatenating and sorting variants "




				##############################SNP filtering, parallel
				
				if [[ -f $OUTDIR/refint01.txt.vcf.gz_Qualtemp.vcf.gz ]]; then
					echo -e "\nALERT: Qualtemp VCFs found; will not repeat quality filtering "
				else
					echo -e "\nALERT: Qualtemp VCFs NOT found; now quality filtering "
					#namelist contains path
					ls $OUTDIR/refint*.txt.vcf > $OUTDIR/namelist
					cat $OUTDIR/namelist | parallel -j ${THREADZb} "gzip {}"

					#namelist contains path
					ls $OUTDIR/refint*.txt.vcf.gz > $OUTDIR/namelist				
					cat $OUTDIR/namelist | parallel -j ${THREADZb} "gunzip -c {} | ${BCFTOOLS} view  -i  'MIN(DP)>'$MINDP' & MIN(QUAL)>'$SNPQ' & MAF[0]> '$MAF' ' -  | gzip > {}_Qualtemp.vcf.gz"

#					if [[ -f $OUTDIR/refint01.txt.vcf.gz_Qualtemp.vcf.gz ]] ; then
#						rm $OUTDIR/refint*.txt.vcf.gz
#					fi
				#end 'make Qualtemp vcf if not already present'
				fi

				if [ ! -f $OUTDIR/vcfheader.txt ]; then
					gunzip -c $OUTDIR/refint01.txt.vcf.gz_Qualtemp.vcf.gz | grep "^#" > $OUTDIR/vcfheader.txt
#					gunzip -c $OUTDIR/refint01.txt.vcf.gz_Qualtemp.vcf.gz | grep "^#" | head -n1 > $OUTDIR/vcfheader.txt
#					gunzip -c $OUTDIR/refint01.txt.vcf.gz_Qualtemp.vcf.gz | grep "^#" | tail -n1 >> $OUTDIR/vcfheader.txt
				fi

				#check for sorted.VCF and that it has header and SNP lines
				if [ -f $OUTDIR/sorted.VCF ]; then
					echo "ALERT: temp concatenated VCF found; checking"
					declare -i sorthead=$(grep "#" $OUTDIR/sorted.VCF | wc -l )
					declare -i sortnonhead=$(grep -v "#" $OUTDIR/sorted.VCF | wc -l )
					if [[ "$sorthead" -gt 0 && "$sortnonhead" -gt 0 ]]; then
						echo "ALERT: temp concatenated VCF contains header and non-header; proceeding"
					else
						echo "ALERT: temp concatenated VCF missing header or non-header; re-concatenating Qualtemp vcfs"
						rm $OUTDIR/sorted.VCF &> /dev/null
						cp $OUTDIR/vcfheader.txt sorted.VCF
						gunzip -c $OUTDIR/*_Qualtemp.vcf.gz | grep -v "^#" | sort -k1,1V -k2,2n -S 50% --parallel ${THREADZb} --temporary-directory=$OUTDIR/tmp >> $OUTDIR/sorted.VCF
					fi
				#end found sorted.VCF, check for header and non-header (SNP) lines
				else
				#begin re-make sorted.VCF from qualtemp VCFs if not found
					echo "ALERT: temp sorted VCF NOT found; concatenating"
					rm $OUTDIR/sorted.VCF &> /dev/null
					cp $OUTDIR/vcfheader.txt $OUTDIR/sorted.VCF
					gunzip -c $OUTDIR/*_Qualtemp.vcf.gz | grep -v "^#" | sort -k1,1V -k2,2n -S 50% --parallel ${THREADZb} --temporary-directory=$OUTDIR/tmp >> $OUTDIR/sorted.VCF
				#end 're-make sorted.VCF from Qualtemp'
				fi

				if [ ! -z "$SCAHEAD" ]; then
					cat $OUTDIR/vcfheader.txt <(grep -v "#" $OUTDIR/sorted.VCF | grep -v "^$SCAHEAD" ) <(grep -v "#" $OUTDIR/sorted.VCF | grep "^$SCAHEAD" ) > $OUTDIR/tmp/resorted.VCF
					mv $OUTDIR/tmp/resorted.VCF $OUTDIR/sorted.VCF
				fi


				declare -i vcflines=$(cat $OUTDIR/sorted.VCF | wc -l ) 
				declare -i head=$(cat $OUTDIR/vcfheader.txt | wc -l )
				declare -i qtlines=$(zcat $OUTDIR/*_Qualtemp.vcf.gz | grep -v "^#" | wc -l)
				declare -i target=$(($qtlines+$head))
				
				#remove temporary files if SNP calling has been successful; else die
				if [[ "$vcflines" -gt "$head" && "$vcflines" -eq "$target" ]] ; then
					mv $OUTDIR/sorted.VCF $OUTDIR/${OUTPOP}.VCF
					echo "ALERT: $qtlines total SNPs retained after SNP calling and filtering"
# 					rm $OUTDIR/bamlist
# 					rm $OUTDIR/bailist
# 					rm $OUTDIR/*reference_intervals.txt
# 					rm $OUTDIR/*reference_intervals_shuf.txt
# 					rm $OUTDIR/namelist
# 					rm $OUTDIR/namelist2
# 					rm $OUTDIR/namelist3
# 					rm $OUTDIR/refint*.txt
# 					rm $OUTDIR/vcfheader.txt
#					rm $OUTDIR/refint*.vcf.gz*
				else
					echo "ALERT: SNP calling, filtering, or sorting failed"
					exit 1
				fi

			#end of 'SNP calling with parallel'	
			fi

	#############################Make variants-only file (chr, pos, qual, depth/class)

			if [[ ! -f $OUTDIR/${OUTPOP}.VCF ]]; then
				echo "ALERT: VCF file missing: SNP calling, filtering, or sorting failed"
				exit 1
			fi

			if [[ ! -s $OUTDIR/${OUTPOP}.VCF ]] ; then
				echo "ERROR: VCF file empty; Check for sufficient memory"
				exit
			fi

			grep -v "^#"  $OUTDIR/${OUTPOP}.VCF | awk '{gsub("=|;", "\t", $0); print;}' \
							| awk '{print $1,$2,$6,$9}'  >  $OUTDIR/${OUTPOP}_variants.txt

			if [[ ! -f $OUTDIR/${OUTPOP}_variants.txt ]]; then
				echo "ALERT: ${OUTPOP}_variants.txt file missing: SNP calling, filtering, or sorting failed"
				exit 1
			fi
		#end of 'no VCF provided; create new, then variants file'
		fi

	#end of 'variants file found' or not
	fi
	
	#count how many variants are SNPs or indels
	awk '!/IDV|,/' $OUTDIR/${OUTPOP}_variants.txt > $OUTDIR/filters/${OUTPOP}_SNP_variants.txt
	awk '/IDV/'  $OUTDIR/${OUTPOP}_variants.txt > $OUTDIR/filters/${OUTPOP}_indel_variants.txt
				numSNP=$(wc -l $OUTDIR/filters/${OUTPOP}_SNP_variants.txt|  cut -f1 -d' ')
				numINDEL=$(wc -l $OUTDIR/filters/${OUTPOP}_indel_variants.txt |  cut -f1 -d' ')
	echo "ALERT: Of the remaining SNPs, there are $numSNP SNPs and $numINDEL INDels"
	echo "ALERT: Variant calling and filtering done at $(date) "
	




	if [[ -f  $OUTDIR/${OUTPOP}.mpileup ]] ; then
		echo "ALERT: ${OUTPOP}.mpileup already exists, skipping mpileup creation "
	else

		echo "ALERT: variant mpileup is being created at  $(date)"
		
#############################Creating mpileup but filtering to variant sites in (ultimately) the VCF tile, runs in background
			gawk 'NR==FNR{a[$1,$2]=$5;next} ($1,$2) in a{print $0, a[$1,$2]}' <(awk '{print $0}' $OUTDIR/${OUTPOP}_variants.txt) <(${SAMTOOLS} mpileup -f ${GENOME} -B $COMBINE)| awk '{gsub(" ","",$0); print;}' > $OUTDIR/${OUTPOP}.mpileup &  PIDIOS=$!

			#wait for stat mpileup and variant mpileup to be created
			wait $PIDIOS
			wait $PIDMIX
			echo "ALERT: Mpileups created at $(date)"
	fi
	
	if [[ ! -s $OUTDIR/${OUTPOP}.mpileup ]]  ; then
		echo "ERROR: Mpileup is empty. Check for sufficient memory during samtools mpileup"
		exit
	fi
		

	#Filter by indels and create sync format
		echo "ALERT: Identifying indel regions and creating sync format at $(date)"
			if [[ -f $OUTDIR/${OUTPOP}.sync ]] ; then
				echo "ALERT: Filtered sync file already exists, skipping this step "
			else
#############################take variant-only mpileup, use indel window, identify indel-masked regions, output as gtf
				perl ${POOL2}/indel_filtering/identify-indel-regions.pl --input $OUTDIR/${OUTPOP}.mpileup --indel-window $INWIN --output ${OUTDIR}/${OUTPOP}.gtffile --min-count 1
#############################take variant-only mpileup, make temporary sync file
				nice -n 19 java -ea -${KMEM} -Djava.io.tmpdir=$OUTDIR/tmp -jar ${POOL2}/mpileup2sync.jar --min-qual 1 --fastq-type sanger --input $OUTDIR/${OUTPOP}.mpileup --output $OUTDIR/${OUTPOP}_temp.sync --threads $THREADZ 
#############################filter temporary sync using indel-masking gtf, output final sync
				perl ${POOL2}/indel_filtering/filter-sync-by-gtf.pl --input ${OUTDIR}/${OUTPOP}_temp.sync --gtf ${OUTDIR}/${OUTPOP}.gtffile --output ${OUTDIR}/${OUTPOP}.sync
			echo "ALERT: Done identifying indel regions and creating sync format at $(date)"
				declare -i before=$(wc -l < $OUTDIR/${OUTPOP}_temp.sync )
				declare -i after=$(wc -l < $OUTDIR/${OUTPOP}.sync)
				declare -i lost="$(($before - $after))"
				lostp=$((100-100*$after/$before))
				echo "ALERT: With an indel window of $INWIN bp you lost $lost SNPs or $lostp % " 
				echo "SNPs remaining: $after"
				rm $OUTDIR/${OUTPOP}_temp.sync 
			fi				
			
		if [[ ! -s $OUTDIR/${OUTPOP}.sync  ]]  ; then
			echo "ERROR: Sync is empty. Check for sufficient memory during popoolation2 pl scripts"
			exit
		fi

#echo "ALERT: Creating allele frequency tables in R"

#REQUIRES R: Frequency and MAF
		echo "ALERT: Rscript called to calculate allele frequencies at $(date) "

#############################create a non-normalized frequency (fz) file from the (variant-only, no-indel) sync file, note 3+ allele sites across  and write to "poly" files
		if [[ -f ${OUTDIR}/${OUTPOP}.fz ]] ; then
					echo "ALERT: Frequency file ${OUTPOP}.fz  exists, skipping this analysis "
		else				
					echo "ALERT: Rscript called to create NON-normalized Frequency file"
					rin=${OUTDIR}/${OUTPOP}.sync
					rout=$OUTDIR/
					rout2=$OUTDIR/filters/
					Rscript $BASEDIR/rscripts/r_frequency.R $rin $rout $rout2 $MAF
					echo "ALERT: Frequency file ${OUTPOP}.fz and its counterparts created at $(date) "
		fi
				
		if [[ ! -s ${OUTDIR}/${OUTPOP}.fz  ]]  ; then
					echo "ERROR: Frequency table is empty. Check for R errors"
					exit
		fi







################normalization below
	
#############################create variant-only mpileup for each population/class with separate individual records, still includes indels, normalization only

# For each population listed in _files_for_pops, combined the filtered bam for individuals in that file as mpileups, while only keeping variants With parrelization
if [[ "$INDCONT" =~(on)$ ]] ; then


		if [ ! -d "$OUTDIR/inds" ]; then
			mkdir $OUTDIR/inds
		fi
		echo "ALERT: Creating individual mpileups contribution in ${#farray[@]} populations at $(date)"


# 	task() {
# 		if [[ -f $OUTDIR/inds/${file}.mpileup ]] ; then
# 			echo "ALERT: ${file}.mpileup already exists, skipping individual mpileup "
# 		else
# 			#Combine all individuals from same population into same mpileup using individual filtered.bam
# 			COMBINE=$(sed 's,^,'$OUTDIR/BAM/', ; s/$/_filtered.bam/' ${OUTDIR}/pops/$file.txt)
# 			gawk 'NR==FNR{a[$1,$2]=$5;next} ($1,$2) in a{print $0, a[$1,$2]}' <(awk '{print $0}' ${OUTDIR}/${OUTPOP}_variants.txt) <(${SAMTOOLS} mpileup -aa -f ${GENOME} -B $COMBINE) | awk '{gsub(" ","",$0); print;}' > $OUTDIR/inds/${file}.mpileup
# 		fi
# 	}
# 
# 		N=$THREADZ
# 		open_sem $N
# 			for file in $(cat $OUTDIR/pops/${OUTPOP}_files_for_pops.txt) ; do
# 				run_with_lock task $i
# 		done ; wait
# 		echo "ALERT: Individual mpileups created at $(date) "

		rm $OUTDIR/inds/mpilist &> /dev/null
		for file in $(cat $OUTDIR/pops/${OUTPOP}_files_for_pops.txt) ; do
			if [[ -f $OUTDIR/inds/${file}.mpileup ]] ; then
				echo "ALERT: ${file}.mpileup already exists, skipping individual mpileup "
			else
				echo ${file} >> $OUTDIR/inds/mpilist
				sed 's,^,'$OUTDIR/BAM/', ; s/$/_filtered.bam/' ${OUTDIR}/pops/${file}.txt > $OUTDIR/inds/"${file}".mpilist
			fi
		done ; 

		if [[ -f $OUTDIR/inds/mpilist ]] ; then

# 			open_sem(){
# 			mkfifo pipe-$$
# 			exec 3<>pipe-$$
# 			rm pipe-$$
# 			local i=$1
# 			for((;i>0;i--)); do
# 				printf %s 000 >&3
# 			done
# 			}
# 			run_with_lock(){
# 			local x
# 			read -u 3 -n 3 x && ((0==x)) || exit $x
# 			(
# 			"$@" 
# 			printf '%.3d' $? >&3
# 			)&
# 			}
# 
# 
# 			task() {
# 					#Combine all individuals from same population into same mpileup using individual filtered.bam
# 					COMBINE=$(sed 's,^,'$OUTDIR/BAM/', ; s/$/_filtered.bam/' ${OUTDIR}/pops/$file.txt)
# 					gawk 'NR==FNR{a[$1,$2]=$5;next} ($1,$2) in a{print $0, a[$1,$2]}' <(awk '{print $0}' ${OUTDIR}/${OUTPOP}_variants.txt) <(${SAMTOOLS} mpileup -aa -f ${GENOME} -B $COMBINE) | awk '{gsub(" ","",$0); print;}' > $OUTDIR/inds/${file}.mpileup
# 			}
# 
# 			N=$THREADZ
# 			open_sem $N
# 				for file in $(cat $OUTDIR/mpilist) ; do
# 					run_with_lock task $i
# 			done ; wait


#			echo "gawk 'NR==FNR{a[\$1,\$2]=\$5;next} (\$1,\$2) in a{print \$0, a[\$1,\$2]}' <(awk '{print \$0}' \${OUTDIR}/\${OUTPOP}_variants.txt) <(\${SAMTOOLS} mpileup -aa -f \${GENOME} -B -b \$OUTDIR/inds/{}.mpilist) | awk '{gsub(\" \",\"\",\$0); print;}' > $OUTDIR/inds/{}.mpileup"

			cat $OUTDIR/inds/mpilist | parallel -j ${THREADZ} "gawk 'NR==FNR{a[\$1,\$2]=\$5;next} (\$1,\$2) in a{print \$0, a[\$1,\$2]}' <(awk '{print \$0}' ${OUTDIR}/${OUTPOP}_variants.txt) <(${SAMTOOLS} mpileup -aa -f ${GENOME} -B -b $OUTDIR/inds/{}.mpilist) | awk '{gsub(\" \",\"\",\$0); print;}' > $OUTDIR/inds/{}.mpileup"

			#gawk 'NR==FNR{a[\$1,\$2]=\$5;next} (\$1,\$2) in a{print \$0, a[\$1,\$2]}' <(awk '{print \$0}' \${OUTDIR}/\${OUTPOP}_variants.txt) <(\${SAMTOOLS} mpileup -aa -f \${GENOME} -B -b \$OUTDIR/inds/{}.mpilist) | awk '{gsub(\" \",\"\",\$0); print;}' > $OUTDIR/inds/{}.mpileup

		fi
		
#		rm $OUTDIR/inds/*mpilist &> /dev/null
		echo "ALERT: Individual mpileups created at $(date) "













#############################create R input for each population/class from pop/class-specific variant-only mpileup created above (still includes indels), normalization only


# 	task() {
# 			if [[ -f $OUTDIR/inds/${file}_snp_stats.txt ]] ; then
# 				echo "ALERT: ${file}_RindSTATs.rin exists, skipping this file "
# 			else
# 
# 			#Get Columns of mpileup and determine population number
# 				declare -i NCOL=$(awk '{ print NF; exit }' $OUTDIR/inds/${file}.mpileup)
# 				declare -i NPOPS=($NCOL-3)/3
# 			# Start at 3, increment by NPOPS until NCOL, i.e, get column of each pop
# 				POPSEQ=$(seq 4 3 $NCOL)
# 				POPNUM=$(seq 1 $NPOPS)
# 			#Select columns to run analyses on 
# 				declare -a arr2=($POPSEQ)
# 			#Subset
# 				ONE=$(echo $POPSEQ | sed -r 's/([^ ]+)/$\1/g')
# 				TWO=$(echo $ONE | awk '$NF=$NF "}" ' OFS=",")
# 				THREE='{print $1,$2,'
# 				FOUR=$(echo ${THREE}${TWO})	
# 				awk "$FOUR" $OUTDIR/inds/${file}.mpileup > $OUTDIR/inds/${file}_RindSTATs.rin 
# 			fi
# 		}
# 
# 		N=$THREADZ
# 		open_sem $N
# 			for file in $(cat $OUTDIR/pops/${OUTPOP}_files_for_pops.txt) ; do
# 				run_with_lock task $i
# 		done ; wait

		rm $OUTDIR/inds/rinlist &> /dev/null
		for file in $(cat $OUTDIR/pops/${OUTPOP}_files_for_pops.txt) ; do
			if [[ -f $OUTDIR/inds/${file}_snp_stats.txt ]] ; then
				echo "ALERT: ${file}_RindSTATs.rin exists, skipping this file "
			else
				echo ${file} >> $OUTDIR/inds/rinlist
			fi
		done


		if [[ -f $OUTDIR/inds/rinlist ]] ; then

# 			open_sem(){
# 			mkfifo pipe-$$
# 			exec 3<>pipe-$$
# 			rm pipe-$$
# 			local i=$1
# 			for((;i>0;i--)); do
# 				printf %s 000 >&3
# 					done
# 			}
# 			run_with_lock(){
# 			local x
# 			read -u 3 -n 3 x && ((0==x)) || exit $x
# 			(
# 			"$@" 
# 			printf '%.3d' $? >&3
# 			)&
# 			}
# 
# 
# 			task() {
# 				#Get Columns of mpileup and determine population number
# 					declare -i NCOL=$(awk '{ print NF; exit }' $OUTDIR/inds/${file}.mpileup)
# 					declare -i NPOPS=($NCOL-3)/3
# 				# Start at 3, increment by NPOPS until NCOL, i.e, get column of each pop
# 					POPSEQ=$(seq 4 3 $NCOL)
# 					POPNUM=$(seq 1 $NPOPS)
# 				#Select columns to run analyses on 
# 					declare -a arr2=($POPSEQ)
# 				#Subset
# 					ONE=$(echo $POPSEQ | sed -r 's/([^ ]+)/$\1/g')
# 					TWO=$(echo $ONE | awk '$NF=$NF "}" ' OFS=",")
# 					THREE='{print $1,$2,'
# 					FOUR=$(echo ${THREE}${TWO})	
# 					awk "$FOUR" $OUTDIR/inds/${file}.mpileup > $OUTDIR/inds/${file}_RindSTATs.rin 
# 			}
# 
# 			N=$THREADZ
# 			open_sem $N
# 				for file in $(cat $OUTDIR/rinlist) ; do
# 					run_with_lock task $i
# 			done ; wait






			rm $OUTDIR/inds/rinlist.commands &> /dev/null
			for file in $(cat $OUTDIR/inds/rinlist) ; do
				#Get Columns of mpileup and determine population number
					declare -i NCOL=$(awk '{ print NF; exit }' $OUTDIR/inds/${file}.mpileup)
					declare -i NPOPS=($NCOL-3)/3
				# Start at 3, increment by NPOPS until NCOL, i.e, get column of each pop
					POPSEQ=$(seq 4 3 $NCOL)
					POPNUM=$(seq 1 $NPOPS)
				#Select columns to run analyses on 
					declare -a arr2=($POPSEQ)
				#Subset
					ONE=$(echo $POPSEQ | sed -r 's/([^ ]+)/\1/g')
					TWO=$(echo $ONE | awk '$NF=$NF' OFS=",")
					THREE='1,2,'
					FOUR=$(echo ${THREE}${TWO})	
					#echo $FOUR
					echo "cut -f"$FOUR" $OUTDIR/inds/${file}.mpileup > $OUTDIR/inds/${file}_RindSTATs.rin" >> $OUTDIR/inds/rinlist.commands
			done

			cat $OUTDIR/inds/rinlist.commands | parallel -j ${THREADZ} "{}"


		fi
		
#		rm $OUTDIR/inds/rinlist*  &> /dev/null
		echo "ALERT: R normalization input created at $(date) "


#############################get SNP (e.g. # individuals scored) and individual (e.g. mean coverage) stats for each population
	#REQUIRES R: Get stats on individual files; these can be used to filter positions
		echo "ALERT: Rscript called to calculate population stats at $(date) "
		for file in $(cat $OUTDIR/pops/${OUTPOP}_files_for_pops.txt) ; do
				if [[ -f $OUTDIR/inds/${file}_snp_stats.txt ]] ; then
					echo "ALERT:${file}_snp_stats.txt exists, skipping this analysis "
				else				
					rin=$OUTDIR/inds/${file}_RindSTATs.rin
					rout=$OUTDIR/inds/
					Rscript $BASEDIR/rscripts/r_ind_stats.R $rin $rout
				fi
		done
#		rm $OUTDIR/inds/*RindSTATs.rin &> /dev/null
		echo "ALERT: Rscript stats done at $(date) "

#############################create sync files for each pop/class using pop/class mpileups (with indels!), normalization only
	#Create individual sync files
		echo "ALERT: Individual syncs being created at $(date) "

		for file in $(cat $OUTDIR/pops/${OUTPOP}_files_for_pops.txt) ; do
				if [[ -f ${OUTDIR}/inds/${file}.sync ]] ; then
					echo "ALERT: ${file}.sync exists, skipping this file "
				else
					syncin=$OUTDIR/inds/${file}.mpileup
					nice -n 19 java -ea -${KMEM} -Djava.io.tmpdir=$OUTDIR/tmp -jar ${POOL2}/mpileup2sync.jar --min-qual 1 --fastq-type sanger --input ${syncin} --output ${OUTDIR}/inds/${file}.sync --threads $THREADZ 
				fi
		done
		echo "ALERT: mpileup2sync.jar finished individual syncs at $(date) "




#############################create sync file normalized by individual for each pop/class, noting loci with individuals with more than 2 alleles
	#REQUIRES R: pop/class-specific syncs

	TEST1=( `ls $OUTDIR/inds/*_norm.sync | wc -l` )
	#echo $TEST1
	TEST2=( `wc -l $OUTDIR/pops/${OUTPOP}_files_for_pops.txt` )
	#echo $TEST2

	if [[ $TEST1 -eq $TEST2 ]] ; then
		echo "ALERT: individual normalized syncs exist, skipping normalization"
	else				
	
		if [ -z $DIVIDE ] ; then 
			declare -i DIVIDE=100
		fi
	
		echo "ALERT: Rscript called to standardize syncs at $(date) ; DIVIDE threshold is $DIVIDE"

		##loop over populations and keep up to two alleles per individual, then tally each population
		for file in $(cat $OUTDIR/pops/${OUTPOP}_files_for_pops.txt) ; do
	
	##check that sync exists
			if [[ -f $OUTDIR/inds/${file}_norm.sync ]] ; then
				echo "ALERT:${file}_norm.sync exists, skipping normalization"
			else
				NINDIV=( `wc -l $OUTDIR/inds/${file}_ind_stats.txt` )
				NINDIV=$(($NINDIV-1))

				echo "ALERT:${file} has $NINDIV individuals"
		
				##reduce and tally all individuals in a pop at once
				if [ $NINDIV -le $DIVIDE ]; then 

					echo "ALERT:making normalized sync for ${file} using pooled normalization script"

					rin=$OUTDIR/inds/${file}.sync
					rout=$OUTDIR/inds/
					Rscript $BASEDIR/rscripts/r_standardize.R $rin $rout

				##reduce and tally subsets of individuals, then tally population together
				elif [ $NINDIV -gt $DIVIDE ]; then 

					declare -i NCHUNKS=$(echo $NINDIV/$DIVIDE | bc | sed 's/\.[0-9]*//g' )

					echo "ALERT:making normalized sync for ${file} using subsetted normalization script; DIVIDE threshold is $DIVIDE, meaning $NCHUNKS+1 subgroups"
					
					###what if NINDIV is an exact multiple of DIVIDE? Script goes through extra empty loop, R throws error.  How to fix it??
					#Does it matter?? It creates an empty file (or no file), and R finds no columns to write to temp file, so dies, then loop proceeds normally
					for ((j = 0; j <= $NCHUNKS; j++));
					do
				
						C=$j
						COLFIRST=$(($C*$DIVIDE+4))
						C=$(($j+1))
						COLLAST=$(($C*$DIVIDE+3))
					
						cut -f1-3,$COLFIRST-$COLLAST $OUTDIR/inds/${file}.sync  > $OUTDIR/inds/temp.sync
					
						rin=$OUTDIR/inds/temp.sync
						rout=$OUTDIR/inds/
						Rscript $BASEDIR/rscripts/r_standardize_partial.R $rin $rout
						#output should be temp_norm.sync
						mv $OUTDIR/inds/temp_norm.sync $OUTDIR/inds/${file}_temp_norm$j.sync
						rm $OUTDIR/inds/temp.sync
				
					done

					paste <(cut -f1-3 $OUTDIR/inds/${file}.sync) $OUTDIR/inds/${file}_temp_norm*.sync > $OUTDIR/inds/${file}_temp_cat.sync
	#				rm $OUTDIR/inds/${file}_temp_norm*.sync

					rin=$OUTDIR/inds/${file}_temp_cat.sync
					rout=$OUTDIR/inds/
					Rscript $BASEDIR/rscripts/r_standardize_reduce.R $rin $rout
	#				rm $OUTDIR/inds/${file}_temp_cat.sync
					mv $OUTDIR/inds/${file}_temp_cat_norm.sync $OUTDIR/inds/${file}_norm.sync
					echo "ALERT:Normalized sync file written for ${file}"
				
				else
					echo "ERROR: Can't figure out how to divide syncs by individuals"
					exit 1
				fi
			fi
		done
		echo "ALERT: Rscript standardization done at $(date) "
	
	fi




#############################spiffy awk code returns unique multi-allelic sites without first pre-sorting [deprecated and superceded by non-normalized version from r_frequency.R]
	#Combine black-listed polymorphic sites into one list of unique positions
#				cat $OUTDIR/inds/*_poly_sites.txt | awk '{$3=$1+$2}1'| awk '!a[$3]++' | awk '{print $1,$2}' | sort -k1 -k2,2n > $OUTDIR/filters/${OUTPOP}_ind_poly_blacklist.txt
#				if [[ -f $OUTDIR/inds/*_poly_sites.txt ]] ; then
#					mv $OUTDIR/inds/*_poly_sites.txt $OUTDIR/filters/
#				fi

#############################combine pop/class normalized syncs into a single sync
	#Combine standardized mpileup into new sync
	echo "ALERT: Normalized sync file order is: "
	ITER=1	
	for file in $OUTDIR/inds/*_norm.sync ; do
		namez=${file##*/}
		echo ${ITER} = ${namez%*_norm.sync}
		ITER="$(($ITER + 1))"
	done	
	if [[ -f $OUTDIR/${OUTPOP}_norm.sync ]] ; then
		echo "ALERT: Normalized sync file already exists; skipping"
	else
		echo "ALERT: Creating normalized sync file"
		for file in $(cat $OUTDIR/pops/${OUTPOP}_files_for_pops.txt| awk 'NR==1{print $1}') ; do
			paste <(awk '{print $1"\t"$2"\t"$3}' $OUTDIR/inds/${file}.mpileup) $OUTDIR/inds/*norm.sync >  $OUTDIR/${OUTPOP}_norm_temp.sync
		done
	fi




	#Filter by indels and create sync format
	echo "ALERT: Now identifying indel regions and filtering sync format at $(date)"



#############################filter combined normalized sync for indels and nearby sites using gtf file
	#Filter indels from new sync by gtffile 
	if [[ -f ${OUTDIR}/${OUTPOP}_norm.sync ]] ; then
		echo "ALERT: Normalized indel file already exists; skipping"
	else
		if [[ -f ${OUTDIR}/${OUTPOP}.gtffile ]] ; then
			echo "ALERT: GTF file exists"
		else
#############################take variant-only mpileup, use indel window, identify indel-masked regions, output as gtf
			echo "ALERT: GTF file does not exist; creating from mpileup"
			perl ${POOL2}/indel_filtering/identify-indel-regions.pl --input $OUTDIR/${OUTPOP}.mpileup --indel-window $INWIN --output ${OUTDIR}/${OUTPOP}.gtffile --min-count 1
		fi
		
		if [[ ! -f ${OUTDIR}/${OUTPOP}.gtffile ]] ; then
			echo "ALERT: GTF file not found or created. Exiting"
			exit 1
		fi
			echo "ALERT: Creating filtered and normalized sync using GTF"
			perl ${POOL2}/indel_filtering/filter-sync-by-gtf.pl --input  $OUTDIR/${OUTPOP}_norm_temp.sync --gtf ${OUTDIR}/${OUTPOP}.gtffile --output ${OUTDIR}/${OUTPOP}_norm.sync
			declare -i before=$(wc -l < $OUTDIR/${OUTPOP}_norm_temp.sync)
			declare -i after=$(wc -l < ${OUTDIR}/${OUTPOP}_norm.sync)
			declare -i lost="$(($before - $after))"
			lostp=$((100-100*$after/$before))
			echo "ALERT: With an indel window of $INWIN bp you lost $lost SNPs or $lostp %" 
			echo "SNPs remaining: $after"
			rm $OUTDIR/${OUTPOP}_norm_temp.sync
	fi

#############################create a normalized frequency (fz) file from the (variant-only) sync file, note 3+ allele sites and write to files
	if [[ -f ${OUTDIR}/${OUTPOP}_norm.fz ]] ; then
		echo "ALERT: Standardized Frequency file ${OUTPOP}_norm.fz exists, skipping this analysis "
	else				
		echo "ALERT: Rscript called to create Standardized Frequency file"
		rin=${OUTDIR}/${OUTPOP}_norm.sync
		rout=$OUTDIR/
		rout2=$OUTDIR/filters/
		Rscript $BASEDIR/rscripts/r_frequency.R $rin $rout $rout2 $MAF $INDCONT
		if [[ ! -s ${OUTDIR}/${OUTPOP}_norm.fz  ]]  ; then
			echo "ERROR: Normalized frequency table is empty. You likely ran out of memory."
			echo "Consider DIVIDE-ing individuals or filter SNPs using higher filter thresholds."
			exit
		fi
		echo "ALERT: Frequency file ${OUTPOP}_norm.fz and its counterparts created at $(date) "
	
	fi
fi

rm -rf $OUTDIR/tmp

echo "ALERT: PPalign completed at $(date)"
echo "Quality reporting may still be running in the background; if PPalign appears stuck, this is usually why." 
echo "Ctrl-C will abort quality reporting"

if [[ -f $OUTDIR/${OUTPOP}.gtffile.params ]] ; then
	rm $OUTDIR/${OUTPOP}.gtffile.params
fi
if [[ -f $OUTDIR/${OUTPOP}.sync.params ]] ; then
	rm $OUTDIR/${OUTPOP}.sync.params
fi
if [[ -f $OUTDIR/${OUTPOP}_norm.sync.params ]] ; then
	rm $OUTDIR/${OUTPOP}_norm.sync.params
fi

exit