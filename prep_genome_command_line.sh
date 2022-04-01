#!/bin/bash

#prepares reference genome in FASTA format for BWA alignment

if [ -z "$2" ]
then
echo "Not enough arguments."
echo "Correct usage: SCRIPT.sh fasta_name picard_jar_path"
echo "Example picard path: ~/miniconda3/envs/poolparty_env/share/picard-2.20.2-0/picard.jar"
exit 1
#else
fi

GENOME=$1
PICARD=$2	

#!/bin/bash
# 1) Index fastA genome file with bwa
	
bwa index -a bwtsw ${GENOME}

# 2) Create fastA file index as well

samtools faidx ${GENOME}

# 3) Create sequence dictionary with picardtools (Java)

#java -Xmx2g -jar ~/miniconda3/envs/poolparty_env/share/picard-2.20.2-0/picard.jar CreateSequenceDictionary REFERENCE= ${GENOME} OUTPUT=${GENOME}.dict
java -Xmx2g -jar ${PICARD} CreateSequenceDictionary REFERENCE= ${GENOME} OUTPUT=${GENOME}.dict

