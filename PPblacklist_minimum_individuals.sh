#!/bin/bash
#v04.01.2022

#collect snps_stats files, awk to find sites below MIN, create temporary list, then sort for non-redundant list

if [ -z "$1" ]
then
echo "Not enough arguments."
echo "Correct usage: SCRIPT.sh MIN"
echo "Output will be called min_MIN_ind_blacklist.txt."
exit 1
#else
fi

min=$1

FILES=( `ls *snp_stats* `)
LENb=( `ls *snp_stats* | wc -l `)

echo "There appear to be "$LENb" snp_stats files to query. Listing sites below "$min" individuals"

LENb=$(($LENb - 1))

rm temp.txt &> /dev/null
touch temp.txt

for ((j = 0; j <= $LENb; j++));
do
		echo "Now querying file "${FILES[j]}
		awk -v pat="$min" '$5<pat' "${FILES[j]}" | awk '{print $1,$2}' >> temp.txt
done

new_out=$(paste -d"_" <(paste -d"_" <(echo "min") <(echo "$min")) <(echo "ind_blacklist.txt"))

echo "Now sorting flagged sites to remove duplicates"
cat temp.txt | sort | uniq > $new_out
rm temp.txt &> /dev/null
