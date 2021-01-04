#!/usr/bin/bash

sample_id=$1
ref_genome=$2
spike_in_genome=$3


for dn in $1/*/hic_results/stats/*;do 
	bn=$(basename "$dn")
	bd=$(echo "$dn" | cut -d "/" -f2)
	genome=$ref_genome
	if [[ $bd == "HiCproOUTPUT.${spike_in_genome}" ]];then
		genome=$spike_in_genome
	fi
	awk -F"\t" 'BEGIN {OFS = FS} {print;
		if ($1=="valid_interaction") all=$2; if ($1=="valid_interaction_rmdup") print "percent_valid_interaction_rmdup",$2/all*100 "%" 
		if ($1=="trans_interaction") trans=$2; if ($1=="cis_interaction") print "percent_cis_interaction",$2/(trans+$2)*100 "%" 
		if ($1=="cis_shortRange") short=$2; if ($1=="cis_longRange") print "percent_longRange_interaction",$2/(short+$2)*100 "%"}' $sample_id/$bd/hic_results/stats/$bn/${bn}_allValidPairs.mergestat > $sample_id/$bd/hic_results/stats/$bn/${bn}_allValidPairs.mergestat_with_percent
	printf "\t$bn.$genome" >> $sample_id/header.txt
	if [ ! -f $sample_id/mergeStats.txt ];then
		cat $sample_id/$bd/hic_results/stats/$bn/${bn}_allValidPairs.mergestat_with_percent > $sample_id/mergeStats.txt
	else		
		paste $sample_id/mergeStats.txt <(cut -f2 $sample_id/$bd/hic_results/stats/$bn/${bn}_allValidPairs.mergestat_with_percent) > $sample_id/tmp.txt
		mv $sample_id/tmp.txt $sample_id/mergeStats.txt
	fi	
done
echo "" >> $sample_id/header.txt
cat $sample_id/header.txt $sample_id/mergeStats.txt > $sample_id/tmp.txt
rm $sample_id/header.txt
mv $sample_id/tmp.txt $sample_id/mergeStats.txt