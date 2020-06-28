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
	printf "\t$bn.$genome" >> $sample_id/header.txt
	if [ ! -f $sample_id/mergeStats.txt ];then
		cat $sample_id/$bd/hic_results/stats/$bn/${bn}_allValidPairs.mergestat > $sample_id/mergeStats.txt
	else		
		paste $sample_id/mergeStats.txt <(cut -f2 $sample_id/$bd/hic_results/stats/$bn/${bn}_allValidPairs.mergestat) > $sample_id/tmp.txt
		mv $sample_id/tmp.txt $sample_id/mergeStats.txt
	fi	
done
echo "" >> $sample_id/header.txt
cat $sample_id/header.txt $sample_id/mergeStats.txt > $sample_id/tmp.txt
rm $sample_id/header.txt
mv $sample_id/tmp.txt $sample_id/mergeStats.txt