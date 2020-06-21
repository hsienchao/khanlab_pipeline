#!/usr/bin/bash

ref_genome=$1
spike_in_genome=$2

for dn in */hic_results/stats/*;do 
	bn=$(basename "$dn")
	bd=$(echo "$dn" | cut -d "/" -f1)
	genome=$ref_genome
	if [[ $bd == "HiCproAQuAOUTPUT" ]];then
		genome=$spike_in_genome
	fi
	printf "\t$bn.$genome" >> header.txt
	if [ ! -f mergeStats.txt ];then 		
		cat $bd/hic_results/stats/$bn/${bn}_allValidPairs.mergestat > mergeStats.txt
	else		
		paste mergeStats.txt <(cut -f2 $bd/hic_results/stats/$bn/${bn}_allValidPairs.mergestat) > tmp.txt
		mv tmp.txt mergeStats.txt
	fi	
done
echo "" >> header.txt
cat header.txt mergeStats.txt > tmp.txt
rm header.txt
mv tmp.txt mergeStats.txt