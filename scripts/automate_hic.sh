#!/usr/bin/bash

sample_id=$1

d=`dirname $0`
cd $d/..
pipeline_home=`pwd`

master_file=/data/khanlab/projects/HiC/manage_samples/HiC_sample_sheet.xlsx
data_home=/data/khanlab/projects/DATA
processed_data_home=/data/khanlab/projects/HiC/processed_DATA

if [ ! -z $sample_id ];then
	if [ $sample_id == "all" ];then
		samples=`python scripts/getHiCSampleList.py $master_file`
	else
		samples=( $sample_id )
	fi
else
	echo "usage: $0 [sample_id|all]"
	echo
	echo "Example 1: process specific sample"
	echo "   $0 RH4_Ent6_H3K27ac_HiChIP_HH3JVBGX7"
	echo
	echo "Example 2: process all samples"
	echo "   $0 all"
	echo
	exit
fi

module load python/3.6

for sample_id in $samples;do
	if [ -s $data_home/Sample_$sample_id/Sample_${sample_id}_R1.fastq.gz ] && [ -s $data_home/Sample_$sample_id/Sample_${sample_id}_R2.fastq.gz ];then
		if [ ! -f $pipeline_home/$sample_id/successful.txt ];then
			python $pipeline_home/scripts/sampleToYaml.py $master_file $sample_id > $processed_data_home/sample_sheets/$sample_id.hic.yaml
			if [ $? -eq 0 ]
			then
				echo "launching HiC pipeline: $sample_id"
				perl $pipeline_home/launch -t hic -w $processed_data_home/$sample_id -s $processed_data_home/sample_sheets/$sample_id.hic.yaml
			else
				rm $processed_data_home/sample_sheets/$sample_id.hic.yaml
				echo "Sample sheet generation failed: $sample_id" |mutt -s 'Khanlab HiC Pipeline Status' `whoami`@mail.nih.gov
			fi
		fi
	fi
done
