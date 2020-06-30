#!/usr/bin/bash

sample_id=$1

d=`dirname $0`
cd $d/..
pipeline_home=`pwd`
sample_sheet_home=/data/khanlab/projects/pipeline_production/sample_sheets
processed_data_home=/data/khanlab/projects/pipeline_production/processed_DATA

if [ ! -z $sample_id ];then
	samples=( $sample_id )
else
	echo "usage: $0 [sample_id]"
	echo
	echo "Example 1: process specific sample"
	echo "   $0 RH4_Ent6_H3K27ac_HiChIP_HH3JVBGX7"
	echo
	exit
fi

module load python/3.6

type=`python $pipeline_home/scripts/sampleToYaml.py -s $sample_id -o ${sample_id}.yaml`
yaml_file=$sample_sheet_home/$type/${sample_id}.yaml
if [[ "$type" == "hic" || "$type" == "chipseq" || "$type" == "rnaseq" ]];then
	mv ${sample_id}.yaml $yaml_file
	perl $pipeline_home/launch -t $type -w $processed_data_home/$type -s $yaml_file
else
	echo "Error: Sample sheet generation failed: $sample_id. Reason: $type" |mutt -s 'Khanlab Pipeline Status' `whoami`@mail.nih.gov
fi

