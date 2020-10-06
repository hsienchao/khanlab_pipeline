#merge CTR_rep6_H3K27ac_HiChIP_HVHT2BGXB and CTR_rep7_H3K27ac_HiChIP_HVHT2BGXB
merged_sample=CTR_rep6_7_H3K27ac_HiChIP_HVHT2BGXB
hic_home=/data/khanlab/projects/pipeline_dev_hc/processed_DATA/hic
config_file=/data/khanlab/projects/pipeline_dev_hc/config/hicpro/hg19_dpnii.txt
dest_dir=${hic_home}/${merged_sample}/HiCproOUTPUT.hg19/hic_results/data/${merged_sample}
source1_dir=${hic_home}/CTR_rep6_H3K27ac_HiChIP_HVHT2BGXB/HiCproOUTPUT.hg19/hic_results/data/CTR_rep6_H3K27ac_HiChIP_HVHT2BGXB
source2_dir=${hic_home}/CTR_rep7_H3K27ac_HiChIP_HVHT2BGXB/HiCproOUTPUT.hg19/hic_results/data/CTR_rep7_H3K27ac_HiChIP_HVHT2BGXB
mkdir $dest_dir
cp ${source1_dir}/* ${dest_dir}/
cp ${source2_dir}/* ${dest_dir}/
sinteractive --mem=16G --time=10:00:00
module load hicpro
HiC-Pro -i ${hic_home}/${merged_sample}/HiCproOUTPUT.hg19/hic_results/data -o ${hic_home}/${merged_sample}/HiCproOUTPUT.hg19 -c ${config_file} -s quality_checks -s merge_persample -s build_contact_maps -s ice_norm
hicpro2juicebox.sh -i ${dest_dir}/${merged_sample}.allValidPairs -g hg19 -j /usr/local/apps/juicer/juicer-1.5.6/scripts/juicer_tools.jar -o ${merged_sample}

