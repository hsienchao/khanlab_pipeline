#!/usr/bin/python
import re
import requests
import os,subprocess
import pandas as pd
import yaml
import sys
import argparse

def main(args):

    #User input:
    #library_id and run_id
    
    library_id = args.lib.strip()
    run_id = args.run.strip()
    sample_id = library_id + '_' + run_id 
    out_format = args.out_format.strip()
    out_file = "test.out"
    if args.out is not None:
        out_file = args.out.strip()
    masterfile_dir="/data/Clinomics/MasterFiles"
    master_files = ["Sequencing_Tracking_Master_db.txt","ClinOmics_Sequencing_Master_File_db.txt","SequencingMasterFile_OutsidePatients_db.txt"]
    hic_master_file=args.hic.strip()
    chipseq_master_file=args.chip.strip()
    #sample information from Young's master file
    #sample_file = "Sample_" + sample_id
    #master_sample = {"SampleFiles":sample_file}
    #sample = {"SampleFiles":sample_file}
    default_genome="hg19"
    xeno_genome="mm10"
    #columns = ["Type of sequencing","Matched normal","Matched RNA-seq"]
    meta_samples=[]
    lib_col_name = 'Amplified_Sample_Library_Name'
    id_col_name = 'Sample_ID'
    for master_file in master_files:
        file = masterfile_dir + "/" + master_file
        master_df = pd.read_csv(file, sep='\t', encoding = "ISO-8859-1")
        master_df[id_col_name] = master_df['Library ID'] + '_' + master_df['FCID']
        if run_id != '':
            master_df = master_df.loc[master_df[id_col_name] == sample_id]
        else:
            master_df = master_df.loc[master_df['Library ID'] == library_id]
        master_df = master_df.reset_index(drop=True)
        if master_df[id_col_name].count() > 0:
            meta_samples = master_df
            #for column in master_df:
            #    master_sample[column] = str(master_df[column][0])
            break
    if id_col_name not in meta_samples:
        print(library_id + " not found in Khanlab master files\n")
        sys.exit(1)
    lib_type = meta_samples["Type of sequencing"][0]
    samples = {}
    #HiC sample
    if "H-il" in lib_type:
        df = pd.read_excel(hic_master_file)
        df = df.loc[df[lib_col_name] == library_id]
        df = df.reset_index(drop=True)
        if df[lib_col_name].count() > 0:
            for i in range(meta_samples[id_col_name].count()):
                sample = {}
                sample_id = meta_samples[id_col_name][i]
                sample_file = "Sample_" + sample_id
                sample[id_col_name] = sample_id
                sample["SampleFiles"] = sample_file
                for column in df:
                    sample[column] = str(df[column][0])
                samples[sample_id] = sample
            print("hic")
            print(sample["Genome"])
        #else:
        #    sys.stderr.write(sample_id + " not found in Khanlab master files\n")
    #ChIPseq sample
    if "C-il" in lib_type and len(samples) == 0:
        df = pd.read_excel(chipseq_master_file)
        df = df.loc[df[lib_col_name] == library_id]
        df = df.reset_index(drop=True)
        if df[lib_col_name].count() > 0:
            fo = open(out_file,"w")
            for i in range(meta_samples[id_col_name].count()):
                sample = {}
                sample_id = meta_samples[id_col_name][i]
                sample_file = "Sample_" + sample_id
                sample[id_col_name] = sample_id
                sample["SampleFiles"] = sample_file
                for column in df:
                    str_data = str(df[column][0])
                    if column == "PairedRNA_SAMPLE_ID":
                        str_data = str_data.replace("Sample_","")
                    if column == "PairedInput":
                        str_data = str_data.replace("Sample_","")
                    sample[column] = str_data
                if sample["SpikeIn"] == "yes" and not "SpikeInGenome" in sample:
                    sample["SpikeInGenome"] = "dm6"
                sample["SampleFiles"] = sample_file
                samples[sample_id] = sample
            fo.write(yaml.dump({"samples":{sample_id : sample}}))
            fo.close()
            print("chipseq")
            print(sample["Genome"])
        else:
            sys.stderr.write(sample_id + " not found in HiC/Chipseq master files\n")
    if "T-il" in lib_type and len(samples) == 0:
        for i in range(meta_samples[id_col_name].count()):
            sample = {}
            sample_id = meta_samples[id_col_name][i]
            sample_file = "Sample_" + sample_id
            sample[id_col_name] = sample_id
            sample["SampleFiles"] = sample_file
            if "SampleRef" in meta_samples:
                sample["Genome"] = meta_samples["SampleRef"][i]
            else:
                sample["Genome"] = default_genome
            if meta_samples["Type"][i].find("xeno") >=0:
                sample["Xenograft"] = "yes"
                sample["XenograftGenome"] = xeno_genome
            sample["SampleCaptures"] = meta_samples["Enrichment step"][i]            
            for col in meta_samples:
                sample[col] = meta_samples[col][i]
            samples[sample_id] = sample
        print("rnaseq")
        print(sample["Genome"])
    if len(samples) > 0:
        fo = open(out_file,"w")
        if out_format == "yaml":
            fo.write(yaml.dump({"samples":samples}))
        else:
            cols = list(samples[list(samples.keys())[0]].keys())
            fo.write('\t'.join(cols)+'\n')
            for sample_id in samples:
                d = []
                for col in cols:
                    d.append(str(samples[sample_id][col]))
                fo.write('\t'.join(d)+'\n')
        fo.close()
parser = argparse.ArgumentParser(description='Generate sample sheet.')
parser.add_argument("--lib", "-l", metavar="LIBRARY_ID", required=True, help="Library ID")
parser.add_argument("--run", "-r", metavar="RUN_ID", help="Run ID (FCID or SRR ID)", default="")
parser.add_argument("--out_format", "-f", metavar="OUTPUT", help="Output format: table or yaml", default="table")
parser.add_argument("--out", "-o", metavar="OUTPUT", help="Output file")
parser.add_argument("--chip", "-c", metavar="ChIPseq MASTER", help="Chipseq Master file", default="/data/khanlab/projects/ChIP_seq/manage_samples/ChIP_seq_samples.xlsx")
parser.add_argument("--hic", "-i", metavar="ChIPseq MASTER", help="HiC Master file", default="/data/khanlab/projects/HiC/manage_samples/HiC_sample_sheet.xlsx")
args = parser.parse_args()

main(args)