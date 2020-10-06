#!/usr/bin/python
import re
import requests
import os,subprocess
import pandas as pd
import yaml
import sys
import argparse

def main(args):

    library_id = args.lib.strip()
    masterfile_dir="/data/Clinomics/MasterFiles"
    master_files = ["Sequencing_Tracking_Master_db.txt","ClinOmics_Sequencing_Master_File_db.txt","SequencingMasterFile_OutsidePatients_db.txt"]
    lib_col_name = 'Library ID'
    for master_file in master_files:
        file = masterfile_dir + "/" + master_file
        master_df = pd.read_csv(file, sep='\t', encoding = "ISO-8859-1")
        master_df = master_df.loc[master_df[lib_col_name] == library_id]
        master_df = master_df.reset_index(drop=True)
        if master_df[lib_col_name].count() > 0:
            for i in range(master_df[lib_col_name].count()):
                print(master_df['FCID'][i])
                
parser = argparse.ArgumentParser(description='Get Run ID')
parser.add_argument("--lib", "-l", metavar="LIBRARY_ID", required=True, help="Library ID")
args = parser.parse_args()

main(args)