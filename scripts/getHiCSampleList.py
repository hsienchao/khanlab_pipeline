import re
import requests
import os,subprocess
import pandas as pd
import yaml
import sys

def main():

    #User input:
    #1. input Master File
    
    master_file=sys.argv[1]
    #outFile=sys.argv[3]
    #master_file="/data/khanlab/projects/HiC/manage_samples/HiC_sample_sheet.xlsx"
    #library="RH4_D6_H3K27ac_HiChIP_HKJ22BGX7"
    df = pd.read_excel(master_file)
    df['Sample_ID'] = df['Amplified_Sample_Library_Name'] + '_' + df['FlowCell_GSE']
    df = df.loc[df['FlowCell_GSE'] != 'TBD'].loc[df['FlowCell_GSE'] != '.'].loc[df['FlowCell_GSE'] != '']
    df = df.reset_index(drop=True)
    for sid in df['Sample_ID']:
        print(sid)    

main()