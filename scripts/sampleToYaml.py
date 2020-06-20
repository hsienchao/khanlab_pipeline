import re
import requests
import os,subprocess
import pandas as pd
import yaml
import sys

def main():

    #User input:
    #1. input Master File
    #2. library
    #3. output file
    
    master_file=sys.argv[1]
    library=sys.argv[2]
    #outFile=sys.argv[3]
    #master_file="/data/khanlab/projects/HiC/manage_samples/HiC_sample_sheet.xlsx"
    #library="RH4_D6_H3K27ac_HiChIP_HKJ22BGX7"
    df = pd.read_excel(master_file)
    df['Sample_ID'] = df['Amplified_Sample_Library_Name'] + '_' + df['FlowCell_GSE']
    df = df.loc[df['Sample_ID'] == library]
    df = df.reset_index(drop=True)    
    if df['Sample_ID'].count() == 0:
        sys.exit("Sample_ID:" + library + " not found")
    if df['Sample_ID'].count() > 1:
        sys.exit("Multiple sample:" + library + " found")
    sample = {}
    for column in df:
        sample[column] = df[column][0]
    print(yaml.dump({"samples":{sample['Sample_ID'] : sample}}))


main()