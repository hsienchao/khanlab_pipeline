#!/usr/bin/bash

perl launch -t hic -w /data/khanlab/projects/HiC/pipeline_test/SRR400264 -s /data/khanlab/projects/HiC/pipeline_test/sample_sheets/SRR400264.yaml -d /data/khanlab/projects/HiC/pipeline_test/test_data
perl launch -t hic -w /data/khanlab/projects/HiC/pipeline_test/SRR400264_00 -s /data/khanlab/projects/HiC/pipeline_test/sample_sheets/SRR400264_00.yaml -d /data/khanlab/projects/HiC/pipeline_test/test_data
perl launch -t hic -w /data/khanlab/projects/HiC/pipeline_test/SRR400264_01 -s /data/khanlab/projects/HiC/pipeline_test/sample_sheets/SRR400264_01.yaml -d /data/khanlab/projects/HiC/pipeline_test/test_data
