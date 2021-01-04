try:
    #prepare targets    
    HICPRO_SCRIPTS = []    
    MERGED_PAIRS = []    
    MERGE_STATS_SUMMARY = []
    JUICEBOX_OUT = []
    VIEWPOINT_OUT = []
    FASTQCS = []
    FITHICHIP_OUT = []
    PEAKACHU_OUT = []
    PAIRS = {}
    MERGE_STATS = {}
    FASTQS = {}
    

    for sample_id, sample in samples.items():
        sample["Genome"] = config["genome"]
        FASTQS[sample_id] = []
        MERGE_STATS[sample_id] = []
        PAIRS[sample_id] = []
        
        SAMPLES.append(sample_id)
        #add FASTQ targets
        if not "SampleFiles" in sample:
            raise Exception('SampleFiles not found in sample sheet')
        sample_file = sample["SampleFiles"]
        if sample["Trim"] != ".":
            FASTQS[sample_id].append(sample_id + "/DATA/" + sample_id + '/' + sample_id + ".trimmed" + suffix_R1)
            FASTQS[sample_id].append(sample_id + "/DATA/" + sample_id + '/' + sample_id + ".trimmed" + suffix_R2)
        else:
            FASTQS[sample_id].append(sample_id + "/DATA/" + sample_id + '/' + sample_id + suffix_R1)
            FASTQS[sample_id].append(sample_id + "/DATA/" + sample_id + '/' + sample_id + suffix_R2)
        if not os.path.exists(data_dir + "/" + sample_file + "/" + sample_file + suffix_R1):
            raise Exception(data_dir + "/" + sample_file + "/" + sample_file + suffix_R1 + ' not found')
        if not os.path.exists(data_dir + "/" + sample_file + "/" + sample_file + suffix_R2):
            raise Exception(data_dir + "/" + sample_file + "/" + sample_file + suffix_R2 + ' not found')
        FASTQCS.append(sample_id + "/qc/" + sample_id + '_R1_fastqc.html')
        FASTQCS.append(sample_id + "/qc/" + sample_id + '_R2_fastqc.html')
        #add HiC-Pro targets
        if not "Genome" in sample:
            raise Exception('Genome not found in sample sheet')
        ref_genome = sample["Genome"]
        if "SpikeIn" in sample and (sample["SpikeIn"] == "yes" or sample["SpikeIn"] == True):
            print("We have SpikeIn\n")
            #only add script targets once
            if len(HICPRO_SCRIPTS) == 0:
                if not "SpikeInGenome" in sample:
                    raise Exception('SpikeInGenome not found in sample sheet')
                spike_in_genome = sample["SpikeInGenome"]
                HICPRO_SCRIPTS.append(sample_id + "/HiCproOUTPUT." + spike_in_genome + "/HiCPro_step1_" + spike_in_genome + "_HiCpro.sh")
                HICPRO_SCRIPTS.append(sample_id + "/HiCproOUTPUT." + spike_in_genome + "/HiCPro_step2_" + spike_in_genome + "_HiCpro.sh")
            PAIRS[sample_id].append(sample_id + "//HiCproOUTPUT." + spike_in_genome + "/hic_results/data/" + sample_id + "/" + sample_id + "_" + config[spike_in_genome]["bowtie2_index"] + ".bwt2pairs.validPairs")
            MERGED_PAIRS.append(sample_id + "/HiCproOUTPUT." + spike_in_genome + "/hic_results/data/" + sample_id + "/" + sample_id + ".allValidPairs")
            MERGE_STATS[sample_id].append(sample_id + "/HiCproOUTPUT." + spike_in_genome + "/hic_results/stats/" + sample_id + "/" + sample_id + "_allValidPairs.mergestat")
        else:
            print("We have no SpikeIn\n")
        HICPRO_SCRIPTS.append(sample_id + "/HiCproOUTPUT." + ref_genome + "/HiCPro_step1_" + ref_genome + "_HiCpro.sh")
        HICPRO_SCRIPTS.append(sample_id + "/HiCproOUTPUT." + ref_genome + "/HiCPro_step2_" + ref_genome + "_HiCpro.sh")
        PAIRS[sample_id].append(sample_id + "/HiCproOUTPUT." + ref_genome + "/hic_results/data/" + sample_id + "/" + sample_id + "_" + config[ref_genome]["bowtie2_index"] + ".bwt2pairs.validPairs")
        MERGED_PAIRS.append(sample_id + "/HiCproOUTPUT." + ref_genome + "/hic_results/data/" + sample_id + "/" + sample_id + ".allValidPairs")    
        MERGE_STATS[sample_id].append(sample_id + "/HiCproOUTPUT." + ref_genome + "/hic_results/stats/" + sample_id + "/" + sample_id + "_allValidPairs.mergestat")    
        JUICEBOX_OUT.append(sample_id + "/" + sample_id + ".allValidPairs.hic")
        FITHICHIP_OUT.append(sample_id + "/FitHiChIP_Out/peaks/MACS2_input.bed")
        FITHICHIP_OUT.append(sample_id + "/FitHiChIP_Out/Summary_results_FitHiChIP.html")
        PEAKACHU_OUT.append(sample_id + "/peakachu_Out/interactions.bed")
        MERGE_STATS_SUMMARY.append(sample_id + "/mergeStats.txt")
        #VIEWPOINT_OUT.append(sample_id + "/HiCproOUTPUT." + ref_genome + "/hic_results/data/" + sample_id + "/" + sample_id + ".bedgraph")
except Exception as err:
    exc_type, exc_value, exc_traceback = sys.exc_info()
    output = io.StringIO()
    traceback.print_exception(exc_type, exc_value, exc_traceback, file=output)
    contents = output.getvalue()
    output.close()
    print(contents)    
    shell("echo 'HiC pipeline has exception: reason " + contents + ". Working Dir:  {work_dir}' |mutt -e 'my_hdr From:chouh@nih.gov' -s 'Khanlab HiC Pipeline Status' `whoami`@mail.nih.gov {emails} ")
    sys.exit()    

TARGETS = MERGED_PAIRS+JUICEBOX_OUT+MERGE_STATS_SUMMARY+FASTQCS+FITHICHIP_OUT+PEAKACHU_OUT

localrules: HiCPro, prepareFASTQ, HiC_pipeline, mergeStats
#print(HICPRO_SCRIPTS,MERGED_PAIRS,JUICEBOX_OUT,MERGE_STATS_SUMMARY)

include: "hicpro.smk"

rule peakachu:
    input:
            "{sample}/{sample}.allValidPairs.hic"
    output: 
            "{sample}/peakachu_Out/interactions.bed"
    version:
            config["version_common"]["sigularity"]
    params:
            work_dir = config["work_dir"],
            pipeline_home = config["pipeline_home"],
            conda_path = config["peakachu"]["conda_path"],
            peakachu_env = config["peakachu"]["env"],
            peakachu_script = config["peakachu"]["script"],
            prob_cutoff = config["peakachu"]["prob_cutoff"],
            model = config["peakachu"]["model"],
            rulename = "peakachu",
            log_dir = lambda wildcards: wildcards.sample + '/log',
            genome = lambda wildcards: samples[wildcards.sample]["Genome"],
            batch    = config["cluster_common"]["medium"],
    benchmark:
            "{sample}/benchmark/peakachu.benchmark.txt"
    shell:
            """            
            source {params.conda_path}/etc/profile.d/conda.sh
            conda activate {params.peakachu_env}
            rm -f {output}
            python {params.peakachu_script} score_genome --balance -p {input} --output {wildcards.sample}/peakachu_Out -m {params.pipeline_home}/{params.model}
            if [ "$(ls -A {wildcards.sample}/peakachu_Out)" ];then
                for i in {wildcards.sample}/peakachu_Out/*; do python {params.peakachu_script} pool -i ${{i}} -t {params.prob_cutoff} >> {output}; done
            else
                # no significant interactions
                touch {output}
            fi
            """
            
rule FitHiChIP:
    input:
            "{sample}/FitHiChIP_Out/peaks/MACS2_input.bed"
    output: 
            "{sample}/FitHiChIP_Out/Summary_results_FitHiChIP.html"
    version:
            config["version_common"]["sigularity"]
    params:
            work_dir = config["work_dir"],
            fithichip_instance = config["app_home"] + '/' + config["fithichip_instance"],
            rulename = "FitHiChIP",
            log_dir = lambda wildcards: wildcards.sample + '/log',
            genome = lambda wildcards: samples[wildcards.sample]["Genome"],
            chr_file= lambda wildcards: config[samples[wildcards.sample]["Genome"]]["chr_size"],
            batch    = config["cluster_common"]["medium"]
    benchmark:
            "{sample}/benchmark/FitHiChIP.benchmark.txt"
    shell:
            """
            module load singularity/{version}
            echo "ValidPairs=/hicpro_data/hic_results/data/{wildcards.sample}/{wildcards.sample}.allValidPairs" > {wildcards.sample}/FitHiChIP_Out/fithichip.conf
            echo "PeakFile=/FitHiChIP_Out/peaks/MACS2_input.bed" >> {wildcards.sample}/FitHiChIP_Out/fithichip.conf
            echo "ChrSizeFile={params.chr_file}" >> {wildcards.sample}/FitHiChIP_Out/fithichip.conf
            echo "OutDir=/FitHiChIP_Out" >> {wildcards.sample}/FitHiChIP_Out/fithichip.conf
            echo "HiCProBasedir=/HiC-Pro-2.11.1/" >> {wildcards.sample}/FitHiChIP_Out/fithichip.conf
            singularity exec -B /data/khanlab/projects/HiC/reference_files/,{params.work_dir}/{wildcards.sample}/FitHiChIP_Out:/FitHiChIP_Out,{params.work_dir}/{wildcards.sample}/HiCproOUTPUT.{params.genome}:/hicpro_data {params.fithichip_instance} bash /FitHiChIP/FitHiChIP_HiCPro.sh -C /FitHiChIP_Out/fithichip.conf
            """
rule FitHiChIPCallPeak:
    input:            
            lambda wildcards: wildcards.sample + "/HiCproOUTPUT." + samples[wildcards.sample]["Genome"] + "/hic_results/data/" + wildcards.sample + "/" + wildcards.sample + ".allValidPairs"
    output: 
            "{sample}/FitHiChIP_Out/peaks/MACS2_input.bed"
    version:
            config["version_common"]["sigularity"]
    params:
            work_dir = config["work_dir"],
            pipeline_home = config["pipeline_home"],
            fithichip_instance = config["app_home"] + '/' + config["fithichip_instance"],
            rulename = "FitHiChIPCallPeak",
            log_dir = lambda wildcards: wildcards.sample + '/log',
            genome = lambda wildcards: samples[wildcards.sample]["Genome"],
            batch    = config["cluster_common"]["medium"],
            suffix_R1 = suffix_R1
    benchmark:
            "{sample}/benchmark/FitHiChIPCallPeak.benchmark.txt"
    shell:
            """            
            module load singularity/{version}
            mkdir -p {wildcards.sample}/FitHiChIP_Out
            mkdir -p {wildcards.sample}/FitHiChIP_Out/peaks
            read_len=`bash {params.pipeline_home}/scripts/getReadLength.sh {wildcards.sample}/DATA/{wildcards.sample}/{wildcards.sample}{params.suffix_R1}`
            echo "${{read_len}}"
            singularity exec -B {params.work_dir}/{wildcards.sample}/FitHiChIP_Out:/FitHiChIP_Out,{params.work_dir}/{wildcards.sample}/HiCproOUTPUT.{params.genome}:/hicpro_data {params.fithichip_instance} bash /FitHiChIP/Imp_Scripts/PeakInferHiChIP.sh -H /hicpro_data -D /FitHiChIP_Out/peaks -L ${{read_len}}
            """
rule cutadapt:
    input:
            lambda wildcards: wildcards.sample + "/DATA/" + wildcards.sample + "/" + wildcards.sample + "_" + wildcards.suffix
    output:
            "{sample}/DATA/{sample}/{sample}.trimmed_{suffix}"
    params:
            adapter_f = lambda wildcards: samples[wildcards.sample]["Trim"],
            adapter_r = lambda wildcards: reverse_complement(samples[wildcards.sample]["Trim"]),
            work_dir = config["work_dir"],
            rulename = "cutadapt",
            log_dir = lambda wildcards: wildcards.sample + '/log',
            batch    = config["cluster_common"]["medium"]
    benchmark:
            "{sample}/benchmark/cutadapt.{suffix}.benchmark.txt"
    version:
            config["version"]["cutadapt"]
    shell:
            """
            module load cutadapt/{version}
            cutadapt -a {params.adapter_f} -a {params.adapter_r} -o {output} -j ${{THREADS}} {input}
            """

rule fastqc:
    input:
            lambda wildcards: wildcards.sample + "/DATA/" + wildcards.sample + "/" + wildcards.sample + "_" + wildcards.suffix + ".fastq.gz"
    output:
            "{sample}/qc/{sample}_{suffix}_fastqc.html"
    params:
            work_dir = config["work_dir"],
            rulename = "fastqc",
            log_dir = lambda wildcards: wildcards.sample + '/log',
            batch    = config["cluster_common"]["medium"]
    benchmark:
            "{sample}/benchmark/fastqc.{suffix}.benchmark.txt"
    version:
            config["version_common"]["fastqc"]
    shell:
            """
            module load fastqc/{version}
            mkdir -p {wildcards.sample}/qc
            fastqc -t ${{THREADS}} -o {wildcards.sample}/qc {input}
            """
            
rule prepareFASTQ:
    input: 
            lambda wildcards: data_dir + "/" + samples[wildcards.sample]["SampleFiles"] + "/" + samples[wildcards.sample]["SampleFiles"] + "_" + wildcards.suffix
    output: 
            "{sample}/DATA/{sample}/{sample}_{suffix}"
    shell:
            """
            mkdir -p {wildcards.sample}
            mkdir -p {wildcards.sample}/log
            mkdir -p {wildcards.sample}/DATA            
            mkdir -p {wildcards.sample}/DATA/{wildcards.sample}
            ln -s {input} {output}
            """