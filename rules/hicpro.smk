rule makeViewPoint:
    input:  
            "{sample}/HiCproOUTPUT.{genome}/hic_results/data/{sample}/{sample}.allValidPairs"
    output: 
            "{sample}/HiCproOUTPUT.{genome}/hic_results/data/{sample}/{sample}.bedgraph"
    params:
            fragment_bed=lambda wildcards: config[wildcards.genome]["fragment_bed"],
            capture_bed=lambda wildcards: config[wildcards.genome]["capture_bed"],
            rulename = "makeViewPoint",
            log_dir = lambda wildcards: wildcards.sample + '/log',
            batch    = config["cluster_common"]["medium"]
    shell:
            """
            module load hicpro/{version}
            make_viewpoints.py -i {input} -f {params.fragment_bed} -t {params.capture_bed} -e 1000 -v -o {output}
            """

rule HiCpro2Juicebox:
    input: 
            lambda wildcards: wildcards.sample + "/HiCproOUTPUT." + samples[wildcards.sample]["Genome"] + "/hic_results/data/" + wildcards.sample + "/" + wildcards.sample + ".allValidPairs"
    output: 
            "{sample}/{sample}.allValidPairs.hic"
    version:
            config["version"]["hicpro"]
    params:
            juicer_jar=config["juicer_jar"],
            juicer_genome=lambda wildcards: config[samples[wildcards.sample]["Genome"]]["juicer_genome"],
            rulename = "HiCpro2Juicebox",
            log_dir = lambda wildcards: wildcards.sample + '/log',
            batch    = config["cluster_common"]["medium"]
    benchmark:
            "{sample}/benchmark/HiCpro2Juicebox.benchmark.txt"
    shell:
            """            
            module load hicpro/{version}
            hicpro2juicebox.sh -i {input} -g {params.juicer_genome} -j {params.juicer_jar} -o {wildcards.sample}
            """

rule mergeStats:
    input:
            lambda wildcards: MERGE_STATS[wildcards.sample]
    output:
            "{sample}/mergeStats.txt"
    params:
            pipeline_home = config["pipeline_home"],
            ref_genome = lambda wildcards: samples[wildcards.sample]["Genome"],
            spike_in_genome = lambda wildcards: samples[wildcards.sample]["SpikeInGenome"] if "SpikeInGenome" in samples[wildcards.sample] else "",
    shell:
            """
            bash {params.pipeline_home}/scripts/mergeStats.sh {wildcards.sample} {params.ref_genome} {params.spike_in_genome}
            """
            
rule HiCProStep2:
    input:
            lambda wildcards: wildcards.sample+"/HiCproOUTPUT." + wildcards.genome + "/hic_results/data/" + wildcards.sample + "/" + wildcards.sample + "_" + config[wildcards.genome]["bowtie2_index"] + ".bwt2pairs.validPairs"
    output: 
            "{sample}/HiCproOUTPUT.{genome}/hic_results/data/{sample}/{sample}.allValidPairs",
            "{sample}/HiCproOUTPUT.{genome}/hic_results/stats/{sample}/{sample}_allValidPairs.mergestat"
    version:
            config["version"]["hicpro"]
    params:
            work_dir = config["work_dir"],
            rulename = "HiCProStep2",
            log_dir = lambda wildcards: wildcards.sample + '/log',
            batch    = config["cluster_common"]["medium"]
    benchmark:
            "{sample}/benchmark/HiCProStep2.{genome}.benchmark.txt"
    shell:
            """            
            module load hicpro/{version}
            sed -i 's/$SLURM_SUBMIT_DIR/{wildcards.sample}\/HiCproOUTPUT\.{wildcards.genome}/' {wildcards.sample}/HiCproOUTPUT.{wildcards.genome}/HiCPro_step2_{wildcards.genome}_HiCpro.sh
            bash {wildcards.sample}/HiCproOUTPUT.{wildcards.genome}/HiCPro_step2_{wildcards.genome}_HiCpro.sh
            rm -rf {wildcards.sample}/HiCproOUTPUT.{wildcards.genome}/bowtie_results
            #rm -f {wildcards.sample}/DATA/{wildcards.sample}/*trimmed*
            cd ..
            """

rule HiCProStep1:
    input: 
            "{sample}/HiCproOUTPUT.{genome}/HiCPro_step1_{genome}_HiCpro.sh"
    output: 
            "{sample}/HiCproOUTPUT.{genome}/hic_results/data/{sample}/{sample}_{bowtie_idx}.bwt2pairs.validPairs"
    version:
            config["version"]["hicpro"]
    params:
            work_dir = config["work_dir"],
            rulename = "HiCProStep1",
            log_dir = "{sample}/log",
            batch    = config["cluster"]["hicpro"]
    benchmark:
            "{sample}/benchmark/HiCProStep1.{genome}.{bowtie_idx}.benchmark.txt"
    shell:
            """
            module load hicpro/{version}
            cd {wildcards.sample}/HiCproOUTPUT.{wildcards.genome}
            sed -i 's/$SLURM_SUBMIT_DIR\///' HiCPro_step1_{wildcards.genome}_HiCpro.sh
            bash HiCPro_step1_{wildcards.genome}_HiCpro.sh
            """

rule HiCPro:
    input:
            lambda wildcards: FASTQS[wildcards.sample]
    output: 
            "{sample}/HiCproOUTPUT.{genome}/HiCPro_step1_{genome}_HiCpro.sh",
            "{sample}/HiCproOUTPUT.{genome}/HiCPro_step2_{genome}_HiCpro.sh"
    version:
            config["version"]["hicpro"]
    params:
            out_dir = lambda wildcards: wildcards.sample + '/HiCproOUTPUT.' + wildcards.genome,
            config_file = lambda wildcards: config["pipeline_home"] + '/config/hicpro/' + wildcards.genome + '_' + samples[wildcards.sample]["Digest"] + '.txt'
    shell:
            """
            module load hicpro/{version}
            rm -rf {params.out_dir}            
            HiC-Pro -i {wildcards.sample}/DATA/ -o {params.out_dir}/ -c {params.config_file} -p
            """
