rule bwa_mem:
    input:
        R1 = "results/demultiplexed/{sample}_R1.fastq",
        R2 = "results/demultiplexed/{sample}_R2.fastq",
        index_flag = "results/flags/genome_indexed",
        demultiplex_flag = "results/flags/demultiplexed_done"
    output:
        "results/mapped_reads/{sample}.sam"
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        index_prefix=f"{data_dir}/{genome}"
    shell:
        """
        (
        bwa mem -t4 {params.index_prefix} {input.R1} {input.R2} > '{output}'
        ) > {log} 2>&1
        """

rule bwa_index:
    input:
        genome_file=data_dir + "/" + genome
    output:
        index_flag = touch("results/flags/genome_indexed")
    # log:
    #     "logs/bwa_index.log"
    message: 
        "Indexing the genome"
    shell:
        """
        (
        if [ "$(ls {data_dir} | grep "{genome}" | wc -l)" -lt 6 ]; then
            bwa index {input.genome_file}
        fi
        ) > {log} 2>&1
        """
