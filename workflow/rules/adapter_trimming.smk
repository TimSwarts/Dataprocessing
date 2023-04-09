# Set rule variables
trim_minimum = config["consider_already_trimmed"]
length = config["trim_galore_min_length"]

rule trim_galore:
    input:
        R1 = f"{results_dir}/demultiplexed/{{sample}}_R1.fastq",
        R2 = f"{results_dir}/demultiplexed/{{sample}}_R2.fastq",
        demultiplex_flag = "results/flags/demultiplexed_done"
    output:
        R1 = f"{results_dir}/trimmed/{{sample}}_R1.fastq",
        R2 = f"{results_dir}/trimmed/{{sample}}_R2.fastq"
    log:
        "logs/trim_galore/{sample}.log"
    params:
        output_dir = f"{results_dir}/trimmed",
        trim_minimum = trim_minimum,
        length = length
    run:
        try:
            shell("""
            (
            trim_galore \
                --paired \
                --dont_gzip \
                --consider_already_trimmed {params.trim_minimum} \
                --length {params.length} \
                --no_report_file \
                --suppress_warn \
                --output_dir {params.output_dir} \
                --basename {wildcards.sample} \
                {input.R1} \
                {input.R2} \
            ) > {log} 2>&1
            """)
        except:
            shell("cp {input.R1} > {output.R1}; cp {input.R2} > {output.R2}")