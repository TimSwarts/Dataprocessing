def get_demultiplex_output_files():
    return [f"{results_dir}/demultiplexed/{sample}_R{read}.fastq" for read in [1, 2] for sample in get_sample_names()]
    

rule demultiplex_sabre:
    input:
        R1 = f"{data_dir}/{read_name}_R1.fastq",
        R2 = f"{data_dir}/{read_name}_R2.fastq",
        barcodes = f"{data_dir}/{barcode_file}"
    output:
        demultiplex_flag = touch("results/flags/demultiplexed_done"),
        output_files = get_demultiplex_output_files()
    log:
        f"logs/demultiplex_sabre/{read_name}.log"
    params:
        output_dir = f"{results_dir}/demultiplexed"
    shell:
        """
        (
        mkdir -p {params.output_dir}
        cd {params.output_dir}
        sabre pe -f {input.R1} -r {input.R2} -b {input.barcodes} -u unkown_sample_R1.fastq -w unkown_sample_R2.fastq -c
        ) > {log} 2>&1
        """