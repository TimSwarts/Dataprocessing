def get_demultiplex_output_files():
    return [f"results/demultiplexed/{sample}_R{read}.fastq" for read in [1, 2] for sample in get_sample_names()]

# def get_demultiplex_output_files():
#     output_files = list()
#     for sample in get_sample_names():
#         output_files.append(f"results/demultiplexed/{sample}_R1.fastq")
#         output_files.append(f"results/demultiplexed/{sample}_R2.fastq")
#     return output_files
    

rule demultiplex_sabre:
    input:
        R1 = f"{data_dir}/{read_name}_R1.fastq",
        R2 = f"{data_dir}/{read_name}_R2.fastq",
        barcodes = f"{data_dir}/{barcode_file}"
    output:
        #directory("results/demultiplexed"),
        demultiplex_flag = touch("results/flags/demultiplexed_done"),
        output_files = get_demultiplex_output_files()
    log:
        f"logs/demultiplexing/{read_name}.log"
    params:
        output_dir = f"{results_dir}/demultiplexed"
    shell:
        """
        mkdir -p {params.output_dir}
        cd {params.output_dir}
        sabre pe -f {input.R1} -r {input.R2} -b {input.barcodes} -u unkown_sample_R1.fastq -w unkown_sample_R2.fastq -c
        """