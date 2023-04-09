rule call_variants:
    input:
        bam_files=expand("results/mapped_reads/{sample}.bam", sample=get_sample_names()),
        ref_genome=f"{data_dir}/{genome}"
    output:
        "results/variants.vcf"
    log:
        "logs/platypus.log"
    shell:
        "platypus callVariants --bamFiles {input.bam_files} "
        "--refFile {input.ref_genome} "
        "--output {output} "
        "--logFileName platypus.log"
