rule plot_variants:
    input:
        vcf_file="results/variants.vcf"
    output:
        plot_file="results/variants_scatterplot.png"
    log:
        "logs/plot_variants.log"
    shell:
        """
        Rscript scripts/plot_variants.R {input.vcf_file} {output.plot_file} > {log} 2>&1
        """
