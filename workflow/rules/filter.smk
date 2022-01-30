
rule filter_calls:
    input:
        "results/calling/calls/annotated/{group}.bcf",
    output:
        calls="results/calling/calls/filtered/{group}.bcf",
        by_exons=temp("results/calling/calls/filtered/{group}.at_least_one_exon.bcf"),
    log:
        "logs/filter-calls/{group}.log",
    benchmark:
        "benchmarks/filter-calls/{group}.txt"
    conda:
        "../envs/vembrane.yaml"
    threads: 1
    params:
        fdr=config["calling"]["filter"]["fdr-control"].get("threshold", 0.05),
        mode=varlociraptor_filtering_mode,
        filter_expression="True",  # 'INFO["NUM_EXONS"] != 0',
    shell:
        """
        vembrane filter '{params.filter_expression}' {input} -O bcf > {output.by_exons} 2> {log}
        varlociraptor filter-calls control-fdr {params.mode} --events PRESENT --var BND --fdr {params.fdr} {output.by_exons} | varlociraptor decode-phred | bcftools sort -m 4G -Ob > {output.calls} 2>> {log}
        """
