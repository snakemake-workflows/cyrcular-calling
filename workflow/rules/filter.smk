
rule filter_overview_table:
    input:
        table="results/calling/tables/{group}/{group}_overview.tsv",
    output:
        coding="results/calling/tables/{group}/{group}_overview.coding.tsv",
        regulatory="results/calling/tables/{group}/{group}_overview.regulatory.tsv",
        intronic="results/calling/tables/{group}/{group}_overview.intronic.tsv",
        other="results/calling/tables/{group}/{group}_overview.other.tsv",
        categorized="results/calling/tables/{group}/{group}_categorized_overview.tsv",
    log:
        "logs/filter_overview_table/{group}.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/filter_overview_table.py"


rule filter_varlociraptor:
    input:
        calls="results/calling/calls/merged/{group}.bcf",
    output:
        fdr_calls="results/calling/calls/filtered_fdr/{group}.bcf",
    log:
        "logs/filter-varlociraptor/{group}.log",
    benchmark:
        "benchmarks/filter-varlociraptor/{group}.txt"
    conda:
        "../envs/varlociraptor.yaml"
    threads: 1
    params:
        fdr=config["filter"]["fdr-control"].get("threshold", 0.05),
        mode=varlociraptor_filtering_mode,
    shell:
        """
        varlociraptor filter-calls control-fdr --mode {params.mode} --events PRESENT --var BND --fdr {params.fdr} {input.calls} | varlociraptor decode-phred | bcftools sort -m 4G -Ob > {output.fdr_calls} 2> {log}
        """
