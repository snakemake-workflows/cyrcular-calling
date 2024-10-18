
rule filter_overview_table:
    input:
        table="results/circle_tables/{group}/{group}_overview.tsv",
    output:
        coding="results/circle_tables/{group}/{group}_overview.coding.tsv",
        regulatory="results/circle_tables/{group}/{group}_overview.regulatory.tsv",
        intronic="results/circle_tables/{group}/{group}_overview.intronic.tsv",
        other="results/circle_tables/{group}/{group}_overview.other.tsv",
        categorized="results/circle_tables/{group}/{group}_categorized_overview.tsv",
    log:
        "logs/filter_overview_table/{group}.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/filter_overview_table.py"


rule varlociraptor_filter_fdr:
    input:
        calls="results/circle_calls_initial/{group}.bcf",
    output:
        fdr_calls="results/circle_calls_fdr_filtered/{group}.bcf",
    log:
        "logs/varlociraptor/filter_fdr/filter_fdr.{group}.log",
    benchmark:
        "benchmarks/varlociraptor/filter_fdr/filter_fdr.{group}.txt"
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
