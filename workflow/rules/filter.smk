
rule filter_overview_table:
    input:
        table="results/calling/tables/{group}/{group}_overview.tsv",
    output:
        coding="results/calling/tables/{group}/{group}_overview.coding.tsv",
        regulatory="results/calling/tables/{group}/{group}_overview.regulatory.tsv",
        intronic="results/calling/tables/{group}/{group}_overview.intronic.tsv",
        discarded="results/calling/tables/{group}/{group}_overview.discarded.tsv",
    log:
        "logs/filter_overview_table/{group}.log",
    conda:
        "../envs/csvtk.yaml"
    # params:
    #     min_length=f"$circle_length >= {config['calling'].get('filter', {}).get('circles', {}).get('min-length', 0)}",
    #     max_length=f"$circle_length <= {config['calling'].get('filter', {}).get('circles', {}).get('max-length', 18446744073709551616)}",
    shell:
        """
        csvtk tab2csv {input.table} | csvtk filter2 -f '$num_exons > 0' | csvtk csv2tab > {output.coding} 2> {log}
        csvtk tab2csv {input.table} | csvtk filter2 -f '$num_exons <= 0 && $regulatory_features != ""' | csvtk csv2tab > {output.regulatory} 2>> {log}
        csvtk tab2csv {input.table} | csvtk filter2 -f '$num_exons <= 0 && $regulatory_features == "" && $gene_names != ""' | csvtk csv2tab > {output.intronic} 2>> {log}
        csvtk tab2csv {input.table} | csvtk filter2 -f '$num_exons <= 0 && $regulatory_features == "" && $gene_names == ""' | csvtk csv2tab > {output.discarded} 2>> {log}
        """


rule filter_varlociraptor:
    input:
        calls="results/calling/calls/annotated/{group}.bcf",
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
        fdr=config["calling"]["filter"]["fdr-control"].get("threshold", 0.05),
        mode=varlociraptor_filtering_mode,
    shell:
        """
        varlociraptor filter-calls control-fdr {params.mode} --events PRESENT --var BND --fdr {params.fdr} {input.calls} | varlociraptor decode-phred | bcftools sort -m 4G -Ob > {output.fdr_calls} 2> {log}
        """
