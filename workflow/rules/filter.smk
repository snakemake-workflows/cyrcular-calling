
rule filter_overview_table:
    input:
        table="results/calling/tables/{group}/{group}_overview.tsv",
    output:
        coding="results/calling/tables/{group}/{group}_overview.coding.tsv",
        regulatory="results/calling/tables/{group}/{group}_overview.regulatory.tsv",
        intronic="results/calling/tables/{group}/{group}_overview.intronic.tsv",
        discarded="results/calling/tables/{group}/{group}_overview.discarded.tsv",
        categorized="results/calling/tables/{group}/{group}_categorized_overview.tsv",
    log:
        "logs/filter_overview_table/{group}.log",
    run:
        df = pd.read_csv(input.table, sep="\t")
        df.loc[:, ["gene_names", "gene_ids", "regulatory_features"]] = df.loc[
            :, ["gene_names", "gene_ids", "regulatory_features"]
        ].fillna("")
        df["category"] = df[["num_exons", "regulatory_features", "gene_names"]].apply(
            lambda r: "coding"
            if r["num_exons"] > 0
            else (
                "regulatory"
                if r["regulatory_features"]
                else ("intronic" if r["gene_names"] else "discarded")
            ),
            axis=1,
        )
        for kind in ["coding", "regulatory", "intronic", "discarded"]:
            part = df.query(f"category == '{kind}'")
            part.to_csv(getattr(output, kind), sep="\t", index=False)
        df.to_csv(output.categorized, sep="\t", index=False)


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
