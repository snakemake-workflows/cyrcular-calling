## annotate genes, exons and breakpoint sequences; produce one overview table containing circles, and one table with details of each segment for each circle


rule cyrcular_generate_tables:
    input:
        reference=rules.genome_remove_chromosomes.output.fastx,
        graph="results/circle_graphs/{group}.annotated.graph",
        bcf="results/circle_calls_fdr_filtered/reheader/{group}.{event}.bcf",
    output:
        overview="results/circle_tables/{group}/{group}.{event}_overview.tsv",
        details=directory("results/circle_tables/{group}/{group}.{event}_details/"),
    log:
        "logs/cyrcular/generate_tables/generate_tables.{group}.{event}.log",
    benchmark:
        "benchmarks/cyrcular/generate_tables/{group}.{event}.txt"
    conda:
        "../envs/cyrcular.yaml"
    params:
        event_names=lambda wc: ",".join(lookup(dpath=f"filter/fdr-control/events/{wc.event}/varlociraptor", within=config)),
    threads: 1
    shell:
        "( cyrcular graph table "
        "    {input.graph} "
        "    {input.bcf} "
        "    --reference {input.reference} "
        "    --event-names {params.event_names} "
        "    --joint-event-name {wildcards.event} "
        "    --circle-table {output.overview} "
        "    --segment-tables {output.details} "
        ") 2> {log}"""


rule cyrcular_annotate_graph:
    input:
        reference=rules.genome_remove_chromosomes.output.fastx,
        graph="results/circle_graphs/{group}.graph",
        gene_annotation="resources/gene_annotation.gff3.gz",
        regulatory_annotation="resources/regulatory_annotation.gff3.gz",
        repeat_annotation=lambda wc: (
            "resources/repeat_masker.fa.out.gz"
            if config["reference"].get("repeat_masker_download_link", "")
            else ""
        ),
    output:
        annotated="results/circle_graphs/{group}.annotated.graph",
    threads: 1
    log:
        "logs/cyrcular_annotate_graph/{group}.log",
    benchmark:
        "benchmarks/cyrcular_annotate_graph/{group}.txt"
    conda:
        "../envs/cyrcular.yaml"
    params:
        repeat_annotation=lambda wc, input: (
            f"  --repeat-annotation {input.repeat_annotation} "
            if config["reference"].get("repeat_masker_download_link", "")
            else ""
        ),
    resources:
        mem_mb=lambda wc, input: 20 * input.size_mb
    shell:
        "cyrcular graph annotate "
        "  --reference {input.reference} "
        "  --gene-annotation {input.gene_annotation} "
        "  --regulatory-annotation {input.regulatory_annotation} "
        "  --output {output.annotated} "
        "  {input.graph} "
        "2> {log} "


rule reheader_filtered_bcf:
    input:
        bcf="results/circle_calls_fdr_filtered/{group}.{event}.bcf",
        sorted_header="results/circle_calls_fdr_filtered/reheader/{group}.{event}.header.sorted.txt",
    output:
        bcf="results/circle_calls_fdr_filtered/reheader/{group}.{event}.bcf",
    log:
        "logs/reheader_filtered_bcf/{group}.{event}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        ## bcftools re-header seems to re-order entries
        # bcftools reheader --header {input.sorted_header} --output {output.bcf} {input.bcf}
        ## so we have to re-header ourselves
        ## TODO: reinvestigate to find a cleaner solution
        """
        cat {input.sorted_header} <(bcftools view -H {input.bcf}) | bcftools view -Ob > {output.bcf} 2> {log}
        """


rule sort_bcf_header:
    input:
        bcf="results/circle_calls_fdr_filtered/{group}.{event}.bcf",
        header="results/circle_calls_fdr_filtered/{group}.{event}.header.txt",
    output:
        sorted_header="results/circle_calls_fdr_filtered/reheader/{group}.{event}.header.sorted.txt",
    log:
        "logs/sort_bcf_header/{group}.{event}.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/sort_bcf_header.py"


rule get_bcf_header:
    input:
        bcf="{file}.bcf",
    output:
        header="{file}.header.txt",
    log:
        "logs/get_bcf_header/{file}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools view -h {input.bcf} > {output.header} 2> {log}
        """


rule extract_vcf_header_lines_for_bcftools_annotate:
    input:
        vcf="results/candidates/{sample}.sorted.bcf",
    output:
        header=temp("results/calling/annotation/{sample}.header_lines.txt"),
    params:
        fields="|".join(CYRCULAR_INFO_FIELDS),
    conda:
        "../envs/vcf_annotate.yaml"
    log:
        "logs/re-annotate/header_{sample}.log",
    shell:
        """
        bcftools view -h {input.vcf} | rg {params.fields:q} > {output.header} 2> {log}
        """
