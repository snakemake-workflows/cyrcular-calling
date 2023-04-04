## annotate genes, exons and breakpoint sequences; produce one overview table containing circles, and one table with details of each segment for each circle


rule cyrcular_generate_tables:
    input:
        reference="resources/genome.fasta",
        graph="results/calling/graphs/{group}.annotated.graph",
        bcf="results/calling/calls/filtered_fdr/reheader/{group}.bcf",
    output:
        overview="results/calling/tables/{group}/{group}_overview.tsv",
        details=directory("results/calling/tables/{group}/{group}_details/"),
    threads: 1
    log:
        "logs/cyrcular_generate_tables/{group}.log",
    benchmark:
        "benchmarks/cyrcular_generate_tables/{group}.txt"
    conda:
        "../envs/cyrcular.yaml"
    shell:
        """cyrcular graph table {input.graph} {input.bcf} --reference {input.reference} --circle-table {output.overview} --segment-tables {output.details} 2> {log}"""


rule cyrcular_annotate_graph:
    input:
        reference="resources/genome.fasta",
        graph="results/calling/graphs/{group}.graph",
        gene_annotation="resources/gene_annotation.gff3.gz",
        regulatory_annotation="resources/regulatory_annotation.gff3.gz",
        repeat_annotation="resources/repeat_masker.hg38.fa.out.gz",
    output:
        annotated="results/calling/graphs/{group}.annotated.graph",
    threads: 1
    log:
        "logs/cyrcular_annotate_graph/{group}.log",
    benchmark:
        "benchmarks/cyrcular_annotate_graph/{group}.txt"
    conda:
        "../envs/cyrcular.yaml"
    shell:
        """
        cyrcular graph annotate --reference {input.reference} --gene-annotation {input.gene_annotation} --regulatory-annotation {input.regulatory_annotation} --repeat-annotation {input.repeat_annotation} --output {output.annotated} {input.graph} 2> {log}
        """


rule reheader_filtered_bcf:
    input:
        bcf="results/calling/calls/filtered_fdr/{group}.bcf",
        sorted_header="results/calling/calls/filtered_fdr/reheader/{group}.header.sorted.txt",
    output:
        bcf="results/calling/calls/filtered_fdr/reheader/{group}.bcf",
    log:
        "logs/reheader_filtered_bcf/{group}.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        ## bcftools re-header seems to re-order entries
        # bcftools reheader --header {input.sorted_header} --output {output.bcf} {input.bcf}
        ## so we have to re-header ourselves
        """
        cat {input.sorted_header} <(bcftools view -H {input.bcf}) | bcftools view -Ob > {output.bcf} 2> {log}
        """


rule sort_bcf_header:
    input:
        bcf="results/calling/calls/filtered_fdr/{group}.bcf",
        header="results/calling/calls/filtered_fdr/{group}.header.txt",
    output:
        sorted_header="results/calling/calls/filtered_fdr/reheader/{group}.header.sorted.txt",
    log:
        "logs/sort_bcf_header/{group}.log",
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
        vcf="results/calling/candidates/{sample}.sorted.bcf",
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


rule download_regulatory_annotation:
    output:
        "resources/regulatory_annotation.gff3.gz",
    log:
        "logs/download_regulatory_annotation.log",
    params:
        release=get_annotation_release,
    benchmark:
        "benchmarks/download_regulatory_annotation.txt"
    cache: True
    conda:
        "../envs/wget.yaml"
    shell:
        """wget https://ftp.ensembl.org/pub/release-{params.release}/regulation/homo_sapiens/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20220201.gff.gz --no-check-certificate -O {output} 2> {log}"""


rule download_repeatmasker_annotation:
    output:
        "resources/repeat_masker.hg38.fa.out.gz",
    log:
        "logs/download_repeatmasker_annotation.log",
    params:
        release=get_annotation_release,
    benchmark:
        "benchmarks/download_repeatmasker_annotation.txt"
    cache: True
    conda:
        "../envs/wget.yaml"
    shell:
        """wget https://repeatmasker.org/genomes/hg38/RepeatMasker-rm405-db20140131/hg38.fa.out.gz --no-check-certificate -O {output} 2> {log}"""


rule download_gene_annotation:
    output:
        "resources/gene_annotation.gff3.gz",
    log:
        "logs/download_gene_annotation.log",
    params:
        release=get_annotation_release,
    benchmark:
        "benchmarks/download_gene_annotation.txt"
    cache: True
    conda:
        "../envs/wget.yaml"
    shell:
        """wget https://ftp.ensembl.org/pub/release-{params.release}/gff3/homo_sapiens/Homo_sapiens.GRCh38.{params.release}.gff3.gz --no-check-certificate -O {output} 2> {log}"""
