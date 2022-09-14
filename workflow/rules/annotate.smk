CYRCULAR_INFO_FIELDS = ["CircleLength", "CircleSegmentCount", "SplitReads", "Support"]


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


rule copy_annotation_from_cyrcular:
    input:
        variants="results/calling/calls/pre_annotated/{sample}.bcf",
        variants_index="results/calling/calls/pre_annotated/{sample}.bcf.csi",
        candidates_with_annotation="results/calling/candidates/{sample}.sorted.bcf",
        candidates_with_annotation_index="results/calling/candidates/{sample}.sorted.bcf.csi",
        header_lines="results/calling/annotation/{sample}.header_lines.txt",
    output:
        variants="results/calling/calls/annotated/{sample}.bcf",
        annotation=temp("results/calling/candidates/{sample}.sorted.bcf.tab"),
        annotation_bgzip=temp("results/calling/candidates/{sample}.sorted.bcf.tab.bgz"),
        annotation_bgzip_tabix=temp(
            "results/calling/candidates/{sample}.sorted.bcf.tab.bgz.tbi"
        ),
    log:
        "logs/re-annotate/{sample}.log",
    benchmark:
        "benchmarks/re-annotate/{sample}.txt"
    conda:
        "../envs/vcf_annotate.yaml"
    params:
        header="CHROM,POS,ID,REF,ALT," + ",".join(CYRCULAR_INFO_FIELDS),
        table_expr="CHROM,POS,ID,REF,ALT," + ",".join(map(lambda s: "INFO['" + s + "']", CYRCULAR_INFO_FIELDS)),
        columns="CHROM,POS,~ID,REF,ALT," + ",".join(CYRCULAR_INFO_FIELDS),
    shell:
        """
        vembrane table --header {params.header:q} {params.table_expr:q} {input.candidates_with_annotation} > {output.annotation} 2> {log}
        bgzip -c {output.annotation} > {output.annotation_bgzip} 2>> {log}
        tabix -p vcf --zero-based -S 1 -f {output.annotation_bgzip} 2>> {log}
        bcftools annotate --header-lines {input.header_lines} --annotations {output.annotation_bgzip} --columns {params.columns} --output-type b --output {output.variants} {input.variants}  2>> {log}
        """


rule annotate_genes:
    input:
        annotation="resources/gencode.v38.annotation.sorted.gff3.gz",
        variants="results/calling/calls/merged/{sample}.bcf",
    output:
        variants="results/calling/calls/pre_annotated/{sample}.bcf",
    threads: 1
    log:
        "logs/annotate_genes/{sample}.log",
    benchmark:
        "benchmarks/annotate_genes/{sample}.txt"
    conda:
        "../envs/gff.yaml"
    script:
        "../scripts/gff_annotate.py"


rule sort_annotation:
    input:
        "resources/gencode.v38.annotation.gff3.gz",
    output:
        gff="resources/gencode.v38.annotation.sorted.gff3.gz",
        tbi="resources/gencode.v38.annotation.sorted.gff3.gz.tbi",
        csi="resources/gencode.v38.annotation.sorted.gff3.gz.csi",
    conda:
        "../envs/gff.yaml"
    log:
        "logs/sort_annotation/gencode.v38.annotation.gff3.log",
    benchmark:
        "benchmarks/sort_annotation/gencode.v38.annotation.gff3.txt"
    threads: 48
    shell:
        """
        gff3sort.pl <(pigz -dc {input}) | bgzip -@ {threads} -c > {output.gff} 2> {log}
        tabix {output.gff} 2>> {log}
        tabix --csi {output.gff} 2>> {log}
        """


rule download_annotation:
    output:
        "resources/gencode.v38.annotation.gff3.gz",
    log:
        "logs/download_annotation/gencode.v38.annotation.gff3.log",
    benchmark:
        "benchmarks/download_annotation/gencode.v38.annotation.gff3.txt"
    cache: True
    conda:
        "../envs/wget.yaml"
    shell:
        """wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gff3.gz --no-check-certificate -O {output} 2> {log}"""
