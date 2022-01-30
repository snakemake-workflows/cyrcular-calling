rule copy_annotation_from_cyrcular:
    input:
        variants="results/calling/calls/pre_annotated/{sample}.bcf",
        variants_index="results/calling/calls/pre_annotated/{sample}.bcf.csi",
        candidates_with_annotation="results/calling/candidates/{sample}.sorted.bcf",
        candidates_with_annotation_index="results/calling/candidates/{sample}.sorted.bcf.csi",
    output:
        variants="results/calling/calls/annotated/{sample}.bcf",
    log:
        "logs/re-annotate/{sample}.log",
    benchmark:
        "benchmarks/re-annotate/{sample}.txt"
    conda:
        "../envs/bcftools.yaml"
    params:
        columns=",".join(
            [
                "INFO/CircleLength",
                "INFO/CircleSegmentCount",
                "INFO/SplitReads",
                "INFO/Support",
            ]
        ),
    shell:
        """
        bcftools annotate --annotations {input.candidates_with_annotation} --columns {params.columns} --output-type b --output {output.variants} {input.variants}
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
