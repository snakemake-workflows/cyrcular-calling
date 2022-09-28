## annotate genes, exons and breakpoint sequences; produce one overview table containing circles, and one table with details of each segment for each circle


rule cyrcular_annotate:
    input:
        reference=config["calling"]["reference"]["path"],
        graph="results/calling/graphs/{group}.graph",
        bcf="results/calling/calls/filtered/reheader/{group}.bcf",
        annotation="resources/gencode.annotation.sorted.gff3.gz",
    output:
        overview="results/calling/tables/{group}/{group}_overview.tsv",
        details=directory("results/calling/tables/{group}/{group}_details/"),
    threads: 1
    log:
        "logs/cyrcular_annotate/{group}.log",
    benchmark:
        "benchmarks/cyrcular_annotate/{group}.txt"
    conda:
        "../envs/cyrcular.yaml"
    shell:
        """cyrcular graph annotate {input.graph} {input.bcf} --reference {input.reference} --annotation {input.annotation} --circle-table {output.overview} --segment-tables {output.details} 2> {log}"""


rule reheader_filtered_bcf:
    input:
        bcf="results/calling/calls/filtered/{group}.bcf",
        sorted_header="results/calling/calls/filtered/reheader/{group}.header.sorted.txt",
    output:
        bcf="results/calling/calls/filtered/reheader/{group}.bcf",
    conda:
        "../envs/bcftools.yaml"
    shell:
        ## bcftools re-header seems to re-order entries
        # bcftools reheader --header {input.sorted_header} --output {output.bcf} {input.bcf}
        ## so we have to re-header ourselves
        """
        cat {input.sorted_header} <(bcftools view -H {input.bcf}) | bcftools view -Ob > {output.bcf}
        """


rule sort_bcf_header:
    input:
        bcf="results/calling/calls/filtered/{group}.bcf",
        header="results/calling/calls/filtered/{group}.header.txt",
    output:
        sorted_header="results/calling/calls/filtered/reheader/{group}.header.sorted.txt",
    run:
        with open(input.header, "rt") as header_file:
            header = [l.strip() for l in header_file.readlines()]
            file_format_line = header[0]
            chrom_line = header[-1]
            other_lines = header[1:-1]
            kinds = [
                "fileDate",
                "source",
                "reference",
                "contig",
                "phasing",
                "FILTER",
                "INFO",
                "FORMAT",
                "ALT",
                "assembly",
                "META",
                "SAMPLE",
                "PEDIGREE",
                "pedigreeDB",
            ]
            categories = {kind: [] for kind in kinds}
            others = []
            for line in other_lines:
                if "=" in line:
                    kind = line.split("=")[0].lstrip("#")
                    group = categories.get(kind, others)
                else:
                    group = others
                group.append(line)

            with open(output.sorted_header, "wt") as out:
                print(file_format_line, file=out)
                for kind in kinds:
                    lines = categories[kind]
                    for line in lines:
                        print(line, file=out)
                for line in others:
                    print(line, file=out)
                print(chrom_line, file=out)


rule get_bcf_header:
    input:
        bcf="{file}.bcf",
    output:
        header="{file}.header.txt",
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools view -h {input.bcf} > {output.header}
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


rule copy_annotation_from_cyrcular:
    input:
        variants="results/calling/calls/merged/{sample}.bcf",
        variants_index="results/calling/calls/merged/{sample}.bcf.csi",
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
        header=copy_annotation_vembrane_header_expr(),
        table_expr=copy_annotation_table_expr(),
        columns=copy_annotation_bcftools_annotate_columns(),
    shell:
        """
        vembrane table --header {params.header:q} {params.table_expr:q} {input.candidates_with_annotation} > {output.annotation} 2> {log}
        bgzip -c {output.annotation} > {output.annotation_bgzip} 2>> {log}
        tabix -p vcf --zero-based -S 1 -f {output.annotation_bgzip} 2>> {log}
        bcftools annotate --header-lines {input.header_lines} --annotations {output.annotation_bgzip} --columns {params.columns} --output-type b --output {output.variants} {input.variants}  2>> {log}
        """


rule sort_annotation:
    input:
        gff="resources/gencode.annotation.gff3.gz",
    output:
        gff="resources/gencode.annotation.sorted.gff3.gz",
        tbi="resources/gencode.annotation.sorted.gff3.gz.tbi",
        csi="resources/gencode.annotation.sorted.gff3.gz.csi",
    conda:
        "../envs/gff.yaml"
    log:
        "logs/sort_annotation/gencode.annotation.gff3.log",
    benchmark:
        "benchmarks/sort_annotation/gencode.annotation.gff3.txt"
    threads: 48
    shell:
        """
        gff3sort.pl <(pigz -dc {input.gff}) | bgzip -@ {threads} -c > {output.gff} 2> {log}
        tabix {output.gff} 2>> {log}
        tabix --csi {output.gff} 2>> {log}
        """


rule download_annotation:
    output:
        "resources/gencode.annotation.gff3.gz",
    log:
        "logs/download_annotation.log",
    params:
        release=get_annotation_release,
    benchmark:
        "benchmarks/download_annotation.txt"
    cache: True
    conda:
        "../envs/wget.yaml"
    shell:
        """wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{params.release}/gencode.v{params.release}.annotation.gff3.gz --no-check-certificate -O {output} 2> {log}"""
