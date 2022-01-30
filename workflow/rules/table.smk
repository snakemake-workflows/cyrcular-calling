from pathlib import Path


rule annotate_breakpoint_sequence:
    input:
        csv="results/calling/tables/{group}.sorted.csv",
        reference=config["calling"]["reference"]["path"],
    output:
        csv="results/calling/tables/{group}.sorted.annotated.csv",
    log:
        "logs/sort_and_tidy_tsv/{group}.log",
    benchmark:
        "benchmarks/sort_and_tidy_tsv/{group}.txt"
    threads: 8
    conda:
        "../envs/table.yaml"
    script:
        "../scripts/annotate_breakpoint_sequence.py"


rule sort_and_tidy_tsv:
    input:
        tsv="results/calling/tables/{group}.tsv",
    output:
        csv="results/calling/tables/{group}.sorted.csv",
    log:
        "logs/sort_and_tidy_tsv/{group}.log",
    benchmark:
        "benchmarks/sort_and_tidy_tsv/{group}.txt"
    conda:
        "../envs/table.yaml"
    script:
        "../scripts/sort_and_tidy_tsv.py"


rule vembrane_table:
    input:
        vcf="results/calling/calls/filtered/{group}.bcf",
    output:
        tsv="results/calling/tables/{group}.tsv",
    threads: 1
    params:
        expression=", ".join(
            [
                'INFO["EVENT"]',
                "ID",
                "CHROM",
                "POS",
                "ALT",
                'INFO["CircleLength"]',
                'INFO["CircleSegmentCount"]',
                'INFO["SplitReads"]',
                'INFO["NUM_EXONS"]',
                '";".join(INFO["GENES"] or [])',
                'for_each_sample(lambda s: FORMAT["AF"][s])',
                'INFO["PROB_PRESENT"]',
                'INFO["PROB_ABSENT"]',
                'INFO["PROB_ARTIFACT"]',
            ]
        ),
        extra=lambda wc: "--header \"EVENT, ID, CHROM, POS, ALT, circle_length, num_segments, split_reads, NUM_EXONS, GENES, for_each_sample(lambda s: 'AF_' + s), PROB_PRESENT, PROB_ABSENT, PROB_ARTIFACT\"",
    log:
        "logs/vembrane_table/{group}.log",
    benchmark:
        "benchmarks/vembrane_table/{group}.txt"
    conda:
        "../envs/vembrane.yaml"
    shell:
        """vembrane table {params.extra} '{params.expression}' {input.vcf} > {output.tsv}"""
