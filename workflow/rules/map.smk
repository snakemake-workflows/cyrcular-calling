rule merge_fastqs:
    input:
        fastqs=get_fastqs,
    output:
        "results/merged_fastqs/{sample}_{read}.fastq.gz",
    log:
        "logs/merge-fastqs/{sample}_{read}.log",
    wildcard_constraints:
        read="single|R1|R2",
    params:
        cmd=lambda wc: (
            "pigz -dc"
            if (any(map(lambda f: f.endswith(".gz"), get_fastqs(wc))))
            else "cat"
        ),
    conda:
        "../envs/pigz.yaml"
    shell:
        """{params.cmd} {input.fastqs} | pigz -c > {output} 2> {log}"""


rule minimap2_bam:
    input:
        target=rules.minimap2_index.output.index,  # can be either genome index or genome fasta
        query=get_minimap2_input,
    output:
        "results/mapped/{sample}.bam",
    log:
        "logs/minimap2/{sample}.log",
    benchmark:
        "benchmarks/minimap2/{sample}.txt"
    params:
        extra=get_minimap2_mapping_params,  # optional
        sorting="coordinate",  # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra=lambda wc, threads: f"-@ {min(threads , 4)}",  # optional: extra arguments for samtools/picard
    resources:
        mem_mb=lambda wc, input: input.size_mb * 2,
    threads: 8
    wrapper:
        "v1.25.0/bio/minimap2/aligner"


rule samtools_merge_nanopore_per_group:
    input:
        expand(
            "results/mapped/{sample}.bam",
            sample=lookup(query="group == '{group}'", within=samples, cols="sample_name"),
        ),
    output:
        bam="results/mapped/{group}.nanopore_samples.sorted.bam",
    log:
        "logs/calling/mapping/{group}.nanopore_samples.merging.log",
    params:
        extra="",  # optional additional parameters as string
    threads: 4
    wrapper:
        "v4.7.2/bio/samtools/merge"


rule samtools_index:
    input:
        "results/mapped/{file}.bam",
    output:
        "results/mapped/{file}.bam.bai",
    log:
        "logs/mapping/samtools_index/{file}.log",
    benchmark:
        "benchmarks/mapping/samtools_index/{file}.txt"
    params:
        "",  # optional params string
    # Samtools takes additional threads through its option -@
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v1.25.0/bio/samtools/index"


rule samtools_faidx:
    input:
        "{sample}.fasta",
    output:
        "{sample}.fasta.fai",
    log:
        "logs/samtools_index/{sample}.log",
    benchmark:
        "benchmarks/samtools_index/{sample}.txt"
    params:
        "",  # optional params string
    # Samtools takes additional threads through its option -@
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v1.25.0/bio/samtools/faidx"
