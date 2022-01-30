
rule minimap2_bam:
    input:
        target=f"results/calling/index/{REFERENCE}.mmi",  # can be either genome index or genome fasta
        query=get_minimap2_input,
    output:
        "results/calling/mapping/{sample}.bam",
    log:
        "logs/minimap2/{sample}.log",
    benchmark:
        "benchmarks/minimap2/{sample}.txt"
    params:
        extra=get_minimap2_mapping_params,  # optional
        sorting="coordinate",  # optional: Enable sorting. Possible values: 'none', 'queryname' or 'coordinate'
        sort_extra=lambda wc: f"-@ {workflow.cores * 0.5}",  # optional: extra arguments for samtools/picard
    threads: workflow.cores * 0.5
    wrapper:
        "v1.0.0/bio/minimap2/aligner"


rule minimap2_index:
    input:
        target=config["calling"]["reference"]["path"],
    output:
        f"results/calling/index/{REFERENCE}.mmi",
    log:
        f"logs/minimap2_index/{REFERENCE}.log",
    benchmark:
        f"benchmarks/minimap2_index/{REFERENCE}.txt"
    params:
        extra="",  # optional additional args
    cache: True
    threads: workflow.cores
    wrapper:
        "v1.0.0/bio/minimap2/index"


rule merge_fastqs:
    input:
        fastqs=get_fastqs,
    output:
        "results/calling/merged/{sample}_{read}.fastq.gz",
    log:
        "logs/merge-fastqs/{sample}_{read}.log",
    wildcard_constraints:
        read="single|R1|R2",
    run:
        gzipped = any(map(lambda f: f.endswith(".gz"), input.fastqs))
        if len(input.fastqs) == 1 and gzipped:
            shell("ln {input.fastqs[0]} {output}")
        else:
            if gzipped:
                shell("pigz -dc {input.fastqs} | pigz -c > {output} 2> {log}")
            else:
                shell("cat {input.fastqs} | pigz -c > {output} 2> {log}")


rule samtools_index:
    input:
        "results/calling/mapping/{sample}.bam",
    output:
        "results/calling/mapping/{sample}.bam.bai",
    log:
        "logs/samtools_index/{sample}.log",
    benchmark:
        "benchmarks/samtools_index/{sample}.txt"
    params:
        "",  # optional params string
    # Samtools takes additional threads through its option -@
    threads: 12  # This value - 1 will be sent to -@
    wrapper:
        "v1.0.0/bio/samtools/index"


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
    threads: 12  # This value - 1 will be sent to -@
    wrapper:
        "v1.0.0/bio/samtools/faidx"
