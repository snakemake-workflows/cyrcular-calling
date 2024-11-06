rule mosdepth_coverage:
    input:
        bam="results/mapped/{sample}.bam",
        bai="results/mapped/{sample}.bam.bai",
    output:
        "results/coverage/{sample}.mosdepth.global.dist.txt",
        summary="results/coverage/{sample}.mosdepth.summary.txt",  # this named output is required for prefix parsing
    log:
        "logs/coverage/mosdepth_{sample}.log",
    params:
        extra="--no-per-base",  # optional
    resources:
        mem_mb=18000
    # additional decompression threads through `--threads`
    threads: 4  # This value - 1 will be sent to `--threads`
    wrapper:
        "v5.0.1/bio/mosdepth"