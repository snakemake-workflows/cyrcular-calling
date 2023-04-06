
rule bcf_index:
    input:
        "{prefix}.bcf",
    output:
        "{prefix}.bcf.csi",
    threads: 2
    log:
        "logs/bcftools/index/{prefix}.log",
    benchmark:
        "benchmarks/bcftools/index/{prefix}.txt"
    wrapper:
        "v1.25.0/bio/bcftools/index"


rule bcftools_concat:
    input:
        calls=gather.calling(
            "results/calling/calls/initial_sorted/{{group}}.{scatteritem}.bcf"
        ),
        indexes=gather.calling(
            "results/calling/calls/initial_sorted/{{group}}.{scatteritem}.bcf.csi"
        ),
    output:
        "results/calling/calls/merged/{group}.bcf",
    log:
        "logs/concat_calls/{group}.log",
    benchmark:
        "benchmarks/concat_calls/{group}.txt"
    threads: 2
    params:
        uncompressed_bcf=False,
        extra="-a",  # optional parameters for bcftools concat (except -o)
    wrapper:
        "v1.25.0/bio/bcftools/concat"


rule bcftools_sort:
    input:
        "results/calling/calls/initial/{group}.{scatteritem}.bcf",
    output:
        "results/calling/calls/initial_sorted/{group}.{scatteritem}.bcf",
    log:
        "logs/sort_calls/{group}.{scatteritem}.log",
    benchmark:
        "benchmarks/sort_calls/{group}.{scatteritem}.txt"
    wrapper:
        "v1.25.0/bio/bcftools/sort"


rule varlociraptor_call:
    input:
        obs=get_observations,
        scenario=get_scenario,
    output:
        "results/calling/calls/initial/{group}.{scatteritem}.bcf",
    log:
        "logs/varlociraptor/call/{group}.{scatteritem}.log",
    conda:
        "../envs/varlociraptor.yaml"
    threads: 1
    benchmark:
        "benchmarks/varlociraptor/call/{group}.{scatteritem}.txt"
    params:
        obs=observation_string,
    shell:
        "varlociraptor "
        "call variants generic --obs {params.obs} "
        "--scenario {input.scenario} > {output} 2> {log}"


rule varlociraptor_alignment_properties:
    input:
        ref=rules.get_genome.output.genome,
        ref_idx=rules.genome_faidx.output.index,
        bam="results/calling/mapping/{sample}.bam",
    output:
        "results/calling/alignment-properties/{sample}.json",
    log:
        "logs/varlociraptor/estimate-alignment-properties/{sample}.log",
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor estimate alignment-properties {input.ref} --bam {input.bam} > {output} 2> {log}"


rule varlociraptor_preprocess:
    input:
        ref=rules.get_genome.output.genome,
        ref_idx=rules.genome_faidx.output.index,
        candidates=get_group_candidates,
        bam="results/calling/mapping/{sample}.bam",
        bai="results/calling/mapping/{sample}.bam.bai",
        alignment_props="results/calling/alignment-properties/{sample}.json",
    output:
        "results/calling/calls/observations/{sample}.{scatteritem}.bcf",
    params:
        mode=pairhmm_mode,
    log:
        "logs/varlociraptor/preprocess/{sample}.{scatteritem}.log",
    benchmark:
        "benchmarks/varlociraptor/preprocess/{sample}.{scatteritem}.txt"
    threads: 1
    benchmark:
        "benchmarks/varlociraptor/preprocess/{sample}.{scatteritem}.tsv"
    conda:
        "../envs/varlociraptor.yaml"
    shell:
        "varlociraptor preprocess variants {input.ref} "
        "--candidates {input.candidates} "
        "--max-depth 200 "
        "--alignment-properties {input.alignment_props} "
        "--pairhmm-mode {params.mode} "
        "--bam {input.bam} --output {output} 2> {log}"


rule scatter_candidates:
    input:
        "results/calling/candidates/{sample}.sorted.bcf",
    output:
        scatter.calling("results/calling/candidate-calls/{{sample}}.{scatteritem}.bcf"),
    log:
        "logs/scatter_candidates/{sample}.log",
    benchmark:
        "benchmarks/scatter_candidates/{sample}.txt"
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-split {input} {output}"


rule sort_bnd_bcfs:
    input:
        "results/calling/candidates/{sample}.bcf",
    output:
        "results/calling/candidates/{sample}.sorted.bcf",
    log:
        "logs/sort_bnds/{sample}.log",
    benchmark:
        "benchmarks/sort_bnds/{sample}.txt"
    wrapper:
        "v1.25.0/bio/bcftools/sort"


rule circle_bnds:
    input:
        bam="results/calling/mapping/{sample}.bam",
        bai="results/calling/mapping/{sample}.bam.bai",
        ref=rules.get_genome.output.genome,
        ref_index=rules.genome_faidx.output.index,
    output:
        bnds="results/calling/candidates/{sample}.bcf",
        graph="results/calling/graphs/{sample}.graph",
        dot=directory("results/calling/graphs/{sample}/"),
    log:
        "logs/cyrcular/{sample}.log",
    benchmark:
        "benchmarks/cyrcular/{sample}.txt"
    threads: 4
    params:
        min_read_depth=config["cyrcular"]["min_read_depth"],  # 2
        min_split_reads=config["cyrcular"]["min_split_reads"],  # 5
        max_paths_per_component=config["cyrcular"]["max_paths_per_component"],  # 15
        max_deletion_length=config["cyrcular"]["max_deletion_length"],  # 10000,
    conda:
        "../envs/cyrcular.yaml"
    shell:
        """cyrcular\
        graph\
        breakends\
        {input.bam}\
        --reference {input.ref}\
        --min-read-depth {params.min_read_depth}\
        --min-split-reads {params.min_split_reads}\
        --max-paths-per-component {params.max_paths_per_component}\
        --max-deletion-length {params.max_deletion_length}\
        -t {threads}\
        --output {output.bnds}\
        --graph {output.graph}\
        --dot {output.dot}\
        2> {log}"""
