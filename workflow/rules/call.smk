
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
            "results/circle_calls_initial/scattered/{{group}}.{scatteritem}.sorted.bcf"
        ),
        indexes=gather.calling(
            "results/circle_calls_initial/scattered/{{group}}.{scatteritem}.sorted.bcf.csi"
        ),
    output:
        "results/circle_calls_initial/{group}.bcf",
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
        "results/circle_calls_initial/scattered/{group}.{scatteritem}.unsorted.bcf",
    output:
        "results/circle_calls_initial/scattered/{group}.{scatteritem}.sorted.bcf",
    log:
        "logs/sort_calls/scattered/{group}.{scatteritem}.log",
    benchmark:
        "benchmarks/sort_calls/scattered/{group}.{scatteritem}.txt"
    wrapper:
        "v1.25.0/bio/bcftools/sort"


rule varlociraptor_call:
    input:
        obs=expand(
            "results/observations/{sample}.{{scatteritem}}.bcf",
            sample=lookup(query="group == '{group}'", within=samples, cols="sample_name")
        ),
        scenario="results/scenarios/{group}.yaml",
    output:
        temp("results/circle_calls_initial/scattered/{group}.{scatteritem}.unsorted.bcf"),
    log:
        "logs/varlociraptor/calls/scattered/{group}.{scatteritem}.log",
    params:
        obs=get_varlociraptor_obs_args,
    conda:
        "../envs/varlociraptor.yaml"
    benchmark:
        "benchmarks/varlociraptor/calls/scattered/{group}.{scatteritem}.tsv"
    threads: 1
    shell:
        "( varlociraptor call variants generic "
        "   --obs {params.obs} "
        "   --scenario {input.scenario} "
        "   {output} "
        ") 2> {log}"


rule render_scenario:
    input:
        template=lookup(dpath="calling_scenario", within=config),
    output:
        report(
            "results/scenarios/{group}.yaml",
            caption="../report/scenario.rst",
            category="Circle calling scenarios",
            labels={"sample group": "{group}"},
        ),
    log:
        "logs/render-scenario/{group}.log",
    conda:
        None
    template_engine:
        "yte"


rule varlociraptor_alignment_properties:
    input:
        ref=rules.get_genome.output.genome,
        ref_idx=rules.genome_faidx.output.index,
        bam="results/mapped/{sample}.bam",
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
        candidates=expand(
            "results/candidates/scattered/{group}.{{scatteritem}}.bcf",
            group=lookup(query="sample_name == '{sample}'", within=samples, cols="group"),
        ),
        bam="results/mapped/{sample}.bam",
        bai="results/mapped/{sample}.bam.bai",
        alignment_props="results/calling/alignment-properties/{sample}.json",
    output:
        "results/observations/{sample}.{scatteritem}.bcf",
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
        "results/candidates/{group}.sorted.bcf",
    output:
        scatter.calling("results/candidates/scattered/{{group}}.{scatteritem}.bcf"),
    log:
        "logs/candidates/scattered/scatter.{group}.log",
    benchmark:
        "benchmarks/candidates/scattered/scatter.{group}.txt"
    conda:
        "../envs/rbt.yaml"
    shell:
        "rbt vcf-split {input} {output}"


rule sort_bnd_bcfs:
    input:
        "results/candidates/{group}.bcf",
    output:
        "results/candidates/{group}.sorted.bcf",
    log:
        "logs/candidates/sort_bnds.{group}.log",
    benchmark:
        "benchmarks/candidates/sort_bnds.{group}.txt"
    wrapper:
        "v1.25.0/bio/bcftools/sort"


rule cyrcular_call_circle_bnds:
    input:
        bam="results/mapped/{group}.nanopore_samples.sorted.bam",
        bai="results/mapped/{group}.nanopore_samples.sorted.bam.bai",
        ref=rules.get_genome.output.genome,
        ref_index=rules.genome_faidx.output.index,
    output:
        bnds="results/candidates/{group}.bcf",
        graph="results/circle_graphs/{group}.graph",
        dot=directory("results/circle_graphs/{group}/"),
    log:
        "logs/cyrcular/call_circle_bnds/call_circle_bnds.{group}.log",
    benchmark:
        "benchmarks/cyrcular/call_circle_bnds/call_circle_bnds.{group}.txt"
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
