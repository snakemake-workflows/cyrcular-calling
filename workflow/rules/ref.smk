rule get_genome:
    output:
        genome=expand(
            "resources/{species}.{build}.{release}.fasta",
            species=config["reference"]["species"],
            build=config["reference"]["build"],
            release=config["reference"]["release"],
        ),
    log:
        "logs/get-genome.log",
    params:
        species=lookup(dpath="reference/species", within=config),
        datatype="dna",
        build=lookup(dpath="reference/build", within=config),
        release=lookup(dpath="reference/release", within=config),
        chromosome=lookup(dpath="reference/chromosome", within=config, default=[]),  # optional: restrict to one or multiple chromosomes
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v4.7.1/bio/reference/ensembl-sequence"


rule genome_faidx:
    input:
        rules.get_genome.output.genome,
    output:
        index=expand(
            "{genome}.fai",
            genome=rules.get_genome.output.genome,
        ),
    log:
        "logs/genome-faidx.log",
    cache: True
    wrapper:
        "v1.25.0/bio/samtools/faidx"


rule minimap2_index:
    input:
        target=rules.get_genome.output.genome,
    output:
        index=expand(
            "resources/{species}.{build}.{release}.mmi",
            species=config["reference"]["species"],
            build=config["reference"]["build"],
            release=config["reference"]["release"],
        ),
    log:
        "logs/minimap2_index/genome.log",
    params:
        extra="",  # optional additional args
    cache: True
    # Minimap2 uses at most three threads when indexing target sequences:
    # https://lh3.github.io/minimap2/minimap2.html
    threads: 3
    wrapper:
        "v1.25.0/bio/minimap2/index"


rule download_regulatory_annotation:
    output:
        "resources/regulatory_annotation.gff3.gz",
    params:
        species=lookup(dpath="reference/species", within=config),
        build=lookup(dpath="reference/build", within=config),
        release=lookup(dpath="reference/release", within=config),
    log:
        "logs/download_regulatory_annotation.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs); for data downloads, software does not affect the resulting data
    wrapper:
        "v4.7.2/bio/reference/ensembl-regulation"



rule download_repeatmasker_annotation:
    output:
        "resources/repeat_masker.fa.out.gz",
    log:
        "logs/download_repeatmasker_annotation.log",
    params:
        download_link=config["reference"].get("repeat_masker_download_link", ""),
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    conda:
        "../envs/wget.yaml"
    shell:
        """wget {params.download_link} --no-check-certificate -O {output} 2> {log}"""


rule download_gene_annotation:
    output:
        "resources/gene_annotation.gff3.gz",
    params:
        species=config["reference"]["species"],
        build=config["reference"]["build"],
        release=config["reference"]["release"],
        flavor="",  # optional, e.g. chr_patch_hapl_scaff, see Ensembl FTP.
        branch="",  # optional: specify branch
    log:
        "logs/download_gene_annotation.log",
    cache: "omit-software"  # save space and time with between workflow caching (see docs)
    wrapper:
        "v1.25.0/bio/reference/ensembl-annotation"
