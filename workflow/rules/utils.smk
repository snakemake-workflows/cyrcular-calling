rule get_genome:
    output:
        "resources/genome.fasta",
    log:
        "logs/get-genome.log",
    params:
        species=config["calling"]["reference"]["species"],
        datatype="dna",
        build=config["calling"]["reference"]["build"],
        release=config["calling"]["reference"]["release"],
    cache: True
    wrapper:
        "v1.0.0/bio/reference/ensembl-sequence"


rule genome_faidx:
    input:
        f"resources/{config['calling']['reference']['name']}.fasta",
    output:
        f"resources/{config['calling']['reference']['name']}.fasta.fai",
    log:
        "logs/genome-faidx.log",
    cache: True
    wrapper:
        "v1.0.0/bio/samtools/faidx"
