rule get_genome:
    output:
        "resources/genome.fasta",
    log:
        "logs/get-genome.log",
    params:
        species=config["reference"]["species"],
        datatype="dna",
        build=config["reference"]["build"],
        release=config["reference"]["release"],
    cache: True
    wrapper:
        "v1.25.0/bio/reference/ensembl-sequence"


rule genome_faidx:
    input:
        f"resources/{config["reference"]["name"]}.fasta",
    output:
        f"resources/{config["reference"]["name"]}.fasta.fai",
    log:
        "logs/genome-faidx.log",
    cache: True
    wrapper:
        "v1.25.0/bio/samtools/faidx"
