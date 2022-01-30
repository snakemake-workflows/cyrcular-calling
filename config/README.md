
# General settings
To configure this workflow, modify ``config/config.yaml`` according to your needs, following the explanations provided in the file.

# Sample sheet

Add samples to `config/samples.tsv`. For each sample, the columns `sample_name`, `platform`, and `group` have to be defined. 
* Samples within the same `group` will be called jointly. 
* The `platform` column needs to contain the used sequencing plaform ('illumina','nanopore'). Note that Illumina data is supplementary, i.e. circles will only be called for sample groups that have at at least 'nanopore' as their platform.

Missing values can be specified by empty columns or by writing `NA`. Lines can be commented out with `#`.

# Unit sheet

For each sample, add one or more sequencing units (runs, lanes or replicates) to the unit sheet `config/units.tsv`.
* Each unit has a `unit` name, which can be e.g. a running number, or an actual run, lane or replicate id.
* Each unit has a `sample` name, which associates it with the biological sample it comes from.
* For each unit, define either one (column `fq1`) or two (columns `fq1`, `fq2`) FASTQ files (these can point to anywhere in your system). 
* Alternatively, you can define an SRA (sequence read archive) accession (starting with e.g. ERR or SRR) by using a column `sra`. In the latter case, the pipeline will automatically download the corresponding paired end reads from SRA. If both local files and SRA accession are available, the local files will be preferred.

Missing values can be specified by empty columns or by writing `NA`. Lines can be commented out with `#`.
