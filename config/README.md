# general settings

To configure this workflow, modify the `config/config.yaml` files according to your needs. 
Comments in the example file provide explanations of the possible entries.

# samples sheet

Add samples to `config/samples.tsv` (or whichever file you specify under `samples:` in the `config/config.yaml` file).
For each sample, the columns `sample_name`, `group`, `alias` and `platform`:
* The `sample_name` needs to be a unique identifier of a biological sample.
* The `group` name groups samples together to be called jointly, for example all sample from a particular cell line or a particular patient.
* The `alias` defines the type of a sample within a group, for example a type of treatment or whether it is a Nanopore or an Illumina sample.
  Each group should contain the same set of aliases and the aliases should be used / reference in the [calling scenario](#calling-scenario).
* The `platform` column needs to contain the used sequencing plaform ('ILLUMINA','NANOPORE').
  Note that Illumina data is supplementary, i.e. circles will only be called for sample groups that have at at least 'nanopore' as their platform.

Missing values can be specified by empty columns or by writing `NA`. Lines can be commented out with `#`.

# units sheet

For each sample, add one or more sequencing units (runs or lanes or replicates) to the unit sheet `config/units.tsv` (or whichever file you specify under `units:` in the `config/config.yaml` file).
* Each unit is associared with one `sample_name` that needs to also be present in the [samples sheet](#samples-sheet).
* Each unit has a `unit` name, which can be the ID of a run, lane or replicate id.
  If you don't have any ID, you can just make one up.
  Unit IDs just need to be distinct for all the units of a particular sample.
* For each unit, define either one (column `fq1`) or two (columns `fq1`, `fq2`) FASTQ files (these can point to anywhere in your system). 
* Alternatively, you can define an SRA (sequence read archive) accession (starting with e.g. ERR or SRR) by using a column `sra`.
  In the latter case, the pipeline will automatically download the corresponding paired end reads from SRA.
  If both local files and SRA accession are available, the local files will be preferred.

Missing values can be specified by empty columns or by writing `NA`. Lines can be commented out with `#`.

# calling scenario

You need to specify a [varlociraptor](https://varlociraptor.github.io) [variant calling scenario](https://varlociraptor.github.io/docs/calling), to define which circle frequencies in which type of sample (`alias`) are of interest to you.
For details regarding the scenario format, see the [variant calling scenario documentation](https://varlociraptor.github.io/docs/calling).
To get an idea of what is possible with such scenarios, you can have a look at the [varlociraptor scenario catalog](https://varlociraptor.github.io/varlociraptor-scenarios/landing/).
