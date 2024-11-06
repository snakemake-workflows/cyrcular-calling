import pandas as pd
import re
from snakemake.utils import validate

validate(config, schema="../schemas/config.schema.yaml")


def read_samples():
    samples = (
        pd.read_csv(
            config["samples"],
            sep="\t",
            dtype={"sample_name": str, "group": str},
            comment="#",
        )
        .set_index("sample_name", drop=False)
        .sort_index()
    )

    def _group_or_sample(row):
        group = row.get("group", None)
        if pd.isnull(group):
            return row["sample_name"]
        return group

    samples["group"] = [_group_or_sample(row) for _, row in samples.iterrows()]
    validate(samples, schema="../schemas/samples.schema.yaml")
    return samples


samples = read_samples()
SAMPLES = list(sorted(set(samples["sample_name"])))
GROUPS = list(sorted(set(samples["group"])))
CATEGORIES = ["coding", "regulatory", "intronic", "other"]
EVENTS=lookup(dpath="filter/fdr-control/events", within=config)

wildcard_constraints:
    sample="|".join(SAMPLES),
    group="|".join(GROUPS),
    event="|".join(EVENTS),


def read_units():
    units = (
        pd.read_csv(
            config["units"],
            sep="\t",
            dtype={"sample_name": str, "unit_name": str},
            comment="#",
        )
        .set_index(["sample_name", "unit_name"], drop=False)
        .sort_index()
    )
    validate(units, schema="../schemas/units.schema.yaml")
    return units


units = read_units()


def get_all_input(wildcards):
    targets = []
    targets += expand(
        "results/coverage/{sample}.mosdepth.summary.txt", 
        sample=SAMPLES,
    )
    targets += expand(
        "results/datavzrd-report/{group}.{event}.fdr-controlled",
        group=GROUPS,
        event=EVENTS,
    )
    targets += expand(
        "results/tmp/{group}.{event}.{category}.graph_plots.marker",
        group=GROUPS,
        event=EVENTS,
        category=CATEGORIES,
    )
    for group in GROUPS:
        targets += expand(
            "results/tmp/{group}.{event}.{sample}.{category}.qc_plots.marker",
            group=group,
            event=EVENTS,
            sample=samples.loc[samples["group"] == group, "sample_name"],
            category=CATEGORIES,
        )
    return targets


def pairhmm_mode(wildcards):
    if samples.loc[wildcards.sample]["platform"].lower() == "nanopore":
        mode = "homopolymer"
    else:
        mode = "exact"
    return mode


def get_varlociraptor_obs_args(wildcards, input):
    aliases = samples.loc[samples["group"] == wildcards.group, "alias"]
    return [
        "{}={}".format(s, f)
        for s, f in zip(aliases, input.obs)
    ]


def get_minimap2_mapping_params(wildcards):
    if samples.loc[wildcards.sample]["platform"].lower() == "nanopore":
        return "-x map-ont"
    elif samples.loc[wildcards.sample]["platform"].lower() == "illumina":
        return "-x sr"
    else:
        return ""


def get_minimap2_input(wildcards):
    if units.loc[wildcards.sample]["fq2"].any():
        return [
            "results/merged_fastqs/{sample}_R1.fastq.gz",
            "results/merged_fastqs/{sample}_R2.fastq.gz",
        ]
    else:
        return ["results/merged_fastqs/{sample}_single.fastq.gz"]


def get_fastqs(wildcards):
    if wildcards.read == "single":
        fq1 = units.loc[wildcards.sample]["fq1"]
        return list(fq1)
    elif wildcards.read == "R1":
        fq1 = units.loc[wildcards.sample]["fq1"]
        return list(fq1)
    elif wildcards.read == "R2":
        fq2 = units.loc[wildcards.sample]["fq2"]
        return list(fq2)


# black wasn't happy about the inline version of this in the params section
def varlociraptor_filtering_mode(wildcards):
    return config["filter"]["fdr-control"].get("mode", "local-smart")


class Region:
    def __init__(self, target: str, start: int, end: int):
        self.target = target
        self.start = start
        self.end = end

    @classmethod
    def from_str(cls, s: str) -> "Region":
        target, pos = s.split(":")
        start, end = pos.split("-")
        start, end = int(start), int(end)
        assert end >= start
        return cls(target, start, end)

    def __str__(self) -> str:
        return f"{self.target}:{self.start}-{self.end}"

    def __repr__(self) -> str:
        return str(self)

    def __len__(self) -> int:
        return self.end - self.start

    def overlap(self, other: "Region") -> int:
        return min(self.end, other.end) - max(self.start, other.start)

    def merge(self, other: "Region", min_overlap: int = 0) -> "Region":
        if self.target != other.target:
            raise ValueError("targets do not match")
        if self.overlap(other) < min_overlap:
            raise ValueError("ranges do not overlap")
        return Region(
            self.target, min(self.start, other.start), max(self.end, other.end)
        )


BND_RE = re.compile(r""".*([]\[])((?P<seqname>.+):(?P<position>[0-9]+))([]\[])?.*""")


def parse_bnd_alt(s: str):
    return BND_RE.search(s)["seqname"], int(BND_RE.search(s)["position"])


## rules/annotation.smk specific
CYRCULAR_INFO_FIELDS = ["CircleLength", "CircleSegmentCount", "SplitReads", "Support"]


def get_detail_tables_group_circle_path_for_report(wildcards, input):
    from pathlib import Path

    group = wildcards.group
    res = []
    folder = input.detail_tables
    group_tsv = pd.read_csv(input.categorized_overview_table, sep="\t")
    keep_event_ids = set(group_tsv["event_id"])
    detail_table_files = [f for f in os.listdir(folder) if f.endswith(".tsv")]
    event_ids = [
        f"{m.group(1)}-{m.group(2)}"
        for m in [
            re.match("graph_([^_]+)_circle_([^_]+)_segments.tsv", path)
            for path in detail_table_files
        ]
    ]
    for event_id, detail_file in zip(event_ids, detail_table_files):
        if event_id not in keep_event_ids:
            continue
        f = pd.read_csv(f"{folder}/{detail_file}", sep="\t")
        if f.empty:
            continue
        res.append((group, event_id, f"{folder}/{detail_file}"))
    return res
