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


wildcard_constraints:
    sample="|".join(SAMPLES),
    group="|".join(GROUPS),


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
    targets += expand("results/datavzrd-report/{group}.fdr-controlled", group=GROUPS)
    targets += expand(
        "results/tmp/{group}.{category}.qc_plots.marker",
        group=GROUPS,
        category=CATEGORIES,
    )
    targets += expand(
        "results/tmp/{group}.{category}.graph_plots.marker",
        group=GROUPS,
        category=CATEGORIES,
    )
    return targets


def pairhmm_mode(wildcards):
    if samples.loc[wildcards.sample]["platform"] == "nanopore":
        mode = "homopolymer"
    else:
        mode = "exact"
    return mode


def get_group_candidates(wildcards):
    sample = wildcards.sample
    group = samples.loc[sample]["group"]
    wildcards.group = group
    scenario = scenario_name(wildcards)
    if scenario == "nanopore_only":
        sample = list(
            samples.query(f"group == '{group}' & platform == 'nanopore'")["sample_name"]
        )[0]
        return f"results/calling/candidate-calls/{sample}.{{scatteritem}}.bcf"
    elif scenario == "illumina_only":
        sample = list(
            samples.query(f"group == '{group}' & platform == 'illumina'")["sample_name"]
        )[0]
        return f"results/calling/candidate-calls/{sample}.{{scatteritem}}.bcf"
    elif scenario == "nanopore_with_illumina_support":
        sample = list(
            samples.query(f"group == '{group}' & platform == 'nanopore'")["sample_name"]
        )[0]
        return f"results/calling/candidate-calls/{sample}.{{scatteritem}}.bcf"
    else:
        raise ValueError(f"Unknown scenario: {scenario}")


def observation_string(wildcards, input):
    group_size = len(samples.query(f"group == '{wildcards.group}'"))
    scenario = scenario_name(wildcards)
    if group_size == 1:
        if scenario == "nanopore_only":
            return f"nanopore={input.obs[0]}"
        elif scenario == "illumina_only":
            return f"illumina={input.obs[0]}"
    elif group_size == 2:
        return f"nanopore={input.obs[0]} illumina={input.obs[1]}"
    else:
        raise ValueError(
            f"Too many samples ({group_size}) for this group ({wildcards.group}) and scenario ({scenario})"
        )


def get_observations(wildcards):
    s = samples.query(f"group == '{wildcards.group}'")
    if len(s) == 0:
        raise ValueError(f"No samples for group {wildcards.group}")

    observations = []

    has_nanopore = len(s.query("platform == 'nanopore'")["sample_name"]) > 0
    has_illumina = len(s.query("platform == 'illumina'")["sample_name"]) > 0

    if has_nanopore:
        for sample_nanopore in list(s.query("platform == 'nanopore'")["sample_name"]):
            observations.append(
                f"results/calling/calls/observations/{sample_nanopore}.{{scatteritem}}.bcf"
            )

    if has_illumina:
        for sample_illumina in list(s.query("platform == 'illumina'")["sample_name"]):
            observations.append(
                f"results/calling/calls/observations/{sample_illumina}.{{scatteritem}}.bcf"
            )

    return observations


def scenario_name(wildcards):
    s = samples.query(f"group == '{wildcards.group}'")
    num_samples_in_group = len(s)
    if num_samples_in_group == 1:
        if "illumina" in set(s["platform"]):
            return "illumina_only"
        elif "nanopore" in set(s["platform"]):
            return "nanopore_only"
        else:
            platforms = ", ".join(set(s["platform"]))
            raise ValueError(
                f"Single sample scenario not defined for platforms {platforms}"
            )
    elif num_samples_in_group == 2:
        if len(set(s["platform"]) - {"illumina", "nanopore"}) == 0:
            return "nanopore_with_illumina_support"
        else:
            raise ValueError(
                "Need both illumina and nanopore samples for this scenario"
            )
    else:
        raise ValueError(
            "Scenarios with more than two samples per group currently not supported"
        )


def get_scenario(wildcards):
    scenario = scenario_name(wildcards)
    if scenario == "nanopore_only":
        return workflow.source_path(
            "../resources/scenarios/nanopore_circle_scenario.yaml"
        )
    elif scenario == "illumina_only":
        return workflow.source_path(
            "../resources/scenarios/illumina_circle_scenario.yaml"
        )
    elif scenario == "nanopore_with_illumina_support":
        return workflow.source_path(
            "../resources/scenarios/nanopore_illumina_joint_circle_scenario.yaml"
        )
    else:
        raise ValueError(f"Unknown scenario: {scenario}")


def get_minimap2_mapping_params(wildcards):
    if samples.loc[wildcards.sample]["platform"] == "nanopore":
        return "-x map-ont"
    elif samples.loc[wildcards.sample]["platform"] == "illumina":
        return "-x sr"
    else:
        return ""


def get_minimap2_input(wildcards):
    if units.loc[wildcards.sample]["fq2"].any():
        return [
            "results/calling/merged/{sample}_R1.fastq.gz",
            "results/calling/merged/{sample}_R2.fastq.gz",
        ]
    else:
        return ["results/calling/merged/{sample}_single.fastq.gz"]


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
