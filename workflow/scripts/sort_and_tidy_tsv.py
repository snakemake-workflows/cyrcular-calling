import pandas as pd
from math import ceil
import re


BND_RE = re.compile(r""".*([]\[])((?P<seqname>.+):(?P<position>[0-9]+))([]\[])?.*""")


def parse_bnd_alt(s: str):
    return BND_RE.search(s)["seqname"], int(BND_RE.search(s)["position"])


def sort_table(df: pd.DataFrame):
    ordered = (
        df.groupby(["EVENT"])[["PROB_PRESENT", "length", "split_reads"]]
        .max()
        .sort_values(["PROB_PRESENT", "length", "split_reads"], ascending=False)
        .index
    )
    categories = list(ordered)
    df["EVENT"] = pd.Categorical(df["EVENT"], categories=categories)
    df.sort_values("EVENT", inplace=True)
    return df


def main(snakemake):
    df = pd.read_csv(snakemake.input.tsv, sep="\t")
    columns = [
        "graph",
        "EVENT",
        "CHROM",
        "start",
        "OTHER_CHROM",
        "stop",
        "ALT",
        "direction",
        "circle_length",
        "num_segments",
        "split_reads",
        "AF_nanopore",
        "PROB_PRESENT",
        "PROB_ABSENT",
        "PROB_ARTIFACT",
        "length",
        "NUM_EXONS",
        "GENES",
    ]
    # if there are no entries, make sure an empty table with the correct header is written anyways
    if df.empty:
        df = pd.DataFrame(columns=columns)
        df.to_csv(snakemake.output.csv, index=False)
        return

    # the list of genes may be too long for excel spreadsheets, hence we split the column into multiple ones
    # from https://stackoverflow.com/questions/58645152/splitting-a-long-string-in-pandas-cell-near-the-n-th-character-position-into-mul
    limit = 32767 - 1
    sep = ";"
    longest = df["GENES"].apply(lambda s: len(s) if isinstance(s, str) else 0).max()
    num_seps = (
        df["GENES"].apply(lambda s: s.count(sep) if isinstance(s, str) else 0).max()
    )
    num_cols = max(1, ceil(longest / (limit - num_seps)))

    for index, row in df.iterrows():
        gene_names = row["GENES"]
        num_chars = len(gene_names) if isinstance(gene_names, str) else 0
        if num_chars == 0:
            continue
        genes = gene_names.split(sep)
        num_genes = len(genes)
        col_index = ceil(len(genes) / num_cols)

        for _ in range(num_cols - 1):
            genes.insert(col_index, "|")
            col_index += col_index
        new = sep.join(genes)
        df.at[index, "GENES"] = new

    df[["OTHER_CHROM", "stop"]] = pd.DataFrame(
        df["ALT"].apply(parse_bnd_alt).tolist(), index=df.index
    )
    df["direction"] = df["ALT"].apply(
        lambda a: "to" if "[" in a else ("from" if "]" in a else "*")
    )
    df.rename(columns=dict(POS="start"), inplace=True)

    df["length"] = df["circle_length"]
    df.query("CHROM != OTHER_CHROM")["length"] = float("nan")
    df.drop(columns=["ID"], inplace=True)

    new_cols = [f"GENES_{i + 1}" for i in range(num_cols)]
    df[new_cols] = df["GENES"].astype(dtype=str).str.split("|", expand=True)
    df[new_cols] = df[new_cols].apply(lambda s: s.str.strip(sep))
    del df["GENES"]
    df.rename(columns=dict(GENES_1="GENES"), inplace=True)
    gene_cols = [c for c in df.columns if c.startswith("GENES")]

    df = sort_table(df)

    df.drop_duplicates(
        subset=["EVENT", "CHROM", "start", "OTHER_CHROM", "stop"], inplace=True
    )

    df["graph"] = [re.sub(r"_circle_.*", "", d) for d in df["EVENT"]]

    # reorder columns
    df = df[columns[:-1] + gene_cols]

    df.to_csv(snakemake.output.csv, index=False, sep=",")


if __name__ == "__main__":
    main(snakemake)
