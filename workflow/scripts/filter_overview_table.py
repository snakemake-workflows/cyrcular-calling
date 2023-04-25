from contextlib import redirect_stderr

with open(snakemake.log[0], "w") as logfile:
    with redirect_stderr(logfile):
        import pandas as pd

        df = pd.read_csv(snakemake.input.table, sep="\t")
        df.loc[:, ["gene_names", "gene_ids", "regulatory_features"]] = df.loc[
            :, ["gene_names", "gene_ids", "regulatory_features"]
        ].fillna("")
        df["category"] = df[["num_exons", "regulatory_features", "gene_names"]].apply(
            lambda r: "coding"
            if r["num_exons"] > 0
            else (
                "regulatory"
                if r["regulatory_features"]
                else ("intronic" if r["gene_names"] else "other")
            ),
            axis=1,
        )
        for kind in ["coding", "regulatory", "intronic", "other"]:
            part = df.query(f"category == '{kind}'")
            part.to_csv(getattr(snakemake.output, kind), sep="\t", index=False)
        df.to_csv(snakemake.output.categorized, sep="\t", index=False)
