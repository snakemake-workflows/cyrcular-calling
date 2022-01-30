import pandas as pd
from pandarallel import pandarallel
from dinopy import FastaReader

pandarallel.initialize(nb_workers=snakemake.threads, progress_bar=True)


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
    df = pd.read_csv(snakemake.input.csv)
    if df.empty:
        df["breakpoint_seq"] = None
        df.to_csv(snakemake.output.csv, index=False, sep=",")
        exit(0)

    far = FastaReader(snakemake.input.reference, write_fai=True)

    # TODO: needs dir information as well.
    def get_seq(chrom1, start, chrom2, stop, direction, up_down=1000):
        if direction == "from":
            chrom_to, pos_to = chrom1, start - 1
            chrom_from, pos_from = chrom2, stop - 1
        elif direction == "to":
            chrom_to, pos_to = chrom2, stop - 1
            chrom_from, pos_from = chrom1, start - 1
        else:
            print("foo")
            raise ValueError("unknown direction")

        entry_from = next(far[chrom_from])
        entry_to = next(far[chrom_to])
        chrom_from_len = entry_from.length
        chrom_to_len = entry_to.length

        from_start = max(0, pos_from - up_down)
        from_stop = pos_from

        to_start = pos_to
        to_stop = min(pos_to + up_down, chrom_to_len - 1)

        seq_from = entry_from.sequence[from_start:from_stop].decode()
        seq_to = entry_to.sequence[to_start:to_stop].decode()
        return seq_from + seq_to

    df["breakpoint_seq"] = df.parallel_apply(
        lambda row: get_seq(
            row["CHROM"],
            row["start"],
            row["OTHER_CHROM"],
            row["stop"],
            row["direction"],
        ),
        axis=1,
    )
    df = sort_table(df)
    df.to_csv(snakemake.output.csv, index=False, sep=",")


if __name__ == "__main__":
    main(snakemake)
