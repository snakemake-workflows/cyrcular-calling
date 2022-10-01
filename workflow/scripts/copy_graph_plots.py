from contextlib import redirect_stderr

with open(snakemake.log[0], "w") as logfile:
    with redirect_stderr(logfile):
        logfile.write("A")
        import os
        import shutil
        import pandas as pd
        from pathlib import Path
        logfile.write("B")

        os.makedirs(snakemake.params.output_dir, exist_ok=True)
        logfile.write("C")
        overview = pd.read_csv(snakemake.input.overview, sep="\t")
        graph_ids = set(overview["graph_id"])
        logfile.write("D")
        for graph_id in graph_ids:
            shutil.copy(
                f"{snakemake.input.plots}/graph_{graph_id}.pdf", snakemake.params.output_dir
            )
        logfile.write("E")
        Path(snakemake.output.marker).touch()
        logfile.write("F")