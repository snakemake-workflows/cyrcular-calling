from contextlib import redirect_stderr

with open(snakemake.log[0], "w") as logfile:
    with redirect_stderr(logfile):
        import os
        import shutil
        import pandas as pd
        from pathlib import Path

        os.makedirs(snakemake.params.output_dir, exist_ok=True)
        overview = pd.read_csv(snakemake.input.overview, sep="\t")
        event_ids = {s.replace("-", "_") for s in overview["event_id"]}
        for event_id in event_ids:
            shutil.copy(
                f"{snakemake.input.plots}/graph_{event_id}.html",
                snakemake.params.output_dir,
            )
        Path(snakemake.output.marker).touch()
