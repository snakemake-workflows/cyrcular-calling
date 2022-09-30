from pathlib import Path


rule render_datavzrd_config:
    input:
        template=workflow.source_path(
            "../resources/datavzrd/circle-table-template.datavzrd.yaml"
        ),
        overview_tables=expand(
            "results/calling/tables/{group}/{group}_overview.{kind}.tsv",
            group=GROUPS,
            allow_missing=True,
        ),
        detail_tables=expand(
            "results/calling/tables/{group}/{group}_details/", group=GROUPS
        ),
    output:
        "results/datavzrd/{kind}.datavzrd.yaml",
    params:
        groups=lambda wc: GROUPS,
        overview_tables=lambda wc: [
            (group, f"results/calling/tables/{group}/{group}_overview.{wc.kind}.tsv")
            for group in GROUPS
        ],
        detail_tables=get_detail_tables_for_report,
    log:
        "logs/datavzrd_render/{kind}.log",
    template_engine:
        "yte"


rule copy_qc_plots_for_datavzrd:
    input:
        plots="results/calling/coverage_graphs/{group}",
        overview="results/calling/tables/{group}/{group}_overview.{kind}.tsv",
        report="results/datavzrd-report/{kind}.fdr-controlled",
    output:
        marker="results/tmp/{group}.{kind}.qc_plots.marker",
    params:
        output_dir=lambda wc: directory(
            f"results/datavzrd-report/{wc.kind}.fdr-controlled/circles-{wc.group}/qc_plots"
        ),
    log:
        "logs/datavzrd/copy_qc_plots/{group}.{kind}.log",
    run:
        import os
        import shutil
        from pathlib import Path

        os.makedirs(params.output_dir)
        overview = pd.read_csv(input.overview, sep="\t")
        event_ids = {s.replace("-", "_") for s in overview["event_id"]}
        for event_id in event_ids:
            shutil.copy(
                Path(input.plots).join(f"graph_{event_id}.html"),
                Path(params.output_dir),
            )
        Path(output.marker).touch()


rule copy_graph_plots_for_datavzrd:
    input:
        plots="results/calling/graphs/rendered/{group}",
        overview="results/calling/tables/{group}/{group}_overview.{kind}.tsv",
        report="results/datavzrd-report/{kind}.fdr-controlled",
    output:
        marker="results/tmp/{group}.{kind}.graph_plots.marker",
    params:
        output_dir=lambda wc: directory(
            f"results/datavzrd-report/{wc.kind}.fdr-controlled/circles-{wc.group}/graphs"
        ),
    log:
        "logs/datavzrd/copy_graph_plots/{group}.{kind}.log",
    run:
        import os
        import shutil
        from pathlib import Path

        os.makedirs(params.output_dir)
        overview = pd.read_csv(input.overview, sep="\t")
        graph_ids = set(overview["graph_id"])
        for graph_id in graph_ids:
            shutil.copy(
                Path(input.plots).join(f"graph_{graph_id}.pdf"), Path(params.output_dir)
            )
        Path(output.marker).touch()


rule datavzrd_circle_calls:
    input:
        config="results/datavzrd/{kind}.datavzrd.yaml",
        overview_tables=expand(
            "results/calling/tables/{group}/{group}_overview.{kind}.tsv",
            group=GROUPS,
            allow_missing=True,
        ),
        detail_tables=expand(
            "results/calling/tables/{group}/{group}_details/", group=GROUPS
        ),
    output:
        report(
            directory("results/datavzrd-report/{kind}.fdr-controlled"),
            htmlindex="index.html",
            category="Circle calls",
        ),
    conda:
        "../envs/datavzrd.yaml"
    log:
        "logs/datavzrd_report/{kind}.log",
    shell:
        "datavzrd {input.config} --output {output} &> {log}"
