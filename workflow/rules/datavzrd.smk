from pathlib import Path


rule render_datavzrd_config:
    input:
        template=workflow.source_path(
            "../resources/datavzrd/circle-table-template.datavzrd.yaml"
        ),
        categorized_overview_table="results/calling/tables/{group}/{group}_categorized_overview.tsv",
        overview_tables=expand(
            "results/calling/tables/{{group}}/{{group}}_overview.{category}.tsv",
            category=CATEGORIES,
        ),
        detail_tables="results/calling/tables/{group}/{group}_details/",
    output:
        "results/datavzrd/{group}.datavzrd.yaml",
    params:
        categories=CATEGORIES,
        overview_tables=lambda wc, input: [
                (file.split("_overview.")[1].replace(".tsv",""), file) for file in list(input.overview_tables)
            ],
        group=lambda wc: wc.group,
        detail_tables=get_detail_tables_for_report,
        summary_spec=workflow.source_path("../resources/datavzrd/summary_plot.json"),
        species=config["reference"]["species"],
    log:
        "logs/datavzrd_render/{group}.log",
    template_engine:
        "yte"


rule copy_qc_plots_for_datavzrd:
    input:
        plots="results/calling/coverage_graphs/{group}",
        overview="results/calling/tables/{group}/{group}_overview.{category}.tsv",
        report="results/datavzrd-report/{group}.fdr-controlled",
    output:
        marker="results/tmp/{group}.{category}.qc_plots.marker",
    params:
        output_dir=lambda wc: directory(
            f"results/datavzrd-report/{wc.group}.fdr-controlled/circles-{wc.category}/qc_plots"
        ),
    log:
        "logs/datavzrd/copy_qc_plots/{group}.{category}.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/copy_qc_plots.py"


rule copy_graph_plots_for_datavzrd:
    input:
        plots="results/calling/graphs/rendered/{group}",
        overview="results/calling/tables/{group}/{group}_overview.{category}.tsv",
        report="results/datavzrd-report/{group}.fdr-controlled",
    output:
        marker="results/tmp/{group}.{category}.graph_plots.marker",
    params:
        output_dir=lambda wc: directory(
            f"results/datavzrd-report/{wc.group}.fdr-controlled/circles-{wc.category}/graphs"
        ),
    log:
        "logs/datavzrd/copy_graph_plots/{group}.{category}.log",
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/copy_graph_plots.py"


rule datavzrd_circle_calls:
    input:
        config="results/datavzrd/{group}.datavzrd.yaml",
        overview_tables=expand(
            "results/calling/tables/{group}/{group}_overview.{category}.tsv",
            category=CATEGORIES,
            allow_missing=True,
        ),
        detail_tables="results/calling/tables/{group}/{group}_details/",
    output:
        report(
            directory("results/datavzrd-report/{group}.fdr-controlled"),
            htmlindex="index.html",
            category="Circle calls",
        ),
    conda:
        "../envs/datavzrd.yaml"
    log:
        "logs/datavzrd_report/{group}.log",
    shell:
        "datavzrd {input.config} --output {output} &> {log}"
