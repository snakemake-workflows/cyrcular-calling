from pathlib import Path


# TODO: rethink overview_tables and detail_tables
#       I am not really sure how to do this well, but I think we could have
#       just one circle-centered and one segment-centered table each, with
#       proper filtering (probably also in the linkouts from the overview plot_
#       and with current plot linkouts for graph and circle as collapsed /
#       hidden plots right there in the table. We'll need to discuss.
rule render_datavzrd_config:
    input:
        template=workflow.source_path(
            "../resources/datavzrd/circle-table-template.datavzrd.yaml"
        ),
        summary_spec=workflow.source_path("../resources/datavzrd/summary_plot.json"),
        categorized_overview_table="results/calling/tables/{group}/{group}_categorized_overview.tsv",
        overview_tables=expand(
            "results/calling/tables/{{group}}/{{group}}_overview.{category}.tsv",
            category=CATEGORIES,
        ),
        detail_tables="results/calling/tables/{group}/{group}_details/",
        circle_qc_plot_link_formatter=workflow.source_path(
            "../scripts/circle_qc_plot_link_formatter.js"
        ),
        graph_link_formatter=workflow.source_path("../scripts/graph_link_formatter.js"),
        gene_card_link_formatter=workflow.source_path(
            "../scripts/gene_card_link_formatter.js"
        ),
    output:
        "results/datavzrd/{group}.datavzrd.yaml",
    params:
        categories=CATEGORIES,
        overview_tables=lambda wc, input: [
            (file.split("_overview.")[1].replace(".tsv", ""), file)
            for file in list(input.overview_tables)
        ],
        detail_tables=lambda wc, input: get_detail_tables_group_circle_path_for_report(
            wc, input
        ),
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
            f"results/datavzrd-report/{wc.group}.fdr-controlled/qc_plots"
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
            f"results/datavzrd-report/{wc.group}.fdr-controlled/graphs"
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
