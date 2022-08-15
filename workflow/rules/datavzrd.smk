rule render_datavzrd_config:
    input:
        template=workflow.source_path(
            "../resources/datavzrd/circle-table-template.datavzrd.yaml"
        ),
    output:
        "resources/datavzrd/all.datavzrd.yaml",
    params:
        groups=lambda wc: GROUPS,
        calls=lambda wc: [
            f"results/calling/tables/{group}.sorted.annotated.csv" for group in GROUPS
        ],
    log:
        "logs/datavzrd_render/all.log",
    template_engine:
        "yte"


rule copy_qc_plots_for_datavzrd:
    input:
        plots="results/calling/coverage_graphs/{group}",
        report="results/datavzrd-report/all.fdr-controlled",
    output:
        marker="results/tmp/{group}.qc_plots.marker",
    params:
        output_dir=lambda wc: directory(
            f"results/datavzrd-report/all.fdr-controlled/circles-{wc.group}/qc_plots"
        ),
    shell:
        """
        mkdir -p {params.output_dir}
        for f in `find {input.plots} -regex '.*/graph_[0-9]+_[0-9]+\.html'`; do ( cp "$f" {params.output_dir} ); done
        touch {output.marker}
        """


rule copy_graph_plots_for_datavzrd:
    input:
        plots="results/calling/graphs/rendered/{group}",
        report="results/datavzrd-report/all.fdr-controlled",
    output:
        marker="results/tmp/{group}.graph_plots.marker",
    params:
        output_dir=lambda wc: directory(
            f"results/datavzrd-report/all.fdr-controlled/circles-{wc.group}/graphs"
        ),
    shell:
        """
        mkdir -p {params.output_dir}
        for f in `find {input.plots} -regex '.*/graph_[0-9]+\.pdf'`; do ( cp "$f" {params.output_dir} ); done
        touch {output.marker}
        """


rule datavzrd_circle_calls:
    input:
        config="resources/datavzrd/all.datavzrd.yaml",
        calls=expand(
            "results/calling/tables/{group}.sorted.annotated.csv", group=GROUPS
        ),
    output:
        report(
            directory("results/datavzrd-report/all.fdr-controlled"),
            htmlindex="index.html",
            category="Circle calls",
        ),
    conda:
        "../envs/datavzrd.yaml"
    log:
        "logs/datavzrd_report/all.log",
    shell:
        "datavzrd {input.config} --output {output} &> {log}"
