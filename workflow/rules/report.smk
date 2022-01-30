
rule csv_report:
    input:
        csv="results/calling/tables/{group}.sorted.annotated.csv",
        plots="results/calling/coverage_graphs/{group}",
    output:
        report(
            directory("results/calling/report/tables/{group}"),
            category="Circle calls",
            htmlindex="index.html",
        ),
    params:
        extra="--formatter resources/report-table-formatter.js",
    log:
        "logs/rbt-csv-report/{group}.log",
    benchmark:
        "benchmarks/rbt-csv-report/{group}.txt"
    conda:
        "../envs/rbt.yaml"
    shell:
        """
        mkdir -p {output}/qc_plots
        for f in {input.plots}/*.html; do ( cp "$f" {output}/qc_plots/ ); done
        rbt csv-report {input.csv} {output} {params.extra} &> {log}
        """
