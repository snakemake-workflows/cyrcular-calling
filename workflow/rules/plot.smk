import os


rule circle_coverage_plot:
    input:
        bam="results/calling/mapping/{sample}.bam",
        bai="results/calling/mapping/{sample}.bam.bai",
        graph="results/calling/graphs/{sample}.graph",
    output:
        plots=directory("results/calling/coverage_graphs/{sample}/"),
    benchmark:
        "benchmarks/cyrcular_plots/{sample}.txt"
    log:
        "logs/cyrcular_plots/{sample}.log",
    conda:
        "../envs/cyrcular.yaml"
    threads: 1
    shell:
        """cyrcular\
        graph\
        plot\
        {input.bam}\
        --graph {input.graph}\
        --output {output.plots}\
        2> {log}"""


rule circle_graph_plot:
    input:
        graph="results/calling/graphs/{sample}",
    output:
        pdf_dir=directory("results/calling/graphs/pdfs/{sample}"),
        pdf=report(
            "results/calling/graphs/{sample}.pdf",
            category="Graphs",
            subcategory="{sample}",
        ),
    log:
        "logs/graphviz/{sample}.log",
    conda:
        "../envs/graphviz.yaml"
    threads: 1
    shell:
        """
        mkdir -p {output.pdf_dir}
        count=`ls -1 {input.graph}/ 2>{log} | wc -l`
        if [ $count != 0 ]; then
            for f in {input.graph}/*.dot; do (dot $f -Tpdf > "{output.pdf_dir}/${{f##*/}}.pdf" 2>>{log}); done
            convert {output.pdf_dir}/*.pdf {output.pdf} 2>> {log}
        else
            convert xc:none -page 2x2 {output.pdf} 2>> {log}
        fi
        """
