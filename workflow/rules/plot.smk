import os


rule circle_graphs_coverage_plot:
    input:
        bam="results/mapped/{sample}.bam",
        bai="results/mapped/{sample}.bam.bai",
        graph="results/circle_graphs/{group}.graph",
    output:
        plots=directory("results/circle_graphs/{group}/coverage/{sample}/"),
    benchmark:
        "benchmarks/cyrcular_plots/{group}/{sample}.txt"
    log:
        "logs/cyrcular_plots/{group}/{sample}.log",
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


rule circle_graphs_rendering:
    input:
        graph="results/circle_graphs/{group}",
    output:
        pdf_dir=directory("results/circle_graphs/{group}/rendered/"),
    log:
        "logs/circle_graphs_rendering/{group}.log",
    conda:
        "../envs/graphviz.yaml"
    threads: 1
    shell:
        """
        mkdir -p {output.pdf_dir}
        count=`ls -1 {input.graph}/ 2>{log} | wc -l`
        for f in {input.graph}/*.dot; do (dot $f -Tpdf > "{output.pdf_dir}/graph_$(basename ${{f}} .dot).pdf" 2>>{log}); done
        """
