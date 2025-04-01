rule run_clustering:
    input:
        data="data/Pollen2014.txt",
        labels="data/SupplementaryLabels.txt"
    output:
        plot="results/arandi_plot.png",
        rds="results/arandi_output.rds"
    shell:
        """
        Rscript scripts/Main.R {input.data} {input.labels} {output.plot} {output.rds}

        """
