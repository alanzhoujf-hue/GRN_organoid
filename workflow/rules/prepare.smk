rule prepare_loom:
    input:
        h5ad=config["input_h5ad"]

    output:
        loom=f"{RESULTS}/{{drug}}/input.loom"

    log:
        f"{LOGS}/prepare_loom/{{drug}}.log"

    benchmark:
        f"{BENCHMARKS}/prepare_loom/{{drug}}.tsv"

    threads:
        20

    resources:
        mem_mb=128000,
        runtime=480,
        io_slots=1

    conda:
        "../envs/preprocess.yaml"

    params:
        timepoint=config["timepoint_value"],
        drug_column=config["drug_column"],
        timepoint_column=config["timepoint_column"],
        celltype_column=config["celltype_column"],
        counts_layer=config["counts_layer"]

    shell:
        """
        python workflow/scripts/prepare_loom.py \
          --input {input.h5ad:q} \
          --output {output.loom:q} \
          --drug {wildcards.drug:q} \
          --timepoint {params.timepoint:q} \
          --drug-column {params.drug_column:q} \
          --timepoint-column {params.timepoint_column:q} \
          --celltype-column {params.celltype_column:q} \
          --counts-layer {params.counts_layer:q} \
          > {log:q} 2>&1
        """