rule grnboost2:
    input:
        loom=f"{RESULTS}/{{drug}}/input.loom",
        tf_list=config["resources"]["tf_list"],
        image=config["container"]["pyscenic"]

    output:
        adjacency=f"{RESULTS}/{{drug}}/adjacencies.tsv"

    log:
        f"{LOGS}/grnboost2/{{drug}}.log"

    benchmark:
        f"{BENCHMARKS}/grnboost2/{{drug}}.tsv"

    threads: 28

    resources:
        mem_mb=192000,
        runtime=1440

    params:
        bind=f"{config['project_dir']}:{config['project_dir']}",
        method="grnboost2",
        workers=8

    shell:
        """
        apptainer exec \
          --cleanenv \
          -B {params.bind:q} \
          {input.image:q} \
          pyscenic grn \
            {input.loom:q} \
            {input.tf_list:q} \
            --method {params.method:q} \
            --num_workers {params.workers} \
            --output {output.adjacency:q} \
          > {log:q} 2>&1
        """


rule cistarget:
    input:
        adjacency=f"{RESULTS}/{{drug}}/adjacencies.tsv",
        loom=f"{RESULTS}/{{drug}}/input.loom",
        databases=config["resources"]["ranking_databases"],
        annotations=config["resources"]["motif_annotations"],
        image=config["container"]["pyscenic"]

    output:
        regulons=f"{RESULTS}/{{drug}}/regulons.csv"

    log:
        f"{LOGS}/cistarget/{{drug}}.log"

    benchmark:
        f"{BENCHMARKS}/cistarget/{{drug}}.tsv"

    threads:
        config["pyscenic"]["num_workers"]

    resources:
        mem_mb=128000,
        runtime=1440
    
    params:
        bind=f"{config['project_dir']}:{config['project_dir']}"

    shell:
        """
        apptainer exec \
          --cleanenv \
          -B {params.bind:q} \
          {input.image:q} \
          pyscenic ctx \
            {input.adjacency:q} \
            {input.databases:q} \
            --annotations_fname {input.annotations:q} \
            --expression_mtx_fname {input.loom:q} \
            --mode custom_multiprocessing \
            --num_workers {threads} \
            --output {output.regulons:q} \
          > {log:q} 2>&1
        """

rule aucell:
    input:
        loom=f"{RESULTS}/{{drug}}/input.loom",
        regulons=f"{RESULTS}/{{drug}}/regulons.csv",
        image=config["container"]["pyscenic"]

    output:
        auc=f"{RESULTS}/{{drug}}/auc_mtx.csv"

    log:
        f"{LOGS}/aucell/{{drug}}.log"

    benchmark:
        f"{BENCHMARKS}/aucell/{{drug}}.tsv"

    threads:
        config["pyscenic"]["num_workers"]

    resources:
        mem_mb=64000,
        runtime=720

    params:
        bind=f"{config['project_dir']}:{config['project_dir']}",
        auc_threshold=config["pyscenic"]["auc_threshold"],
        seed=config["random_seed"]

    shell:
        """
        apptainer exec \
          --cleanenv \
          -B {params.bind:q} \
          {input.image:q} \
          pyscenic aucell \
            {input.loom:q} \
            {input.regulons:q} \
            --num_workers {threads} \
            --auc_threshold {params.auc_threshold} \
            --seed {params.seed} \
            --output {output.auc:q} \
          > {log:q} 2>&1
        """

rule aucell_dmso_ref:
    input:
        loom=f"{RESULTS}/{{drug}}/input.loom",
        regulons=f"{RESULTS}/dmso/regulons.csv",
        image=config["container"]["pyscenic"]

    output:
        auc=f"{RESULTS}/{{drug}}/auc_mtx_dmso_ref.csv"

    log:
        f"{LOGS}/aucell/{{drug}}_dmso_ref.log"

    benchmark:
        f"{BENCHMARKS}/aucell/{{drug}}_dmso_ref.tsv"

    threads:
        config["pyscenic"]["num_workers"]

    resources:
        mem_mb=64000,
        runtime=720

    params:
        bind=f"{config['project_dir']}:{config['project_dir']}",
        auc_threshold=config["pyscenic"]["auc_threshold"],
        seed=config["random_seed"]

    shell:
        """
        apptainer exec \
          --cleanenv \
          -B {params.bind:q} \
          {input.image:q} \
          pyscenic aucell \
            {input.loom:q} \
            {input.regulons:q} \
            --num_workers {threads} \
            --auc_threshold {params.auc_threshold} \
            --seed {params.seed} \
            --output {output.auc:q} \
          > {log:q} 2>&1
        """