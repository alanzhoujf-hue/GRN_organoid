rule prepare_pooled_loom:
    input:
        expand(f"{RESULTS}/{{drug}}/input.loom", drug=DRUGS) #a Path list
    output:
        f"{POOLED}/pooled_allcells.loom" #a Path list need [0]
    conda:
        "../envs/preprocess.yaml"
    threads: 10
    resources:
        mem_mb=64000,
        runtime=480,
        io_slots=1
    script:
        "../scripts/prepare_pooled_loom.py"


rule pooled_grnboost2:
    input:
        loom=f"{POOLED}/pooled_allcells.loom",
        tf_list=config["resources"]["tf_list"],
        image=config["container"]["pyscenic"]
    output:
        f"{POOLED}/adjacencies.tsv"
    threads: 40
    resources:
        mem_mb=256000,
        runtime=2880
    params:
        bind=f"{config['project_dir']}:{config['project_dir']}",
        workers=8
    shell:
        """
        apptainer exec --cleanenv -B {params.bind:q} {input.image:q} \
          pyscenic grn {input.loom:q} {input.tf_list:q} \
          --method grnboost2 --num_workers {params.workers} \
          --output {output:q}
        """


rule pooled_cistarget:
    input:
        adjacency=f"{POOLED}/adjacencies.tsv",
        loom=f"{POOLED}/pooled_allcells.loom",
        databases=config["resources"]["ranking_databases"],
        annotations=config["resources"]["motif_annotations"],
        image=config["container"]["pyscenic"]
    output:
        f"{POOLED}/pooled_regulons.csv"
    threads: 20
    resources:
        mem_mb=128000,
        runtime=1440
    params:
        bind=f"{config['project_dir']}:{config['project_dir']}",
        workers=8
    shell:
        """
        apptainer exec --cleanenv -B {params.bind:q} {input.image:q} \
          pyscenic ctx {input.adjacency:q} {input.databases:q} \
          --annotations_fname {input.annotations:q} \
          --expression_mtx_fname {input.loom:q} \
          --mode custom_multiprocessing \
          --num_workers {params.workers} \
          --output {output:q}
        """


rule aucell_pooled_ref:
    input:
        loom=f"{RESULTS}/{{drug}}/input.loom",
        regulons=f"{POOLED}/pooled_regulons.csv",
        image=config["container"]["pyscenic"]
    output:
        f"{RESULTS}/{{drug}}/auc_mtx_pooled_allcells_ref.csv"
    threads: 10
    resources:
        mem_mb=64000,
        runtime=720
    params:
        bind=f"{config['project_dir']}:{config['project_dir']}",
        workers=8
    shell:
        """
        apptainer exec --cleanenv -B {params.bind:q} {input.image:q} \
          pyscenic aucell {input.loom:q} {input.regulons:q} \
          --num_workers {params.workers} --output {output:q}
        """


rule mean_auc_celltype_pooled_ref:
    input:
        auc=f"{RESULTS}/{{drug}}/auc_mtx_pooled_allcells_ref.csv",
        loom=f"{RESULTS}/{{drug}}/input.loom"
    output:
        f"{RESULTS}/{{drug}}/mean_auc_celltype_pooled_allcells_ref.csv"
    conda:
        "../envs/preprocess.yaml"
    script:
        "../scripts/mean_auc_celltype.py"
