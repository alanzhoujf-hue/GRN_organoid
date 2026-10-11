rule rss_celltype:
    input:
        auc=f"{RESULTS}/{{drug}}/auc_mtx.csv",
        loom=f"{RESULTS}/{{drug}}/input.loom",
        image=config["container"]["pyscenic"]
    
    output:
        rss=f"{RESULTS}/{{drug}}/rss_cellType.csv"
    
    log:
        f"{LOGS}/rss_celltype/{{drug}}.log"

    benchmark:
        f"{BENCHMARKS}/rss_celltype/{{drug}}.tsv"

    threads:
        4

    resources:
        mem_mb=16000,
        runtime=240
    
    params:
        bind=f"{config['project_dir']}:{config['project_dir']}"
    
    
    shell:
        """
        apptainer exec \
          --cleanenv \
          -B {params.bind:q} \
          {input.image:q} \
          python workflow/scripts/calculate_rss.py \
            --auc {input.auc:q} \
            --loom {input.loom:q} \
            --output {output.rss:q} \
          > {log:q} 2>&1
        """


rule rss_celltype_dmso_ref:
    input:
        auc=f"{RESULTS}/{{drug}}/auc_mtx_dmso_ref.csv",
        loom=f"{RESULTS}/{{drug}}/input.loom",
        image=config["container"]["pyscenic"]
    
    output:
        rss=f"{RESULTS}/{{drug}}/rss_celltype_dmso_ref.csv"
    
    log:
        f"{LOGS}/rss_celltype/{{drug}}_dmso_ref.log"

    benchmark:
        f"{BENCHMARKS}/rss_celltype/{{drug}}_dmso_ref.tsv"

    threads:
        4

    resources:
        mem_mb=16000,
        runtime=240
    
    params:
        bind=f"{config['project_dir']}:{config['project_dir']}"
    
    
    shell:
        """
        apptainer exec \
          --cleanenv \
          -B {params.bind:q} \
          {input.image:q} \
          python workflow/scripts/calculate_rss.py \
            --auc {input.auc:q} \
            --loom {input.loom:q} \
            --output {output.rss:q} \
          > {log:q} 2>&1
        """


