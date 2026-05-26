pixi run python src/manuscript_blood/comparative_analysis/jaccard_marker_genes.py \
    --adata_raw  /data_nfs/og86asub/netmap/netmap-evaluation/data/blood/reprocessed/bd-rhap-rep1/bd-rhap-rep1.h5ad \
    --results_dir /data_nfs/og86asub/netmap/netmap-evaluation/netmap/case_studies/blood \
    --analysis_dir /data_nfs/og86asub/netmap/netmap-evaluation/results/case_studies/blood/ \
    --experiment_name bd-rhap-rep1-X \
    --cluster_column 'leiden_remap' 

pixi run python src/manuscript_blood/comparative_analysis/jaccard_marker_genes.py \
    --adata_raw  /data_nfs/og86asub/netmap/netmap-evaluation/data/blood/reprocessed/bd-rhap-rep2/bd-rhap-rep2.h5ad \
    --results_dir /data_nfs/og86asub/netmap/netmap-evaluation/netmap/case_studies/blood \
    --analysis_dir /data_nfs/og86asub/netmap/netmap-evaluation/results/case_studies/blood/ \
    --experiment_name bd-rhap-rep2-X \
    --cluster_column 'leiden_remap' 


pixi run python src/manuscript_blood/comparative_analysis/jaccard_marker_genes.py \
    --adata_raw  /data_nfs/og86asub/netmap/netmap-evaluation/data/blood/reprocessed/10x-rep1-kallisto-cellbender/10x-rep1-kallisto-cellbender.h5ad \
    --results_dir /data_nfs/og86asub/netmap/netmap-evaluation/netmap/case_studies/blood \
    --analysis_dir /data_nfs/og86asub/netmap/netmap-evaluation/results/case_studies/blood/ \
    --experiment_name 10x-rep1-kallisto-cellbender-X \
    --cluster_column 'leiden_remap' 

pixi run python src/manuscript_blood/comparative_analysis/jaccard_marker_genes.py \
    --adata_raw  /data_nfs/og86asub/netmap/netmap-evaluation/data/blood/reprocessed/10x-rep2-kallisto-cellbender/10x-rep2-kallisto-cellbender.h5ad \
    --results_dir /data_nfs/og86asub/netmap/netmap-evaluation/netmap/case_studies/blood \
    --analysis_dir /data_nfs/og86asub/netmap/netmap-evaluation/results/case_studies/blood/ \
    --experiment_name blood-10x-rep2-kallisto-cellbender-X \
    --cluster_column 'leiden_remap' 