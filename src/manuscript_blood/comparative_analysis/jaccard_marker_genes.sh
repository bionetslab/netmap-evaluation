pixi run python src/manuscript_blood/comparative_analysis/jaccard_marker_genes.py \
    --adata_raw  /data_nfs/og86asub/netmap/netmap-evaluation/data/blood/reprocessed/bd-rhap-rep1/bd-rhap-rep1_with_markers.h5ad \
    --results_dir /data_nfs/og86asub/netmap/netmap-evaluation/netmap/case_studies/blood \
    --analysis_dir /data_nfs/og86asub/netmap/netmap-evaluation/results/case_studies/blood/ \
    --experiment_name bd-rhap-rep1-markers \
    --cluster_column 'leiden_remap' 

pixi run python src/manuscript_blood/comparative_analysis/jaccard_marker_genes.py \
    --adata_raw  /data_nfs/og86asub/netmap/netmap-evaluation/data/blood/reprocessed/bd-rhap-rep2/bd-rhap-rep2_with_markers.h5ad \
    --results_dir /data_nfs/og86asub/netmap/netmap-evaluation/netmap/case_studies/blood \
    --analysis_dir /data_nfs/og86asub/netmap/netmap-evaluation/results/case_studies/blood/ \
    --experiment_name bd-rhap-rep2-markers \
    --cluster_column 'leiden_remap' 


pixi run python src/manuscript_blood/comparative_analysis/jaccard_marker_genes.py \
    --adata_raw  /data_nfs/og86asub/netmap/netmap-evaluation/data/blood/reprocessed/10x-rep1-kallisto-cellbender/10x-rep1-kallisto-cellbender_with_markers.h5ad \
    --results_dir /data_nfs/og86asub/netmap/netmap-evaluation/netmap/case_studies/blood \
    --analysis_dir /data_nfs/og86asub/netmap/netmap-evaluation/results/case_studies/blood/ \
    --experiment_name 10x-rep1-kallisto-cellbender-markers \
    --cluster_column 'leiden_remap' 

pixi run python src/manuscript_blood/comparative_analysis/jaccard_marker_genes.py \
    --adata_raw  /data_nfs/og86asub/netmap/netmap-evaluation/data/blood/reprocessed/10x-rep2-kallisto-cellbender/10x-rep2-kallisto-cellbender_with_markers.h5ad \
    --results_dir /data_nfs/og86asub/netmap/netmap-evaluation/netmap/case_studies/blood \
    --analysis_dir /data_nfs/og86asub/netmap/netmap-evaluation/results/case_studies/blood/ \
    --experiment_name 10x-rep2-kallisto-cellbender-markers \
    --cluster_column 'leiden_remap' 