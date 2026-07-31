# for i in {0..9}; do
#   EXP_NAME="10x-rep1-kallisto-cellbender-X-$i"
  
#   echo "------------------------------------------------"
#   echo "Running analysis for: $EXP_NAME"
#   echo "------------------------------------------------"

#   python src/manuscript_blood/netmap_downstream/batch_analyse_output.py \
#     --adata_raw '/data_nfs/og86asub/netmap/netmap-evaluation/data/blood/reprocessed/10x-rep1-kallisto-cellbender/10x-rep1-kallisto-cellbender.h5ad' \
#     --collectri '/data_nfs/og86asub/netmap/netmap-evaluation/data/input_network/collectri.tsv' \
#     --results_dir /data_nfs/og86asub/netmap/netmap-evaluation/netmap/case_studies/blood/ \
#     --experiment_name "$EXP_NAME" \
#     --min_score 0.85 \
#     --cluster_column 'leiden_remap' \
#     --analysis_dir '/data_nfs/og86asub/netmap/netmap-evaluation/results/case_studies/blood/'

# done


# for i in {0..9}; do
#   EXP_NAME="bd-rhap-rep2-X-$i"
  
#   echo "------------------------------------------------"
#   echo "Running analysis for: $EXP_NAME"
#   echo "------------------------------------------------"

#   python src/manuscript_blood/netmap_downstream/batch_analyse_output.py \
#     --adata_raw '/data_nfs/og86asub/netmap/netmap-evaluation/data/blood/reprocessed/bd-rhap-rep2/bd-rhap-rep2.h5ad' \
#     --collectri '/data_nfs/og86asub/netmap/netmap-evaluation/data/input_network/collectri.tsv' \
#     --results_dir /data_nfs/og86asub/netmap/netmap-evaluation/netmap/case_studies/blood/ \
#     --experiment_name "$EXP_NAME" \
#     --min_score 0.85 \
#     --cluster_column 'leiden_remap' \
#     --analysis_dir '/data_nfs/og86asub/netmap/netmap-evaluation/results/case_studies/blood/'

# done


# EXP_NAME="bd-rhap-rep2-X"
# python src/manuscript_blood/netmap_downstream/batch_analyse_output.py \
#   --adata_raw '/data_nfs/og86asub/netmap/netmap-evaluation/data/blood/reprocessed/bd-rhap-rep2/bd-rhap-rep2.h5ad' \
#   --collectri '/data_nfs/og86asub/netmap/netmap-evaluation/data/input_network/collectri.tsv' \
#   --results_dir /data_nfs/og86asub/netmap/netmap-evaluation/netmap/case_studies/blood/ \
#   --experiment_name "$EXP_NAME" \
#   --min_score 0.85 \
#   --cluster_column 'leiden_remap' \
#   --analysis_dir '/data_nfs/og86asub/netmap/netmap-evaluation/results/case_studies/blood/'\
#   --save_objects \
#   --rerun


# EXP_NAME="bd-rhap-rep1-X"
# python src/manuscript_blood/netmap_downstream/batch_analyse_output.py \
#   --adata_raw '/data_nfs/og86asub/netmap/netmap-evaluation/data/blood/reprocessed/bd-rhap-rep1/bd-rhap-rep1.h5ad' \
#   --collectri '/data_nfs/og86asub/netmap/netmap-evaluation/data/input_network/collectri.tsv' \
#   --results_dir /data_nfs/og86asub/netmap/netmap-evaluation/netmap/case_studies/blood/ \
#   --experiment_name "$EXP_NAME" \
#   --min_score 0.85 \
#   --cluster_column 'leiden_remap' \
#   --analysis_dir '/data_nfs/og86asub/netmap/netmap-evaluation/results/case_studies/blood/' \
#   --save_objects \
#   --rerun


# EXP_NAME="10x-rep1-kallisto-cellbender-X"
# python src/manuscript_blood/netmap_downstream/batch_analyse_output.py \
#   --adata_raw '/data_nfs/og86asub/netmap/netmap-evaluation/data/blood/reprocessed/10x-rep1-kallisto-cellbender/10x-rep1-kallisto-cellbender.h5ad' \
#   --collectri '/data_nfs/og86asub/netmap/netmap-evaluation/data/input_network/collectri.tsv' \
#   --results_dir /data_nfs/og86asub/netmap/netmap-evaluation/netmap/case_studies/blood/ \
#   --experiment_name "$EXP_NAME" \
#   --min_score 0.85 \
#   --cluster_column 'leiden_remap' \
#   --analysis_dir '/data_nfs/og86asub/netmap/netmap-evaluation/results/case_studies/blood/'\
#   --save_objects \
#   --rerun

# EXP_NAME="blood-10x-rep2-kallisto-cellbender-X"
# python src/manuscript_blood/netmap_downstream/batch_analyse_output.py \
#   --adata_raw '/data_nfs/og86asub/netmap/netmap-evaluation/data/blood/reprocessed/10x-rep2-kallisto-cellbender/10x-rep2-kallisto-cellbender.h5ad' \
#   --collectri '/data_nfs/og86asub/netmap/netmap-evaluation/data/input_network/collectri.tsv' \
#   --results_dir /data_nfs/og86asub/netmap/netmap-evaluation/netmap/case_studies/blood/ \
#   --experiment_name "$EXP_NAME" \
#   --min_score 0.85 \
#   --cluster_column 'leiden_remap' \
#   --analysis_dir '/data_nfs/og86asub/netmap/netmap-evaluation/results/case_studies/blood/' \
#   --save_objects \
#   --rerun


EXP_NAME="10x-rep2-kallisto-cellbender-markers"
python src/manuscript_blood/netmap_downstream/batch_analyse_output.py \
  --adata_raw '/data_nfs/og86asub/netmap/netmap-evaluation/data/blood/reprocessed/10x-rep2-kallisto-cellbender/10x-rep2-kallisto-cellbender_with_markers.h5ad' \
  --collectri '/data_nfs/og86asub/netmap/netmap-evaluation/data/input_network/collectri.tsv' \
  --results_dir /data_nfs/og86asub/netmap/netmap-evaluation/netmap/case_studies/blood/ \
  --experiment_name "$EXP_NAME" \
  --min_score 0.85 \
  --cluster_column 'leiden_remap' \
  --analysis_dir '/data_nfs/og86asub/netmap/netmap-evaluation/results/case_studies/blood/' \
  --save_objects \
  --rerun

EXP_NAME="10x-rep1-kallisto-cellbender-markers"
python src/manuscript_blood/netmap_downstream/batch_analyse_output.py \
  --adata_raw '/data_nfs/og86asub/netmap/netmap-evaluation/data/blood/reprocessed/10x-rep1-kallisto-cellbender/10x-rep1-kallisto-cellbender_with_markers.h5ad' \
  --collectri '/data_nfs/og86asub/netmap/netmap-evaluation/data/input_network/collectri.tsv' \
  --results_dir /data_nfs/og86asub/netmap/netmap-evaluation/netmap/case_studies/blood/ \
  --experiment_name "$EXP_NAME" \
  --min_score 0.85 \
  --cluster_column 'leiden_remap' \
  --analysis_dir '/data_nfs/og86asub/netmap/netmap-evaluation/results/case_studies/blood/' \
  --save_objects \
  --rerun


EXP_NAME="bd-rhap-rep1-markers"
python src/manuscript_blood/netmap_downstream/batch_analyse_output.py \
  --adata_raw '/data_nfs/og86asub/netmap/netmap-evaluation/data/blood/reprocessed/bd-rhap-rep1/bd-rhap-rep1_with_markers.h5ad' \
  --collectri '/data_nfs/og86asub/netmap/netmap-evaluation/data/input_network/collectri.tsv' \
  --results_dir /data_nfs/og86asub/netmap/netmap-evaluation/netmap/case_studies/blood/ \
  --experiment_name "$EXP_NAME" \
  --min_score 0.85 \
  --cluster_column 'leiden_remap' \
  --analysis_dir '/data_nfs/og86asub/netmap/netmap-evaluation/results/case_studies/blood/' \
  --save_objects \
  --rerun

EXP_NAME="bd-rhap-rep2-markers"
python src/manuscript_blood/netmap_downstream/batch_analyse_output.py \
  --adata_raw '/data_nfs/og86asub/netmap/netmap-evaluation/data/blood/reprocessed/bd-rhap-rep2/bd-rhap-rep2_with_markers.h5ad' \
  --collectri '/data_nfs/og86asub/netmap/netmap-evaluation/data/input_network/collectri.tsv' \
  --results_dir /data_nfs/og86asub/netmap/netmap-evaluation/netmap/case_studies/blood/ \
  --experiment_name "$EXP_NAME" \
  --min_score 0.85 \
  --cluster_column 'leiden_remap' \
  --analysis_dir '/data_nfs/og86asub/netmap/netmap-evaluation/results/case_studies/blood/' \
  --save_objects \
  --rerun
