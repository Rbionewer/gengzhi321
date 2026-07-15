Generated NetAct input files

Files:
- expression_fpkm.tsv: gene-by-sample FPKM matrix (required expression input).
- expression_counts.tsv: gene-by-sample read count matrix.
- sample_metadata.tsv: sample annotations (WT/Toluene groups).
- gene_ranking_toluene_vs_wt.tsv: log2 fold-change ranking for enrichment workflows.
- tf_candidates_in_dataset.tsv: TFs present in expression matrix, based on pySCENIC human TF list.
- tf_expression_fpkm.tsv: expression matrix restricted to TF candidates.
- reference/hs_hgnc_tfs.txt: downloaded human TF reference list.

Gene symbols were extracted from the NT annotation column by parsing symbols inside parentheses, then aggregated by mean FPKM / summed read counts across contigs.
