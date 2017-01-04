if(exists("MAIN_DIR")) setwd(MAIN_DIR)
source("parse_config.R")
parse_config("example_config.txt")

setwd(script_dir)
source("parse_gtf.R")
source("heatmap.3-split.R")
source("heatmap.3-kmeans_wrapper.R")

setwd(data_dir)
setup_clustering_by_correlation(ref_file = REF_FILE, 
                                mRNA_counts_file = MRNA_COUNT_FILE, 
                                linc_counts_file = LNC_COUNT_FILE)

setwd(output_dir)
cluster_file = "cluster_info.save"
if(file.exists(cluster_file)){
  warning("loading previously calculated cluster info and omitting large heatmap.")
  load(cluster_file)
  
}else{
  cluster_info = plot_2d_kmeans(hmap_data = cors, 
                                nclust1 = LNC_NCLUST, 
                                nclust2 = MRNA_NCLUST, 
                                kmean_seed = KMEAN_SEED)
  save(cluster_info, file = cluster_file)
}

plot_summary_heatmap(as_plotted = cluster_info$as_plotted, col_clust = cluster_info$col_clust, row_clust = cluster_info$row_clust)
