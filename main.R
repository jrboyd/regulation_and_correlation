#this script sets the configuration to use and loads it.
#loads helper functions and executes pipelines.

#run install packages if you're missing libraries
library("RCurl")
library("httr")
library("png")
library("pbapply")
library("gplots")
library("gtools")
library("xlsx")
library("GenomicRanges")

CFG_FILE = "example_config.txt"
if(exists("MAIN_DIR")) setwd(MAIN_DIR)
source("parse_config.R")

#procure example data if needed
if(CFG_FILE == "example_config.txt"){
  if(!dir.exists("example_data")){
    print("downloading example data...")
    library(httr)
    set_config(config(ssl_verifypeer = 0L))
    download.file(url = "https://galaxy.med.uvm.edu/static/example_data.tar.gz", destfile = "example_data.tar.gz", extra = "--no-check-certificate")
    print("untar...")
    untar("example_data.tar.gz")
  }
}

#load configuration variables
parse_config(CFG_FILE)

#for input linc, return mRNA within certain distance, 
#check if mRNA are DE and also if their expression correlates with that of the linc.


#load functions
setwd(script_dir)
source("parse_gtf.R")
source("parse_exon_gtf_and_merge.R")
source('fetch_ucsc_image.R')
source("functions_correlate_cis_regulatory.R")
source("heatmap.3-split.R")
source("heatmap.3-kmeans_wrapper.R")
source("functions_clustering_by_correlation.R")

#setup reference
setwd(data_dir)
if(!exists("ensg_ref")){
  print(paste("loading ensg_ref from:", REF_FILE))
  ensg_ref = parse_gtf(REF_FILE, rownames_attrib = "gene_id", feature_type = "gene", additional_attrib = "gene_type")
  if(exists("ADD_EXONS_REF_FILE") && ADD_EXONS_REF_FILE != ""){
    to_add = parse_exon_gtf_and_merge(ADD_EXONS_REF_FILE)
    ensg_ref = rbind(ensg_ref, to_add)
  }
}else{
  warning("using previously loaded ensg_ref.")
}

genetype = ensg_ref[, c("gene_name", "gene_type")]
ensg2sym = as.character(genetype[,1])
names(ensg2sym) = rownames(genetype)



# run the pipelines
setwd(MAIN_DIR)
source("run_clustering_by_correlation.R")
setwd(MAIN_DIR)
source("run_correlate_cis_regulatory.R")
setwd(MAIN_DIR)