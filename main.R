CFG_FILE = "example_config.txt"
if(exists("MAIN_DIR")) setwd(MAIN_DIR)
source("parse_config.R")

#procure example data if needed
if(CFG_FILE == "example_config.txt"){
  if(!dir.exists("example_data")){
    print("downloading example data...")
    library(httr)
    set_config(config(ssl_verifypeer = 0L))
    download.file(url = "https://galaxy.med.uvm.edu/static/example_data.tar.gz", destfile = "example_data.tar.gz")
    print("untar...")
    untar("example_data.tar.gz")
  }
}

#load configuration variables
parse_config(CFG_FILE)

#for input linc, return mRNA within certain distance, 
#check if mRNA are DE and also if their expression correlates with that of the linc.
library(GenomicRanges)
library(pbapply)

#load functions
setwd(script_dir)
source("parse_gtf.R")
source('fetch_ucsc_image.R')
source("functions_correlate_cis_regulatory.R")
source("heatmap.3-split.R")
source("heatmap.3-kmeans_wrapper.R")
source("functions_clustering_by_correlation.R")

# run the pipelines
source("run_clustering_by_correlation.R")
setwd(MAIN_DIR)
source("run_correlate_cis_regulatory.R")
