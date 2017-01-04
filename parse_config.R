#parse a config file to assign variables that will control analysis
parse_config = function(config_file = "example_config.txt"){
  cfg = read.table(config_file, comment.char = "#", sep = "=", stringsAsFactors = F)  
  for(i in 1:nrow(cfg)){
    assign(cfg[i,1], cfg[i,2], envir = .GlobalEnv)
  }
  
  setwd(MAIN_DIR)
  script_dir <<- normalizePath(SCRIPT_DIR)
  dir.create(OUTPUT_DIR, showWarnings = F)
  output_dir <<- normalizePath(OUTPUT_DIR)
  data_dir <<- normalizePath(DATA_DIR)
}




