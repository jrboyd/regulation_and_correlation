#parse a config file to assign variables that will control analysis
parse_config = function(config_file = "example_config.txt"){
  cfg = read.table(config_file, comment.char = "#", sep = "=", stringsAsFactors = F)  
  for(i in 1:nrow(cfg)){
    val = cfg[i,2]
    # print(val)
    #convert to integer/numeric if possible
    suppressWarnings(
    if(!is.na(as.integer(val))){
      if(as.character(as.integer(val)) == val){
        val = as.integer(val)
      }else if(as.character(as.numeric(val)) == val){
        val = as.integer(val)
      }
    })
    # print(val)
    assign(cfg[i,1], val, envir = .GlobalEnv)
  }
  
  setwd(MAIN_DIR)
  script_dir <<- normalizePath(SCRIPT_DIR)
  dir.create(OUTPUT_DIR, showWarnings = F)
  output_dir <<- normalizePath(OUTPUT_DIR)
  data_dir <<- normalizePath(DATA_DIR)
}




