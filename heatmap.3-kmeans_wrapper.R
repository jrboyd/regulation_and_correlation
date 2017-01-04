library(gplots)
library(xlsx)
library(RColorBrewer)
# source('H:/R_workspace/jrb_R_scripts/heatmap.3-split.R')
# source('H:/R_workspace/jrb_R_scripts/heatmap.ngsplots_kmeans_with_sideplot.R')
#a wrapper for heatmap.2, primarily to add kmeans clustering(much faster for larger datasets)
#changes default colors
#outputs clusters to excel file
example_heatmap.3_kmeans_wrapper = function(n_row = 2000, n_col = 6, n_spikes = 4, n_clust = 5){
  #creates ranomized data with dimensions n_row and n_col
  #n_spikes is the number of simulated clusters added to dataset
  #n_clust is the number of clusters used by kmeans (should capture the n_spikes clusters)
  MIN = 95
  MAX = 105
  set.seed(1)
  dat = matrix(runif(n_row * n_col, min = MIN, max = MAX), nrow = n_row, ncol = n_col)
  rand =  function(n, of){
    order(runif(of))[1:n]
  }
  for(i in 1:n_spikes){
    n_choose = (1 / n_spikes)*n_row
    rows_key = rand(n_choose, n_row )
    cols_key = rand(2, of = n_col)
    bump = runif(1, min = 10, max = 30)
    print(round(bump))
    if(runif(1) > .5){
      bump = -bump
    }
    dat[rows_key, cols_key] = dat[rows_key, cols_key] + bump
  }
  heatmap.3_kmeans_wrapper(dat = dat, nclust = n_clust)
}

heatmap.3_kmeans_wrapper = function(dat, nclust = 4, col = c("blue", 'white', 'red'), xlsxname = NULL, skip_plot = F, seed = 1, ...){
  set.seed(seed)
  kclust = kmeans(dat, centers = nclust)
  o = order(kclust$cluster)
  kclust$cluster = kclust$cluster[o]
  dat = dat[o,]
  
  for(i in 1:nclust){
    keep = kclust$cluster == i
    subset = dat[keep,]
    o = order(apply(subset, 1, function(x){
      return(max(x)-min(x))
    }), decreasing = T)
    o = order(apply(subset, 1, max), decreasing = T)
    rownames(dat)[keep] = rownames(subset)[o]
    dat[keep,] = subset[o,]
  }
  
  o = order(rowSums(kclust$centers), decreasing = T)
  kclust$centers = kclust$centers[o,]
  plot_dat = matrix(0, nrow = 0, ncol = ncol(dat))
  new_cluster = integer()
  colnames(plot_dat) = colnames(dat)
  j = 1
  for(i in o){
    keep = kclust$cluster == i
    subset = dat[keep,]
    plot_dat = rbind(plot_dat, subset)
    new_cluster = c(new_cluster, rep(j, sum(keep)))
    j = j + 1
  }
  
  clust_sizes = sapply(o, function(x){
    sum(kclust$cluster == x)
  })
  
  cr = colorRamp(col)
  colors = rgb(cr(0:100/100)/255)
  # pdf("diff_mRNA_mMSCs_heatmap.pdf", width = 6, height = 8)
  seps = cumsum(sapply(1:nclust, function(x)sum(x == new_cluster)))
  
  if(!skip_plot){
    override_o = cbind(1:nclust, clust_sizes)
    res = heatmap.3(plot_dat, trace = 'n', Rowv = F, Colv = F, scale = 'n',  cexCol = 1.6, cexRow = .4, col = colors, density.info = 'n', key.title = "", labRow = "", override_o = override_o, ...)
    res = res[1:3]
    if(!is.null(xlsxname)){
      if(is.null(rownames(plot_dat))){
        rownames(plot_dat) = paste0("row_", 1:nrow(plot_dat))
      }  
      if(is.null(colnames(plot_dat))){
        colnames(plot_dat) = paste0("column_", 1:ncol(plot_dat))
      } 
      
      clust_memb_wb = createWorkbook()
      sheets = list()
      for(i in 1:nclust){
        new_sheet = createSheet(clust_memb_wb, paste("cluster", i))
        keep = new_cluster == i
        members = rownames(plot_dat)[keep]
        df = as.data.frame(plot_dat[members,])
        # df = cbind(rownames(my_fpkm[members,]), df)
        addDataFrame(df, new_sheet)
        
      }
      bad = F
      
      newxlsx = xlsxname
      i = 1
      while(file.exists(newxlsx)){
        bad = T
        newxlsx = sub(".xlsx", paste0("(", i ,").xlsx"), xlsxname)
        i = i + 1
      }
      xlsxname = newxlsx
      saveWorkbook(clust_memb_wb, xlsxname)
    }
  }else{
    res = list(clust_sizes, 1:nrow(plot_dat), plot_dat)
  }
  res = res[1:3]
  names(res) = c("clust_sizes", "o", "dat")
  return(res)
}

if(F) example_heatmap.3_kmeans_wrapper()
