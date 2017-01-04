

#sets a series of global values after loading data
setup_clustering_by_correlation = function(ref_file, mRNA_counts_file, linc_counts_file){
  if(!exists("ensg_ref")){
    print(paste("loading ensg_ref from:", ref_file))
    ensg_ref = parse_gtf(ref_file, rownames_attrib = "gene_id", feature_type = "gene", additional_attrib = "gene_type")
  }else{
    warning("using previously loaded ensg_ref.")
  }
  genetype = ensg_ref[, c("gene_name", "gene_type")]
  ensg2sym = as.character(genetype[,1])
  names(ensg2sym) = rownames(genetype)
  norm_mRNA = read.table(mRNA_counts_file)
  norm_lincs = read.table(linc_counts_file)
  #remove sense intronic lincs
  is_sense_intronic = genetype[rownames(norm_lincs),]$gene_type == "sense_intronic"
  norm_lincs = norm_lincs[!is_sense_intronic,]
  
  #filter down to columns that contain count data
  is_num = grepl(NUMERIC_COLUMN_PATTERN, colnames(norm_mRNA))
  #correlate lincs to mRNA after log scaling
  cors = cor(t(log10(norm_lincs[,is_num]+1)), t(log10(norm_mRNA[,is_num]+1)))#columns are protein_coding, rows are lincs
  
  #aggregate reps
  agg = function(to_agg){
    is_num = grepl(NUMERIC_COLUMN_PATTERN, colnames(to_agg))
    to_agg = to_agg[, is_num]
    col_ids = sapply(strsplit(colnames(to_agg), COLNAME_SEP), function(x)x[POOL_BY])
    done_agg = matrix(0, nrow = nrow(to_agg), ncol = length(unique(col_ids)))
    rownames(done_agg) = rownames(to_agg)
    colnames(done_agg) = unique(col_ids)
    for(id in unique(col_ids)){
      k = col_ids == id
      done_agg[, id] = rowMeans(to_agg[,k])
    }
    return(done_agg)
  }
  
  agg_mRNA = agg(norm_mRNA)
  agg_linc = agg(norm_lincs)
  
  ensg_ref <<- ensg_ref
  ensg2sym <<- ensg2sym
  norm_mRNA <<- norm_mRNA
  norm_lincs <<- norm_lincs
  cors <<- cors
  agg_mRNA <<- agg_mRNA
  agg_linc <<- agg_linc
}

plot_2d_kmeans = function(hmap_data, nclust1, nclust2, kmean_seed = 1){
  # png("cor_hmaps1.png", width = 2400, height = 2400)
  print("clustering lncs...")
  hmap_res = heatmap.3_kmeans_wrapper(t(hmap_data), nclust = nclust1, skip_plot = T, seed = kmean_seed)
  # dev.off()
  print("clustering mrnas...")
  setwd(output_dir)
  png("mRNA_to_linc_correlation_heatmap.png", width = 2400, height = 2400)
  setwd(MAIN_DIR)
  hmap_res2 = heatmap.3_kmeans_wrapper(t(hmap_res$dat), nclust = nclust2, skip_plot = F, seed = kmean_seed)
  dev.off()
  print("done!")
  as_plotted = hmap_res2$dat
  #plot profiles
  hmap_res2clust_infor = function(res){
    sizes = res$clust_sizes
    ends = cumsum(sizes)
    starts = c(1, ends[-length(ends)]+1)
    return(list(sizes = sizes, starts = starts, ends = ends))
  }
  col_clust = hmap_res2clust_infor(hmap_res)
  row_clust = hmap_res2clust_infor(hmap_res2)
  return(list(col_clust = col_clust, row_clust = row_clust, as_plotted = as_plotted))
}

plot_summary_heatmap = function(as_plotted, col_clust, row_clust){
  # pb <- txtProgressBar(min = 0, max = sum(col_clust$sizes) * sum(row_clust$sizes), style = 3)
  sum_mat = t(pbsapply(1:length(col_clust$sizes), function(x){
    sapply(1:length(row_clust$sizes), function(y){
      # print(paste(x,y))
      x_start = col_clust$starts[x]
      x_end = col_clust$ends[x]
      y_start = row_clust$starts[y]
      y_end = row_clust$ends[y]
      ensg_lincs = rownames(as_plotted)[y_start:y_end]
      ensg_mRNA = colnames(as_plotted)[x_start:x_end]
      return(mean(as_plotted[ensg_lincs, ensg_mRNA]))
    })
  }))
  colors = rgb(colorRamp(c("darkblue", "white", "darkred"))(0:21/21)/255)
  # setwd(output_dir)
  pdf("avg correlation between clusters.pdf")
  heatmap.2(sum_mat, Colv = F, Rowv = F, trace = "n", col = colors, dendrogram = "n",
            srtRow = 90, labRow = c(rep("", 33), "lincRNA clusters"), 
            srtCol = 0, labCol = c(rep("", 33), "mRNA clusters"), 
            main = "Average correlation between clusters")
  dev.off()
  # setwd(MAIN_DIR)
}

pos_corr = list()
neg_corr = list()
sum_mat = t(pbsapply(1:length(col_clust$sizes), function(x){
  sapply(1:length(row_clust$sizes), function(y){
    x_start = col_clust$starts[x]
    x_end = col_clust$ends[x]
    y_start = row_clust$starts[y]
    y_end = row_clust$ends[y]
    ensg_lincs = rownames(as_plotted)[y_start:y_end]
    ensg_mRNA = colnames(as_plotted)[x_start:x_end]
    avg = mean(as_plotted[ensg_lincs, ensg_mRNA])
    if(avg > .7){
      pos_corr[[paste(y, x, sep = ',')]] <<- list(ensg_mRNA = ensg_mRNA, ensg_lincs = ensg_lincs)
    }else if(avg < -.7){
      neg_corr[[paste(y, x, sep = ',')]] <<- list(ensg_mRNA = ensg_mRNA, ensg_lincs = ensg_lincs)
    }
    return(mean(as_plotted[ensg_lincs, ensg_mRNA]))
  })
}))





plot_corr_blocks = function(corr_list, pdf_name){
  pdf(pdf_name)
  for(i in 1:length(corr_list)){
    nam = names(corr_list)[i]
    block = corr_list[[i]]
    corr_vals = as_plotted[block$ensg_lincs, block$ensg_mRNA]
    keep_linc = apply(corr_vals, 1, function(x)max(abs(x))) > .85
    keep_mRNA = apply(corr_vals, 2, function(x)max(abs(x))) > .85
    mRNA_dat = agg_mRNA[block$ensg_mRNA[keep_mRNA], ]
    mRNA_dat = mRNA_dat - mRNA_dat[,1]
    linc_dat = agg_linc[block$ensg_linc[keep_linc], ]
    linc_dat = linc_dat - linc_dat[,1]
    comb_dat = rbind(mRNA_dat, linc_dat)
    colors = c(rep("darkorange", nrow(mRNA_dat)), rep("darkblue", nrow(linc_dat)))
    ylim = range(comb_dat)
    plot(0, xlim = c(1,4), ylim = ylim, type = "n", axes = F, ylab = "log2 FC from day 0", xlab = "day")
    axis(side = 1, at = 1:4, labels = colnames(comb_dat))
    axis(side = 2)
    for(j in 1:nrow(comb_dat)){
      lines(comb_dat[j,], col = colors[j], lwd = 2)
    }
    y = as.integer(strsplit(nam, ",")[[1]][1])
    x = as.integer(strsplit(nam, ",")[[1]][2])
    title(paste0("linc cluster ", y, "\n", "mRNA cluster ", x, "\n", "r2 = ", round(sum_mat[x, y], 2)))
    legend("left", legend = c("lincRNAs", "mRNAs"), fill = c("darkblue", "darkorange"), bg = "white")
  }
  dev.off()
}
plot_corr_blocks(pos_corr, "positively_corr_plots.pdf")
plot_corr_blocks(neg_corr, "negatively_corr_plots.pdf")

write_corr_lists = function(corr_list, xlsx_name){
  wb = createWorkbook()
  sheet_mRNA_sym = createSheet(wb, "mRNA_sym")
  sheet_mRNA_ensg = createSheet(wb, "mRNA_ensg")
  sheet_linc_sym = createSheet(wb, "linc_sym")
  sheet_linc_ensg = createSheet(wb, "linc_ensg")
  for(i in 1:length(corr_list)){
    nam = names(corr_list)[i]
    block = corr_list[[i]]
    corr_vals = as_plotted[block$ensg_lincs, block$ensg_mRNA]
    keep_linc = apply(corr_vals, 1, function(x)max(abs(x))) > .85
    keep_mRNA = apply(corr_vals, 2, function(x)max(abs(x))) > .85
    ensg_mRNA = block$ensg_mRNA[keep_mRNA]
    ensg_lincs = block$ensg_lincs[keep_linc]
    sym_mRNA = unique(ensg2sym[ensg_mRNA])
    sym_lincs = unique(ensg2sym[ensg_lincs])
    y = as.integer(strsplit(nam, ",")[[1]][1])
    x = as.integer(strsplit(nam, ",")[[1]][2])
    cname = paste0("linc ", y, ", mRNA ", x)
    
    CB = CellBlock(sheet = sheet_mRNA_sym, startRow = 1, startColumn = i, noRows = length(sym_mRNA) + 1, noColumns = 1)
    CB.setColData(CB, c(cname, sym_mRNA), colIndex = 1)
    
    CB = CellBlock(sheet = sheet_mRNA_ensg, startRow = 1, startColumn = i, noRows = length(ensg_mRNA) + 1, noColumns = 1)
    CB.setColData(CB, c(cname, ensg_mRNA), colIndex = 1)
    
    CB = CellBlock(sheet = sheet_linc_sym, startRow = 1, startColumn = i, noRows = length(sym_lincs) + 1, noColumns = 1)
    CB.setColData(CB, c(cname, sym_lincs), colIndex = 1)
    
    CB = CellBlock(sheet = sheet_linc_ensg, startRow = 1, startColumn = i, noRows = length(ensg_lincs) + 1, noColumns = 1)
    CB.setColData(CB, c(cname, ensg_lincs), colIndex = 1)
  }
  saveWorkbook(wb, file = xlsx_name)
}

write_corr_lists(corr_list = pos_corr, xlsx_name = "positively_corr_lists.xlsx")
write_corr_lists(corr_list = neg_corr, xlsx_name = "negatively_corr_lists.xlsx")
