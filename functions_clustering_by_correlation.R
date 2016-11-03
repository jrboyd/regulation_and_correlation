if(exists("script_dir")) setwd(script_dir)
source("parse_gtf.R")
script_dir = getwd()
data_dir = "Z:/Coralee/From_Joe/Human_Guilt-by-association/"
setwd(data_dir)

ensg_ref = parse_gtf("gencode.v25.annotation_geneonly.gff", rownames_attrib = "gene_id", feature_type = "gene", additional_attrib = "gene_type")
genetype = ensg_ref[, c("gene_name", "gene_type")]
ensg2sym = as.character(genetype[,1])
names(ensg2sym) = rownames(genetype)
norm_mRNA = read.table("DESeq2_hMSC_normalized_counts_mRNA_v25.txt")
norm_lincs = read.table("DESeq2_hMSC_normalized_counts_lncRNA_v25.txt")
is_sense_intronic = genetype[rownames(norm_lincs),]$gene_type == "sense_intronic"
norm_lincs = norm_lincs[!is_sense_intronic,]
is_num = grepl("day", colnames(norm_mRNA))
cors = cor(t(log10(norm_lincs[,is_num]+1)), t(log10(norm_mRNA[,is_num]+1)))#columns are protein_coding, rows are lincs
library(gplots)
source("H:/R_workspace/jrb_R_scripts/heatmap.3-split.R")
source("H:/R_workspace/jrb_R_scripts/heatmap.3-kmeans_wrapper.R")

png("cor_hmaps1.png", width = 2400, height = 2400)
hmap_res = heatmap.3_kmeans_wrapper(t(cors), nclust = 50, skip_plot = F)
dev.off()

png("cor_hmaps2.png", width = 2400, height = 2400)
hmap_res2 = heatmap.3_kmeans_wrapper(t(hmap_res[[3]]), nclust = 60, skip_plot = F)
dev.off()

res2clust_infor = function(res){
  sizes = res[[1]]
  ends = cumsum(sizes)
  starts = c(1, ends[-length(ends)]+1)
  return(list(sizes = sizes, starts = starts, ends = ends))
}

as_plotted = hmap_res2[[3]]
#plot profiles
col_clust = res2clust_infor(hmap_res)
row_clust = res2clust_infor(hmap_res2)

sum_mat = t(sapply(1:length(col_clust$sizes), function(x){
  sapply(1:length(row_clust$sizes), function(y){
    print(paste(x,y))
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
pdf("avg correlation between clusters.pdf")
heatmap.2(sum_mat, Colv = F, Rowv = F, trace = "n", col = colors, 
          srtRow = 90, labRow = c(rep("", 33), "lincRNA clusters"), 
          srtCol = 0, labCol = c(rep("", 33), "mRNA clusters"), 
          main = "Average correlation between clusters")
dev.off()

pos_corr = list()
neg_corr = list()
sum_mat = t(sapply(1:length(col_clust$sizes), function(x){
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


agg_mRNA = t(apply(norm_mRNA[,is_num], 1, function(x){
  x = c(mean(x[1:3]), mean(x[4:6]), mean(x[7:9]), mean(x[10:12]))
  return(log2(x + 8)-3)
})); colnames(agg_mRNA) = paste("day", 0:3*7)
agg_linc = t(apply(norm_lincs[,is_num], 1, function(x){
  x = c(mean(x[1:3]), mean(x[4:6]), mean(x[7:9]), mean(x[10:12]))
  return(log2(x + 8)-3)
})); colnames(agg_linc) = paste("day", 0:3*7)


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
