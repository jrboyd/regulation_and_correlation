#for input linc, return mRNA within certain distance, 
#check if mRNA are DE and also if their expression correlates with that of the linc.
if(exists("script_dir")) setwd(script_dir)
script_dir = getwd()
data_dir = "Z:/Coralee/From_Joe/Human_Guilt-by-association/"
setwd(data_dir)

source("features_within_dist.R")
source("n_closest_features.R")
source("parse_gtf.R")
source('fetch_ucsc_image.R')


#load reference with gene features that covers mRNA and lncRNA - make sure versions match!
ensg_ref = parse_gtf("gencode.v25.annotation_geneonly.gff", rownames_attrib = "gene_id", feature_type = "gene")

#lncRNA DE result files
linc_DE_files = dir(pattern = "_Significant_lncRNAs_v25.txt", full.names = T)
all_DE_lincs = character()
for(f in linc_DE_files){
  res = rownames(read.table(f))
  all_DE_lincs = union(all_DE_lincs, res)
}

#mRNA DE result files
mRNA_DE_files = dir(pattern = "_Significant_mRNAs_v25.txt", full.names = T)
all_DE_mRNA = character()
for(f in mRNA_DE_files){
  res = rownames(read.table(f))
  all_DE_mRNA = union(all_DE_mRNA, res)
}

#find mRNA nearby DE lincs
mRNA_nearby = n_closest_features(test_res = all_DE_lincs, ref_dict = ensg_ref, n_closest = 4, strand = "either")

#filter nearby mRNA by whether they are DE
mRNA_nearby_DE = lapply(mRNA_nearby, function(x){
  intersect(x, all_DE_mRNA)
})
is_DE = sapply(mRNA_nearby_DE, length) > 0
n_pos = sum(is_DE)
perc_pos = round((n_pos / length(all_DE_lincs)) * 100, 2)
print(paste0(perc_pos, "% of DE lincs near DE mRNA"))

#this file must contain counts for all mRNA and lncs in DE results
norm_counts = read.table("DESeq2_hMSC_normalized_counts_v25.txt")

mRNA_nearby_corr = lapply(names(mRNA_nearby_DE)[is_DE], function(x){
  mRNA = mRNA_nearby_DE[[x]]
  dat = norm_counts[c(x, mRNA),]
  dat_cor = cor(t(dat))[1,-1]
  names(dat_cor) = mRNA
  return(dat_cor)
})
names(mRNA_nearby_corr) = names(mRNA_nearby_DE)[is_DE]




plot_all_nearby = function(DE_lists){
  MAX = max(sapply(DE_lists, length)) + 1
  l_mat = matrix(c(rep(1, MAX), 1:MAX + 1), ncol = MAX, byrow = T)
  l_hei = c(1,3)
  n = 1
  mRNA_nearby_plots = lapply(names(DE_lists), function(name){
    print(n / length(DE_lists)); n <<- n + 1
    layout(l_mat, heights = l_hei)
    mRNA = DE_lists[[name]]
    mRNA_o = order(mRNA_nearby_corr[[name]][mRNA], decreasing = T)
    mRNA = mRNA[mRNA_o]
    ref = ensg_ref[c(name, mRNA),]
    s = min(ref$start)
    e = max(ref$end)
    chrm = ref$chrm[1]
    ucsc_img = fetch_ucsc_image(chrm, s, e)
    sH = dev.size()[2] / 4 #height and width on screen
    sW = dev.size()[1]
    iH = nrow(ucsc_img) #height and width of image
    iW = ncol(ucsc_img)
    pH = max(iH, 0)
    pW = max(iW, 0)
    par(mai = rep(0, 4))
    
    ratio = (sH / sW) / (iH / iW)
    if(ratio <= 1){
      pW = pW / ratio
      xmin = .5 * pW - .5 * iW
      xmax = .5 * pW + .5 * iW
      ymin = 0
      ymax = pH
    }else{
      pH = pH * ratio
      xmin = 0
      xmax = pW
      ymin = .5 * pH - .5 * iH
      ymax = .5 * pH + .5 * iH
    }
    plot(c(0, pW), c(0, pH), xlab = "", ylab = "", axes = F, type = "n")
    plot_ucsc_image(ucsc_img, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, lab_width = 0)
    
    dat = norm_counts[c(name, mRNA),]
    # grps = sapply(strsplit(colnames(dat), "_"), function(x){
    #   paste(x[1:2], collapse = "_")
    # })
    day = as.integer(sub("day", "", sapply(strsplit(colnames(dat), "_"), function(x){
      x[3]
    })))
    
    for(i in 1:nrow(dat)){
      x_wiggle = rnorm(ncol(dat), sd = .3)
      # y_wiggle = rnorm(ncol(dat), sd = .01 * max(dat[i,]))
      
      par(mai = rep(.8, 4))
      ylab = ""; if(i == 1) ylab = "Normalized HTSeq-Counts"
      plot_title = paste0("(", ref$strand[i], ") ", ref$gene_name[i])
      if(i > 1) plot_title = paste(plot_title, round(mRNA_nearby_corr[[name]][rownames(ref)[i]], 2), "r2")
      if(i == 1) plot_title = paste(plot_title, "\n", name)
      plot(x = day + x_wiggle, y = dat[i, ], axes = F, xlab = "Day", ylab = ylab)
      title(plot_title, cex.main = 1)
      axis(side = 2)
      axis(side = 1, at = unique(day))
      y_means = sapply(unique(day), function(x){
        keep = day == x
        rowMeans(dat[i, keep])
      })
      lines(unique(day), y_means)
    }
  })
}

DE_mRNA_nearby_DE = mRNA_nearby_DE[is_DE]

pdf("DE_lincs_with_DE_mRNA_full.human.V25.pdf", width = 2.5 * 4, height = 2.5 * 4 / 3)
o = order(sapply(mRNA_nearby_corr, function(x)max(abs(x))), decreasing = T)
plot_all_nearby(DE_mRNA_nearby_DE[o])
dev.off()

pdf("DE_lincs_with_DE_mRNA_top20.human.V25.pdf", width = 2.5 * 4, height = 2.5 * 4 / 3)
o = order(sapply(mRNA_nearby_corr, function(x)max((x))), decreasing = T)
plot_all_nearby(DE_mRNA_nearby_DE[o][1:20])
dev.off()

pdf("DE_lincs_with_DE_mRNA_bottom20.human.V25.pdf", width = 2.5 * 4, height = 2.5 * 4 / 3)
o = order(sapply(mRNA_nearby_corr, function(x)min((x))), decreasing = F)
plot_all_nearby(DE_mRNA_nearby_DE[o][1:20])
dev.off()
