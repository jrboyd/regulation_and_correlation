setup_correlation_cis_regulatory = function(ref_file, lnc_de_pattern, mrna_de_pattern){
  # setwd(data_dir)
  #load reference with gene features that covers mRNA and lncRNA - make sure versions match!
  if(!exists("ensg_ref")){
    print(paste("loading ensg_ref from:", ensg_ref_file))
    ensg_ref = parse_gtf(ensg_ref_file, rownames_attrib = "gene_id", feature_type = "gene", additional_attrib = "gene_type")
  }else{
    warning("using previously loaded ensg_ref.")
  }
  #lncRNA DE result files
  linc_DE_files = dir(pattern = lnc_de_pattern, full.names = T)
  all_DE_lincs = character()
  for(f in linc_DE_files){
    res = rownames(read.table(f))
    all_DE_lincs = union(all_DE_lincs, res)
  }
  #mRNA DE result files
  mRNA_DE_files = dir(pattern = mrna_de_pattern, full.names = T)
  all_DE_mRNA = character()
  for(f in mRNA_DE_files){
    res = rownames(read.table(f))
    all_DE_mRNA = union(all_DE_mRNA, res)
  }
  
  ensg_ref <<- ensg_ref
  all_DE_lincs <<- all_DE_lincs
  all_DE_mRNA <<- all_DE_mRNA
  #this file must contain counts for all mRNA and lncs in DE results
  norm_counts <<- read.table(COUNTS_FILE)
}

#find mRNA nearby DE lincs
#source("Z://Coralee/scripts/guilt-by-assoication-and-correlation/n_closest_features.R")
n_closest_features = function(test_res, ref_dict, n_closest, max_distance = 10^5, strandedness = c("either", "same", "opposite")[1]){
  #test_res must be array of gene_id or whatever ref_dict uses as rownames or GenomicRanges
  #convert to GRanges
  ref_dict$seqnames = ref_dict$chrm
  
  #older version of GRanges aren't clever enough to create GRanges from data.frame based on colnames
  smart_gr_from_df = function(df){
    if(is.null(df$seqnames)){
      df$seqnames = df$chrm
    }
    df$chrm  = NULL
    # gr = makeGRangesFromDataFrame(df, keep.extra.columns = T)
    gr = GRanges(seqnames = df$seqnames, IRanges(df$start, df$end))
    if(!is.null(df$strand)) strand(gr) = df$strand
    for(cn in setdiff(colnames(df), c("seqnames", "start", "end", "strand", "width"))){
      mcols(gr)[[cn]] = df[[cn]]
    }
    names(gr) = rownames(df)
    return(gr)
  }
  
  ref_gr = smart_gr_from_df(ref_dict)
  if(class(test_res) == "character"){
    test_res = ref_gr[test_res,]
  }
  if(class(test_res) != "GRanges"){
    stop(paste("test_res must be character array of rownames of ref_dict or a GRanges object. was", class(test_res)))
  }
  # setTxtProgressBar(pb, i)
  n_nearest = pblapply(test_res, function(x){
    m = merge(data.frame(seqnames = seqnames(x)), as.data.frame(ref_gr), by = "seqnames")
    chrm_gr = smart_gr_from_df(m)
    keep = chrm_gr$gene_name != x$gene_name
    chrm_gr = chrm_gr[keep]
    if(strandedness == "opposite"){
      if(as.data.frame(x)$strand == "+"){
        strand(x) = "-"
      }else if(as.data.frame(x)$strand == "-"){
        strand(x) = "+"
      }
    }
    d = distance(x, chrm_gr, ignore.strand = strandedness == "either")
    chrm_gr$distance = d
    o = order(d, decreasing = F, na.last = T)
    chrm_gr = chrm_gr[o][1:n_closest]
    keep = chrm_gr$distance <= max_distance
    chrm_gr = chrm_gr[keep]
    chrm_gr$gene_id
  })
  n_nearest
}

#each array in list_to_filter is interesected with filter_to_apply
filter_lists = function(list_to_filter, filter_to_apply){
  filtered_list = lapply(list_to_filter, function(x){
    intersect(x, filter_to_apply)
  })
  return(filtered_list)
}

correlate_nearby = function(){
  mRNA_nearby_corr = lapply(names(mRNA_nearby_DE)[is_DE], function(x){
    mRNA = mRNA_nearby_DE[[x]]
    dat = norm_counts[c(x, mRNA),]
    dat_cor = cor(t(dat))[1,-1]
    names(dat_cor) = mRNA
    return(dat_cor)
  })
  names(mRNA_nearby_corr) = names(mRNA_nearby_DE)[is_DE]
  return(mRNA_nearby_corr)
}

plot_nearby = function(DE_lists){
  MAX = max(sapply(DE_lists, length)) + 1
  l_mat = matrix(c(rep(1, MAX), 1:MAX + 1), ncol = MAX, byrow = T)
  l_hei = c(1,3)
  n = 1
  mRNA_nearby_plots = pblapply(names(DE_lists), function(name){
    # print(n / length(DE_lists)); n <<- n + 1
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



