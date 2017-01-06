
parse_exon_gtf_and_merge = function(exon_gtf_file, gene_type = "novel_lincRNA"){
  exon_gtf = read.table(exon_gtf_file, sep = "\t", stringsAsFactors = F)
  
  all_attribs = strsplit(exon_gtf[,9], "; ")
  
  get_attrib = function(key) {
    out = lapply(all_attribs, function(x) {
      keep = grepl(key, x)
      str = "NA"
      if (sum(keep) > 0) 
        str = strsplit(x[keep], " ")[[1]][2]
      return(str)
    })
    return(unlist(out))
  }
  
  exon_gtf$gene_id = get_attrib("gene_id")
  
  novel_ref = ensg_ref[0,]
  novel_ref$gene_type = factor(levels = gene_type)
  for(g in unique(exon_gtf$gene_id)){
    k = exon_gtf$gene_id == g
    chrm = exon_gtf[1,1]
    start = min(exon_gtf[k,4])
    end = max(exon_gtf[k,5])
    strand = exon_gtf[1,7]
    novel_ref[g, c(1:3,6) ] = c(g, g, chrm,  strand)
    novel_ref[g, c(4:5) ] = c(start, end)
    novel_ref[g,7] = gene_type
  }
  return(novel_ref)
}

