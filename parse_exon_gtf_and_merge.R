exon_gtf = read.table("C:/Users/jrboyd/Downloads/Galaxy156-[hMSC_BMP2_novel_lncRNAs_gtf_no_22.txt].gff", sep = "\t", stringsAsFactors = F)


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
novel_ref$gene_type = factor(levels = "novel_lincRNA")
for(g in unique(exon_gtf$gene_id)){
  k = exon_gtf$gene_id == g
  chrm = exon_gtf[1,1]
  start = min(exon_gtf[k,4])
  end = max(exon_gtf[k,5])
  strand = exon_gtf[1,7]
  novel_ref[g, ] = c(g, g, chrm, start, end, strand, "novel_lincRNA")
}
