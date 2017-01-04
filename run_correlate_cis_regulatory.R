source("parse_config.R")
parse_config()
#for input linc, return mRNA within certain distance, 
#check if mRNA are DE and also if their expression correlates with that of the linc.
library(GenomicRanges)
library(pbapply)

setwd(script_dir)
source("parse_gtf.R")
source('fetch_ucsc_image.R')
source("functions_correlate_cis_regulatory.R")

setwd(data_dir)
load_correlation_cis_regulatory( 
                                ref_file = REF_FILE, 
                                lnc_de_pattern = LNC_DE_PATTERN, 
                                mrna_de_pattern = MRNA_DE_PATTERN)

mRNA_nearby = n_closest_features(test_res = all_DE_lincs, 
                                 ref_dict = ensg_ref, 
                                 n_closest = 4, 
                                 max_distance = MAX_DIST, 
                                 strandedness = "either")

#filter nearby mRNA by whether they are DE
mRNA_nearby_DE = lapply(mRNA_nearby, function(x){
  intersect(x, all_DE_mRNA)
})
is_DE = sapply(mRNA_nearby_DE, length) > 0
n_pos = sum(is_DE)
perc_pos = round((n_pos / length(all_DE_lincs)) * 100, 2)
print(paste0(perc_pos, "% of DE lincs near DE mRNA"))

mRNA_nearby_corr = correlate_nearby()

DE_mRNA_nearby_DE = mRNA_nearby_DE[is_DE]

setwd(output_dir)

pdf("DE_lincs_with_DE_mRNA_top20.human.V25.pdf", width = 2.5 * 4, height = 2.5 * 4 / 3)
o = order(sapply(mRNA_nearby_corr, function(x)max((x))), decreasing = T)
plot_nearby(DE_mRNA_nearby_DE[o][1:20])
dev.off()

pdf("DE_lincs_with_DE_mRNA_bottom20.human.V25.pdf", width = 2.5 * 4, height = 2.5 * 4 / 3)
o = order(sapply(mRNA_nearby_corr, function(x)min((x))), decreasing = F)
plot_nearby(DE_mRNA_nearby_DE[o][1:20])
dev.off()

setwd(MAIN_DIR)