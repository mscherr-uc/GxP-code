# Generate the input for mashr (pvalues,SEs)
# based on SCAIP-genetic/mashr_eQTL/mashr_prep.R
# 9/16/25 MS

library(ashr)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

if (length(args)>0){
      dataset <- args[1]
    }

eQTL_dir <- "/rs/rs_grp_gxp/RNAseq_analysis/GxP_20250730/eQTL_mapping/tensor/tensorqtl_output_cis-nominal_SV15_100kb/combined_nominal_results/"

# read in, grab effect estimate and its standard error
m <- fread(paste0(eQTL_dir, dataset, "_all_nominal_pairs.txt.gz"),sep="\t")
colnames(m) <- c("phenotype_id","variant_id", "start_distance", "end_distance", "af", "ma_samples", "ma_count", "pval_nominal", "slope", "slope_se")

# add unique SNP identifier:
m$uniqID <- paste0(m$phenotype_id,"_",m$variant_id)

# calculate lfsr:
ashr=ash(m$slope,m$slope_se)
m$lfsr <- ashr$result$lfsr

# save the slope, SE and pvals in separate files:
slope <- m[,c("uniqID","slope")]
SE <- m[,c("uniqID","slope_se")]  # Use slope_se
pvals <- m[,c("uniqID","pval_nominal")]  # Use pval_nominal

colnames(slope)[2] <- dataset
colnames(SE)[2] <- dataset
colnames(pvals)[2] <- dataset

write.table(slope, paste0("input/", dataset, "_slope.txt"), sep="\t", col.names=T, row.names=F, quote=F)
write.table(SE, paste0("input/", dataset, "_SE.txt"), sep="\t", col.names=T, row.names=F, quote=F)
write.table(pvals, paste0("input/", dataset, "_pvalue.txt"), sep="\t", col.names=T, row.names=F, quote=F)

# save lfsr to use instead of pvalues
lfsr <- m[,c("uniqID","lfsr")]
colnames(lfsr)[2] <- dataset
write.table(lfsr, paste0("input/", dataset, "_lfsr.txt"), sep="\t", col.names=T, row.names=F, quote=F)

# gzip all:
system(paste0("for file in input/", dataset,"*txt ; do gzip $file; done"))

### END
