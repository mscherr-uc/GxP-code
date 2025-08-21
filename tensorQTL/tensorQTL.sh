module load plink/2.0
mkdir -p /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf_plink

#run one of the 2 vcf to plink options
#cluster="C3"
treat="CTRL"
for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt `;do
for chr in {1..22};do
  if [ ! -f /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf_plink/ref.ac1.$cluster.$treat.${chr}.psam ];then
echo running cluster $cluster chr $chr
plink2 --vcf /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf/ref.ac1.$cluster.$treat.filtered.vcf.gz --chr chr${chr} --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf_plink/ref.ac1.$cluster.$treat.${chr}
fi
done
done

  if [ ! -f /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf_plink/ref.ac1.psam ];then
echo running
plink2 --vcf /rs/rs_grp_scaloft/genotypes_liftOver2hg38/ref.ac1.reheader.gsubRI.vcf --out /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/vcf_plink/ref.ac1
fi

#my residuals were not designed for tensorqtl, this reformats it
R
library(data.table)
library(plyr)
library("AnnotationHub")
ah <- AnnotationHub()
if(length(ah["AH98047"]) == 0) {
  edb <- ah[["AH75011"]]
} else {
  edb <- ah[["AH98047"]]
}
geneIDs <- genes(edb) %>%
  as.data.frame() %>% 
  setDT(keep.rownames = "ensembl_gene_id") %>%
  .[, c("ensembl_gene_id","entrezid","symbol","seqnames","start","end","strand","gene_biotype", "description")]
names(geneIDs)[c(1,2,4,8)] <- c("ensgene","entrez","chr","biotype")
geneIDs <- subset(geneIDs, biotype=="protein_coding")
geneIDs <- transform(geneIDs, chr=as.character(chr),strand=as.character(strand))
#this 1bp positioning for bed file is based on examples from qtltools and tensorqtl
geneIDs <- transform(geneIDs, strand_start=ifelse(strand=="+",start,end),strand_end=ifelse(strand=="+",start+1,end+1))
#geneIDs <- transform(geneIDs, strand_start=start,strand_end=start+1)
geneIDs <- subset(geneIDs, chr %in% c(1:22))

base="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/"
method="demux"
project="ALL"
resset=0.2
dimset=50
combatrun="income_PCs_sex_age_and_treats_adjusted"
baseoutFolder=paste0(base,method,"_pseudobulk_ctrl/lessfilt/cell20filt/")
opfn <- paste0(base,method,"_pseudobulk_ctrl/lessfilt/",project,".",resset,".",dimset,".DESeq_countlists_wavefilt.icfilt.RData")
#opfn <- paste0(base,method,"_pseudobulk_ctrl/lessfilt/",project,".",resset,".",dimset,".DESeq_countlists_wavefilt.bticfilt.RData")
load(opfn)
outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/"
if (!file.exists(outFolder)) dir.create(outFolder, showWarnings=F)
treat="CTRL"
method="voom"

lapply(names(counts_ls),function(cluster){
    cat("running ",cluster,treat,"\n")
## load normalized data / or residuals:
  Res <- fread(paste0(outFolder,cluster,".",treat,".residuals_",method,".txt"))
  all_counts_bed <- merge(geneIDs[,c("chr","strand_start","strand_end","ensgene")],Res,by.x="ensgene",by.y="V1")
  all_counts_bed <- all_counts_bed %>% relocate(ensgene,.after =strand_end)  # have to use ensgene as there was multi gene symbols
  #all_counts_bed <- transform(all_counts_bed, chr=paste0("chr",chr))
  names(all_counts_bed)[c(2:4)] <- c("start","end","gene_id")
  #all_counts_bed <- transform(all_counts_bed, chr=paste0("chr",chr))
  fwrite(all_counts_bed, sep='\t', quote=F, row.names=F, col.names=T, file=paste0(outFolder,"phenotypes.",cluster,".",treat,".residuals_",method,".bed"))
})

#treat="CTRL"
#method="voom"
#fwrite(cvs, sep='\t', quote=F, row.names=F, col.names=T, file=paste0(outFolder,"covariates/PCcovariates_",method,".",cluster,".",treat,".txt"))
#for i in $(seq 1 20); do 
#  head -n $(($i+1)) /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/covariates/PCcovariates_${method}.$cluster.$treat.txt > /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/covariates/$cluster.$treat.PC1-$i.covariates_${method}.txt
#done

#add # to header and sort and gzip
module swap gnu9 gnu7/7.3.0
module load bedtools/2.25.0
method="voom"
#ALOFT
treat="CTRL"
for cluster in C0 C1 C10 C11 C2 C3 C4 C5 C6 C7 C8 C9;do
i=`ls -1 /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/residuals/phenotypes.$cluster.$treat.residuals_${method}.bed | grep -v 'sort'`
  echo "running " $i
  #less $i | awk 'NR == 1{print "#"$0;next}; NR > 1 {print $0 | "sortBed -i"}' > ${i%.*}.sort.bed
  less $i | awk 'NR == 1{print "#"$0;next}; NR > 1 {print $0 | "sort -k1,1 -k2,2n "}' > ${i%.*}.sort.bed
  bgzip ${i%.*}.sort.bed && tabix -p bed ${i%.*}.sort.bed.gz
done
module swap gnu7/7.3.0 gnu9

#####################################################
#this is all part of me getting the enivronment right
cd /wsu/home/groups/piquelab/cindy/jaxqtl
pip3 install tensorqtl

  WARNING: The scripts f2py, f2py3 and f2py3.7 are installed in '/wsu/home/fh/fh85/fh8591/.local/bin' which is not on PATH.
  Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.
mlxtend 0.18.0 requires scipy>=1.2.1, but you have scipy 1.1.0 which is incompatible.
multiqc 1.7 requires matplotlib<3.0.0,>=2.1.1, but you have matplotlib 3.0.3 which is incompatible.
mxnet 1.5.0 requires requests<3,>=2.20.0, but you have requests 2.19.1 which is incompatible.

cat >> ~/.bashrc
export PATH="/wsu/home/fh/fh85/fh8591/.local/bin:$PATH"
source ~/.bashrc

pip3 install scipy -U
pip3 install matplotlib -U
pip3 install multiqc -U
pip3 install requests -U

pip3 install tensorqtl
pip3 install importlib-metadata -U

conda create -n tensorqtl_p3.11 python=3.11
conda activate tensorqtl_p3.11
pip3 install tensorqtl
python3 -m tensorqtl --help
#ModuleNotFoundError: No module named 'pandas_plink'
pip3 install pandas_plink

python3 -m tensorqtl --help
Warning: 'rfunc' cannot be imported. R with the 'qvalue' library and the 'rpy2' Python package are needed to compute q-values.
pip3 install 'rpy2<=3.5.12'

conda config --add channels conda-forge
conda config --set channel_priority strict 
conda install r-essentials

conda install -c bioconda bioconductor-qvalue

python3
from rpy2.robjects.packages import importr

OSError: cannot load library '/wsu/el7/groups/piquelab/R/4.3.2/lib64/R/lib/libR.so': libRblas.so: cannot open shared object file: No such file or directory

export LD_LIBRARY_PATH=/wsu/el7/groups/piquelab/R/4.3.2/lib64/R/lib:$LD_LIBRARY_PATH

python3
from rpy2.robjects.packages import importr
#works!

python3 -m tensorqtl --help

conda env export > /wsu/home/groups/piquelab/cindy/tensorqtl.environment.yaml
##################################################
#make and check environment is working
conda env create --name tensorqtl_p3.11_env --file=/wsu/home/groups/piquelab/cindy/tensorqtl.environment.yaml
conda activate tensorqtl_p3.11_env
export LD_LIBRARY_PATH=/wsu/el7/groups/piquelab/R/4.3.2/lib64/R/lib:$LD_LIBRARY_PATH
python3 -m tensorqtl --help | head
##################################################

#interactive test run
data_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis"
out_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/"
mkdir -p ${out_path}
mkdir -p ${data_path}/logs

cluster="C4"
treat="CTRL"
chr=1
PC=2
#pheno="${data_path}/counts/phenotypes.${cluster}.${treat}.sort.bed.gz"
geno="${data_path}/residuals/vcf_plink/ref.ac1.${cluster}.${treat}.${chr}" # prefix for plink triplet files
# prefix for output file
out_prefix="${cluster}_${treat}_${chr}_PC${PC}_tensorqtl"
pheno=${data_path}/residuals/phenotypes.$cluster.$treat.PC1-$PC.residuals_voom_PCregress2step.sort.bed.gz

echo $cluster $treat $chr $PC
echo $geno
echo $pheno
echo $out_prefix

#test
python3 -m tensorqtl ${geno} ${pheno} ${out_prefix} \
    --mode cis \
    --fdr 0.1 \
    -o ${out_path}
#    --chunk_size 'chr' \

#test with non-pc regressed residualscluster="C4"
treat="CTRL"
chr=1
PC=2
for chr in {2..22};do
  echo running chr $chr
geno="${data_path}/residuals/vcf_plink/ref.ac1.${cluster}.${treat}.${chr}" # prefix for plink triplet files
# prefix for output file
out_prefix="${cluster}_${treat}_${chr}_PC${PC}_tensorqtlr"
pheno=${data_path}/residuals/phenotypes.$cluster.$treat.residuals_voom.sort.bed.gz
covar="${data_path}/residuals/covariates/${cluster}.${treat}.PC1-$PC.covariates_voom.txt"

echo $cluster $treat $chr $PC
echo $geno
echo $pheno
echo $covar
echo $out_prefix

python3 -m tensorqtl ${geno} ${pheno} ${out_prefix} \
    --mode cis \
    --covariates ${covar} \
    --fdr 0.1 \
    -o ${out_path}
done

#HPC run
data_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis"
out_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/"
cluster="C0"
treat="CTRL"
PC=2
mkdir -p ${out_path}
for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
#for cluster in C0 C4; do
for chr in {1..22};do
  for PC in {1..10}; do
  njobs=`squeue -u fh8591 -r| wc -l`
  maxjobs=500
  while [ "$njobs" -gt "$maxjobs" ];do 
  echo waiting for jobspace, sleeping ...
  sleep 300 
  njobs=`squeue -u fh8591 -r| wc -l`
  done #end while
      if [ -f "${out_path}/${cluster}_${treat}_${chr}_PC${PC}_tensorqtlr.cis_qtl.txt.gz" ] ; then
      echo " Output already exists. Skipping..."
    else
    echo submitting cluster $cluster chr $chr PC $PC
    sbatch --export=cluster="${cluster}",treat="${treat}",chr="${chr}",PC="$PC",data_path="$data_path",out_path="$out_path" ${data_path}/tensorQTL/src/run_tensorqtl.sh 
    sleep 1
  fi
  done #PC end
done #chr end
done #cluster end

#check output has all cols (19)
for i in tensorQTL/output/*.txt.gz; do  
  numcol=`less $i | head -n1 | awk '{print NF}'`; 
  echo $i has $numcol; done

#check what has run
treat="CTRL"
rm ${out_path}/not_completed.txt ${out_path}/completed.txt
rm ${out_path}/not_completed.txt ${out_path}/not_completed.txt
for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt `;do
for chr in {1..22};do
  for PC in {1..10}; do
    out_prefix="${cluster}_${treat}_${chr}_PC${PC}_tensorqtlr" # prefix for output file
    out="${out_path}/${out_prefix}.cis_qtl.txt.gz"
    #echo $out
    if [ -f "$out" ] ; then
      #echo " Output already exists."
      echo "${cluster}_${treat}_${chr}_${PC}" >> ${out_path}/completed.txt
    else
      echo $out not run
      echo "${cluster}_${treat}_${chr}_${PC}" >> ${out_path}/not_completed.txt
    fi
done
    #sleep 1
done
done
#25812 - 10305
less ${out_path}/completed.txt | wc -l #440
less ${out_path}/not_completed.txt | wc -l #0


#################### INTERACTION
#The interaction term is a tab-delimited text file or dataframe mapping sample ID to interaction value(s) 
#(if multiple interactions are used, the file must include a header with variable names).
#each row is sample id, each col is variable

#make each int file in R
R
library(data.table)
library(plyr);library(dplyr)

data_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/"
outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/"
treat="CTRL"
pheno=${data_path}/residuals/phenotypes.$cluster.$treat.residuals_voom.sort.bed.gz
covar="${data_path}/residuals/covariates/${cluster}.${treat}.PC1-$PC.covariates_voom.txt"
base="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/"
method="demux"
project="ALL"
resset=0.2
dimset=50
clusters <- fread("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt")
variables <- c("pedu", "pincme", 
                "psesl", 
                "pnsi", "cddstf",
                "cdres", "cpwm",  "cpeqcm", "criskf", "cditsm", 
                "IL5_co", "IL13_co", "IFNG_co", "IL5_hc", "IL13_hc", "IFNG_hc",
                "baso_av", "eosi_av", "lymp_av", "mono_av", "neut_av",
                "aBPFAM", "aBPFPM", "FEVPP", "FVCPP", "FFPP",
                "cdatot", "csasg", "ctasfq", "ctasev", 
                "csnuma",  "cssdh", "csslpq", 
                "genPC1", "genPC2", "genPC3",  "Sex", "cage1", "ceth1", "cwght1", "chght1", "csex1",
                "cgpd5", "cgpd", "cbpd"
                )
 
variable_names <- c("Parental Education", "Parental Income",
                "Subjective SES", 
                "Neighborhood Stress", "Self-disclosure", 
                "Perceived responsiveness", "Parental Warmth", "YR Parent Child Conflict",  "Risky Family", "Youth Depression", 
                "Stimulated IL-5", "Stimulated IL-13", "Stimulated IFNy", "Stimulated IL-5 cortisol trt", "Stimulated IL-13 cortisol trt", "Stimulated IFNy cortisol trt",    
                "Basophils", "Esosinophils", "Lymphocytes", "Monocytes", "Neutrophils",
                "Peak flow AM", "Peak flow PM", "FEV1 Percent Predicted", "FVC Precent Predicted", "FEV1/FVC Precent Predicted", 
                "DD Asthma symptoms", "Nightly Asthma", "Asthma Frequency", "Asthma Severity", 
                "Nightly Awakenings",  "Sleep Duration", "Sleep quality",
                "New Genotype PC1", "New Genotype PC2", "New Genotype PC3","Sex_alph", "Age", "Ethnicity", "Weight", "Height", "sex",
                "Female Menarche Status", "Female Puberty Score", "Male puberty score"
                )
variables_df <- data.frame(variable=variables, description=variable_names)
variables_dftorun <- subset(variables_df, !variable %in% c("csex1","Sex","cage1","Wave","genPC1","genPC2","genPC3"))
fwrite(variables_dftorun, file=paste0(data_path,"variables_dftorun.txt"),sep="\t",col.names=F,row.names=F, quote=F)

opfn <- paste0(base,method,"_pseudobulk_ctrl/lessfilt/",project,".",resset,".",dimset,".DESeq_countlists_wavefilt.icfilt.RData")
load(opfn)

filenames_res <- list.files(paste0(data_path,"residuals/")) #file list from directory
filenames_cov <- list.files(paste0(data_path,"residuals/covariates/")) #file list from directory
best_df <- fread(file=paste0(outFolder,"bestPCs_table.txt"))

for(c in clusters$cluster){
  for(v in variables_dftorun$variable){
    best.PCs <- subset(best_df, cluster==c)$PCs
    pheno <- fread(paste0(data_path,"residuals/",filenames_res[grep(paste0("phenotypes.",c,".",treat,".residuals_voom.sort.bed.gz$"),filenames_res)]))
    covar <- fread(paste0(data_path,"residuals/covariates/",filenames_cov[grep(paste0(c,".",treat,".PC1-",best.PCs,".covariates_voom.txt"),filenames_cov)]))
    pc_signif_pairs <- subset(fread(paste0(outFolder,"results/",c,".",treat,".best_", best.PCs, ".GEPCs.txt")),qval<0.1)
    phenosub <- subset(pheno, gene_id %in% pc_signif_pairs$phenotype_id)
    fwrite(phenosub,paste0(data_path,"residuals/","phenotypes.",c,".",treat,".residuals_voom.eQTLonly.sort.bed.gz"),sep="\t",col.names=T,row.names=F, quote=F)

    phenoind <- colnames(phenosub[,-c(1:4)])
    covarind <- colnames(covar)[-1]
    cluster_metadata_sce <- metadata_ls[[c]]
    cluster_metadata <- data.frame(cluster_metadata_sce)
    cluster_metadata_t <- subset(cluster_metadata, treats==treat & Sample_ID %in% phenoind)
    covvar <- as.data.frame(cluster_metadata_t[,v])
    notna <- complete.cases(covvar)
    inds <- cluster_metadata_t$Sample_ID[notna]

    covvar <- data.frame(Sample_ID=inds,variable=covvar[notna, ])
    covvarv <- covvar[,-1, drop=F]
    rownames(covvarv) <- inds #NEED ROWNAMES
    names(covvarv)[1] <- v
    fwrite(covvarv,file=paste0(data_path,"residuals/covariates/",c,".",treat,".",v,".for_tensorqtl_int.txt"),sep="\t",col.names=T,row.names=T, quote=F)
    phenosubcov <- subset(phenosub, select = c(colnames(phenosub)[1:4],inds))
    fwrite(phenosubcov,paste0(data_path,"residuals/","phenotypes.",c,".",treat,".",v,".residuals_voom.eQTLonly.sort.bed.gz"),sep="\t",col.names=T,row.names=F, quote=F)
    covarsubcov <- subset(covar, select = c(colnames(covar)[1],inds))
    fwrite(covarsubcov,paste0(data_path,"residuals/covariates/",c,".",treat,".",v,".PC1-",best.PCs,".covariates_voom.txt"),sep="\t",col.names=T,row.names=F, quote=F)

  }
}


conda init bash
source ~/.bashrc
module swap gnu7/7.3.0 gnu9

conda activate tensorqtl_p3.11_env
export LD_LIBRARY_PATH=/wsu/el7/groups/piquelab/R/4.3.2/lib64/R/lib:$LD_LIBRARY_PATH
module load R/4.3.2 #old: module load r/4.2.0

data_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis"
outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output"
mkdir -p ${outFolder}/interaction
out_path="${outFolder}/interaction"
treat="CTRL"
chr=1
var="cditsm"
cluster="C9"
PC=`grep $cluster ${outFolder}/bestPCs_table.txt | cut -f2`
geno="${data_path}/residuals/vcf_plink/ref.ac1.${cluster}.${treat}.${chr}" # prefix for plink triplet files
out_prefix="${cluster}_${treat}_${chr}_${var}_tensorqtlint" # prefix for output file
pheno=${data_path}/residuals/phenotypes.$cluster.$treat.${var}.residuals_voom.eQTLonly.sort.bed.gz
covar="${data_path}/residuals/covariates/${cluster}.${treat}.${var}.PC1-$PC.covariates_voom.txt"
int="${data_path}/residuals/covariates/${cluster}.${treat}.${var}.for_tensorqtl_int.txt"

if [ -s "${out_path}/${out_prefix}.cis_qtl_top_assoc.txt.gz" ] ; then
    echo "[$(date)] Output already exists. Skipping..."
    exit 0
fi

### -------- LOGGING -------- ###
echo "[$(date)] Starting tensorQTL interaction on node: $(hostname)"
echo "Variable=$var, CLUSTER=$cluster, TREAT=$treat, PC=$PC"
echo "pheno file: $pheno"
echo "geno file: $geno"
echo "covar file: $covar"
echo "interactions file: $int"
echo "Output file: ${out_path}/${out_prefix}"
echo "settings: mode=cis_nominal;best_only;fdr=0.1"

#testing using one single plink not chr for gxe
geno="${data_path}/residuals/vcf_plink/ref.ac1" # prefix for plink triplet files


python3 -m tensorqtl ${geno} ${pheno} ${out_prefix} \
    --covariates ${covar} \
    --interaction ${int} \
    --best_only \
    --fdr 0.1 \
    --mode cis_nominal \
    -o ${out_path}
 
echo end on "[$(date)]"

#HPC test
data_path="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis"
outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output"
mkdir -p ${outFolder}/interaction
out_path="${outFolder}/interaction"

#cluster="C0"
treat="CTRL"

for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt`;do
PC=`grep $cluster ${outFolder}/bestPCs_table.txt | cut -f2`
for var in `cut -f1 ${data_path}/variables_dftorun.txt`; do
#for chr in {1..22};do
  chr=NA #can run all chr at once
  njobs=`squeue -u fh8591 -r| wc -l`
  maxjobs=500
  while [ "$njobs" -gt "$maxjobs" ];do 
  echo waiting for jobspace, sleeping ...
  sleep 300 
  njobs=`squeue -u fh8591 -r| wc -l`
  done #end while
      if [ -s "${out_path}/${cluster}_${treat}_${chr}_${var}_tensorqtlint.cis_qtl_top_assoc.txt.gz" ] ; then
      echo " Output already exists. Skipping..."
    else
    echo submitting cluster $cluster chr $chr variable $var 
    sbatch --export=cluster="${cluster}",treat="${treat}",chr="${chr}",var="$var",PC="${PC}",data_path="$data_path",out_path="$out_path" ${data_path}/tensorQTL/src/run_tensorqtl_int.sh 
    sleep 1
  fi
#done #chr end
done #var end
done #cluster end
mkdir ${out_path}/logs
mv ${out_path}/*.log ${out_path}/logs/

#combine
outFolder="/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/interaction/"
clusters <- fread("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt")

myDir <- outFolder #directory to load from
filenames <- list.files(myDir) #file list from directory
filenames <- filenames[grep("tensorqtlint.cis_qtl_top_assoc.txt.gz", filenames)] #pick specific files from list
int <- ldply(lapply(clusters$cluster,function(cluster){
  vl <- ldply(lapply(variables_dftorun$variable,function(v){
    cat("running",cluster,v,"\n")
    filenamesc <- filenames[grepl(paste0(cluster,".",treat,".",v), filenames)] #pick specific files from list
    if(length(filenamesc)>0){
    data_names <- gsub("_tensorqtlint.cis_qtl_top_assoc.txt.gz", "", filenamesc) #remove file ending
    shortnames <- gsub("[.].*", "", data_names)
    for(i in 1:length(filenamesc)) assign(data_names[i], fread(file.path(myDir, filenamesc[i]),header = T)[,analysis:=shortnames[i]]) #read in specific files and set the df object names. can dro0p unwanted columns
    #combining data
    all_chrs <- lapply(data_names, function(x) get(x)) #grab data from list of df names
    names(all_chrs) <- data_names
    all_chrsdf <- ldply(all_chrs, data.frame)
    singelchrrun <- all_chrsdf[grep(".NA",all_chrsdf$.id),] #after testing decided to run all chr together for interaction test
    #table(singelchrrun$analysis,singelchrrun$pval_adj_bh<0.1)
    #pval_emt = pval_gi * tests_emt and pval_adj_bh=p.adjust(pval_emt,method="BH")
    singelchrrun <- transform(singelchrrun, genotype_pval_emt=pval_g*tests_emt, variable_pval_emt=pval_i*tests_emt)
    singelchrrun <- transform(singelchrrun, genotype_padj=p.adjust(genotype_pval_emt,method="BH"), variable_padj=p.adjust(variable_pval_emt,method="BH"),
      cluster=cluster,treat=treat,variable=v)
    singelchrrun <- singelchrrun %>% relocate(c(cluster,treat,variable),.after =.id)  # have to use ensgene as there was multi gene symbols
    #table(singelchrrun$analysis,singelchrrun$genotype_padj<0.1)
    #table(singelchrrun$analysis,singelchrrun$variable_padj<0.1)
    return(singelchrrun)
    } 
  }),data.frame)
}),data.frame)
#counts_ls[sapply(counts_ls, is.null)] <- NULL

fwrite(int, paste0(outFolder,treat,".GxE_abundance_perclus.txt"), sep='\t', quote=F, row.names=F)

tested <- as.data.frame(table(int$cluster,int$variable))
intsig <- subset(int,pval_adj_bh<0.1 )
tablesig <- as.data.frame(table(intsig$cluster,intsig$variable))
outtable <- merge(tested,tablesig,by=c("Var1","Var2"))
colnames(outtable) <- c("cluster","variable","tested","int_FDR10")
outtabledf <- merge(outtable,variables_dftorun)
fwrite(outtabledf[,c("cluster","variable","description","tested","int_FDR10")], paste0(outFolder,treat,".GxE_summary_perclus.txt"), sep='\t', quote=F, row.names=F)

outtablem <- melt(outtabledf)
names(outtablem)[4] <- "count"
p <- ggplot(outtablem, aes(fill=count, y=value, x=cluster)) + 
    geom_bar(position="stack", stat="identity")+
    facet_wrap(.~description,ncol=4)
    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/figures/tensorint_GxE_summary_bar.png")
png(width = 12, height = 10, file=figfn, pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()

intsig <- transform(intsig, direction=if_else(b_gi<0,"down","up"))
intsigc <- plyr::count(intsig, c("cluster","variable","direction"))
intsigcdf <- merge(intsigc, intsig)
intsigcdf <- merge(intsigcdf,variables_dftorun)
intsigcdf <- transform(intsigcdf, DEG_direction= if_else(direction=="down",-(freq),freq))
p <- ggplot(unique(intsigcdf[,c("description","cluster","DEG_direction","direction")]), aes(fill=direction, y=DEG_direction, x=cluster)) + 
    geom_bar(position="stack", stat="identity")+
    facet_wrap(.~description,ncol=4,scales="free_y")+
    labs(y="# Interaction eGenes")
    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/figures/tensorint_GxE_summary_bar_degonly.png")
png(width = 13, height = 10, file=figfn, pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p)
dev.off()

library(ggrastr)
i="CTRL"

lapply(split(int,int$variable),function(v){
    
    v <- transform(v, cluster=as.factor(cluster))
    v <- v %>%
   group_by(cluster)%>%
   arrange(pval_gi) %>%
   mutate(observed=-log10(pval_gi), expected=-log10(ppoints(length(pval_gi))))

p0 <- ggplot(v, aes(x=expected, y=observed, color=cluster))+
    geom_point()+
    geom_abline(color="grey")+
    xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
    ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
    ggtitle(paste0(unique(v$variable)))+
    theme_bw()
    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/figures/tensorint_GxE_abundance_pvalues_",unique(v$variable),".pcl_qqplot.png")
png(width = 8, height = 8, file=figfn, pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p0)
dev.off()
p0 <- ggplot(v, aes(x=expected, y=observed, color=cluster))+
    geom_point()+
    geom_abline(color="grey")+
    #scale_color_manual(values=c("C0"="#F8766D", "C1"="#D39200", "C2"="#93AA00", "C3"="#00BA38",
    #    "C4"="#00C19F", "C5"="#00B9E3", "C6"="#619CFF", "C7"="#DB72FB", "C8"="","C9"="#FF61C3"),
    #    guide=guide_legend(override.aes=list(size=3)))+
    facet_grid(.~cluster, scales="free_y")+
    xlab(bquote("Expected"~-log[10]~"("~italic(p)~")"))+
    ylab(bquote("observed"~-log[10]~"("~italic(p)~")"))+
    ggtitle(paste0(unique(v$variable)))+
    theme_bw()

    figfn <- paste0("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/tensorQTL/output/figures/tensorint_GxE_abundance_pvalues_",unique(v$variable),".facetcluster_pcl_qqplot.png")
png(width = 12, height = 5, file=figfn, pointsize=12, 
      bg = "transparent", canvas = "white", units = "in", res = 600)
print(p0)
dev.off()

})

####################
#on reflection, I don't htink it will work for cluster or treatment since that is multiple readings for the same individual
R
library(data.table)
library(plyr);library(dplyr)
library(rlist)

args <- c("/rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/","/rs/rs_grp_scaloft/scALOFT_2024/covariates/ALOFT_covariate_issues_fixed_uniq-n265_psesl-a2_fixed_12-19-2024.txt","ALL","demux","/rs/rs_grp_scaloft/scALOFT_2024/covariates/scALOFT_samples_batch2.txt",0.2) 
base <- args[1]
cov_file=fread(args[2]) #this is the psych cov file
project=args[3]
method=args[4]
sample_batch <- args[5]
resset <- args[6]
dimset=50
outFolder=paste0(base,method,"_pseudobulk_ctrl/")
#for resolution 0.2 and including V2 chem
outFolder=paste0(outFolder,"lessfilt/")
opfn <- paste0(outFolder,project,".",resset,".",dimset,".DESeq_countlists_wavefilt.icfilt.RData")
load(opfn)
treatments=c("CTRL","LPS","LPS-DEX","PHA","PHA-DEX")
treatmentsfirst=c("CTRL","PHA")
contrastdf <- data.frame(control=c("CTRL","CTRL","PHA"),treatment=c("LPS","PHA","PHA-DEX")) #no lps vs lps-dex as too few ind

metadata <- list.rbind(metadata_ls)
lapply(1:length(contrastdf$control),function(x){
    con=contrastdf[x,]
    contrast=paste0(con$treatment,"_vs_",con$control)
    meta_t <- subset(metadata, treats %in% c(con$control,con$treatment))

c <- plyr::count(unique(meta_t[,c("Sample_ID","treats")]),c("Sample_ID"))

############################################################
for cluster in `awk 'NR>1{print $1}' /rs/rs_grp_scaloft/scALOFT_2024/cindy_analysis/demux_pseudobulk_ctrl/lessfilt/ALL.0.2.50.cluster_celltype.txt | grep -v "C4"`;do
for chr in {1..22};do
  for PC in {1..20}; do

python3 -m tensorqtl ${geno} ${pheno} ${out} \
    --mode cis \
    --chunk_size 'chr' \
    --fdr 0.1
    -o ${out_path}

    --covariates ${covar} \
--chunk_size CHUNK_SIZE
                        For cis-QTL mapping, load genotypes into CPU memory in chunks of chunk_size variants, or by chromosome if chunk_size is
                        'chr'.
  --fdr FDR             FDR for cis-QTLs
-o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory

