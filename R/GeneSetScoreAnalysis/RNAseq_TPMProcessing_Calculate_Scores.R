############################################################
# PROCESS INSPIRE RNASEQ DATA 
# TPM Quantile Normalize -> 
# Remove QC failed samples ->
# log2 scale ->
# Combat batch correction ->
# Annotate gene names/ collapse duplicates ->
# Calculate IFN/Ayers/CYT scores ->
# ssGSEA for ESTIMATE and immune gene sets ->
# write matrices for downstream data analyses
# April 9, 2019 Cindy Yang
# #########################################################
LOCAL="TRUE"
if (LOCAL){
  path.head = "/Volumes/Samwise"
}else{
  path.head = "/mnt/work1/users/pughlab"
  .libPaths(c(paste(path.head,"bin/r_lib/library/3.2.2", sep = "/"), .libPaths()));
}

library("RColorBrewer")
library ("pheatmap")
library("sva")
library("ggfortify")

workdir <- paste( path.head, "projects/INSPIRE/Review/RNAseq/data", sep = "/")
outdir <- paste( path.head, "projects/INSPIRE/Review/RNAseq/results", sep = "/")
setwd(workdir)
##################################################
# Quantile normalize data
# module load perl
# Normalize separately BASELINE and TREATED
# quan_file="/mnt/work1/users/pughlab/projects/INSPIRE/Review/RNAseq/data/INSPIRE_RNASEQ_12032019CY_RSEM_tpm_renamed_BASELINE.txt"
# outfile="/mnt/work1/users/pughlab/projects/INSPIRE/Review/RNAseq/data/INSPIRE_RNASEQ_12032019CY_RSEM_tpm_renamed_BASELINE_QuantileNormalized.txt.txt"
# perl ~/pughlab/projects/INSPIRE/Review/RNAseq/ubu/perl/quartile_norm.pl -c -1 -q 75 -t 1000 -s 1 -o $outfile $quan_file

# quan_file="/mnt/work1/users/pughlab/projects/INSPIRE/Review/RNAseq/data/INSPIRE_RNASEQ_12032019CY_RSEM_tpm_renamed_TREATED.txt"
# outfile="/mnt/work1/users/pughlab/projects/INSPIRE/Review/RNAseq/data/INSPIRE_RNASEQ_12032019CY_RSEM_tpm_renamed_TREATED_QuantileNormalized.txt.txt"
# perl ~/pughlab/projects/INSPIRE/Review/RNAseq/ubu/perl/quartile_norm.pl -c -1 -q 75 -t 1000 -s 1 -o $outfile $quan_file

# Merge BASELINE and TREATED samples to normalize together
# quan_file="/mnt/work1/users/pughlab/projects/INSPIRE/Review/RNAseq/data/INSPIRE_RNASEQ_12032019CY_RSEM_tpm_renamed_BASELINE.txt"
# quan_file1="/mnt/work1/users/pughlab/projects/INSPIRE/Review/RNAseq/data/INSPIRE_RNASEQ_12032019CY_RSEM_tpm_renamed_TREATED.txt"
# quan_file2="/mnt/work1/users/pughlab/projects/INSPIRE/Review/RNAseq/data/INSPIRE_RNASEQ_12032019CY_RSEM_tpm_renamed_BASELINE_TREATED.txt"
# paste $quan_file <(cut -f2- $quan_file1) > $quan_file2
# outfile="/mnt/work1/users/pughlab/projects/INSPIRE/Review/RNAseq/data/INSPIRE_RNASEQ_12032019CY_RSEM_tpm_renamed_BASELINE_TREATED_QuantileNormalized.txt"
# perl ~/pughlab/projects/INSPIRE/Review/RNAseq/ubu/perl/quartile_norm.pl -c -1 -q 75 -t 1000 -s 1 -o $outfile $quan_file2
############
# FUNCTION #
############
RNAseqNormCombat <- function(workdir, outdir, suffix, group, log.v = TRUE, plot.stat = FALSE) {
  unit = "tpm"
  raw <- paste(workdir,"/INSPIRE_RNASEQ_12032019CY_RSEM_tpm_renamed_", group,".txt", sep = "")
  normalized <- paste(workdir,"/INSPIRE_RNASEQ_12032019CY_RSEM_tpm_renamed_", group,suffix, sep = "")
  inventory <- paste(workdir,"/RNASeQCMetricsInventory_032019.csv", sep = "")
  
  raw.df <- read.csv(raw, sep = "\t", header = T, stringsAsFactors = F, check.names = F)
  normalized.df <- read.csv(normalized, sep = "\t", header = T, stringsAsFactors = F, check.names = F)
  inventory.df <- read.csv(inventory, header = T,  stringsAsFactors = F, row.names = 1 ,check.names = F)
  
  ############
  # Visualize raw and quantile-norm data PCA to identify outliers
  ###########
  # log scale applied for better data visualization
  pcaraw <- prcomp(t(log2(raw.df[,-1]+1)))
  pcanorm <- prcomp(t(log2(normalized.df[,-1]+1)))
  plotinv <- inventory.df[row.names(pcaraw$x),]
  colnames(plotinv)[colnames(plotinv) %in% "QC Flags"] <- "QC.Flags"
  
  if (plot.stat){
    # raw count data
    # https://www.rdocumentation.org/packages/ggfortify/versions/0.4.5/topics/autoplot.pca_common
    pdf(paste(outdir, "/","PCA.raw.", unit, suffix, ".pdf", sep = ""))
    autoplot(pcaraw, data = plotinv, colour = "Centre", x = 1, y = 2)
    autoplot(pcaraw, data = plotinv, colour = "Cohort", x = 1, y = 2)
    autoplot(pcaraw, data = plotinv, colour = "QC.Flags", x = 1, y = 2)
    dev.off()
    
    # quantile normalized data
    pdf(paste(outdir, "/","PCA.qnorm.", unit, suffix, ".pdf", sep = ""))
    autoplot(pcanorm, data = plotinv, colour = "Centre", x = 1, y = 2)
    autoplot(pcanorm, data = plotinv, colour = "Cohort", x = 1, y = 2)
    autoplot(pcanorm, data = plotinv, colour = "QC.Flags", x = 1, y = 2, label = TRUE, label.size = 2.5)
    # kmeans on PCs
    autoplot(kmeans(pcanorm$x,2), data = pcanorm$x, label = TRUE, label.size = 2.5, frame = TRUE)
    # pam on PCs
    #library(cluster)
    #autoplot(pam(pcanorm$x,3), label = TRUE, label.size = 2.5, frame = TRUE, frame.type = "norm")
    dev.off()
  }
  
  ############
  # Remove outlier samples prior to ComBat batch correction
  # Version 1. with High stranded bias and <25k genes detected (n = 109)
  v1.normalized <- normalized.df[,c("ENSEMBLID", row.names(plotinv[plotinv$QC.Flags %in% c(""),]))] 
  v1.normalized[,-1] <- log2(v1.normalized[,-1]+1)
  v1.plotinv <- plotinv[plotinv$QC.Flags %in% c(""),]                                                                  
  
  #############
  # Apply ComBat 
  # Combat function
  
  applyCombat <- function (dat, inv, batch, parametric = TRUE){
    coldata <- data.frame(
      "batch" = inv[,batch], 
      row.names = row.names(inv)
    )
    print(coldata)
    modcombat <- model.matrix(~1, data=coldata)
    batch <- coldata[,"batch"]
    # removing all genes that are 0 across all samples
    n.genes.no.expr <- sum(rowSums(v1.normalized[,-1])==0)
    print (paste("Number of transcripts not expressed across all samples:", n.genes.no.expr))
    genes.removed.ids <- dat[rowSums(dat[,-1])==0,1]
    dat <- dat[rowSums(dat[,-1])>0,]
    combat_edata <- ComBat(dat=dat[,-1], batch=batch, mod=modcombat, 
                           par.prior=parametric, prior.plots=FALSE) 
    row.names(combat_edata) <- dat[,1]
    return(list(combat_edata = combat_edata,
                genes.removed = genes.removed.ids))
  }
  # Apply to both datasets
  v1.normalized.batchcorr <- applyCombat (dat = v1.normalized, inv = v1.plotinv, batch = "Centre")
  
  # Visualize any differences in PCs
  if (plot.stat){
    pcav1 <- prcomp(t(v1.normalized.batchcorr$combat_edata))
    pdf(paste(outdir, "/","PCA.ComBat.", unit, suffix, ".pdf", sep = ""))
    autoplot(pcav1, data = v1.plotinv, colour = "Centre", x = 1, y = 2)
    autoplot(pcav1, data = v1.plotinv, colour = "Cohort", x = 1, y = 2)
    autoplot(pcav1, data = v1.plotinv, colour = "QC.Flags", x = 1, y = 2)
    autoplot(pcav1, data = v1.plotinv, colour = "Timepoint", x = 1, y = 2)
    autoplot(pcav1, data = v1.plotinv, colour = "Cohort", x = 1, y = 2, label = T, label_size = 2.5)
    dev.off()
    
    # Plot sample distributions
    geneids <- row.names(v1.normalized.batchcorr$combat_edata)
    sampleids <- colnames(v1.normalized.batchcorr$combat_edata)
    
    raw.dat <- raw.df[,sampleids]
    row.names(raw.dat) <- raw.df[,1]
    raw.dat <- raw.dat[geneids,]
    
    quant.norm.dat <- normalized.df[,sampleids]
    row.names(quant.norm.dat) <- normalized.df[,1]
    quant.norm.dat <- quant.norm.dat[geneids,]
    
    combat.dat <- v1.normalized.batchcorr$combat_edata
    
    if (log.v){
      raw.dat <- log2(raw.dat+1)
      quant.norm.dat <- log2(quant.norm.dat +1)
      combat.dat[combat.dat<0] <- 0
      combat.dat <- log2(combat.dat +1)
    }
    
    pdf (paste(outdir,"/","Expr.distribution.profiles.comparison.normalization.",unit,suffix,".pdf", sep = ""), width = 14, height = 6)
    par(mfrow=c(1,3))
    plot(density(as.numeric(raw.dat[,1])), ylim = c(0, 3), xlim = c(-2, 20), main = "Uncorrected")
    for (i in 2:ncol(raw.dat)){
      lines(density(raw.dat[,i]), col = rainbow (100)[i], add = TRUE)}
    # Quantile normalized log2 TPM + 1
    plot(density(quant.norm.dat[,1]), ylim = c(0, 0.6), xlim = c(-2, 25), main = "Quantile normalized")
    for (i in 2:ncol(quant.norm.dat)){
      lines(density(quant.norm.dat[,i]), col = rainbow (100)[i], add = TRUE)}
    # Combat corrected quantile normalized log2 TPM + 1
    plot(density(combat.dat[,1]), ylim = c(0, 0.50), xlim = c(-2, 25), main = "Combat corrected")
    for (i in 2:ncol(combat.dat)){
      lines(density(combat.dat[,i]), col = rainbow (100)[i], add = TRUE)}
    dev.off()
  }

  return (v1.normalized.batchcorr)
}

annotateGeneNames <- function(dat, rnames = "hugo_names", merge.dup = TRUE ,genome = "hg38"){
  
  if (genome == "hg38"){
    gencode.annot <- paste(path.head,"references/gencode/GRCh38/gencode.v26.ensg_annotation.txt", sep = "/")
  }else{
    print ("This function currently only supports hg38 annotations")
  }
  
  gencode.annot <- read.delim(gencode.annot, sep = "\t", header = F, stringsAsFactors = F, check.names = F)[,1:8]
  gencode.annot$stable_ids <- sapply(gencode.annot$V1, function(x) return(unlist(strsplit(x = x, split = ".", fixed = T))[1]))
  
  # create lookup table for gene stable id and gene names
  dat.names <- data.frame(
    stable_ids = sapply(row.names(dat), function(x) return(unlist(strsplit(x = x, split = ".", fixed = T))[1])),
    hugo_names = gencode.annot$V8[match (sapply(row.names(dat), function(x) return(unlist(strsplit(x = x, split = ".", fixed = T))[1])),
                                         gencode.annot$stable_ids)],
    row.names = row.names(dat)
  )
  dat.names$category <- gencode.annot$V7[match(dat.names$stable_ids, gencode.annot$stable_ids)]
  
  print(head(dat.names))
  
  # Take care of all duplicated data (return either stable_id or hugo_symbol as row.names)
  # TODO: return stats for number of duplicates collapsed and what they are?
  if (merge.dup){
    dat.names$duplicate <- duplicated(dat.names[,rnames])
    if (sum(dat.names$duplicate)>0){
      collapseDups <- function (genename, dat, dat.names,rnames) {
        dups <- colSums(dat[row.names(dat.names[dat.names[,rnames] %in% genename,]),])
      }
      new.means <- do.call("rbind", lapply(unique(dat.names[dat.names$duplicate, rnames]), collapseDups, dat = dat, dat.names = dat.names))
      row.names(new.means) <- unique(dat.names[dat.names$duplicate, rnames])
      
      new.dat <- dat[!dat.names[,rnames] %in% row.names(new.means),]
      row.names(new.dat) <- dat.names[row.names(new.dat), rnames]
      new.dat <- rbind(new.dat,new.means)
      rm (new.means)
    }else{
      new.dat <- dat
      row.names(new.dat) <- dat.names[row.names(new.dat), rnames]
    }
    
 }else {
    print ("Currently only duplicate merging by sums is implemented.")
  }
  
  return(list (expr=new.dat, gene.annotation = dat.names))
}

########
# MAIN # 
########
# STEP 1: COMBAT SEQUENCING BATCH CORRECTION 
##############
group <- "BASELINE_TREATED"
suffix <- "_QuantileNormalized.txt"

expr.files <- c("Combat.corrected.log2quantile.norm.TPM_BASELINE_TREATED.txt",
                "Combat.corrected.log2quantile.norm.TPM_BASELINE_TREATED.HugoGenes_deduped.txt",
                "Combat.corrected.log2quantile.norm.TPM_BASELINE_TREATED.EnsemblIDs_deduped.txt",
                "Unlogged.Combat.corrected.log2quantile.norm.TPM_BASELINE_TREATED.HugoGenes_deduped.txt",
                "Unlogged.Combat.corrected.log2quantile.norm.TPM_BASELINE_TREATED.EnsemblIDs_deduped.txt")
expr.files <- paste(workdir, expr.files , sep = "/")

combat.expr <- RNAseqNormCombat(workdir, outdir, suffix, group, log.v = TRUE, plot.stat = FALSE)
# write original data files with rownames annotated with ENSEMBL ids
write.table(combat.expr$combat_edata, expr.files[1], sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

##############
# STEP 2: ANNOTATE expr matrix with gene names
##############
# write matrix of log2 transformed values annotated with hugo gene names (for visualization and score calculations)
# write matrix of unlogged values annotated with hugo gene names (for CIBERSORT)
new.dat <- annotateGeneNames(dat= combat.expr$combat_edata, rnames = "hugo_names", merge.dup = TRUE, genome = "hg38")
write.table(new.dat$expr, expr.files[2], sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(2^(new.dat$expr), expr.files[4], sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

# write matrix of log2 transformed values annotated with ensembl stable ids (for ssGSEA)
# write matrix of unlogged values annotated with ensembl stable ids
new.dat.ensembl <- annotateGeneNames(dat= combat.expr$combat_edata, rnames = "stable_ids", merge.dup = TRUE, genome = "hg38")
write.table(new.dat.ensembl$expr, expr.files[3], sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
write.table(2^(new.dat.ensembl$expr), expr.files[5], sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

##############
# STEP 3: CALCULATE SCORES CYT, IFNg, Ayers
##############
# Using a limited gene set of immune related genes
# Genelists# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5531419/
# Methods from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5531419/
# After performance of quantile normalization, a log10 transformation was applied, 
# and signature scores were calculated by averaging of the included genes for the IFN-γ (6-gene) and expanded immune (18-gene) signatures
AYERS.18 <- c("CXCR6", "TIGIT", "CD27", "CD274", "PDCD1LG2", "LAG3", "NKG7", "PSMB10", "CMKLR1", "CD8A", 
                    "IDO1", "CCL5", "CXCL9", "HLA-DQA1", "CD276", "HLA-DRB1", "STAT1", "HLA-E")
RIBAS.IFN6 <- c("IDO1", "CXCL10", "CXCL9", "HLA-DRA", "STAT1", "IFNG")
ROONEY.CYT <- c("GZMA", "PRF1")
#ayers.10.genes <-c("IFNG", "STAT1", "CCR5", "CXCL9", "CXCL10", "CXCL11", "IDO1", "PRF1", "GZMA", "HLA-DRA")
#ribas.exp.imm.signature <- c("CD3D", "IDO1", "CIITA", "CD3E", "CCL5", "GZMK", "CD2", "HLA-DRA", "CXCL13", "IL2RG",
#                             "NKG7", "HLA-G","CXCR6", "LAG3", "TAGAP", "CXCL10", "STAT1", "GZMB")
#brooks.ifn.signature <- c("BST2", "EIF2AK2", "IFI27", "IFI6", "IFIT1", "IFIT2", "ISG15", "MX1", "OAS1", "OASL", "RNASEL", "USP18", "XAF1")

ayers.expr.18 <- new.dat$expr[AYERS.18,]
ribas.ifng.expr.6 <- new.dat$expr[RIBAS.IFN6,]
rooney.cyt <- new.dat$expr[ROONEY.CYT,]

GEP.scores <- data.frame(
  Ayers.Score = apply(ayers.expr.18, 2, mean),
  IFNgamma.Score = apply(ribas.ifng.expr.6, 2, mean),
  CYT.Score = apply(rooney.cyt, 2, mean)
)

write.table(GEP.scores, file = paste(outdir,"/IFNg",".scores.tsv", sep = ""), sep = "\t", quote = F)
write.table(ayers.expr.18, file = paste(outdir,"/Ayers18Genes.log2TPM.tsv", sep = ""), sep = "\t", quote = F)
write.table(ribas.ifng.expr.6, file = paste(outdir,"/Ribas6IFNg.log2TPM.tsv", sep = ""), sep = "\t", quote = F)

##############
# STEP 4: CALCULATE ssGSEA scores
##############
library ("GSEABase")
library ("GSVA");
load(paste(path.head,"src/R_wrappers/ssGSEA/data/YoshiharaGeneSetCollection.RData", sep = "/"))
load(paste(path.head,"src/R_wrappers/ssGSEA/data/Angelova.gsc.RData", sep = "/"))
load(paste(path.head,"src/R_wrappers/ssGSEA/data/Bindea.gsc.RData", sep = "/"))
load(paste(path.head,"src/R_wrappers/ssGSEA/data/Bowtell.gsc.RData", sep = "/"))
load(paste(path.head,"src/R_wrappers/ssGSEA/data/hallmark.gsc.RData", sep = "/"))
data.set <- "INSPIRE_normalized_ComBat";
genesets <- list(yoshihara = Yoshihara.gsc, bowtell = Bowtell.gsc, hallmark = hallmark.gsc, bindea = Bindea.gsc, angelova = Angelova.gsc)

for (i in 1:5){
  Scores <- gsva(as.matrix(new.dat.ensembl$expr),
                 genesets[[i]],
                 method = "ssgsea",
                 min.sz = 2,
                 max.sz = 500,
                 ssgsea.norm = T)
  
  write.table(t(Scores), file = paste(outdir,"/", data.set, ".",names(genesets)[i],".scores.tsv", sep = ""), sep = "\t", quote = F);
}
rm(Scores)



