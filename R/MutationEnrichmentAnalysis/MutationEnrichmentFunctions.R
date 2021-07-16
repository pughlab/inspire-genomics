####################################
# MutationEnrichmentFunctions.R
#####################################
library (maftools)
library(plyr)
library(dplyr)

splitData <- function (x, df){
  tmp <- df[,c("Study.ID", x)]
  colnames(tmp)[2] <- "gene"
  return (tmp)
}

## Split into annotation and score
splitScore <- function (x, type){
  switch (type,
          category = unlist(strsplit (x, split = "(", fixed = TRUE))[1],
          score = as.numeric(gsub (".*(\\d.\\d+).*","\\1",x)))
}

calcBkgMutRate <- function(gene, perc.dat, tcga.summary ,tcga.dat) {
  df <- tcga.summary[[gene]]
  print(head(df))
  idx.match <- match(row.names(perc.dat[perc.dat$tcga.match ==1,]), as.character(tcga.dat$Cancer.type))
  print(idx.match)
  df.match <- df[idx.match,]
  
  idx.others <- which(!tcga.dat$matching.groups & as.character(tcga.dat$Tissue) == "SOLID")
  df.others <- df[idx.others,]
  df.others$ratio.samples <- df.others$n.sample/sum(df.others$n.sample)
  df.others$scaled.mut.rate <- df.others$perc.mutated*df.others$ratio.samples
  print(df.others)
  print(sum(df.others$scaled.mut.rate))
  df.match <- rbind(df.match,
                    c("OTHERS",
                      sum(df.others$n.sample), 
                      sum(df.others$n.mutated), 
                      round(sum(df.others$scaled.mut.rate), 2)))
  print(df.match)
  mut.rate <- data.frame(
    Study.ID = df.match$Study.ID,
    GR1.samples = round(as.numeric(df.match$n.sample)* as.numeric(perc.dat[,"perc.grp1"])/100),
    GR2.samples = round(as.numeric(df.match$n.sample)* as.numeric(perc.dat[,"perc.grp2"])/100),
    GR1.perc.mut = as.numeric(df.match$perc.mutated)* as.numeric(perc.dat[,"perc.grp1"])/100,
    GR2.perc.mut = as.numeric(df.match$perc.mutated) * as.numeric(perc.dat[,"perc.grp2"])/100
  )
  
  print(mut.rate)
  return(mut.rate)
}
testGeneEnrichmentBkg <- function (gene, enrichment.counts){
  bkg.mut.tcga <- calcBkgMutRate(gene, group.type.counts.merged, pancan.muts.summary ,tcga.cancer.types)
  group.counts <- enrichment.counts[enrichment.counts$Hugo_Symbol %in% gene & enrichment.counts$keep,]
  group.counts <- as.numeric(unlist(strsplit(group.counts[,"n_mutated_group1"],split = " ", fixed = T))[c(1,3)])
  print(group.counts)
  response.group <- paste(enrichment.counts[enrichment.counts$Hugo_Symbol %in% gene & enrichment.counts$keep,"Group1"], "perc.mut", sep = ".")
  response.group.n <- paste(enrichment.counts[enrichment.counts$Hugo_Symbol %in% gene & enrichment.counts$keep ,"Group1"], "samples", sep = ".")
  bkg.mut.perc <- sum(bkg.mut.tcga[,response.group])
  bkg.n <- sum(bkg.mut.tcga[,response.group.n])
  bkg.mut.n <- round(as.numeric(bkg.mut.perc)/100*as.numeric(bkg.n),0)
  
  mat <- matrix(as.numeric(c(group.counts[1], group.counts[2]-group.counts[1], bkg.mut.n, bkg.n - bkg.mut.n)), nrow = 2, byrow = TRUE)
  print(mat)
  test <- fisher.test(mat,alternative="two.sided")
  rec <- enrichment.counts[enrichment.counts$Hugo_Symbol %in% gene & enrichment.counts$keep,]
  stats <- c(
    gene = gene, 
    Group1 = rec$Group1,
    Group2 = rec$Group2,
    n_mutated_group1 = rec$n_mutated_group1,
    n_mutated_group2 = rec$n_mutated_group2,
    n_mutated_BkgPanCan = bkg.mut.n,
    PanCan_n = bkg.n,
    Pval = test$p.value,
    OR = test$estimate,
    Conf.int = test$conf.int)
  return(stats)
}

mutationEnrichment <- function (maf, clin.dat, group.feature = "Pembro_CB_groups", genes_of_interest = NA){
  MAF <- read.maf (maf = maf,clinicalData = clin.dat)
  response.ce <- clinicalEnrichment(maf = MAF, clinicalFeature = group.feature)
  enrichment.counts <- data.frame(response.ce$groupwise_comparision)
  enrichment.counts$group1.perc <-  round(100*as.integer(gsub("(.*).of.*", "\\1", enrichment.counts$n_mutated_group1)) /
                                            as.integer(gsub(".*of.(.*$)", "\\1", enrichment.counts$n_mutated_group1)))
  enrichment.counts$group2.perc <- round(100*as.integer(gsub("(.*).of.*", "\\1", enrichment.counts$n_mutated_group2)) /
                                           as.integer(gsub(".*of.(.*$)", "\\1", enrichment.counts$n_mutated_group2)))
  # selection criteria: at least 25% in one group and at most 5% in the other
  enrichment.counts$keep <- enrichment.counts$group1.perc >=24 & enrichment.counts$group2.perc <= 4 
  return (results <- list (enrichment.counts, response.ce))                                       
}