###############################################
# INSPIRE mutation enrichment analysis
# Original analysis without filtering
# Feb 2021 - updated with filtering variants by SIFT and PolyPhen
###############################################
library (maftools)
library(plyr)
library(dplyr)

setwd ("/Users/yangc/PughLab/Repos/inspire-genomics/")
source (paste(getwd(),"R","MutationEnrichmentFunctions.R", sep = "/"))

# Step 1 - load source mutation and clinical annotation data
load (paste(getwd(),"data","PanCanMut.data.RData", sep = "/"))
ctDNA.annotation <- read.csv(paste(getwd(),"data","ctDNA.annotations.csv", sep = "/"),
                             stringsAsFactors = F, header = T, row.names = 1)

# make corrections on cancer typ annotations
pt.dat$CANCER_TYPE[pt.dat$CANCER_TYPE == "COADRED"] <- "COADREAD"
pt.dat$tcga.tumor.types <- 
  factor(pt.dat$CANCER_TYPE, 
         levels = unique(pt.dat$CANCER_TYPE),
         labels = c("HNSC", "HNSC", "HNSC", "HNSC", "HNSC", "BRCA",
                    "OV", "SKCM", "UVM", "SKCM", "VMM", "SARC", 
                    "SARC", "CESC", "NPC", "UCEC", "SARC", "SCUP",
                    "MCC", "SARC", "CENE", "GEJ", "UCS", "PSCC",
                    "ANSC", "ACYC", "COADREAD", "SKCM", "CHS", "UCEC",
                    "SARC", "CHOL", "PANET"))

# Step 2 - prepare group comparisons annotation for maftools object
clin.dat <- data.frame(
  Tumor_Sample_Barcode = unique(maf$Tumor_Sample_Barcode),
  ptid = gsub("(INS-.-...).*", "\\1",unique(maf$Tumor_Sample_Barcode)),
  stringsAsFactors = FALSE
)
clin.dat <- cbind(clin.dat, pt.dat[as.character(clin.dat$ptid),c("CANCER_TYPE","tcga.tumor.types")])
clin.dat <- cbind(clin.dat, ctDNA.annotation[as.character(clin.dat$ptid), c("COHORT","CB","Pembro_CB_groups")])
clin.dat$Pembro_CB_groups <- gsub(" ", "", clin.dat$Pembro_CB_groups)

# Step 3 - Process SIFT annotation (deleterious) - added for rebuttal
maf$SIFT_category <- sapply(maf$SIFT, splitScore, "category")
maf$SIFT_score <- sapply(maf$SIFT, splitScore, "score")

# Step 4 - Process PolyPhen annotation (damaging) - added for rebuttal
maf$PolyPhen_category <- sapply(maf$PolyPhen, splitScore, "category")
maf$PolyPhen_score <- sapply(maf$PolyPhen, splitScore, "score")

# BRCA2 mutation categories
brca2 <- maf[maf$Hugo_Symbol == "BRCA2",]
table (brca2$SIFT_category, brca2$PolyPhen_category) # only 2/7 mut are damaging and deleterious, 1/7 deleterious+benign

# Step 5 - Enrichment analysis and mutation counts
# mutation counts using all non-synonymous mutations in the MAF file
all.ns.mut.enrichment <- mutationEnrichment (maf = maf, clin.dat = clin.dat)

# mutation counts using mutations that are found to be deterious and damaging by SIFT and PolyPhen only
mut.enrichment.stringent <- mutationEnrichment (
  maf = maf[maf$SIFT_category %in% c("deleterious") & maf$PolyPhen_category %in% c("probably_damaging"),], 
  clin.dat = clin.dat)

# mutation counts using mutations that are found to be either deterious or damaing by SIFT and PolyPhen
mut.enrichment.filtered <- mutationEnrichment (
  maf = maf[maf$SIFT_category %in% c("deleterious") | maf$PolyPhen_category %in% c("probably_damaging"),], 
  clin.dat = clin.dat)

# Step 6 - Calculate TCGA background mutation rate for selected genes (genes short listed in the enrichment analysis above, downloaded through cbioportal.org)
pancan.muts <- "/mnt/work1/users/pughlab/projects/INSPIRE/Review/maftools/PanCancer/mutation_enrichment/alterations_across_samples_new.tsv"
pancan.muts <- read.csv(pancan.muts, sep = "\t", header = T, stringsAsFactors = F)
pancan.muts <- pancan.muts[,c(1:3, grep ("MUT", colnames(pancan.muts)))] 
colnames(pancan.muts) <- gsub ("..MUT", ".MUT", colnames(pancan.muts))

data.list <- lapply(colnames(pancan.muts)[grep("MUT", colnames(pancan.muts))], 
                    splitData, df = pancan.muts)
names(data.list) <- gsub(".MUT", "", colnames(pancan.muts)[grep("MUT", colnames(pancan.muts))])
pancan.muts.summary <- lapply(data.list, ddply, .variables = ~Study.ID, 
                              .fun = summarize, 
                              n.sample = sum(!gene %in% "not profiled"), 
                              n.mutated = sum(!gene %in% c("no alteration", "not profiled")),
                              perc.mutated = 100*round(sum(!(gene %in% c("no alteration", "not profiled")))/length(Study.ID), digits =4))

tcga.cancer.types <- data.frame(
  Study.ID = unique(pancan.muts[,1]),
  Cancer.type = toupper(gsub("_.*", "", unique(pancan.muts[,1])))
)

tcga.cancer.types$Tissue <- factor(tcga.cancer.types$Cancer.type %in% c("LAML", "DLBC"),
                                   levels = c("TRUE", "FALSE"),
                                   labels = c("BLOOD", "SOLID"))

# Get the proportion of tumor types in Responder group and non responder group
group.type.counts <- as.data.frame.matrix (table(clin.dat$tcga.tumor.types, clin.dat$Pembro_CB_groups))
group.type.counts$perc.grp1 <- 100*round(group.type.counts[,2]/sum(group.type.counts[,2]), digit = 2)
group.type.counts$perc.grp2 <- 100*round(group.type.counts[,3]/sum(group.type.counts[,3]), digit = 2)
group.type.counts$tcga.match <- row.names(group.type.counts) %in% tcga.cancer.types$Cancer.type

# Calculate background mutation rate of theoretical cohorts from tcga
group.type.counts.merged <- rbind(group.type.counts[group.type.counts$tcga.match,],
                                  colSums(group.type.counts[!group.type.counts$tcga.match,]))
row.names(group.type.counts.merged)[12] <- "OTHER"
tcga.cancer.types$matching.groups <- tcga.cancer.types$Cancer.type %in% row.names(group.type.counts.merged)


