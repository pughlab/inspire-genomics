########
# Setup library path 
.libPaths(c("/cluster/projects/pughlab/src/r_lib/library/3.6", .libPaths()))
library ("tximport")
library ("DESeq2")
library ("pheatmap")

########
# FUNCTIONS
########
# return complete nested model matrix design
createModelMatrix <- function (formula.text, design.matrix){
  design.formula <- as.formula (formula.text)
  mm <- model.matrix (design.formula, design.matrix)
  mm <- mm[,colSums(mm)>0]
  return (mm)
}

runDESeq <- function (txi, model.matrix, design.matrix, formula.text){
  design.formula <- as.formula (formula.text)
  txi$length[txi$length == 0] <- 1
  dd <- DESeqDataSetFromTximport(
    txi = txi,
    colData = design.matrix[row.names(model.matrix),],
    design = design.formula
  )
  dd <- dd[rowSums(counts(dd))>10,]
  dd <- DESeq (dd, full = model.matrix, betaPrior = FALSE)
  return (dd)
}

########
# VARIABLES
#########
# Protein coding genes/ GTF annotation
GENCODE_v26 <- read.csv("/cluster/projects/pughlab/references/gencode/gencode.v26.ensg_annotation.txt",
                               sep = "\t", header = FALSE, stringsAsFactors = FALSE)
GENCODE_v26_PC <- GENCODE_v26[GENCODE_v26[,7]=="protein_coding",c(1,8)] 
print(head(GENCODE_v26_PC, n = 5))

# Step 0 - Verify input argument information
args <- c(
  "/cluster/projects/pughlab/projects/INSPIRE_ms1/RNAseq/inventory/INSPIRE_RNASEQ_bamfiles_12032019CY.avail.csv",
  "default",
  "/cluster/projects/pughlab/projects/INSPIRE_ms1/RNAseq/DESeq2",
  "INSPIRE_NCB"
)

#args <- commandArgs(trailingOnly = TRUE)

if (length(args)<4){
  print ("Usage: Rscript ProcessRSEM.R <INVENTORY_FILE> <OUT_DIR> <PROJECT_CODE>\n Output: Tab-deliminated matrix of RSEM counts by gene")
}else{
  args <- c(args, paste0(args[3], '/', args[4], "_", Sys.Date() ,"_txi.rsem.counts.tsv"))
  names(args) <- c("INVENTORY","SAMPLE_TABLE","OUT_DIR","PCODE","MATRIX_FILE")
}

# Print the input arguments to confirm filepaths to user
print (c(args))

########
# MAIN
########
# Step 1 - Read the inventory file or load default design matrix

if (args["SAMPLE_TABLE"] == "default"){
  #load ("/cluster/projects/pughlab/projects/INSPIRE_ms1/RNAseq/R/DEseq2.Paired.Design.Matrix.RData")
  load ("/cluster/projects/pughlab/projects/INSPIRE_ms1/RNAseq/R/DEseq2.Paired.Design.Matrix.ctDNA.RData")
  sampleTable <- DE.design.full[as.character(DE.design.full$Response)=="N",]
  rm(DE.design.full)
  print(head(sampleTable, n=5))
  
}else{
  INVENTORY <- read.csv (args["SAMPLE_FILE"], sep = ",",stringsAsFactors = FALSE, header = FALSE)
  print(dim(INVENTORY))
  colnames(INVENTORY) <- c("PATIENT_ID",
                           "SAMPLE_ID",
                           "TIMEGROUP",
                           "TIMEPOINT",
                           "COHORT",
                           "SEQ_BATCH",
                           "CLINICAL_BENEFIT",
                           "CTDNA_CHANGE",
                           "ICB_SENSITIVITY")
  print(head(INVENTORY, n = 5))
  
  # Step 3 - Create the design table (Factorize)
  sampleTable <- data.frame(
    nPatient.ID = ,
    Treatment = ,
    Response = ,
    stringsAsFactors = FALSE
  )
  row.names(sampleTable)
}

# Step 2 - Import RSEM count data as indicated by sampleTable
INVENTORY <- read.csv (args["INVENTORY"], sep = ",",stringsAsFactors = FALSE, header = FALSE)
row.names(INVENTORY) <- make.names(INVENTORY[,2])
print(dim(INVENTORY))
colnames(INVENTORY) <- c("PATIENT_ID",
                         "SAMPLE_ID",
                         "TIMEGROUP",
                         "TIMEPOINT",
                         "COHORT",
                         "FILTER",
                         "KEEP",
                         "BAM_DIR",
                         "BAM_FILE",
                         "RSEM_PATH")
print(head(INVENTORY, n = 5))
files <- INVENTORY[row.names(sampleTable),"RSEM_PATH"]
names (files) <- INVENTORY[row.names(sampleTable),"SAMPLE_ID"]
# Step 2 - use tximport to read the output files and create count matrix
# all samples
txi.rsem <- tximport(files, type = "rsem", txIn = TRUE, txOut = TRUE)
# Step 3 - write countfile to output_dir
# optional - checkif outputput dir exist, if not create new and print out warning
write.table (txi.rsem$counts, sep = "\t", quote = FALSE, file = args["MATRIX_FILE"])

# Step 4 - DESeq2 Analysis 
# if analyzing single timepoint, model design including seqeuncing batch effect (PMGC vs TGL)
# if analyzing paired samples, batch is nested within patient when samples are paired
# Models 
# "~Batch+Response"
# "~Batch+Cohort+Response"
# "~Patient.ID+Treatment"
# "~Response + Response:nPatient.ID + Response:Treatment" 
#DE_MODEL <- "~Response + Response:nPatient.ID + Response:Treatment"
DE_MODEL <- "~Patient.ID + Treatment"
print ("Created Model Matrix")
sampleTable[,"Patient.ID"]<- factor (as.character(sampleTable[,"Patient.ID"]))
print(str(sampleTable))
model.matrix <- createModelMatrix (formula.text=DE_MODEL, design.matrix=sampleTable[,c("Patient.ID","Treatment")])
print (head (sampleTable, n = 5))
dimnames(txi.rsem$counts)[[2]] <- make.names(dimnames(txi.rsem$counts)[[2]])
dimnames(txi.rsem$abundance)[[2]] <- make.names(dimnames(txi.rsem$abundance)[[2]])
dimnames(txi.rsem$length)[[2]] <- make.names(dimnames(txi.rsem$length)[[2]])

print ("Running DESeq2")
dd <- runDESeq (txi=txi.rsem, model.matrix=model.matrix, design.matrix=sampleTable[,c("Patient.ID", "Treatment")], formula.text="~Treatment")
res <- results(dd, alpha = 0.10, lfcThreshold = 0)

print ("Write Results Table")
save(dd, res , file=paste0(args["OUT_DIR"],"/", args["PCODE"],"_",Sys.Date(),"_DESeq2.RData"))
df <- as.data.frame(res[order(res$pvalue),])
df <- cbind(df, "hugo.name" = GENCODE_v26[match(row.names(df), GENCODE_v26[,1]), 8])
write.csv(df,quote = F, file=paste0(args["OUT_DIR"],"/", args["PCODE"],"_",Sys.Date(),"_response_treated_results.csv"))

print ("Save ranked values for GSEA GOBP")
df$fcsign <- sign(df$log2FoldChange)
df$logP <- -log10(df$pvalue)
df$metric <- df$logP/df$fcsign
genelist.rk <- setNames (df$metric,df$hugo.name)
save (genelist.rk, df, GENCODE_v26, file = paste0(args["OUT_DIR"],"/",args["PCODE"],"_DESeq2.DE.genelist.rnk.sep.Rdata")
