############
# Comut plot Figure 3 Manuscript
library (XLConnect)
################
# Categorical data to plot:
# ICB Category
# Best Response RECIST
# COHORT
# Cancer type
# MSI status
################
#####
# VARIABLES
######
x.limit <- c(0,90)
main.df <- read.csv (file = main.csv, nrows = 106, header = T, check.names = F, stringsAsFactors = F)[,1:314]
sequenza.df <- readWorksheetFromFile(file = "~/Desktop/INSPIRE_ms2/tables/INSPIRE_Baseline_WES_Enrichment_Table (version 1).xlsx", 
                                     sheet = "Sequenza", header = T ,rownames = 1, check.names = F,
                                     startRow = 2, endRow = 211,startCol = 1, endCol = 14)
mut.enriched <- read.csv ("/Volumes/Samwise/projects/INSPIRE/Review/maftools/data/enriched.genes.selected.maf", 
                          sep = "\t", header = T, stringsAsFactors = F)
gene.list <- read.csv ("/Volumes/Samwise/projects/INSPIRE/Review/maftools/data/enriched.gene.list.txt", 
                       sep = "\t", header = T, stringsAsFactors = FALSE)
imm.gene.list <- read.csv ("/Volumes/Samwise/projects/INSPIRE/Review/maftools/data/immune.genes.list.txt", 
                           sep = "\t", header = T, stringsAsFactors = FALSE)
imm.mut <- read.csv ("/Volumes/Samwise/projects/INSPIRE/Review/maftools/data/immune.genes.selected.maf",
                     sep = "\t", header = T, stringsAsFactors = F)
COLS_MUTATIONS <- setNames(
  c("#2e8b57", "#ffb300", "#fd5e53", "#ca2c92", "#ca2c92", "#7851a9", "#7851a9", "#a50026", "#967117","white"),
  nm = c("Missense_Mutation", "Nonsense_Mutation", "Splice_Site", 
         "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", 
         "Translation_Start_Site","Multiple" ,"none")
)
variant.nm <- c("Missense", "Nonsense", "Splice site", "Frameshift indel", "In frame indel", "Translation start site")

# load CN data, BL.seg.LOH, BL.seg.tCN
load("~/Desktop/INSPIRE_ms2/data/BL.CN.by.Gene.RData")

#############
# FUNCTIONS #
#############
plotHistogram <- function (data.df, sample.order, metric,
                           xlimit, ylimit, target = "", 
                           ylabel, color = "#734f96"){
  #par(mar=c(6,5,3,1))
  #par(mar=c(0,6,0,3))
  par(mar=c(0,5,0,3))
  plot(1, type="n", axes=F, xlab="", ylab="", xlim = xlimit, ylim = ylimit,xaxs = "i")
  mp <- barplot( data.df[sample.order,metric], 
                 xlim = xlimit, ylim = ylimit, cex.names = 0.8, border=NA, xaxs = "i",
                 xaxt = "n", yaxt = "n",  add = T, col = color, lwd = 2)
  axis (2, at = c(0,1,2,2.5), labels = c(0, 10, 100, 316), cex.axis = 0.5, las = 2, col.axis = "black")
  title (ylab = gsub("_", " ", ylabel), font.lab = 2, cex.lab = 0.50, line = 2)
  return (mp)
}

plotHistogramStacked <- function (data.df, sample.order, metric,
                           xlimit, ylimit, target = "", 
                           ylabel){
  #par(mar=c(6,5,3,1))
  #par(mar=c(0,3,0,3))
  par(mar=c(0,5,0,3))
  plot(1, type="n", axes=F, xlab="", ylab="", xlim = xlimit, ylim = ylimit,xaxs = "i")
  data.df <- data.df[sample.order,metric]
  data.df <- do.call ("cbind", lapply (seq (1, length(metric), by =1), 
                                       function(x) as.numeric(gsub ("%", "", data.df[,x]))))
  print(head(data.df))
  data.df <- cbind(data.df, 100-rowSums(data.df))
  print(head(data.df))
  mp <- barplot( t(data.df), 
                 xlim = xlimit, ylim = ylimit, cex.names = 0.8, border=NA, xaxs = "i",
                 xaxt = "n", yaxt = "n",  add = T, 
                 col = c("#0038a8","#73c2fb","#0067a5","#e0b0ff","#734f96","#ffdf00", "white"), 
                 lwd = 2, beside = FALSE)
  axis (2, at = target, labels = target, cex.axis = 0.5, las = 2, col.axis = "black")
  #title (ylab = gsub("_", " ", ylabel), line = 1, font.lab = 2, cex.lab = 0.75)
  return (mp)
}

plotCovariate <- function (data.df, sample.order, metric, x.limit,colmap, mp, label.name, add.germline = NA, font = 1) {
  par(mar=c(0,5,0,3))
  y.limit = c(0,1.5)
  plot(x.limit, y.limit, type = "n", xlab = "", ylab = "",
       xaxs = "i",xaxt = 'n', yaxt = 'n',bty = 'n')
  rect ( xleft = mp-0.7, 
         ybottom = as.vector(rep (max(y.limit), length(mp))),
         xright = mp+0.7,
         ytop= as.vector(rep(min(y.limit), length (mp))),
         col = colmap[data.df[sample.order, metric]],
         border = "white")
  if (class(add.germline)=="matrix"){
    gl.colmap <- c("germline" = "black", "none" = adjustcolor("white", alpha.f = 0))
    rect ( xleft = mp-0.45, 
           ybottom = as.vector(rep (max(y.limit)*0.60, length(mp))),
           xright = mp+0.45,
           ytop= as.vector(rep(max(y.limit)*0.30, length (mp))), border = NA,
           col = gl.colmap[add.germline[sample.order, metric]])
  }
  #mtext(side = 2, at = 0.75, labels = label.name, cex.axis = 0.5, tick = FALSE, las = 2, line = 0.5)
  mtext (text = label.name, side = 2, at = 0.75, line = 0.2, cex = 0.50, las =2,font = font)
}

plotCovContinuous <- function (data.df, sample.order, metric, x.limit, 
                               colRamp = colorRampPalette(brewer.pal(5,"Blues")), mp, label.name){
  par(mar=c(0,5,0,3))
  y.limit = c(0,1.5)
  plot(x.limit, y.limit, type = "n", xlab = "", ylab = "",
       xaxs = "i",xaxt = 'n', yaxt = 'n',bty = 'n')
  
  #rank varaible for color assignment
  data.df$order <- findInterval(data.df[,metric], sort(data.df[,metric]))
  data.df$cols <- colRamp(nrow(data.df))[data.df$order]
  data.df$cols[is.na(data.df$order)] <- "grey"
  rect ( xleft = mp-0.7, 
         ybottom = as.vector(rep (max(y.limit), length(mp))),
         xright = mp+0.7,
         ytop= as.vector(rep(min(y.limit), length (mp))),
         col = data.df[sample.order,"cols"],
         border = "white")
  mtext (text = label.name, side = 2, at = 0.75, line = 0.2, cex = 0.50, las =2)
}
plotTrgl <- function (side = "right", i, xleft, xright, y.trgl, col = NA, border = NA) {
  if (side == "left"){
    x.trgl <- c(xleft[i], xright[i]-0.05, xright[i]-0.05)
  }else{
    x.trgl <- c(xleft[i]+0.05, xleft[i]+0.05, xright[i])
  }
  polygon (
    x = x.trgl, 
    y = y.trgl, 
    density = NULL, angle = 45,
    border = border[i], col = col[i], lty = par("lty"))
}



########
# MAIN #
########
#####
# Process Data
######
source ("~/Desktop/INSPIRE_ms2/src/SetUp.Rdata")
source("~/Desktop/INSPIRE_ms2/src/Color.Palettes.R")
setwd (work.dir)

sequenza.df$pid <- as.character(gsub ("(INS-.-...)-.*", "\\1", row.names(sequenza.df)))
sequenza.bl.df <- sequenza.df[grep ("ST",as.character(sequenza.df$Time)),]

main.df$logTMB <- log10(main.df$TMB+1)
data.df <- main.df[main.df$`WES_BASELINE_AVAIL (N=72)`==1,]
sample.order <- order(data.df[,"logTMB"], decreasing = T)
sample.order <-order(data.df[,"ICB_SENSITIVTY_CAT"], data.df[,"CR_PR_SD>6cycles_MARCO"], decreasing = T)
sample.order <-order(data.df[,"CR_PR_SD>6cycles_MARCO"],data.df[,"ICB_SENSITIVTY_CAT"], decreasing = T)

wes.df <- data.frame(
  sid = sort(unique(mut.enriched$Tumor_Sample_Barcode)), 
  pid = gsub ("(INS-.-...).*", "\\1", sort(unique(mut.enriched$Tumor_Sample_Barcode))),
  timept = gsub ("INS-.-...-(.*)", "\\1", sort(unique(mut.enriched$Tumor_Sample_Barcode))), 
  stringsAsFactors = FALSE
)

wes.df$germline <- ifelse (wes.df$timept=="SB", yes = "germline", no = "somatic")
wes.df$baseline <- ifelse (wes.df$timept %in% c("ST", "ST-S"), yes = "baseline", no =  "")

filtered.mut.enriched <- mut.enriched [mut.enriched$Tumor_Sample_Barcode %in% wes.df$sid[wes.df$baseline == "baseline"] |
                                         mut.enriched$Tumor_Sample_Barcode %in% wes.df$sid[wes.df$germline == "germline"], ]
filtered.mut.enriched <- filtered.mut.enriched [filtered.mut.enriched$Variant_Classification %in% names(COLS_MUTATIONS)[1:8],]
filtered.mut.enriched <- filtered.mut.enriched[filtered.mut.enriched$Hugo_Symbol %in% gene.list[,1],]
filtered.mut.enriched$pid <- gsub ("(INS-.-...).*", "\\1", filtered.mut.enriched$Tumor_Sample_Barcode)
filtered.mut.enriched <- filtered.mut.enriched[filtered.mut.enriched$pid %in% data.df$PATIENT_ID[sample.order],]
filtered.mut.enriched$pid <- factor (filtered.mut.enriched$pid, levels = data.df$PATIENT_ID[sample.order])
filtered.mut.enriched$Variant_Classification <- factor (filtered.mut.enriched$Variant_Classification, levels = names(COLS_MUTATIONS)[1:8])

mut.matrix <- setNames(lapply (gene.list[,1], function (x) 
  table (subset(filtered.mut.enriched, Hugo_Symbol == x)$pid,
         subset(filtered.mut.enriched, Hugo_Symbol == x)$Variant_Classification)), nm = gene.list[,1])

collpased.mut.matrix <- do.call ("cbind", lapply (mut.matrix, function(x) 
  ifelse (rowSums (x) > 1, yes = "Multiple", no = ifelse (
    x[,1] > 0, yes = "Missense_Mutation", no = ifelse(
      x[,2] > 0, yes = "Nonsense_Mutation", no = ifelse(
        x[,3] > 0, yes = "Splice_Site", no = ifelse(
          x[,4] > 0, yes = "Frame_Shift_Del", no = ifelse(
            x[,5] > 0, yes = "Frame_Shift_Ins", no = ifelse(
              x[,6] > 0, yes = "In_Frame_Del", no = ifelse(
                x[,7] > 0, yes = "In_Frame_Ins", no = ifelse(
                  x[,8] > 0, yes = "Translation_Start_Site", no = "none")))))))))))

mut.gl.matrix <- setNames(lapply (gene.list[,1], function (x) 
  table (subset(filtered.mut.enriched[filtered.mut.enriched$Tumor_Sample_Barcode %in% wes.df[wes.df$germline=="germline","sid"],], Hugo_Symbol == x)$pid,
         subset(filtered.mut.enriched[filtered.mut.enriched$Tumor_Sample_Barcode %in% wes.df[wes.df$germline=="germline","sid"],], Hugo_Symbol == x)$Variant_Classification)), nm = gene.list[,1])

collpased.gl.mut.matrix <- do.call ("cbind", lapply (mut.gl.matrix, function(x) ifelse (rowSums (x) > 0, yes = "germline", no = "none")))

filtered.mut.imm <- imm.mut [imm.mut$Tumor_Sample_Barcode %in% wes.df$sid[wes.df$baseline == "baseline"] |
                             imm.mut$Tumor_Sample_Barcode %in% wes.df$sid[wes.df$germline == "germline"], ]
filtered.mut.imm <- filtered.mut.imm [filtered.mut.imm$Variant_Classification %in% names(COLS_MUTATIONS)[1:8],]
filtered.mut.imm <- filtered.mut.imm[filtered.mut.imm$Hugo_Symbol %in% imm.gene.list[,1],]
filtered.mut.imm$pid <- gsub ("(INS-.-...).*", "\\1", filtered.mut.imm$Tumor_Sample_Barcode)
filtered.mut.imm <- filtered.mut.imm[filtered.mut.imm$pid %in% data.df$PATIENT_ID[sample.order],]
filtered.mut.imm$pid <- factor (filtered.mut.imm$pid, levels = data.df$PATIENT_ID[sample.order])
filtered.mut.imm$Variant_Classification <- factor (filtered.mut.imm$Variant_Classification, levels = names(COLS_MUTATIONS)[1:8])


# remove germline SNP
table(filtered.mut.imm[filtered.mut.imm$Tumor_Sample_Barcode %in% wes.df$sid[wes.df$germline == "germline"],"dbSNP_RS"])
filtered.mut.imm <- filtered.mut.imm[filtered.mut.imm$dbSNP_RS != "rs397509199,rs1799966",]
imm.mut.matrix <- setNames(lapply (imm.gene.list[,1], function (x) 
  table (subset(filtered.mut.imm, Hugo_Symbol == x)$pid,
         subset(filtered.mut.imm, Hugo_Symbol == x)$Variant_Classification)), nm = imm.gene.list[,1])

collpased.imm.mut.matrix <- do.call ("cbind", lapply (imm.mut.matrix, function(x) 
  ifelse (rowSums (x) > 1, yes = "Multiple", no = ifelse (
    x[,1] > 0, yes = "Missense_Mutation", no = ifelse(
      x[,2] > 0, yes = "Nonsense_Mutation", no = ifelse(
        x[,3] > 0, yes = "Splice_Site", no = ifelse(
          x[,4] > 0, yes = "Frame_Shift_Del", no = ifelse(
            x[,5] > 0, yes = "Frame_Shift_Ins", no = ifelse(
              x[,6] > 0, yes = "In_Frame_Del", no = ifelse(
                x[,7] > 0, yes = "In_Frame_Ins", no = ifelse(
                  x[,8] > 0, yes = "Translation_Start_Site", no = "none")))))))))))

HLA.mut <- main.df[match(row.names(collpased.imm.mut.matrix), as.character(main.df$PATIENT_ID)), 
                   c("HLA-A MUT", "HLA-B MUT", "HLA-C MUT")]

HLA.mut <- do.call ("cbind", lapply (colnames(HLA.mut), function (x) 
  ifelse (HLA.mut[,x] == "none" | HLA.mut[,x] == "", yes = "none", no = ifelse (
    HLA.mut[,x] == "p.S349_splice", yes = "Splice_Site", no = ifelse (
      HLA.mut[,x] == "p.H216Y", yes = "Missense_Mutation", no = "none")))))
#colnames(HLA.mut) <- c("HLA.A", "HLA.B", "HLA.C")
dimnames(HLA.mut) <- list(row.names(collpased.imm.mut.matrix),c("HLA.A", "HLA.B", "HLA.C"))

collpased.imm.mut.matrix <- cbind( collpased.imm.mut.matrix, HLA.mut)

imm.gene.list.collection <- list(
  "Antigen presentation" = c(imm.gene.list[27:29,1], c("HLA.A", "HLA.B","HLA.C")),
  "Interferon signaling" = imm.gene.list[30:41,1],
  "DNA repair" = imm.gene.list[1:10,1], 
  "Immune evasion" = imm.gene.list[11:26,1],
  "Others" = imm.gene.list[42:50,1]
)

collapsed.imm.mut.collection.matrix <- do.call ("cbind", lapply(imm.gene.list.collection, function(x) 
  ifelse (rowSums (collpased.imm.mut.matrix[,x]!= "none") > 0, yes = "altered", no = "unaltered")))
          
imm.mut.gl.matrix <- setNames(lapply (imm.gene.list[,1], function (x) 
  table (subset(filtered.mut.imm[filtered.mut.imm$Tumor_Sample_Barcode %in% wes.df[wes.df$germline=="germline","sid"],], Hugo_Symbol == x)$pid,
         subset(filtered.mut.imm[filtered.mut.imm$Tumor_Sample_Barcode %in% wes.df[wes.df$germline=="germline","sid"],], Hugo_Symbol == x)$Variant_Classification)), nm = imm.gene.list[,1])

collpased.gl.imm.mut.matrix <- cbind(do.call ("cbind", lapply (imm.mut.gl.matrix, 
                                                         function(x) ifelse (rowSums (x) > 0, yes = "germline", no = "none"))),
                                     HLA.A = "none", HLA.B = "none", HLA.C = "none")
collpased.gl.imm.mut.collection.matrix <- do.call ("cbind", lapply(imm.gene.list.collection, function(x) 
  ifelse (rowSums (collpased.gl.imm.mut.matrix[,x]!= "none") > 0, yes = "germline", no = "none")))

imm.mut.matrix <- setNames(lapply (imm.gene.list[,1], function (x) 
  table (subset(filtered.mut.imm, Hugo_Symbol == x)$pid,
         subset(filtered.mut.imm, Hugo_Symbol == x)$Variant_Classification)), nm = imm.gene.list[,1])

### Copy number matrix

colnames(BL.seg.tCN)[grep ("INS", colnames(BL.seg.tCN))] <- gsub ("(INS-.-...)-.*", "\\1", colnames(BL.seg.tCN)[grep ("INS", colnames(BL.seg.tCN))])
colnames(BL.seg.LOH)[grep ("INS", colnames(BL.seg.LOH))] <- gsub ("(INS-.-...)-.*", "\\1", colnames(BL.seg.LOH)[grep ("INS", colnames(BL.seg.LOH))])
sequenza.bl.df["INS-E-018-ST","ploidy"] <- 3
sequenza.bl.df["INS-E-018-ST","cellularity"] <- 0.57
ploidies <- round(as.numeric(sequenza.bl.df[match (colnames(BL.seg.tCN)[grep ("INS", colnames(BL.seg.tCN))], sequenza.bl.df$pid),"ploidy"]))

CN.matrix <- BL.seg.tCN[,grep ("INS", colnames(BL.seg.tCN))]
LOH.matrix <- BL.seg.LOH[,grep ("INS", colnames(BL.seg.LOH))]
CN.change.matrix <- BL.seg.tCN[,grep ("INS", colnames(BL.seg.tCN))] - ploidies

# CN = 0,  gain = CN = ploidy + 2, LOH+loss, nCN+LOH 
sub.gene.list <- c("B2M", "TAP1", "TAP2", "HLA-A", "HLA-B", "HLA-C")
imm.gene.list.collection <- list(
  "Antigen presentation" = c(imm.gene.list[27:29,1], c("HLA-A", "HLA-B","HLA-C")),
  "Interferon signaling" = imm.gene.list[30:41,1],
  "DNA repair" = imm.gene.list[1:10,1], 
  "Immune evasion" = imm.gene.list[11:26,1],
  "Others" = imm.gene.list[42:50,1]
)
sub.gene.list <- unlist(imm.gene.list.collection)
sub.LOH <- LOH.matrix[match(sub.gene.list, BL.seg.LOH$genename),]
sub.tCN <- CN.matrix[match(sub.gene.list, BL.seg.LOH$genename),]
sub.CN.change <- CN.change.matrix[match(sub.gene.list, BL.seg.LOH$genename),]

CN.categories <- ifelse (sub.tCN==0, yes = "tCN0", 
                         no = ifelse (sub.LOH==0 & sub.tCN==1, yes = "LOSS-LOH", 
                                      no = ifelse (sub.LOH==0 & sub.tCN>1 & sub.CN.change > 2, yes = "GAIN-LOH", 
                                                   no = ifelse (sub.LOH==0 & sub.tCN>1 & sub.CN.change <= 2, yes = "CN-LOH",
                                                                no = ifelse(sub.CN.change > 2, yes = "GAIN", no = "no change")))))
COLS_CN <- c("no change" = "white", "no data" = "grey","tCN0" = "#e5b73b" ,"GAIN" = "#a50026",
             "CN-LOH" = "#66ddaa","GAIN-LOH" = "#734f96", "LOSS-LOH" = "#4575b4")

CN.categories <- t(CN.categories)
colnames(CN.categories) <- sub.gene.list

CN.categories <- rbind (CN.categories, 
                        matrix ("no data", nrow = length(setdiff(as.character(data.df$PATIENT_ID)[sample.order], row.names(CN.categories))), 
                                ncol = ncol (CN.categories), 
                                dimnames = list (setdiff(as.character(data.df$PATIENT_ID)[sample.order], row.names(CN.categories)), 
                                                 colnames(CN.categories))))
#######
# Plot comut
#pdf(outfile, width = 12, height = 8)
#pdf(outfile, paper = "a4r")
#######
screen.matrix <-
  rbind(rbind (c(0.01, 0.80, 0.940, 0.990), 
               c(0.01, 0.80, 0.885, 0.935)),
        cbind (rep (0.01, 5), 
               rep (0.80, 5), 
               sort(seq(0.835, length.out = 5, by = 0.012), decreasing = TRUE) - 0.012, 
               sort(seq(0.835, length.out = 5, by = 0.012), decreasing = TRUE)),
        cbind (rep (0.01, 35), 
              rep (0.80, 35), 
              sort(seq(0.515, length.out = 35, by = 0.009), decreasing = TRUE) - 0.008, 
              sort(seq(0.515, length.out = 35, by = 0.009), decreasing = TRUE)),
        cbind (rep (0.01, 58), 
               rep (0.80, 58), 
               sort(seq(0.015, length.out = 58, by = 0.0085), decreasing = TRUE) - 0.008, 
               sort(seq(0.015, length.out = 58, by = 0.0085), decreasing = TRUE))
        )

split.screen (screen.matrix)

data.df$test.groups <- ifelse(data.df$ICB_group=="X1.LS", yes = "1.LS", no = 
                                ifelse (data.df$ICB_group=="X4.HS", yes = "2. HS/CB", no = 
                                          ifelse (data.df$CB_mod == "X1", yes = "2. HS/CB", no = "0.Not included")))
sample.order <- order (data.df$test.groups, data.df$logTMB, decreasing = TRUE)
sample.order <- order (data.df$test.groups, data.df$ICB_group, data.df$logTMB, decreasing = TRUE)
screen(1)             
mp <- plotHistogram (data.df = data.df, 
               sample.order = sample.order, 
               metric = "logTMB",
               xlimit = x.limit, 
               ylimit = c(0,2.5), 
               target = "", 
               ylabel = "TMB\n(ns Mut/Mb)")

screen(2)
plotHistogramStacked (data.df = data.df, 
               sample.order = sample.order, 
               metric = c(
                 "Age_signature",
                 "Alkylating_agents_signature",
                 "APOBEC_signature",
                 "DSB_signature",
                 "MMR_signature",
                 "UV_signature"),
               xlimit = c(0,90), 
               ylimit = c(0,100), 
               target = c(0,50,100), 
               ylabel = "")

screen(3)
main.df$ICB_SENSITIVTY_CAT_mod <- make.names(main.df$ICB_SENSITIVTY_CAT)
data.df$ICB_SENSITIVTY_CAT_mod <- make.names(data.df$ICB_SENSITIVTY_CAT)
main.df$ICB_group <- as.character(factor(
  paste(main.df$CTDNA_DELTA_CATEGORY, main.df$TARGETED_LESION_TM_EARLY>0, sep = ","),
  levels = unique(paste(main.df$CTDNA_DELTA_CATEGORY, main.df$TARGETED_LESION_TM_EARLY>0, sep = ",")),
  labels = c("X1.LS", "X4.HS", "X3.MRPP", "X2.MRSP", "NA", "NA")))
data.df <- main.df[main.df$`WES_BASELINE_AVAIL (N=72)`==1,]

plotCovariate (data.df = data.df,
                   sample.order = sample.order, 
                   metric = "ICB_group", 
                   colmap = setNames(c("gray","#4575b4","#35978f","#dfc27d","#a50026"), 
                                     sort(unique(main.df$ICB_group))), x.limit = x.limit, mp = mp,
                   label.name = "ICB group") 


screen(4)
main.df$CB_mod <- make.names(main.df[,"CR_PR_SD>6cycles_MARCO"])
data.df$CB_mod <- make.names(data.df[,"CR_PR_SD>6cycles_MARCO"])
plotCovariate (data.df = data.df,
               sample.order = sample.order, 
               metric = "CB_mod", 
               colmap = setNames(c("gray","#0038a8","#d73027"), sort(unique(main.df$CB_mod))), 
               x.limit = x.limit, mp = mp,
               label.name = "Clinical benefit")
screen(5)
plotCovariate (data.df = data.df,
               sample.order = sample.order, 
               metric = "BEST_OVERALL_RECIST", 
               colmap = c("CR"= "#d73027", 
                          "PR" = "#d73027",
                          "SD" = "#fdae61", 
                          "PD" = "#74add1", 
                          "NE" = "gray"), x.limit = x.limit, mp = mp,
               label.name = "Best response")

screen(6)
plotCovariate (data.df = data.df,
               sample.order = sample.order, 
               metric = "COHORT", 
               colmap = setNames(COLS_CANCER_COHORT,
                                 sort(unique(data.df$COHORT))), 
               x.limit = x.limit, 
               mp = mp,
               label.name = "Cohort")
screen(7)
#plotCovariate (data.df = data.df,
#               sample.order = sample.order, 
#               metric = "CANCER_CELL_TYPE", 
#               colmap = setNames(sort(COLS_CANCER_CELL_TYPES),
#                                 sort(unique(data.df$CANCER_CELL_TYPE))), x.limit = x.limit, 
#               mp = mp,label.name = "Cell type")

plotCovContinuous(data.df = sequenza.bl.df, 
                  sample.order = match (as.character(data.df$PATIENT_ID[sample.order]), as.character(sequenza.bl.df$pid)), 
                  metric = "cellularity",
                  mp = mp, x.limit = x.limit, label.name = "Cellularity")

#### Plot mutations
for (i in 8:42){
  screen(i)
  plotCovariate (data.df = collpased.mut.matrix,
                 sample.order = as.character(data.df$PATIENT_ID)[sample.order], 
                 metric = gene.list[i-7,1], 
                 colmap = COLS_MUTATIONS, x.limit = x.limit, 
                 mp = mp,label.name = gene.list[i-7,1], 
                 add.germline = collpased.gl.mut.matrix)
}
for (i in 43:47){
  screen (i)
  plotCovariate (data.df = collapsed.imm.mut.collection.matrix,
                 sample.order = as.character(data.df$PATIENT_ID)[sample.order], 
                 metric = names(idx)[i-42], 
                 colmap = c("altered" ="grey", "unaltered"="white"), x.limit = x.limit, 
                 mp = mp,label.name = names(idx)[i-42], 
                 add.germline = collpased.gl.imm.mut.collection.matrix)
}
#### Plot processes
idx  <- cumsum(unlist(lapply(imm.gene.list.collection, length))+1)+1
idx <- c(1, idx[1:4])+42
idx <- setNames(idx, names(imm.gene.list.collection))
for (list.idx in names(idx)){
  screen(idx[list.idx])
  plotCovariate (data.df = collapsed.imm.mut.collection.matrix,
                 sample.order = as.character(data.df$PATIENT_ID)[sample.order], 
                 metric = list.idx, 
                 colmap = c("altered" ="grey", "unaltered"="white"), x.limit = x.limit, 
                 mp = mp,label.name = list.idx, 
                 add.germline = collpased.gl.imm.mut.collection.matrix)
  for (gene.idx in seq(1, length(imm.gene.list.collection[[list.idx]]))){
    screen (idx[list.idx]+gene.idx)
    plotCovariate (data.df = collpased.imm.mut.matrix,
                   sample.order = as.character(data.df$PATIENT_ID)[sample.order], 
                   metric = imm.gene.list.collection[[list.idx]][gene.idx], 
                   colmap = COLS_MUTATIONS, x.limit = x.limit, 
                   mp = mp,label.name = imm.gene.list.collection[[list.idx]][gene.idx], 
                   add.germline = collpased.gl.imm.mut.matrix)
  }
}
 
plot.new()
legend (x = 0.05, y = 1, title = "Mutation signatures",
        legend = c("Aging",
                   "Alkylating",
                   "APOBEC",
                   "DSB",
                   "MMR",
                   "UV",
                   "Other"),
        fill = c("#0038a8","#73c2fb","#0067a5","#e0b0ff","#734f96","#ffdf00", "gray"),
        bty = "n", title.adj = 0)

legend (x = 0.05, y = 0.60, title = "ICB group",
        legend  = c("no delta-ctDNA", "LS","MSER", "MSPP" ,"HS"),
        fill = c("gray","#4575b4","#35978f","#dfc27d","#a50026"),
        bty = "n", title.adj = 0)

legend (x = 0.05, y = 0.30, title = "Clinical benefit",
        legend  = c("no data", "nCB","CB"),
        fill = c("gray","#0038a8","#d73027"),bty = "n",title.adj = 0)

legend (x = 0.35, y = 1, title = "Best response",
        legend  = rev(c("CR/PR", "SD", "PD", "no data")),
        fill = rev(c("#d73027","#fdae61","#74add1","gray")),bty = "n",title.adj = 0)

legend (x = 0.35, y =0.70, title = "Cohort",
        legend = c("HNSCC", "BRCA", "OV", "MM", "MST"),
        fill = COLS_CANCER_COHORT,bty = "n",title.adj = 0)
lg <- rep (NA, 10)
lg[c(1,5,10)] <- rev(c(0,0.5,1))
legend ( x = 0.35, y = 0.4, title = "Cellularity",
         legend = lg,
         fill = rev(colorRampPalette(brewer.pal(5,"Blues"))(10)), 
         border = NA, bty = "n", y.intersp = 0.50, cex = 1,title.adj = 0) 

legend (x = 0.60, y = 1, title = "Mutation type",
        legend = c("Missense", "Nonsense", "Splice site",
                   "Frameshift indel","In frame indel", 
                   "TSS","Complex", "No mutation", "Germline"),
        fill = c(unique(COLS_MUTATIONS), "black"),bty = "n", title.adj = 0)

