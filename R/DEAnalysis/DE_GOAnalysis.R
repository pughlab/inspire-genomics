####################
# GSEA GO BP on INSPIRE DE genes
# Pre-post treatment
###################

############
# VARIABLES
#############
source ("/Volumes/h4h/pughlab/projects/INSPIRE_ms1/RNAseq/R/runGSEA.R")
gmt.dir <- "~/Desktop/INSPIRE_ms2/data/"
gmt.files <- list.files (gmt.dir,pattern = ".gmt",full.names = TRUE)  
##############
# FUNCTIONS
##############
# returns GSEA results object (contains @Results and @Plot)
runGSEA <- function (rnks, gmt.idx = 2, pval = 0.10, 
                     output.file = "~/Desktop/INSPIRE_ms2/data/INSPIRE.DESeq.Paired.Induced.GSEA.results.tsv") {
  # remove any entry that is duplicated
  duplicated.genes <- unique(names(genelist.rk)[duplicated(names(genelist.rk))])
  # Do GSEA analysis on GO Biological Processes
  INSPIRE.gsea <- GSEA (gene_list = genelist.rk[!names(genelist.rk) %in% duplicated.genes], 
                        GO_file = gmt.files[gmt.idx],
                        pval = pval)
  # Write table of results to file
  write.table (INSPIRE.gsea$Results[,c(1:7,9)], sep = "\t", quote = FALSE, row.names = F, 
               file = output.file)
  return (INSPIRE.gsea)
}
# returns data frame of selected genes passing DE threshold
# if labels.df is provided (with modified positions for labels, the locations are plotted as provided)
plotVolcano <- function (df=r, addGeneNames = FALSE, labels.df=NA, threshold = 0,
                         xlimit = c(-35,35), ylimit = c(0,3),
                         main.col = "deepskyblue2", highlight.col = "red", af = 1){
  # Plot volcano plot based on DESeq results with genes labeled
  # x = log2FoldChange
  # y = -log10(P adj)
  colnames (r)[colnames (r)=="hugo.name"] <- "hugoname"
  r$log10Padj <- log10(r$padj)/-1
  plot (log10Padj~log2FoldChange, data = r, las = 1,
        xlim =xlimit , ylim = ylimit,
        xlab = expression ("log"[2]*"(response/progression)"),
        ylab = expression ('-log'[10]*"(adjusted "*italic("P")*" value)"),
        pch = 1, col = adjustcolor(col = main.col, alpha = 0.6))
  abline (v=0)
  #abline (v=c(-2,2), col = "grey")
  #abline (h=1, col = "grey", lty = 2)
  
  # label significant genes in highlighting color
  genes.to.label <- r[!is.na(r$log10Padj) & r$log10Padj>1,c("log2FoldChange", "log10Padj", "hugoname")]
  genes.to.label <- genes.to.label[order(genes.to.label$log10Padj, decreasing = TRUE),]
  if (length(highlight.col)==1){
    points (log10Padj~log2FoldChange, data = genes.to.label, pch = 16, col = highlight.col)
    legend ("topleft", pch = c(16,1), col = c(highlight.col,main.col), legend = c("Differentially Regulated","None"), cex = 0.80, bty = "n")
  }else{
    
    if (length(threshold)==2){
      points (log10Padj~log2FoldChange, data = genes.to.label[genes.to.label$log2FoldChange>= max(threshold),], pch = 16, col = adjustcolor(highlight.col[1], alpha.f = af))
      points (log10Padj~log2FoldChange, data = genes.to.label[genes.to.label$log2FoldChange<= min(threshold),], pch = 16, col = adjustcolor(highlight.col[2], alpha.f = af))
      pt.counts <- c(sum(genes.to.label$log2FoldChange >= max(threshold)), sum(genes.to.label$log2FoldChange <= min(threshold)))
      abline (v = threshold, lty = 2)
    }else{
      points (log10Padj~log2FoldChange, data = genes.to.label[genes.to.label$log2FoldChange>0,], pch = 16, col = adjustcolor(highlight.col[1], alpha.f = af))
      points (log10Padj~log2FoldChange, data = genes.to.label[genes.to.label$log2FoldChange<0,], pch = 16, col = adjustcolor(highlight.col[2], alpha.f = af))
      pt.counts <- c(sum(genes.to.label$log2FoldChange>0), sum(genes.to.label$log2FoldChange<0))
    }
   
 
    legend.labels <- paste0(c("Induced", "Suppressed"), " (", pt.counts, ")")
    legend ("topleft", pch = c(16,16,1), col = c(highlight.col,main.col), legend = c(legend.labels,"None"), cex = 0.80, bty = "n")
  }
  
  
  if (addGeneNames){
    if (is.na(labels.df)){
      # add gene names for DE genes above padj = 0.10 threshold
      genes.to.label$new.x <- ifelse (genes.to.label$log2FoldChange>0, genes.to.label$log2FoldChange + 3, genes.to.label$log2FoldChange -3)
      genes.to.label$new.y <- ifelse (genes.to.label$log10Padj>1.07, genes.to.label$log10Padj + 0.25, genes.to.label$log10Padj - 0.25)
      
      # custom label for "all"
      #genes.to.label[genes.to.label$hugoname=="HMGA2","new.y"] <- genes.to.label[genes.to.label$hugoname=="HMGA2","log10Padj"] + 0.15
      #genes.to.label[genes.to.label$hugoname=="IGSF5","new.y"] <- genes.to.label[genes.to.label$hugoname=="IGSF5","log10Padj"]
    }else{
      genes.to.label <- labels.df
    }
    
    text (x=genes.to.label$new.x, y = genes.to.label$new.y, genes.to.label$hugoname, cex = 0.75, col = "black", font = 3)
    segments(x0=genes.to.label$log2FoldChange,
             y0=genes.to.label$log10Padj,
             x1=genes.to.label$new.x,
             y1=genes.to.label$new.y, col = "grey")
  }
  return (genes.to.label)
}

###########
# MAIN
###########
DE.seq.rk.rds <- list(all = "/Volumes/h4h/pughlab/projects/INSPIRE_ms1/RNAseq/DESeq2/DESeq2.DE.genelist.rnk.sep.RData",
                      cb = "/Volumes/h4h/pughlab/projects/INSPIRE_ms1/RNAseq/DESeq2/INSPIRE_CB_DESeq2.DE.genelist.rnk.sep.Rdata",
                      ncb = "/Volumes/h4h/pughlab/projects/INSPIRE_ms1/RNAseq/DESeq2/INSPIRE_NCB_DESeq2.DE.genelist.rnk.sep.Rdata")

DE.seq.rk.rds <- list(all = "/Volumes/h4h/pughlab/projects/INSPIRE_ms1/RNAseq/DESeq2/INSPIRE_CTDNA__DESeq2.DE.genelist.rnk.sep.Rdata",
                      cb = "/Volumes/h4h/pughlab/projects/INSPIRE_ms1/RNAseq/DESeq2/INSPIRE_CB_CTDNA__DESeq2.DE.genelist.rnk.sep.Rdata",
                      ncb = "/Volumes/h4h/pughlab/projects/INSPIRE_ms1/RNAseq/DESeq2/INSPIRE_LS_CTDNA__DESeq2.DE.genelist.rnk.sep.Rdata")
par(mfrow = c(1,2))
load (DE.seq.rk.rds[[2]])
r <- df
gene.labels1 <-plotVolcano (df = r, addGeneNames = FALSE, labels.df = NA,xlimit = c(-15,15), ylimit = c(0,6),
               main.col = "black", highlight.col = c("red","deepskyblue2"), af = 0.30)
r1 <- r
load (DE.seq.rk.rds[[3]])
r <- df
gene.labels2 <-plotVolcano (df = r, addGeneNames = FALSE, labels.df = NA,xlimit = c(-7,7), ylimit = c(0,6),
               main.col = "black", highlight.col = c("red","deepskyblue2"), af = 0.30)
r2 <- r

load (DE.seq.rk.rds[[3]])
r <- df
gene.labels3 <-plotVolcano (df = r, addGeneNames = FALSE, labels.df = NA, xlimit = c(-7,7), ylimit = c(0,8),
                            main.col = "black", highlight.col = c("red","deepskyblue2"), af = 0.30,
                            threshold = c(-1,1))


# Plot log2(foldChange) significantly regulated in CB and nCB patients
genes.intersect <- intersect(row.names(gene.labels1), row.names(gene.labels2))
genes.all <- unique(c(row.names(gene.labels1), row.names(gene.labels2)))
plot.dat <- data.frame(
  CBlog2FoldChange = r1[genes.all, "log2FoldChange"],
  NCBlog2FoldChange = r2[genes.all, "log2FoldChange"],
  CBpadj = r1[genes.all, "padj"],
  NCBpadj = r2[genes.all, "padj"],
  stringsAsFactors = F,
  row.names = genes.all
)
plot.all.cols <- ifelse (genes.all %in% genes.intersect,
                         yes = "red",
                         no = ifelse( genes.all %in% row.names(gene.labels1), 
                                      yes = "tan1",
                                      no = "white"))
                                 
par(mfrow = c(1,2))
plot (CBlog2FoldChange~NCBlog2FoldChange, plot.dat, 
      pch = 21, col = "gray30", bg = adjustcolor(plot.all.cols, alpha = 0.60), las = 1,
      xlab = expression ("Week 9: log"[2]*" FC LS group"),
      ylab = expression ("Week 9: log"[2]*" FC HS/CB gorup"))
points (CBlog2FoldChange~NCBlog2FoldChange, plot.dat[genes.intersect,],
        pch = 16, col = "red")
abline (h = 0, v = 0, lty = 2, col = "gray50")
abline (a = 0, b = 1, lty = 4, col = "gray30")

legend ("topright", pch = 21, col = "black", pt.bg = c("red", "tan1", "white"), 
        title = "DRGs (# Genes)",
        legend = paste0(c("Both", "HS/CB", "LS"), "(", as.numeric(table(plot.all.cols)), ")"), 
                       cex = 0.80, bty = "n")

# draw rectangle for zoom panel
rect (xleft = -2, 
      ybottom = -5, 
      xright = 3, 
      ytop = 5)

plot (CBlog2FoldChange~NCBlog2FoldChange, plot.dat, 
      pch = 21, col = "gray88", 
      bg = adjustcolor(ifelse(plot.all.cols=="red","red", "white"), alpha = 1), 
      ylim = c(-5,5), xlim = c(-2,3), las = 1,
      xlab = expression ("Week 9: log"[2]*" FC LS group"),
      ylab = expression ("Week 9: log"[2]*" FC HS/CB group"))
points (CBlog2FoldChange~NCBlog2FoldChange, plot.dat[genes.intersect,],
        pch = 16, col = "red")
abline (h = 0, v = 0, lty = 2, col = "gray50")
abline (a = 0, b = 1, lty = 4, col = "gray30")
#abline (a = 3, b = -1, lty = 4, col = "gray30")
# plot gene annotations 
dist_point_line <- function(x1, y1, slope, intercept) {
  a = c(x1,y1)
  b = c(1, intercept+slope)
  c = c(-intercept/slope,0)       
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  print(m)
  return(det(m)/sqrt(sum(v1*v1)))
}

genes.to.label <- plot.dat[genes.intersect,]
genes.to.label$hugo.name <- r1[row.names(genes.to.label), "hugo.name"]
genes.to.label$dist <- sapply (1:nrow(genes.to.label), function (i) dist_point_line(genes.to.label[i,2], genes.to.label[i,1], 1, 0))
genes.to.label$dist2 <- sapply (1:nrow(genes.to.label), function (i) dist_point_line(genes.to.label[i,2], genes.to.label[i,1], -1, 3))
genes.to.label$new.x <- ifelse (genes.to.label$dist>0, genes.to.label[,2] - 1, genes.to.label[,2] + 1)
genes.to.label$new.y <- ifelse (genes.to.label$dist>0, genes.to.label[,1] + 1, genes.to.label[,1] - 1)

segments(x0=genes.to.label[,2],
         y0=genes.to.label[,1],
         x1=genes.to.label$new.x,
         y1=genes.to.label$new.y, col = "gray20")

text (x=genes.to.label$new.x, y = genes.to.label$new.y, genes.to.label$hugo.name, cex = 0.75, col = "black", font = 3)
points (CBlog2FoldChange~NCBlog2FoldChange, plot.dat[genes.intersect,],
        pch = 16, col = "red")

plot (CBlog2FoldChange~NCBlog2FoldChange, plot.dat, 
      pch = 21, col = "gray88", las = 1,
      #ylim = c(-5,5), xlim = c(-2,3), las = 1,
      xlab = expression ("Week 9: log"[2]*" FC No Benefit"),
      ylab = expression ("Week 9: log"[2]*" FC Clinical Benefit"))
abline (h = 0, v = 0, lty = 2, col = "gray50")
abline (a = 0, b = 1, lty = 4, col = "gray30")
#abline (a = 3, b = -1, lty = 4, col = "gray30")
# plot gene annotations 
genes.de <- row.names(r1)[r1$hugo.name %in% c("TNRC6C", "ZFAT", "NCOA2", "LCE2B", "C10orf99", "SPRR2B",
                                  "GSG2", "HMGA2", "IGSF5", "CTAG1B", "GTSE1")]
genes.to.label <- plot.dat[genes.de,]
genes.to.label$hugo.name <- r1[row.names(genes.to.label), "hugo.name"]
genes.to.label$dist <- sapply (1:11, function (i) dist_point_line(genes.to.label[i,2], genes.to.label[i,1], 1, 0))
# only label those outside of the zoomed window
genes.to.label <- genes.to.label[abs(genes.to.label$dist) > 5,] 
genes.to.label$new.x <- ifelse (genes.to.label$dist>0, genes.to.label[,2] + 1, genes.to.label[,2] - 1)
genes.to.label$new.y <- ifelse (genes.to.label$dist>0, genes.to.label[,1] + 1, genes.to.label[,1] - 1)

segments(x0=genes.to.label[,2],
         y0=genes.to.label[,1],
         x1=genes.to.label$new.x,
         y1=genes.to.label$new.y, col = "gray50")

text (x=genes.to.label$new.x, y = genes.to.label$new.y, genes.to.label$hugo.name, cex = 0.75, col = "black", font = 3)
points (CBlog2FoldChange~NCBlog2FoldChange, plot.dat[genes.de,],
        pch = 16, col = "blue")
# draw rectangle for zoom panel
rect (xleft = -2, 
      ybottom = -5, 
      xright = 3, 
      ytop = 5)

plot (CBlog2FoldChange~NCBlog2FoldChange, plot.dat, 
      pch = 21, col = "gray88", 
      ylim = c(-5,5), xlim = c(-2,3), las = 1,
      xlab = expression ("Week 9: log"[2]*" FC No Benefit"),
      ylab = expression ("Week 9: log"[2]*" FC Clinical Benefit"))
abline (h = 0, v = 0, lty = 2, col = "gray50")
abline (a = 0, b = 1, lty = 4, col = "gray30")
#abline (a = 3, b = -1, lty = 4, col = "gray30")
# plot gene annotations 
genes.to.label <- plot.dat[genes.de,]
genes.to.label$hugo.name <- r1[row.names(genes.to.label), "hugo.name"]
genes.to.label$dist <- sapply (1:11, function (i) dist_point_line(genes.to.label[i,2], genes.to.label[i,1], 1, 0))
genes.to.label$new.x <- ifelse (genes.to.label$dist>0, genes.to.label[,2] + 1, genes.to.label[,2] - 1)
genes.to.label$new.y <- ifelse (genes.to.label$dist>0, genes.to.label[,1] + 1, genes.to.label[,1] - 1)

segments(x0=genes.to.label[,2],
         y0=genes.to.label[,1],
         x1=genes.to.label$new.x,
         y1=genes.to.label$new.y, col = "gray50")

text (x=genes.to.label$new.x, y = genes.to.label$new.y, genes.to.label$hugo.name, cex = 0.75, col = "black", font = 3)
points (CBlog2FoldChange~NCBlog2FoldChange, plot.dat[genes.de,],
        pch = 16, col = "blue")

# Analyze genes that are discordantly regulated by pembrolizumab in CB and NCB
plot.dat$hugo.name <- r1[row.names(plot.dat), "hugo.name"]
up.cb <-plot.dat[plot.dat$CBlog2FoldChange>0 & plot.dat$NCBlog2FoldChange<0,]
up.ncb <- plot.dat[plot.dat$CBlog2FoldChange<0 & plot.dat$NCBlog2FoldChange>0,]
# run GSEA GOBP and see if there are any enriched processes
# nothing statistically significant
gsea.up.cb <- GSEA(sort(setNames(as.numeric(up.cb[,1]-up.cb[,2]), up.cb[,3]), decreasing = T),
                    gmt.files[2],
                    1)
gsea.up.ncb <- GSEA(sort(setNames(as.numeric(up.ncb[,2]-up.ncb[,1]), up.ncb[,3]), decreasing = T),
                    gmt.files[2],
                    1)


# Plot genes differentially regulated by anti-PD1 treatment in responder vs non-responders
load ("~/PughLab/INSPIRE/Manuscript2/ImmuneScoreAnalysis/INSPIRE_Immune_Analysis.data.RData")
# Subset data into list of 3 matrices 
expr.dat.plot <- list(
  baseline = expr.sub.dat.full[as.character (genes.to.label$hugoname),row.names(subset(DE.design.full, Treatment == 0))],
  treated = expr.sub.dat.full[as.character (genes.to.label$hugoname),row.names(subset(DE.design.full, Treatment == 1))]
  )
expr.dat.plot[["delta"]] <- expr.dat.plot[["treated"]] - expr.dat.plot[["baseline"]]

# plot each gene 
source("~/Desktop/INSPIRE_ms2/src/Color.Palettes.R")
plotBoxScatter <- function (dat, xcat, ycat, ylimit, grp.order, col.grps, p.val=NA, 
                            xlabel, ylabel, cex.pts  = 0.5, cex.lbs = 0.5){
  
  dat[,xcat] <- factor (dat[,xcat], levels = grp.order)
  print("setting factor")
  print(str(dat))
  # plot box between grps
  box.stats <- boxplot (formula (paste(ycat,xcat,sep="~")), data = dat, 
                        xlab = xlabel, ylab = ylabel, plot = FALSE)
  # Boxplot plot
  boxplot (formula (paste (ycat, xcat, sep = "~")), data = dat, xaxt = "n",
           xlab = xlabel, ylab = "", 
           plot = TRUE, boxwex = 0.50, boxlty = 0, whisklty = 1, 
           staplelwd = 3, medcol = "black", frame.plot = TRUE,
           las = 1, border = "gray", col = col.grps,
           boxlwd = 2, outline = FALSE, ylim = ylimit, 
           cex.axis = cex.lbs)
  axis (1, at = seq(1,length(grp.order)), labels = paste0(box.stats$names,"\n(",box.stats$n,")"), mgp = c(3,1.5,0))
  # plot scatter 
  # fill color (bg) is applied row wise which makes it impossible to plot by group
  # solution is to plot each strip separately for each group
  # for (i in seq(1,length(grp.order))){
  #   stripchart (dat[dat[,xcat] == grp.order[i],ycat] ~ dat[dat[,xcat] == grp.order[i],xcat], 
  #               col = "gray", 
  #               bg = col.grps[i], method = "jitter", jitter = 0.02,
  #               pch = 21, las = 1, vertical = T, cex = 1.2,
  #               add = T)
  # }
  stripchart(
    formula (paste (ycat, xcat, sep = "~")), data = dat,
    col = "black", pch = "o", las = 2, cex = 0.8, vertical = T, add = T
  )

  # Add DESeq2 uncorrected p-val
  if (!is.na(p.val)){
    legend ("topright", legend = bquote(italic("P") == .(formatC(p.val, format = "e", digits = 2))), cex = 0.90, bty = "n")
  }
  
  # Add gene name
    #legend ("bottom", legend = bquote (italic(.(ycat))), inset = c(0,1), xjust = 0.5, xpd = TRUE, horiz = TRUE, bty = "n", cex = 1.2)
    title (main = bquote (bolditalic(.(ycat))), line = 1)
    #title (ylab = ylabel, line = 2)
  }
plotGeneExprGrp <- function (expr.dat, sample.sheet,id.name, xcat, ycat, ylimit, 
                             grp.order, col.grps, p.val, xlabel, ylabel, cex.pts  = 0.5, cex.lbs = 0.5){
  plot.dat <- data.frame(
      V1 = colnames(expr.dat),
      V2 = as.character(sample.sheet[match(colnames(expr.dat), sample.sheet[,id.name]), xcat]),
      V3 = as.numeric(expr.dat[ycat,])
    , 
    stringsAsFactors = FALSE
  )
  
  colnames(plot.dat) <- c(id.name, xcat, ycat)
  
  print (head(plot.dat, n = 5))
  print(str(plot.dat))
  # plot scatter box plot
  plotBoxScatter (dat = plot.dat, 
                  xcat = xcat, ycat = ycat, ylimit = ylimit, grp.order = grp.order, 
                  col.grps = col.grps, p.val=p.val, 
                  xlabel = xlabel, ylabel = ylabel, cex.pts  = cex.pts, cex.lbs = cex.lbs)
  
}
DE.design.full$sid <- row.names(DE.design.full)

par(mfrow = c(2,3))
par(mar = c(3,2.5,2,0.5))
for (gene in as.character(genes.to.label[genes.to.label$log2FoldChange>0,]$hugoname)){
  plotGeneExprGrp(expr.dat = expr.dat.plot[["delta"]], 
                  ycat = gene,
                  sample.sheet = DE.design.full, 
                  id.name = "sid",
                  xcat = "Response",
                  ylimit = 1.05*range(expr.dat.plot[["delta"]][gene,]),
                  grp.order = c("Y","N"),
                  col.grps = c(COLS_CB[1], "gray88"),
                  p.val = r[match(gene, r$hugoname), "pvalue"],
                  xlabel = "Clincal Benefit", 
                  ylabel = "normalized expression",
                  cex.pts = 0.90, 
                  cex.lbs = 0.90)
}

par(mfrow = c(2,3))
par(mar = c(3,2.5,2,0.5))
for (gene in as.character(genes.to.label[genes.to.label$log2FoldChange<0,]$hugoname)){
  plotGeneExprGrp(expr.dat = expr.dat.plot[["delta"]], 
                  ycat = gene,
                  sample.sheet = DE.design.full, 
                  id.name = "sid",
                  xcat = "Response",
                  ylimit = 1.50*range(expr.dat.plot[["delta"]][gene,]),
                  grp.order = c("Y","N"),
                  col.grps = c(COLS_CB[1], "gray88"),
                  p.val = r[match(gene, r$hugoname), "pvalue"],
                  xlabel = "Clincal Benefit", 
                  ylabel = "normalized expression",
                  cex.pts = 0.90, 
                  cex.lbs = 0.90)
}


# plot GSEA enrichment scores (top10 up and top 10 down)
ssgsea.results <- read.csv("~/Desktop/INSPIRE_ms2/data/INSPIRE.DESeq.Paired.Induced.GSEA.results.tsv",
                           sep = "\t", stringsAsFactors = F)
filtRes <- rbind(head(ssgsea.results, n = 10),
                tail(ssgsea.results, n = 10 ))
g <- ggplot(filtRes, aes(reorder(GO_ID, NES), NES)) +
  geom_segment( aes(reorder(GO_ID, NES), xend=GO_ID, y=0, yend=NES)) +
  geom_point( size=5, aes( fill = Enrichment),
              shape=21, stroke=2) +
  scale_fill_manual(values = c("Down-regulated" = "dodgerblue",
                               "Up-regulated" = "firebrick") ) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="GSEA - Top10 Regulated GO Biological Processes") + 
  theme_minimal()
g

# plot overlap of leadingEdgeGenes in Upregulated (Induced) vs Downregulated (Suppressed) BPs
library(pheatmap)
genes.up.unique <- unique(unlist(INSPIRE.gsea$Results[INSPIRE.gsea$Results$Enrichment=="Up-regulated","leadingEdge"]))
# plot smoothed frequency of genes in leading edge of 27 processes that are enriched in responder vs non-responders
plot (density(sort(table(unlist(INSPIRE.gsea$Results[INSPIRE.gsea$Results$Enrichment=="Up-regulated","leadingEdge"])), decreasing = T)),
      main = "Frequency of genes appearing in enriched process", las = 1)
abline (v = 1:20, col = "grey")
genes.up.matrix <- do.call("rbind", lapply (INSPIRE.gsea$Results[INSPIRE.gsea$Results$Enrichment=="Up-regulated","leadingEdge"],
                                            function (x) match(x, genes.up.unique)))
row.names(genes.up.matrix) <- ssgsea.results[ssgsea.results$Enrichment=="Up-regulated","GO_ID"]
colnames(genes.up.matrix) <- genes.up.unique   
pheatmap (genes.up.matrix, border_color = NA, cutree_rows = 5)
# There are 110 genes that 