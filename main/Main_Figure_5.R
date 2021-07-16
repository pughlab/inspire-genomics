###############
# Figure 5 Main
# Figure5.R
###############
setwd ("/Users/yangc/PughLab/Repos/inspire-genomics/")
source (paste(getwd(),"R","PlotFunctions.R", sep = "/"))
#############
# Figure 5A #
#############
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
library (XLConnect)
library (plyr)
data.df <- readWorksheetFromFile (sheet = "Fig. 5A", file = "./data/Source Data/SourceData_Fig5.xlsx")
# See GO enrichment from DE
#############
# Figure 5B #
#############
data.df <- readWorksheetFromFile (sheet = "Fig. 5B", file = "./data/Source Data/SourceData_Fig5.xlsx")
data.df$signature <- gsub ("delta_abs_","",data.df$signature)
data.df$Pembro_CB_groups <- gsub ("^.. ","",data.df$Pembro_CB_groups)
data.df$signature <- factor (data.df$signature, levels = c("B.cells.naive","B.cells.memory","Plasma.cells",
                                                           "T.cells.CD8", "T.cells.CD4.memory.resting",
                                                           "T.cells.CD4.memory.activated","Tfh","Tregs","Tdg",
                                                           "NK.cells.activated","Macrophages.M0",
                                                           "Macrophages.M1","Macrophages.M2"))

data.df$Pembro_CB_groups <- factor (data.df$Pembro_CB_groups, levels = c("HS/CB","LS"))
data.df <- arrange(data.df, signature, Pembro_CB_groups, score)
data.df$x.combined <- paste (data.df$signature, data.df$Pembro_CB_groups, sep=",")
data.df$x.combined <- factor(data.df$x.combined, levels = data.df$x.combined[!duplicated(data.df$x.combined)])
p.vals <- unlist (lapply (levels(data.df$signature), 
                          function (x) wilcox.test (score~Pembro_CB_groups,
                                                    subset(data.df,signature==x), 
                                                    alternative = "two.sided")$p.value))
new.x <- unlist(lapply (seq (from = 1, by =3, length.out = length(levels(data.df$signature))), 
                        function (x) seq (from = x - 0.55, to = x + 0.55, length.out = 2)))
boxplot (score~x.combined, data.df, at = new.x, las = 2, ylim = c(-5,20), lty = 1,
         col=rev(COLS_CB_BOX), whiskcol="black", medcol ="black", border = NA,
         xaxt = "n", outline = F, cex.axis = 0.8, ylab = "Delta CIBERSORT score")
stripchart(score~x.combined, data.df, at = new.x, add = T, pch ="o", cex = 0.5, 
           vertical = T, method = "jitter", jitter = 0.2)
axis (side = 1, las = 2, cex.axis = 0.5, labels = levels(data.df$signature),
      at = seq (from = 1, by =3, length.out = length(levels(data.df$signature))))
axis (side = 3,labels = FALSE,
      at = seq (from = 1, by =3, length.out = length(levels(data.df$signature))))
mtext (side = 3, text = paste("p=",round(p.vals,2)), padj = -2, cex = 0.50,
      at = seq (from = 1, by =3, length.out = length(levels(data.df$signature))),
      col= ifelse(p.vals<0.10,"red","black"))
legend("topleft", fill = rev(COLS_CB_BOX), legend = c("HS/CB (n=11)", "LS (n=11)"), cex = 0.8)
#############
# Figure 5C #
#############
data.df <- readWorksheetFromFile (sheet = "Fig. 5C", file = "./data/Source Data/SourceData_Fig5.xlsx")
data.df$cat <- ifelse(data.df$score>0, yes = "2. increase from baseline", no = "1. decrease from baseline")
data.df$PFSTIME_Months_since_C3 <- as.numeric(data.df$PFSTIME_Months_since_C3)
data.df$PFSevent <- as.numeric(data.df$PFSevent)
survKM (dat.df=data.df, 
        ycat="PFS", 
        ycat.time.label="PFSTIME_Months_since_C3", 
        ycat.event.label="PFSevent",
        grp.cat="cat", 
        grp.levels=c("1. decrease from baseline","2. increase from baseline"), 
        withlegend = TRUE, legend.loc = "bottomleft",
        withPval = TRUE, pval.loc = "topright",
        add.covariate = FALSE, covariate.cat, cov.levels,
        grp.cols = c("black", "red"))

