###############
# Figure 4 Main
# Figure4.R
###############
setwd ("/Users/yangc/PughLab/Repos/inspire-genomics/")
source (paste(getwd(),"R","PlotFunctions.R", sep = "/"))
#############
# Figure 4A #
#############
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
library (XLConnect)
library (plyr)
data.df <- readWorksheetFromFile (sheet = "Fig. 4A", file = "./data/Source Data/SourceData_Fig4/SourceData_Fig4.xlsx")
# See heatmap
#############
# Figure 4B #
#############
data.df <- readWorksheetFromFile (sheet = "Fig. 4B", file = "./data/Source Data/SourceData_Fig4/SourceData_Fig4.xlsx")
# See DE analysis
#############
# Figure 4C #
#############
data.df <- readWorksheetFromFile (sheet = "Fig. 4C", file = "./data/Source Data/SourceData_Fig4/SourceData_Fig4.xlsx")
data.df$groups <- gsub ("_BL","",data.df$groups)
data.df$groups <- factor (data.df$groups , levels = c("HS/CB","LS"))
stat1 <- plotGeneBox (main.df=data.df, GOI = "TIGIT", x.cat.order = levels (data.df$groups), y.limit = c(-1, 6))
stat2 <- plotGeneBox (main.df=data.df, GOI = "PDCD1", x.cat.order = levels (data.df$groups), y.limit = c(-1, 6))
stat3 <-plotGeneBox (main.df=data.df, GOI = "CTLA4", x.cat.order = levels (data.df$groups), y.limit = c(-2, 6))
#############
# Figure 4D #
#############
data.df <- readWorksheetFromFile (sheet = "Fig. 4D", file = "./data/Source Data/SourceData_Fig4/SourceData_Fig4.xlsx")
data.df$groups <- factor (data.df$groups, levels = data.df$groups[!(duplicated(data.df$groups))])
data.df$COHORT <- factor (data.df$COHORT, levels = c("HS/CB","LS"))
plotGeneBoxWithPairs (main.df = data.df, GOI = "PLA2G2D", x.var="groups", 
                      y.var ="log2TPM", y.grp = "COHORT",
                      x.var.name = "", y.var.name = "log2(TMP+1)",
                      x.cat.order=levels(data.df$groups), 
                      group.colors=rev(COLS_CB)[sort(rep(1:2,2))],
                      y.limit = c(0,12),BL="BL",TX="C2_3")
stat4 <- wilcox.test (log2TPM~Timepoint, subset(data.df,COHORT=="HS/CB"), alternative = "two.sided")
stat5 <- wilcox.test (log2TPM~Timepoint, subset(data.df,COHORT=="LS"), alternative = "two.sided")
