###############
# Figure 2 Main
# Figure2.R
###############
setwd ("/Users/yangc/PughLab/Repos/inspire-genomics/")
load (paste(getwd(),"R","PlotFunctions.R", sep = "/"))
#############
# Figure 2A #
#############
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
library (XLConnect)
library (plyr)
data.df <- readWorksheetFromFile (sheet = "Fig. 2A", file = "./data/Source Data/SourceData_Fig2.xlsx")

#############
# Figure 2B #
#############
data.df <- readWorksheetFromFile (sheet = "Fig. 2B", file = "./data/Source Data/SourceData_Fig2.xlsx")
data.df$BRCA2_mutated <- ifelse (data.df[,"BRCA2_MUT_status"] %in% c("germline", "somatic"), no = "none", yes = "mutated")
data.df$BRCA2_mutated <- factor (data.df$BRCA2_mutated, levels = c("none", "mutated"))
data.df <- arrange (data.df, BRCA2_mutated, TMB_n72)

plotGroupsWithScatter (main.df = data.df, 
                       x.var = "BRCA2_mutated", 
                       y.var = "TMB_n72", 
                       x.var.name = "", y.var.name = "TMB (mut/Mb)",
                       x.cat.order = c("none","mutated"), group.colors = COLS_CB, 
                       y.limit = c(0.1, 1000), log.scale = "y", 
                       pch = 16, lwd = 2.5, wex = 0.8, 
                       custom.axis = TRUE, custom.at = c(0.1, 1, 10, 100, 1000), 
                       custom.lbls = c(0.1, 1, 10, 100, 1000))

wilcox.test (TMB_n72 ~ BRCA2_mutated, data.df, alternate = "two.sided")

#############
# Figure 2C #
#############
data.df <- readWorksheetFromFile (sheet = "Fig. 2C", file = "./data/Source Data/SourceData_Fig2.xlsx")
data.df$BRCA2_mutated <- ifelse (data.df[,"BRCA2_MUT_status"] %in% c("germline", "somatic"), no = "none", yes = "mutated")
data.df$BRCA2_mutated <- factor (data.df$BRCA2_mutated, levels = c("none", "mutated"))
data.df <- arrange (data.df, BRCA2_mutated, PDL1_IHC)
data.df <- data.df[!is.na(data.df$PDL1_IHC),] # remove one instance of NA

plotGroupsWithScatter (main.df = data.df, 
                       x.var = "BRCA2_mutated", 
                       y.var = "PDL1_IHC", 
                       x.var.name = "", y.var.name = "PDL1_IHC MPS (%)",
                       x.cat.order = c("none","mutated"), group.colors = COLS_CB, 
                       y.limit = c(0,100), log.scale = "", 
                       pch = 16, lwd = 2.5, wex = 0.8, 
                       custom.axis = FALSE)

wilcox.test (PDL1_IHC ~ BRCA2_mutated, data.df, alternate = "two.sided")

#############
# Figure 2D #
#############
data.df <- readWorksheetFromFile (sheet = "Fig. 2D", file = "./data/Source Data/SourceData_Fig2.xlsx")
count.loh <- table (data.df$Pembrolizumab_sensitivity_groups_n71, data.df$B2M_LOH)
fraction.loh <- prop.table (count.loh, margin = 1)
mp <- barplot (fraction.loh[,2], ylim = c(0,1), las = 1, col = COLS_PEM4, border=NA)
title (ylab = "Faction of samples with B2M LOH", cex.lab = 1) 
axis (side = 3, at = c(mp[1]-0.5,mp), labels = c("n =", rowSums (count.loh)), lty = 0)
text (x = mp, y = fraction.loh[,2]+0.1, labels = count.loh[,2], cex = 0.8)
abline (h = 0, lty = 1.2)  
