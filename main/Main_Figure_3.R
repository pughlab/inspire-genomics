###############
# Figure 3 Main
# Figure3.R
###############
setwd ("/Users/yangc/PughLab/Repos/inspire-genomics/")
load (paste(getwd(),"R","PlotFunctions.R", sep = "/"))
load (paste(getwd(),"R","ImmuneScoreForestPlotFunctions.R", sep = "/"))
load (paste(getwd(),"R","ColorScheme.R", sep = "/"))
#############
# Figure 3A #
#############
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
library (XLConnect)
library(dplyr)
library (tidyr)
data.df <- readWorksheetFromFile (sheet = "Fig. 3A", file = "./data/Source Data/SourceData_Fig3.xlsx")
grp.var <- "Pembrolizumab_sensitivity_groups_n71"

data.df[,grp.var] <- gsub ("^.. ", "", data.df[,grp.var])
data.df[,grp.var] <- factor (data.df[,grp.var], levels = c("LS", "MSER", "MSPP", "HS"))
data.df$signature <- factor (data.df$signature, levels = c(BL,TX))
data.df <- data.df %>% arrange_(c(grp.var, "signature", "score"))
data.df$combined.x <- paste(data.df[,grp.var],data.df$signature, sep = "_")

plotGroupsWithPairs (main.df = data.df, 
                     x.var = "combined.x", 
                     y.var = "score", 
                     x.var.name = "", 
                     y.var.name = "IM Score",
                     y.grp = grp.var,
                     y.limit = c(0.5, 2),
                     x.cat.order = data.df$combined.x[!duplicated(data.df$combined.x)], 
                     group.colors = COLS_PEM4[sort(rep(1:4,2))],
                     BL = "BL_IMS", TX = "TX_IMS")
# Wilcox tests
stats.out <- lapply (levels(data.df[,grp.var]), 
                     function (x) wilcox.test (score~signature,data.df[data.df[,grp.var]==x,], alternate = "two.sided"))
#############
# Figure 3B #
#############
data.df <- readWorksheetFromFile (sheet = "Fig. 3B", file = "./data/Source Data/SourceData_Fig3.xlsx")
grp.var <- "Pembro_CB_groups"

data.df[,grp.var] <- make.names(gsub ("^.. ", "", data.df[,grp.var]))
data.df[,grp.var] <- factor (data.df[,grp.var], levels = c("LS", "HS.CB"))
data.df$signature <- factor (data.df$signature, levels = c(BL,TX))
data.df <- data.df %>% arrange_(c(grp.var, "signature", "score"))
data.df$combined.x <- paste(data.df[,grp.var],data.df$signature, sep = "_")

plotGroupsWithPairs (main.df = data.df, 
                     x.var = "combined.x", 
                     y.var = "score", 
                     x.var.name = "", 
                     y.var.name = "IM Score",
                     y.grp = grp.var,
                     y.limit = c(0.5, 2),
                     x.cat.order = data.df$combined.x[!duplicated(data.df$combined.x)], 
                     group.colors = COLS_CB[sort(rep(1:2,2))],
                     BL = "BL_IMS", TX = "TX_IMS")
# Wilcox tests
stats.out <- lapply (levels(data.df[,grp.var]), 
                     function (x) wilcox.test (score~signature,data.df[data.df[,grp.var]==x,], alternate = "two.sided"))
stats1.out <- lapply (levels(data.df[,"signature"]), 
                     function (x) wilcox.test (as.formula(paste("score",grp.var,sep="~")),data.df[data.df[,"signature"]==x,], alternate = "two.sided"))
#############
# Figure 3C #
#############
data.df <- readWorksheetFromFile (sheet = "Fig. 3C", file = "./data/Source Data/SourceData_Fig3.xlsx")
grp.var <- "Pembro_CB_groups"

data.df[,grp.var] <- make.names(gsub ("^.. ", "", data.df[,grp.var]))
data.df[,grp.var] <- factor (data.df[,grp.var], levels = c("LS", "HS.CB"))
data.df$signature <- factor (data.df$signature, levels = c(BL,TX))
data.df <- data.df %>% arrange_(c(grp.var, "signature", "score"))
data.df$combined.x <- paste(data.df[,grp.var],data.df$signature, sep = "_")

plotGroupsWithPairs (main.df = data.df, 
                     x.var = "combined.x", 
                     y.var = "score", 
                     x.var.name = "", 
                     y.var.name = "IFNG Score",
                     y.grp = grp.var,
                     y.limit = c(0.5, 2),
                     x.cat.order = data.df$combined.x[!duplicated(data.df$combined.x)], 
                     group.colors = COLS_CB[sort(rep(1:2,2))],
                     BL = "BL_IFNG", TX = "TX_IFNG")
# Wilcox tests
stats.out <- lapply (levels(data.df[,grp.var]), 
                     function (x) wilcox.test (score~signature,data.df[data.df[,grp.var]==x,], alternate = "two.sided"))
stats1.out <- lapply (levels(data.df[,"signature"]), 
                      function (x) wilcox.test (as.formula(paste("score",grp.var,sep="~")),data.df[data.df[,"signature"]==x,], alternate = "two.sided"))
#############
# Figure 3D #
#############
data.df <- readWorksheetFromFile (sheet = "Fig. 3D", file = "./data/Source Data/SourceData_Fig3.xlsx")
grp.var <- "Pembro_CB_groups"

data.df[,grp.var] <- make.names(gsub ("^.. ", "", data.df[,grp.var]))
data.df[,grp.var] <- factor (data.df[,grp.var], levels = c("LS", "HS.CB"))
data.df$signature <- factor (data.df$signature, levels = c(BL,TX))
data.df <- data.df %>% arrange_(c(grp.var, "signature", "score"))
data.df$combined.x <- paste(data.df[,grp.var],data.df$signature, sep = "_")

plotGroupsWithPairs (main.df = data.df, 
                     x.var = "combined.x", 
                     y.var = "score", 
                     x.var.name = "", 
                     y.var.name = "CYT Score",
                     y.grp = grp.var,
                     y.limit = c(0.5, 2),
                     x.cat.order = data.df$combined.x[!duplicated(data.df$combined.x)], 
                     group.colors = COLS_CB[sort(rep(1:2,2))],
                     BL = "BL_CYT", TX = "TX_CYT")
# Wilcox tests
stats.out <- lapply (levels(data.df[,grp.var]), 
                     function (x) wilcox.test (score~signature,data.df[data.df[,grp.var]==x,], alternate = "two.sided"))
stats1.out <- lapply (levels(data.df[,"signature"]), 
                      function (x) wilcox.test (as.formula(paste("score",grp.var,sep="~")),data.df[data.df[,"signature"]==x,], alternate = "two.sided"))
#############
# Figure 3E #
#############
data.df <- readWorksheetFromFile (sheet = "Fig. 3E", file = "./data/Source Data/SourceData_Fig3.xlsx")
# subset data to include only baseline scores
data.df <- data.df[grep("BL_",data.df[,"signature"]),]
data.df[,"signature"] <- gsub ("BL_", "", data.df[,"signature"])
data.df[,"signature"] <- gsub ("abs_", "", data.df[,"signature"])
wilcox.results.2 <- do.call (rbind,lapply(unique(data.df[,"signature"]),
                                    wilcox.test.group2, data = data.df, x.var = "Pembro_CB_groups", x.var.levels = c("2. HS/CB", "1. LS")))
summary.tbl <- prepFisherTable (data.df)
tabletext <- cbind (c("Signature", wilcox.results.2[,1]),
                    c("HS/CB mean (range)", paste(unlist(round(summary.tbl[summary.tbl$Pembro_CB_groups=="2. HS/CB",4],digit=2))," (",
                                                  unlist(round(summary.tbl[summary.tbl$Pembro_CB_groups=="2. HS/CB",5],digit=2)),"-",
                                                  unlist(round(summary.tbl[summary.tbl$Pembro_CB_groups=="2. HS/CB",6],digit=2)),")",sep="")),
                    c("LS mean (range)",paste(unlist(round(summary.tbl[summary.tbl$Pembro_CB_groups=="1. LS",4],digit=2))," (",
                                              unlist(round(summary.tbl[summary.tbl$Pembro_CB_groups=="1. LS",5],digit=2)),"-",
                                              unlist(round(summary.tbl[summary.tbl$Pembro_CB_groups=="1. LS",6],digit=2)),")",sep="")),
                    c("Delta-mean (95% CI)", wilcox.results.2[,"CI"]),
                    c("p", wilcox.results.2[,"p"]))
plotForestPlot (tabletext=tabletext, 
                stats.results=wilcox.results.2, 
                title="Wilcoxon ranksum test",
                p.val.cutoff = 0.10) 


