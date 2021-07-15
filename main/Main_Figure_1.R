###############
# Figure 1 Main
# Figure1.R
###############
setwd ("/Users/yangc/PughLab/Repos/inspire-genomics/")
load (paste(getwd(),"R","Setup.R", sep = "/"))
load (paste(getwd(),"R","Plot.Functions.R", sep = "/"))
#############
# Figure 1A #
#############
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))
library (XLConnect)
library (plyr)
data.df <- readWorksheetFromFile (sheet = "Fig. 1A", file = "./data/Source Data/SourceData_Fig1.xlsx")
sample.count <- table (data.df$COHORT, factor (data.df$BEST_OVERALL_RECIST, levels = names (COLS_RECIST_NE)))
num.cohort <- rowSums(sample.count)
num.recist <- colSums(sample.count)
plotWaterfall (main.df = data.df, sortVar = "COHORT", barCol = "COHORT")
#############
# Figure 1B #
#############
data.df <- readWorksheetFromFile (sheet = "Fig. 1B", file = "./data/Source Data/SourceData_Fig1.xlsx")
n.pt <- length(data.df$PATIENT_ID)
n.ctdna <- sum(data.df$Early.TM.avail.=="yes" & data.df$Early.ctDNA.avail=="yes")
n.wes <- sum(data.df$Early.TM.avail.=="yes" & data.df$Early.ctDNA.avail=="yes" & data.df$Baseline.WES.avail=="yes")
pem.sens.group.count <- table (data.df[data.df$Early.TM.avail.=="yes" & data.df$Early.ctDNA.avail=="yes" & data.df$Baseline.WES.avail=="yes",
                                       "Pembrolizumab_sensitivity_groups"])
#############
# Figure 1C #
# Contoured scatter plot of delta-TM and delta-ctDNA
#############
data.df <- readWorksheetFromFile (sheet = "Fig. 1C", file = "./data/Source Data/SourceData_Fig1.xlsx")
fig1c <- plotScatterContour(data.df)
fig1c
#############
# Figure 1D #
#############
data.df <- readWorksheetFromFile (sheet = "Fig. 1D", file = "./data/Source Data/SourceData_Fig1.xlsx")
KM.out <- survKM  (dat.df = subset(data.df, Pembrolizumab_sensitivity_groups!=""), 
                   ycat = "OS", 
                   ycat.time.label = "OSTIME_Months_since_C3", 
                   ycat.event.label = "OSevent",
                   grp.cat = "Pembrolizumab_sensitivity_groups", 
                   grp.levels = names(table(data.df[,"Pembrolizumab_sensitivity_groups"])), 
                   withlegend = FALSE, withPval = FALSE,
                   grp.cols = COLS_PEM4)
#############
# Figure 1E #
#############
data.df <- readWorksheetFromFile (sheet = "Fig. 1E", file = "./data/Source Data/SourceData_Fig1.xlsx")
KM.PFS.out <- survKM  (dat.df = subset(data.df, Pembrolizumab_sensitivity_groups!=""), 
                       ycat = "PFS", 
                       ycat.time.label = "PFSTIME_Months_since_C3", 
                       ycat.event.label = "PFSevent",
                       grp.cat = "Pembrolizumab_sensitivity_groups", 
                       grp.levels = names(table(data.df[,"Pembrolizumab_sensitivity_groups"])), 
                       withlegend = FALSE, withPval = FALSE, legend.loc = "bottomright",
                       grp.cols = COLS_PEM4)
#############
# Figure 1F #
#############
data.df <- readWorksheetFromFile (sheet = "Fig. 1F", file = "./data/Source Data/SourceData_Fig1.xlsx")
library(plyr)
data.df <- arrange (data.df,Pembrolizumab_sensitivity_groups, TMB_n72)
data.df[,"Pembrolizumab_sensitivity_groups"] <- gsub("^.. ","",data.df[,"Pembrolizumab_sensitivity_groups"])
plotGroupsWithScatter (data.df, 
                       x.var="Pembrolizumab_sensitivity_groups", 
                       y.var="TMB_n72", x.var.name="", y.var.name="TMB",
                       x.cat.order=c("LS","MSER","MSPP","HS"), group.colors = COLS_PEM4, 
                       y.limit = c(0.1, 1000), log.scale = "y", 
                       lwd = 2.5, wex = 0.8, pch = 16,
                       custom.axis = TRUE, custom.at = c(0.1, 1, 10, 100, 1000), 
                       custom.lbls = c(0.1, 1, 10, 100, 1000))
  
# Two-sided KS test for enrichment of high TMB in HS group 
ks.test (x = data.df$TMB_n72[data.df$Pembrolizumab_sensitivity_group =="HS"], 
         y = data.df$TMB_n72[data.df$Pembrolizumab_sensitivity_group !="HS"], alternative = "two.sided")

#############
# Figure 1G #
#############
data.df <- readWorksheetFromFile (sheet = "Fig. 1G", file = "./data/Source Data/SourceData_Fig1.xlsx")
data.df$change_ctDNA[is.na(data.df$change_ctDNA)] <- 0 # no ctDNA deteted at baseline but increased after
data.df$TMBcat <- ifelse (as.numeric(data.df$TMB_n72) >= 10, yes = "high-TMB", no = "low-TMB")
data.df$ctDNAcat <- ifelse (as.numeric(data.df$change_ctDNA) >= 0, yes = "increase from baseline", no = "decrease from baseline")
data.df$TMB_ctDNA_category_10 <- paste (data.df$TMBcat,data.df$ctDNAcat, sep = ",")

KM.TMB.OS.out <- survKM  (dat.df = subset(data.df, TMB_ctDNA_category_10!=""), 
                          ycat = "OS", 
                          ycat.time.label = "OSTIME_Months_since_C3", 
                          ycat.event.label = "OSevent",
                          grp.cat = "TMB_ctDNA_category_10", 
                          grp.levels = rev(names(table(data.df[,"TMB_ctDNA_category_10"]))), 
                          withlegend = FALSE, withPval = FALSE,
                          grp.cols = COLS_TMB4)
#############
# Figure 1H #
#############
data.df <- readWorksheetFromFile (sheet = "Fig. 1H", file = "./data/Source Data/SourceData_Fig1.xlsx")
data.df$change_ctDNA[is.na(data.df$change_ctDNA)] <- 0 # no ctDNA deteted at baseline but increased after
data.df$TMBcat <- ifelse (as.numeric(data.df$TMB_n72) >= 10, yes = "high-TMB", no = "low-TMB")
data.df$ctDNAcat <- ifelse (as.numeric(data.df$change_ctDNA) >= 0, yes = "increase from baseline", no = "decrease from baseline")
data.df$TMB_ctDNA_category_10 <- paste (data.df$TMBcat,data.df$ctDNAcat, sep = ",")

KM.TMB.PFS.out <- survKM  (dat.df = subset(data.df, TMB_ctDNA_category_10!=""), 
                           ycat = "PFS", 
                           ycat.time.label = "PFSTIME_Months_since_C3", 
                           ycat.event.label = "PFSevent",
                           grp.cat = "TMB_ctDNA_category_10", 
                           grp.levels = rev(names(table(data.df[,"TMB_ctDNA_category_10"]))), 
                           withlegend = FALSE, withPval = FALSE, legend.loc = "bottomright",
                           grp.cols = COLS_TMB4)