#########
# Immune Score Forest Plots and Analysis
# ImmuneScoreForestPlotFunctions.R
#########
library(plyr)
library(dplyr)
library(tidyr)
library(Hmisc)
library(survival)
library(forestplot)

setwd ("/Users/yangc/PughLab/Repos/inspire-genomics/")
#####
# Do stats test
wilcox.test.group2 <- function (signature.name, data, x.var = "Pembro_CB_groups", x.var.levels = c("2. HS/CB", "1. LS")){
  data <- data[as.character(data$signature) %in% signature.name,]
  data[,x.var] <- factor(data[,x.var], levels = x.var.levels)
  wilcox.result <- wilcox.test(
    as.formula(paste("score",x.var,sep="~")),
    data = data,
    alternative = "two.sided", conf.int = TRUE, conf.level = 0.95)
  estimate <- round(wilcox.result$estimate, digits = 2)
  CI.lower <- round(wilcox.result$conf.int[1], digits = 2)
  CI.upper <- round(wilcox.result$conf.int[2], digits = 2)
  results.string <- paste(estimate," (",CI.lower,",",CI.upper,")",sep="")
  pval <- Hmisc::format.pval(wilcox.result$p.value, digits = 3, eps = 0.001)
  return(c(signature = signature.name,
           hazard.ratio = estimate,
           CI.lower = CI.lower,
           CI.upper = CI.upper,
           CI = results.string,
           p = pval))
}

#####
# Prepare table
prepFisherTable <- function (main.df){
  summary.tbl <- main.df %>% 
    group_by(signature,Pembro_CB_groups) %>% 
    summarise(median = median (score),
              mean = mean (score),
              min = min (score),
              max = max (score),
              N = length (score))
  return(summary.tbl)
}

#####
# Plot data
plotForestPlot <- function (tabletext, stats.results, title, p.val.cutoff = 0.10) {
  forestplot (labeltext = tabletext, 
              mean = c(NA,as.numeric(stats.results[,2])), 
              lower = c(NA,as.numeric(stats.results[,"CI.lower"])), 
              upper = c(NA,as.numeric(stats.results[,"CI.upper"])),
              title = title,
              xlog = F, align = "l", is.summary = c(TRUE, rep(FALSE,25)), 
              col = fpColors(box = "black", lines = "black"), 
              graph.pos = 4, graphwidth = unit(2.5,"cm"),
              colgap=unit(5,"mm"),lwd.ci=1, 
              ci.vertices=TRUE, ci.vertices.height = 0.2, boxsize = 0.40,
              xticks.digits = 1, lineheight = unit(6,"mm") ,
              txt_gp = fpTxtGp(label = gpar(fontfamily = ""), 
                               xlab = gpar(cex = 1), ticks = gpar(cex = 1)),
              fn.ci_norm = local({
                i=0
                p.val.cols <- ifelse (as.numeric(wilcox.results.2[,"p"])<p.val.cutoff,
                                      yes = "red", no = "black")
                function(...,clr.line,clr.marker){
                  i<<-i+1
                  fpDrawNormalCI(...,clr.line=p.val.cols[i],clr.marker=p.val.cols[i])}
              }))
}

