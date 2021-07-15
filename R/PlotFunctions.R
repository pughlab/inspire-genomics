######
# Plot Fuctions
# PlotFunctions.R
######

# WATERFALL PLOT
# In descending best TM, options to sort by Cohort, color by Cohort or RECIST1.1 best response
# Return table of sample count by sort category

plotWaterfall <- function(main.df, sortVar = "COHORT", barCol = "COHORT"){
  library (plyr)
  main.df[,"BEST_RECIST_TM_TARGET"] <- as.numeric(main.df[,"BEST_RECIST_TM_TARGET"])
  
  if (!(sortVar %in% c("COHORT","NONE"))){
    print (paste("ERROR: Invalid variable for ordering", sortVar, sep = ":"))
    print ("Please choose either COHORT or NONE")
    print ("Default order by COHORT will be used" )
    sortVar <- "COHORT"
  }
  
  if (sortVar == "NONE"){
    main.df <- main.df[order(main.df[,"BEST_RECIST_TM_TARGET",], decreasing = TRUE),]
  } else {
    main.df <- arrange (main.df, COHORT, desc (BEST_RECIST_TM_TARGET))
  }
  
  if (barCol == "COHORT"){
    barcols <- as.character(factor(main.df[,barCol], labels = COLS_CANCER_COHORT))
  }else if (barCol == "BEST_OVERALL_RECIST"){
    barcols <- as.character(factor(main.df[,barCol], order = names(COLS_RECIST_NE),  labels = COLS_RECIST_NE))
  }else {
    barcols <- barCol
  }
  
  barplot (main.df[,"BEST_RECIST_TM_TARGET"],
           col = barcols, cex.axis = 0.8,
           ylim = c(-100, 150), border=NA, xaxs = "i",las = 2) 
  
  title (ylab = "% delta TM", cex.lab = 1)
  abline (h = c(20,0,-30), col = c("red", "black", "blue"), lty = c(3,1,3))
  axis (2, at = c(20, -30), labels = c(20,-30), cex.axis = 0.8, las = 2) 
  
  print ("Waterfall plot generated")
  
  # sample count table 
  sample.count <- list (
    cohort = sort(table (main.df[,"COHORT"])),
    recist = sort(table (main.df[,"BEST_OVERALL_RECIST"])))
  #return (sample.count)
}

# SCATTER CONTOUR PLOT
plotScatterContour <- function(main.df){
  library (ggplot2)
  library(scales)
  asinh_breaks <- function(x) {
    br <- function(r) {
      lmin <- round(log10(r[1]))
      lmax <- round(log10(r[2]))
      lbreaks <- seq(lmin, lmax, by = 1)
      breaks <- 10 ^ lbreaks
    }
    p.rng <- range(x[x > 0], na.rm = TRUE)
    breaks <- br(p.rng)
    if (min(x) <= 0) {breaks <- c(0, breaks)}
    if (sum(x < 0) > 1) { #more negative values that expected from expanding scale that includes zero
      n.rng <- -range(x[x < 0], na.rm = TRUE)
      breaks <- c(breaks, -br(n.rng))
    }
    return(sort(breaks))
  }
  asinh_trans <- function() {
    trans_new("asinh",
              transform = asinh,
              inverse   = sinh,
              breaks = asinh_breaks)
  }
  main.df[,"change_ctDNA"] <- as.numeric(main.df[,"change_ctDNA"])*100
  fig2.plot <- ggplot(subset(main.df, Pembrolizumab_sensitivity_groups!="")) + 
    stat_density_2d(aes (x = TARGETED_LESION_TM_EARLY, 
                         y = change_ctDNA, 
                         color = Pembrolizumab_sensitivity_groups, 
                         fill = Pembrolizumab_sensitivity_groups), geom="polygon", bins=4, alpha=.1)+
    geom_point (size = 2, aes (x = TARGETED_LESION_TM_EARLY, 
                               y = change_ctDNA, color = Pembrolizumab_sensitivity_groups)) + 
    scale_fill_manual(values = as.character(COLS_PEM4),
                      name = "Pembrolizumab\nMolecular\nSensitivity\nGroups", 
                      labels = c("LS (n=32)","MSER (n=7)","MSPP (n=16)","HS (n=16)"))+
    scale_color_manual(values = as.character(COLS_PEM4),
                       name = "Pembrolizumab\nMolecular\nSensitivity\nGroups", 
                       labels = c("LS (n=32)","MSER (n=7)","MSPP (n=16)","HS (n=16)"))+
    geom_hline(yintercept = asinh(0), colour="#990000", linetype="dashed") + 
    geom_vline(xintercept = asinh(0), colour="#990000", linetype="dashed") + 
    theme_bw(base_size = 10) + 
    theme(legend.position = "right") +
    xlab("Tumor measurement change (%)") + 
    ylab("ctDNA change (%)") + 
    scale_y_continuous(trans = asinh_trans(), limits = c(-200, 5000),
                       breaks = c(-100,-10,0,10,100,1000),
                       labels = c(-100,-10,0,10,100,1000))+
    scale_x_continuous(trans = asinh_trans(), limits = c(-200, 500),
                       breaks = c(-100,-10,0,10,100,1000),
                       labels = c(-100,-10,0,10,100,1000))   
  return (fig2.plot)
}

# KAPLAN-MEIER PLOT & SURVIVAL ANALYSIS
survKM <- function (dat.df, ycat, ycat.time.label, ycat.event.label,
                    grp.cat, grp.levels, 
                    withlegend = TRUE, legend.loc = "bottomleft",
                    withPval = TRUE, pval.loc = "topright",
                    add.covariate = FALSE, covariate.cat, cov.levels,
                    grp.cols = c("black", "red")){
  library(survival)
  # Set up variables for the survival analysis
  #Time <- ifelse (ycat == "OS", yes = "OSTIME_Months", no = "PFSTIME_Months")
  #Event <- ifelse (ycat == "OS", yes = "OSevent", no = "PFSevent")
  Time <- ycat.time.label
  Event <- ycat.event.label 
  ylabel <- ifelse (ycat == "OS", yes = "Overall Survival", no = "Progression-Free Survival")
  plot.title <- paste ("Stratified by", grp.cat)
  
  # Set levels for comparison
  dat.df[,grp.cat] <- factor (as.character(dat.df[,grp.cat]), levels = grp.levels)
  print(levels(dat.df[,grp.cat]))
  # Build survival formulas for survival fit
  print ("Build Survival Formulas")
  # base model with only main effects
  fit.formula <- as.formula(paste0 ("Surv(",Time,",",Event,")~",grp.cat))
  # Fit coxph model 
  print ("Fit coxph model with main effect")
  print (paste0 ("Surv(",Time,",",Event,")~",grp.cat))
  surv.mod <- survfit (formula = fit.formula,
                       type = "kaplan-meier",
                       conf.type = "log",
                       data = dat.df)
  # Plot KM plot with main effect only
  print ("plot KM Curve")
  plot (surv.mod, main = plot.title, xlab = "Months", ylab = ylabel,lwd = 3.5, col = grp.cols, lty = 1,
        las = 1, cex.axis = 1, cex.lab = 1)
  
  # Add Legend p value of main effect with or without correction of co-variates
  print ("Test KM Fits between groups")
  
  # If sample size legend is to be added to the bottomleft corner
  if (withlegend){
    legend (legend.loc, col = grp.cols, lty = 1, 
            legend = paste0(grp.levels, "(n=", surv.mod$n, ")"), 
            bty = "n", cex = 0.8, lwd = 2.5)
  }
  
  # Stats test if p-value is to be added to the top-right corner
  if(withPval){
    if (add.covariate){
      dat.df[,covariate.cat] <- factor (as.character (dat.df[,covariate.cat]),
                                        levels = cov.levels)
      fit.formula.covariate <- as.formula(paste0 ("Surv(",Time,",",Event,")~",grp.cat,"+",covariate.cat))
      print ("Fit coxph model with covariate")
      print(paste0 ("Surv(",Time,",",Event,")~",grp.cat,"+",covariate.cat))
      print (dat.df[,covariate.cat])
      
      surv.mod.covariate <- coxph  (formula = fit.formula.covariate,
                                    data = dat.df)
      p.val.corrected <- summary(surv.mod.covariate)$coefficients[1,5]
      legend (pval.loc, bty = "n", legend = paste0("p=", round (p.val.corrected, 4)), cex = 1)
    }else {
      surv.mod <- coxph  (formula = fit.formula,data = dat.df)
      p.val.main <- Hmisc::format.pval(summary(surv.mod)$logtest[3], digits = 3, eps = 0.001)
      legend (pval.loc, bty = "n", legend = paste0("p = ", p.val.main), cex = 1)
    }
  }
  if (add.covariate){
    return(surv.mod.covariate)
  }else{
    return(coxph(fit.formula, dat.df))
  }
  # Return fitted model and p-values
  #return (list (surv.mod, test.raw, test.cohort,test.cell, coxph(fit.formula.cohort, dat.df)))   
  #return (list (surv.mod, coxph(fit.formula.cohort, dat.df)))
}

# VIOLIN PLOT BY GROUP
plotGroupsWithScatter <- function (main.df, x.var, y.var, x.var.name, y.var.name,
                                   x.cat.order, group.colors, 
                                   y.limit = c(0.1, 1000), log.scale = "", 
                                   pch = 16, lwd = 2.5, wex = 0.8, 
                                   custom.axis = FALSE, custom.at = c(0.1, 1, 10, 100, 1000), 
                                   custom.lbls = c(0.1, 1, 10, 100, 1000)){
  library(vioplot)
  library (plyr)
  
  plot.formula <- as.formula (paste (y.var, x.var, sep = "~"))
  names(group.colors) <- x.cat.order
  main.df[,x.var] <- factor (main.df[,x.var], levels = x.cat.order)
  
  box.summary <- boxplot (plot.formula, main.df, plot = FALSE)
  new.x <- unlist(lapply (seq (from = 1, by =1, length.out = length(box.summary$names)), 
                          function (x) seq (from = x - 0.1, to = x + 0.1, length.out = box.summary$n[x])))
  
  vioplot (plot.formula, main.df, log = log.scale, border = group.colors, lwd = lwd, wex = wex,
           plotCentre = "line", drawRect = TRUE, col = adjustcolor(group.colors, alpha = 0.10), 
           ylim = y.limit, las = 1, xlab = x.var.name, ylab = y.var.name, yaxt = ifelse (custom.axis, yes = "n", no = "s"))
  
  points(x = new.x, y = main.df[,y.var], log = log.scale, add = T, pch = pch,
         col = as.character(factor(as.character(main.df[,x.var]), levels = x.cat.order, labels = group.colors)))
  
  axis (side = 1, at =  seq(1, length (box.summary$names), by = 1), labels = box.summary$names , las = 1)
  axis (side = 3, at =  c(0.5,seq(1, length(box.summary$n), by = 1)), labels = c("n =", box.summary$n))
  if (custom.axis) axis (side = 2, at = custom.at , labels = custom.lbls , las = 1)
  
  return (box.summary)
}






