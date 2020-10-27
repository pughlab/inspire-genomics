#########
# Custom visulzation functions
#########

#########
# Kaplan-Meier survival analysis 
########
library(survival)
survKM <- function (dat.df, ycat, ycat.time.label, ycat.event.label,
                    grp.cat, grp.levels, 
                    withlegend = TRUE, legend.loc = "bottomleft",
                    withPval = TRUE, pval.loc = "topright",
                    add.covariate = FALSE, covariate.cat, cov.levels,
                    grp.cols = c("black", "red")){
  
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
