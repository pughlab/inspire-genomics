############
# Figure 1 #
############
load ("./R/SetUp.Rdata")
source ("./R/Color.Palettes.R")
source ("./R/Plot.Functions.R")
source ("./R/Helper.functions.R")

# Input data:
main.df <- read.csv ("./data/Fig1.csv", stringsAsFactors = FALSE)

#######################
# Waterfall plot sorted by cancer cohort
# colored by cancer cohort

main.df[,"BEST_RECIST_TM_TARGET"] <- as.numeric(main.df[,"BEST_RECIST_TM_TARGET"])
main.df <- arrange (main.df, COHORT, desc (BEST_RECIST_TM_TARGET))
sample.count <- sort(table (main.df[,"COHORT"]))

barplot (main.df[,"BEST_RECIST_TM_TARGET"], cex.axis = 0.8,
         col = as.character(factor(main.df[,"COHORT"], labels = COLS_CANCER_COHORT)),
         ylim = c(-100, 150), border=NA, xaxs = "i",las = 2) 

title (ylab = "% delta TM", cex.lab = 1)
abline (h = c(20,0,-30), col = c("red", "black", "blue"), lty = c(3,1,3))
axis (2, at = c(20, -30), labels = c(20,-30), cex.axis = 0.8, las = 2)   

###########################
# Contoured scatter plot of delta-TM and delta-ctDNA
main.df[,"change_ctDNA"] <- as.numeric(main.df[,"change_ctDNA"])*100
fig2.plot <- ggplot(subset(main.df, Pembrolizumab_sensitivity_groups_n71!="")) + 
  stat_density_2d(aes (x = TARGETED_LESION_TM_EARLY, 
                       y = change_ctDNA, 
                       color = Pembrolizumab_sensitivity_groups_n71, 
                       fill = Pembrolizumab_sensitivity_groups_n71), geom="polygon", bins=4, alpha=.1)+
  geom_point (size = 2, aes (x = TARGETED_LESION_TM_EARLY, 
                             y = change_ctDNA, color = Pembrolizumab_sensitivity_groups_n71)) + 
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
                       
###########################
# Survival plots by pembro sensitivity
KM.out <- survKM  (dat.df = subset(main.df, Pembrolizumab_sensitivity_groups_n71!=""), 
                   ycat = "OS", 
                   ycat.time.label = "OSTIME_Months_since_C3", 
                   ycat.event.label = "OSevent",
                   grp.cat = "Pembrolizumab_sensitivity_groups_n71", 
                   grp.levels = names(table(main.df[,"Pembrolizumab_sensitivity_groups_n71"]))[-1], 
                   withlegend = FALSE, withPval = FALSE,
                   grp.cols = COLS_PEM4)
                   
KM.PFS.out <- survKM  (dat.df = subset(main.df, Pembrolizumab_sensitivity_groups_n71!=""), 
                   ycat = "PFS", 
                   ycat.time.label = "PFSTIME_Months_since_C3", 
                   ycat.event.label = "PFSevent",
                   grp.cat = "Pembrolizumab_sensitivity_groups_n71", 
                   grp.levels = names(table(main.df[,"Pembrolizumab_sensitivity_groups_n71"]))[-1], 
                   withlegend = FALSE, withPval = FALSE, legend.loc = "bottomright",
                   grp.cols = COLS_PEM4)
                   
###########################
# TMB violin plot by ICB group
main.df <- arrange (subset(main.df, Pembrolizumab_sensitivity_groups_n71!=""), 
                    Pembrolizumab_sensitivity_groups_n71, TMB_n72)
box.summary <- boxplot (TMB_n72 ~ Pembrolizumab_sensitivity_groups_n71, main.df, plot = FALSE)
new.x <- unlist(lapply (seq (from = 1, by =1, length.out = length(box.summary$names)), 
                 function (x) seq (from = x - 0.1, to = x + 0.1, length.out = box.summary$n[x])))
main.df[,"Pembrolizumab_sensitivity_groups_n71"] <- factor (main.df[,"Pembrolizumab_sensitivity_groups_n71"], 
                                                       labels = unique(gsub ("^.. ","",main.df[,"Pembrolizumab_sensitivity_groups_n71"],)))  
                                                                     
vioplot (TMB_n72 ~ Pembrolizumab_sensitivity_groups_n71, main.df, log = "y", border = COLS_PEM4, lwd = 2.5, wex = 0.8,
         plotCentre = "line", drawRect = TRUE, col = adjustcolor(COLS_PEM4, alpha = 0.10), ylim = c(0.1, 1000), las = 1,
         xlab = "", ylab = "TMB", yaxt = "n")
axis (side = 2, at =  c(0.1, 1, 10, 100, 1000), labels = c(0.1, 1, 10, 100, 1000), las = 1)
axis (side = 1, at =  1:4, labels = box.summary$names , las = 1)
axis (side = 3, at =  c(0.5,1:4), labels = c("n =", box.summary$n))
points(x = new.x, y = main.df$TMB_n72[!is.na(main.df$TMB_n72)], log = "y", add = T, pch = 16,
       col = as.character(factor (main.df[,"Pembrolizumab_sensitivity_groups_n71"], labels = COLS_PEM4))[!is.na(main.df$TMB_n72)])
abline (h = quantile (main.df$TMB_n72, na.rm = TRUE, 0.65), lty = 2, col = "red")
# Two-sided KS test for enrichment of high TMB in HS group 
text ( x = 4 , y = 500 , label = paste ("p =", round(       
ks.test (x = main.df$TMB_n72[main.df$Pembrolizumab_sensitivity_groups_n71=="HS"], 
         y = main.df$TMB_n72, alternative = "two.sided")$p.value, 2)))

###########################
# Survival plots by TMB and ctDNA change
KM.TMB.OS.out <- survKM  (dat.df = subset(main.df, TMB_ctDNA_category_n55!=""), 
                   ycat = "OS", 
                   ycat.time.label = "OSTIME_Months_since_C3", 
                   ycat.event.label = "OSevent",
                   grp.cat = "TMB_ctDNA_category_n55", 
                   grp.levels = rev(names(table(main.df[,"TMB_ctDNA_category_n55"]))[-1]), 
                   withlegend = FALSE, withPval = FALSE,
                   grp.cols = COLS_TMB4)
                   
KM.TMB.PFS.out <- survKM  (dat.df = subset(main.df, TMB_ctDNA_category_n55!=""), 
                   ycat = "PFS", 
                   ycat.time.label = "PFSTIME_Months_since_C3", 
                   ycat.event.label = "PFSevent",
                   grp.cat = "TMB_ctDNA_category_n55", 
                   grp.levels = rev(names(table(main.df[,"TMB_ctDNA_category_n55"]))[-1]), 
                   withlegend = FALSE, withPval = FALSE, legend.loc = "bottomright",
                   grp.cols = COLS_TMB4)                   