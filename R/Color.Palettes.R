##################
# Common color palettes used for visualziations
#################

library (RColorBrewer)
# Cancer Cohorts
COLS_CANCER_COHORT <- c("A: Head and Neck" = "#8dd3c7",
                        "B: Breast" = "#fb8072",
                        "C: Ovary" = "#80b1d3",
                        "D: Melanoma" = "#fdb462",
                        "E: Mixed" = "#b3de69")
# RECIST
COLS_RECIST <- c("PR/CR" = "#d73027", 
                 "SD" = "#fdae61", 
                 "PD" = "#74add1")

# Clinical benefit
COLS_CB <- c("nCB" = "black",
             "CB" = "#d73027") 

# Pembrolizumab sensitivity - 4 colors
COLS_PEM4 <- c("low" = "#4575b4",
               "mixed" ="#35978f",
               "psud" = "#dfc27d",
               "high"="#a50026")

# TMB groups - 4 colors
COLS_TMB4 <- c("black", 
               "grey",
               "red", 
               "orange")

# Mutation status - 3 colors
"BRCA2 Mutation Status" <- setNames (c("#ffffbf","#5ab4ac","#d8b365"),c("none", "germline", "somatic"))

