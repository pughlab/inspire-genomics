######################
# Environment setup #
# ColorScheme.R
#####################
# Color schemes used in figures
library (RColorBrewer)

COLS_CANCER_COHORT <- c("A: Head and Neck" = "#8dd3c7",
                        "B: Breast" = "#fb8072",
                        "C: Ovary" = "#80b1d3",
                        "D: Melanoma" = "#fdb462",
                        "E: Mixed" = "#b3de69")

COLS_RECIST_NE <- c("CR" = "#d73027", 
                    "PR" = "#d73027", 
                    "SD" = "#fdae61", 
                    "PD" = "#74add1",
                    "NE" = "grey")

COLS_RECIST <- c("PR/CR" = "#d73027", 
                 "SD" = "#fdae61", 
                 "PD" = "#74add1")

COLS_PEM4 <- c("low" = "#4575b4",
               "mixed" = "#35978f",
               "psud" = "#dfc27d",
               "high"= "#a50026")

COLS_TMB4 <- c("black", "grey","red", "orange")

COLS_CB <- c("nCB" = "black",
             "CB" = "#d73027") 

COLS_CB_BOX <- c("nCB" = "grey",
                 "CB" = "#d73027") 
