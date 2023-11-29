#' Data pre-processing of the enzyme rho-associated protein kinase 1 
#' IC50 bioactivities from the ChEMBL database 


# Workspace and History ----------
save.image("rock1_preprocess.Rdata")
savehistory("rock1_preprocess.Rhistory")


# Libraries ---------- 
library(tidyverse)
library(dplyr)
library(skimr)
library(recipes)
library(janitor)
library(caret)


# Raw dataset --------------------------------------------------
rock1 <- read.csv2("rhoassociatedproteinkinase1_ic50.csv",
                   sep = ";",
                   header = TRUE)  
rock1_new <- rock1
rock1_new[rock1_new == ""] <- NA
sum(is.na(rock1_new)) 
sum(duplicated(rock1_new$Molecule.ChEMBL.ID)) 
screenedData <- distinct(rock1_new, Molecule.ChEMBL.ID,
                         .keep_all = TRUE) 
sum(duplicated(screenedData$Molecule.ChEMBL.ID)) 

# Bioactivity classification ---------------------------------------------

# Analyzing pIC50 values for 
# bioactivity classification ranges.

# What ranges constitute 
# an active and an inactive compound?

# Active compounds --------------------------------------------------
active <- screenedData %>%
  select(Molecule.ChEMBL.ID,
         Smiles,
         Standard.Relation,
         Standard.Value,
         Document.ChEMBL.ID) %>% 
  mutate(Molar = as.numeric(Standard.Value) * 10 ^ -9) %>%
  mutate(pIC50 = -log10(Molar)) %>%
  filter(pIC50 >= 5.4) %>% # 5.4 is the pIC50 of velpatasvir
  arrange(pIC50) %>%
  na.omit()

# Inactive compounds --------------------------------------------------
inactive <- screenedData %>%
  select(Molecule.ChEMBL.ID,
         Smiles,
         Standard.Relation,
         Standard.Value,
         Document.ChEMBL.ID) %>%
  mutate(Molar = as.numeric(Standard.Value) * 10 ^ -9) %>%
  mutate(pIC50 = -log10(Molar)) %>%
  filter(pIC50 < 5.4) %>% 
  arrange(pIC50) %>%
  na.omit()

# Intermediate Compounds --------------------------------------------- 

# What ranges constitute 
# an intermediate compound?

# Analyzing the box and whisker plot 
# of the active and inactive
# using their numerical quartile values.

# Boxplot --------------------------------------------------
boxplot(active$pIC50, inactive$pIC50,
        xlab = "Bioactivity class",
        ylab = "pIC50 values",
        names = c("Active", "Inactive")) 

activeSkim <- skim(active$pIC50) # p0 for intermediate upper limit
inactiveSkim <- skim(inactive$pIC50) # p50 for lower limit

# Intermediate compounds -----------------------------------
intermediate <- screenedData %>%
  select(Molecule.ChEMBL.ID,
         Smiles,
         Standard.Relation,
         Standard.Value,
         Document.ChEMBL.ID) %>%
  mutate(Molar = as.numeric(Standard.Value) * 10 ^ -9) %>%
  mutate(pIC50 = -log10(Molar)) %>%
  filter(between(pIC50, 5, 5.39)) %>% # 5.20 < pIC50 > 5.39
  arrange(pIC50) %>%
  na.omit()

# All classes ------------------------------ 
boxplot(active$pIC50, 
        intermediate$pIC50,
        inactive$pIC50,
        xlab = "Bioactivity class",
        ylab = "pIC50 values",
        names = c("active", "intermediate", "inactive"))


cmpds_all <-
  rbind(active, 
        intermediate, 
        inactive) 

# Duplicate compounds and SMILES -------------------------
sum(duplicated(cmpds_all$Molecule.ChEMBL.ID))
sum(duplicated(cmpds_all$Smiles))

# Final dataset cleaning -----------------------------------
cmpds_all_final <-
  distinct(cmpds_all, 
           Molecule.ChEMBL.ID,
           .keep_all = TRUE)
cmpds_all_final <- distinct(cmpds_all,
                            Smiles,
                            .keep_all = TRUE)

# Classes categorical ----------------------------------------
cmpds_all_final$Bioactivity.Class <-
  as.factor(ifelse(
    cmpds_all_final$pIC50 >= 5.4,
    "active",
    ifelse(cmpds_all_final$pIC50 <= 4.99, "inactive",
           "intermediate")
  ))
summary.factor(cmpds_all_final$Bioactivity.Class) 
final_df <-
  arrange(cmpds_all_final, pIC50, by_group = TRUE) # increasing pIC50


# Data Visualization ---------------------------------------------
final_df %>%
  ggplot(aes(Bioactivity.Class)) +
  geom_bar(fill = "#3983C4") +
  #coord_flip() + # flipping/inversing the coordinate
  theme_bw() +
  labs(x = "Bioactivity Class",
       y = "Frequency",
       title = "Distribution of Bioactivities")
#fct_infreq(Bioactivity.Class) # create/arrange in order

# PNG export ------------------------------
ggsave(
  "boxPlot_bioactDistrib.png",
  width = 10,
  height = 7,
  units = "cm",
  dpi = 300
)

# Final dataset csv export --------------------
write.csv(final_df, file =
            "rock1cmpdsScreened.csv")


#' The rock1 pre-processed data
#' as a SMILES export.

# ChemmineR Toolkit ---------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ChemmineR")
BiocManager::install("ChemmineOB")

library("ChemmineR")
library("ChemmineOB")

library("BiocManager")
library("BiocGenerics")
library("BiocVersion")


# SMILES dataset -------------------------
rock1Smiles <- 
  final_df %>% 
  select(Molecule.ChEMBL.ID, Smiles)

# SMI export ------------------------------

# Write as table.
write.table(rock1Smiles, 
            file = "rock1_sub.txt",
            sep = "\t")
rock1_sub <- read.table("rock1_sub.txt")

# Write as SMI.
rock1set <- rock1_sub
(rock1char <- 
    as.character(
      rock1set$Smiles))
writeLines(rock1char, "rock1.smi") # smiles string
rock1set2 <- read.SMIset("rock1.smi")

chembl_ids <- rock1Smiles$Molecule.ChEMBL.ID
cid(rock1set2) <- chembl_ids

write.SMI(rock1set2, file = "rock1_1217.smi",
          cid = TRUE)



#' PaDEL descriptor 
#' feature extraction 

# Dataset --------------------------------------------------
# PaDEL descriptors for 1D and 2D.
df1 <- read_csv("rock1_1D2D_MACCS.csv")

# Screened compounds.
df2 <- read.csv("rock1cmpdsScreened.csv", sep = ",") %>% 
  select("Molecule.ChEMBL.ID",
         "Bioactivity.Class",
         "pIC50") 

# Export for Matlab analysis.
df3 <- dplyr::bind_cols(df2, df1)

# for checking purposes 
all(df3$Molecule.ChEMBL.ID == df3$Name)

# if you want to subset one of the names
#df3 <- subset(df3, select = -Name or -Molecule.ChEMBL.ID)

write.csv(df3, file = "forMatlab.csv") # for matlab analysis

# Colab/python export
df3_colab <- df3
df3_colab[df3_colab == ""] <- NA
write.csv(df3_colab, file = "forColab.csv") 


# Store fingerprints in a matrix -------------------------
maccs_col <- names(df3)[grepl("^MACCS", names(df3))]
selected_col <- c("Name", maccs_col)
maccs_matrix <- as.matrix(df3[, selected_col])


# Display the MACCS  -----------------------------------
# If you want to find a specific compound
# and display only its MACCS fingerprint
compound_name <- "CHEMBL67352" # compound example   

# Find the row index for the specific compound
compound_row_index <- which(maccs_matrix[, 1] == compound_name)

if (length(compound_row_index) > 0) {
  specific_compound_vector <- as.vector(maccs_matrix[compound_row_index, -1])  
  
  cat(specific_compound_vector, sep = " ")
} else {
  cat("Compound not found")
}


#' Model development

# Data splitting for Machine learning -----------------------------------

# Set seed for reproducibility
set.seed(123)  
indices <- sample(nrow(df3))  

trainSize <- round(0.8 * nrow(df3))
trainingSet <- df3[indices[1:trainSize], ]
testingSet <- df3[indices[(trainSize + 1):nrow(df3)], ]

# Export for MATLAB model development
write.csv(trainingSet, file = "trainingSet.csv")
write.csv(testingSet, file = "testingSet.csv")


