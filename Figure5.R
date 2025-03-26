#***Scripts for [Figure 5] heatmap (Nant cohort), KM plots (TCGA COADREAD), box plots (TCGA COADREAD)***

setwd("C:/Users/dasmohua/Downloads/SERPINE1_TCGA COADREAD data")
rm(list=ls())
dev.off()

library(tidyverse)
library(circlize)
library(ComplexHeatmap)
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)

################
### Box Plots ###
################
#Defining signatures for TCGA-COADREAD
rps <- c("RIF1", "C12orf48", "RAD51", "XRCC5") #Synonym for PARBPB is C12orf48
rcc <- c("FAP", "INHBA", "BGN",  "MKI67", "MYC", "MYBL2","GADD45B")
FrO <- c("CD22","CASP1","CISH","ALCAM")
MDPI <- c("COPE","P4HA1","ATF6","IBTK","PHLDB3")
RESIST_M1 <- c("SERPINE1", "SMARCD3")
RESIST_M2 <- c("CYP51A1", "SC5DL","ACAT2",  "HMGCR", "HMGCS1", "MVD", "FDPS") #synonym for sc5d is sc5dl

`%!in%` <- Negate('%in%')

metadata1 <- read_tsv("fromxena_survival_COADREAD_survival.txt") %>%
  filter(Redaction %!in% "Redacted") %>%
  rename("sample2" = "sample", "sample" = "_PATIENT")

metadata2 <- read_tsv("fromsynapse_clinical_molecular_public_all.txt")

metadata <- metadata1 %>% inner_join(metadata2, by = "sample")

rm(metadata1, metadata2)

data <- read_tsv("HiSeqV2.txt") %>% column_to_rownames("sample")

picked_samples <- intersect(metadata$sample2, colnames(data))
picked_metadata <- metadata %>% filter(sample2 %in% picked_samples)

picked_data <- data[, picked_samples]
picked_data <- t(picked_data) %>% as.data.frame() %>% rownames_to_column("sample2")

merged_data <- inner_join(picked_metadata,picked_data, by = "sample2")

rm(data, picked_data, picked_metadata, picked_samples, metadata)

# Filter CMS4 subtype for boxplot
cms4_data <- subset(merged_data, cms_label != "NOLBL")
#transpose dataframe 
tcgadata <- t(cms4_data [,-1]) #exclude first column and make a matrix first
colnames(tcgadata) <- cms4_data$sample2 #set first col as col names
tcgadata <- as.data.frame(tcgadata)
tcgadata <- subset(tcgadata, rownames(tcgadata) %in% c("msi", "cms_label", "SC5DL","HMGCR","HMGCS1","ACAT2","FDPS",
                                                       "CYP51A1","MVD","SMARCD3","SERPINE1",
                                                       "RIF1","C12orf48","RAD51","XRCC5",           # rps
                                                       "FAP","INHBA","BGN","MKI67","MYC","MYBL2","GADD45B", # rcc
                                                       "CD22","CASP1","CISH","ALCAM",             # FrO
                                                       "COPE","P4HA1","ATF6","IBTK","PHLDB3"))    # MDPI


# Extract the rows for SERPINE1 and SMARCD3
selected_genes <- tcgadata[RESIST_M1, ]

# Transpose the gene expression to align with columns (patients)
expression_data <- as.data.frame(t(selected_genes))
colnames(expression_data) <- RESIST_M1


expression_data[, RESIST_M1] <- lapply(expression_data[, RESIST_M1], as.numeric)


# Calculate the combined gene score (e.g., average or sum of SERPINE1 and SMARCD3)
expression_data$gene_score <- rowMeans(expression_data[, RESIST_M1], na.rm = TRUE)

# Add CMS subtype labels
expression_data$cms_label <- as.vector(t(tcgadata["cms_label", ]))

expression_data$cms_label <- as.factor(expression_data$cms_label)

custom_colors <- c("CMS1" = "#E69F00",  # Yellow
                   "CMS2" = "#0072B2",  # Blue
                   "CMS3" = "#CC79A7",  # Pink
                   "CMS4" = "#009E73")  # Green

# Plot the combined gene score across CMS subtypes
RESIST_M1 <- ggplot(expression_data, aes(x = cms_label, y = gene_score, fill = cms_label)) +
  geom_boxplot(outlier.shape = NA) +
  #geom_jitter(width = 0.2, alpha = 0.7, color = "black") +
  labs(
    #title = "Combined Gene Signature Score across CMS Subtypes",
    x = "CMS",
    y = "RESIST-M1 score"
  ) +
  scale_fill_manual(values=custom_colors) +
  guides(fill = "none")+
  stat_compare_means(method = "wilcox.test", label = "p.signif", size = 6,
                     comparisons = list(
                       c("CMS4", "CMS1"), 
                       c("CMS4", "CMS2"), 
                       c("CMS4", "CMS3")
                     ))+
  theme_gray(base_size = 20) +
  theme(axis.text = element_text(size = 14))


pdf("RESIST_M1_boxplot_25032025.pdf", width=7, height=6)
RESIST_M1 #Figure 5A, for 5A,B,E,F_choose signature accordingly, for 5C,D PETACC-3 data was used by collaborator
dev.off()

################
### KM Plots ###
################
# Filter CMS4 subtype for KM plots
cms4_data <- merged_data
cms4_data$SignatureScore <- rowMeans(cms4_data[, rps])  # Replace with your genes
cms4_data$GeneGroup <- ifelse(cms4_data$SignatureScore > median(cms4_data$SignatureScore), "High", "Low")

# Create survival object
surv_object <- Surv(time = cms4_data$OS.time, event = cms4_data$OS)
#surv_object <- Surv(time = cms4_data$PFI.time, event = cms4_data$PFI)

# Fit KM model
fit <- survfit(surv_object ~ GeneGroup, data = cms4_data)

# Step 4: Fit the Cox Proportional Hazards model to calculate Hazard Ratio (HR) and log-rank p-value
cox_model <- coxph(surv_object ~ GeneGroup, data = cms4_data)

# Step 5: Get Hazard Ratio and 95% Confidence Interval
hr <- exp(coef(cox_model))  # Hazard Ratio
ci <- exp(confint(cox_model))  # Confidence Interval for HR

# Step 6: Get log-rank p-value from the Cox model
logrank_p_value <- summary(cox_model)$logtest[3]

# Plot KM curve


rps <- ggsurvplot(
  fit,
  data = cms4_data,
  pval = F,
  conf.int = F,
  risk.table = TRUE,
  legend.labs = c("High RPS", "Low RPS"),
  title = paste("TCGA COADREAD_Overall Survival", 
                "\nHazard Ratio: ", round(hr, 2),
                " (95% CI: ", round(ci[1], 2), " - ", round(ci[2], 2), ")",
                "\nLog-rank p-value: ", round(logrank_p_value, 4)),
  xlab = "Time (Days)",
  ylab = "Survival Probability",
  palette = c("red", "blue"),
)

pdf("rcc_All CMS.pdf", width = 6, height = 6)
rps #Figure 5M, for 55H-M_choose signature accordingly
dev.off()


#################################
### Heatmap using Nant cohort ###
#################################

rm(list=ls())
### Arrange column in metadata ###
metadata <- read_csv("SG-BULK_patient_clinical_information.csv")%>%
  column_to_rownames("patient_id") %>% arrange(iCMS)

counts <- read_csv("SG-BULK_salmonTPM.csv") 
colnames(counts)[1] <- "ensg"

mart <- read_tsv("features.tsv", col_names = F) %>%
  rename("ensg" = X1,"gene" = X2) %>% select(-X3)

test <- inner_join(mart, counts, by="ensg")
#write_tsv(test, "NANT.txt")

# Trying other paper signatures
rps <- c("RIF1", "PARPBP", "RAD51", "XRCC5")
rcc <- c("FAP", "INHBA", "BGN",  "MKI67", "MYC", "MYBL2", "GADD45B")
rpscc <- c("RIF1", "PARPBP", "RAD51", "XRCC5", "FAP", "INHBA", "BGN",  "MKI67", "MYC", "MYBL2")
FrO <- c("CD22","CASP1","CISH","ALCAM")
MDPI <- c("COPE","P4HA1","ATF6","IBTK","PHLDB3")

sel <- test %>% filter(gene %in% c("SC5D","HMGCR","HMGCS1","ACAT2","FDPS",
                                   "CYP51A1","MVD","SMARCD3","SERPINE1",
                                   "RIF1","PARPBP","RAD51","XRCC5",           # rps
                                   "FAP","INHBA","BGN","MKI67","MYC","MYBL2","GADD45B", # rcc
                                   "CD22","CASP1","CISH","ALCAM",             # FrO
                                   "COPE","P4HA1","ATF6","IBTK","PHLDB3"      # MDPI
                                   )) %>%
  select(-ensg) %>% column_to_rownames("gene")

test2 <- rownames(sel) %>%
  as.data.frame() %>%
  rename("Nautanki" = ".") %>%
  mutate(GG = case_when(
    Nautanki %in% c("SC5D","HMGCR","HMGCS1","ACAT2","FDPS",
                    "CYP51A1","MVD","SMARCD3","SERPINE1") ~ "RESIST-M",
    Nautanki %in% c("RIF1","PARPBP","RAD51","XRCC5") ~ "RPS",
    Nautanki %in% c("FAP","INHBA","BGN","MKI67","MYC","MYBL2","GADD45B") ~ "RCC",
    Nautanki %in% c("CD22","CASP1","CISH","ALCAM") ~ "FrO",
    Nautanki %in% c("COPE","P4HA1","ATF6","IBTK","PHLDB3") ~ "MDPI",
    TRUE ~ NA_character_  # Use NA_character_ since "RESIST-M" is a string
  )) %>%
  column_to_rownames("Nautanki")

my_levels <- c("RESIST-M","FrO","MDPI","RCC","RPS")

test2$GG <- factor(x = test2$GG, levels = my_levels)

### Order heatmap columns as per predecided arranged order ###
sel2 <- sel[, rownames(metadata)]

min(sel2)
max(sel2)

### Calculating z-score
metadata$iCMS %>% unique()
# [1] "iCMS2" "iCMS3" NA  
metadata$CMS %>% unique()
# [1] "CMS2" "CMS4" NA     "CMS3" "CMS1"
metadata$TGFBR2 %>% unique()
# [1] "wt"  "mut"
metadata$group3 %>% unique()
# [1] NA          "iCMS2_MSS" "iCMS3_MSI" "iCMS3_MSS"
metadata$group5 %>% unique()
# [1] NA               "iCMS2_fibrotic" "iCMS2_MSS"      "iCMS3_MSI"      "iCMS3_fibrotic" "iCMS3_MSS" 

ha = HeatmapAnnotation(iCMS = metadata$iCMS,
                       CMS = metadata$CMS,
                       TGFBR2 = metadata$TGFBR2,
                       iM = metadata$group3,
                       iF = metadata$group5,
                       annotation_name_side = "left",
                       col = list(iCMS = c("iCMS2" = "#BCBD22", "iCMS3" = "#17BECF", "indeterminate" = "#9467BD"),
                                  CMS = c("CMS1" = "#D62728", "CMS2" = "#2CA02C", "CMS3" = "#8C564B", "CMS4" = "#9467BD"),
                                  TGFBR2 = c("wt" = "#baffc9", "mut" = "#E890CD"),
                                  iM = c("iCMS2_MSS" = "#1F77B4", "iCMS3_MSI" = "#2CA02C", "iCMS3_MSS" = "#9467BD"),
                                  iF = c("iCMS2_fibrotic" = "#BCBD22","iCMS2_MSS" = "#8C564B","iCMS3_MSI" = "#FF7F0E","iCMS3_fibrotic" = "#baffc9","iCMS3_MSS"= "#D62728")),
                       na_col = "white")

mat_scaled = t(scale(t(sel2)))

min(mat_scaled)
max(mat_scaled)

# Define an order of cluster identities
my_levels <- c("RESIST-M","FrO","MDPI","RCC","RPS")

x <- Heatmap(mat_scaled,
        name = "Z-score",
        show_row_dend = F,
        show_column_dend = F,
        col = colorRamp2(c(-3, 0, 7), c("blue", "white", "red")),
        top_annotation = ha,
        #bottom_annotation_height = unit(4, "mm"),
        cluster_rows = T,
        cluster_columns = T,
        row_names_side = "left",
        #row_split = test2$GG,
        cluster_row_slices = F,
        row_split =  factor(test2$GG, levels = my_levels),
        show_row_names = T,
        show_column_names = F)


pdf("Overall.pdf", width=15, height=7)
x #[Figure 5G]
dev.off()
