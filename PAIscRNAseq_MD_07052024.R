rm(list=ls())
setwd("C:/Users/dasmohua/Downloads/pai1scrna")
dev.off()

library("Seurat")
library("stringr")
library("ggplot2")
library(tidyverse)
library(dplyr)
library(patchwork)
library(ggplot2)
library(sctransform)
library(scDblFinder)
library(DoubletFinder)
library(UCell)
library(CytoTRACE2)

set.seed(1234)
##########################################################################################################################################################
##Read object
HCT116 <- readRDS("C:/Users/dasmohua/Downloads/pai1scrna/HCT116.rds")
hct116_md <- readRDS("C:/Users/dasmohua/Downloads/pai1scrna/hct116_md.rds")
#Leiden color palette extended not checked
colors <- c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b',
                     '#e377c2', '#b5bd61', '#17becf', '#aec7e8', '#ffbb78', '#98df8a',
                     '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#dbdb8d', '#9edae5',
                     '#ad494a', '#8c6d31', '#FF34FF', '#763568', '#63FFAC', '#B79762',
                     '#004D43', '#002dff', '#997D87', '#5A0007', '#809693', '#6A3A4C',
                     '#1B4400', '#4FC601', '#3B5DFF', '#4A3B53', '#FF2F80', '#61615A',
                     '#BA0900', '#6B7900', '#00C2A0', '#FFAA92', '#FF90C9', '#B903AA',
                     '#D16100', '#DDEFFF', '#000035', '#7B4F4B', '#A1C299', '#300018',
                     '#0AA6D8', '#013349', '#00846F', '#372101', '#FFB500', '#C2FFED',
                     '#A079BF', '#CC0744', 'yellow', 'green', 'magenta', 'red')  # added 4 more colors for the 'individuals' plot
                     
##########################################################################################################################################################
Idents(hct116_md) <- "Dose"
Idents(HCT116) <- "seurat_clusters"

png("featuresofobjectfromNir.png", width = 1200, height = 800)
VlnPlot(HCT116, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
dev.off()

##########################################################################################################################################################
##For publication figures::::

png("bbUMAPclusterssub2.png", width = 1400, height = 1200)
DimPlot(bb, reduction = "umap", cols = colors,label = TRUE, label.size = 10, pt.size=2) +
  guides(colour = guide_legend(override.aes = list(size=25))) +
  theme(legend.text=element_text(face="bold", size=35))
dev.off()


hct116_md$Dose <- factor(x = hct116_md$Dose, levels = c("PAR", "HCT116-LD", "HCT116-MD","HCT116-HD"))

png("splitUMAPsubtype2.png", width = 3000, height = 1000)
DimPlot(bb, reduction = "umap", cols = colors,split.by = "Dose", label = TRUE, label.size = 10, pt.size=2) +
  guides(colour = guide_legend(override.aes = list(size=25))) +
  theme(legend.text=element_text(face="bold", size=35))
dev.off()
##################################################################################################
#Finding markers
hct116markers<-FindAllMarkers(HCT116,min.pct = 0.25,logfc.threshold = 0.25)
write.csv(hct116markers,"./hct116markers.csv",row.names = T)
hct116markers <- read_csv("hct116markers.csv")

newsig <- c("OASL","IFIT3","IFIT2","RARRES3","RSAD2","IFI44","CCL5","DDX58","IFIT1","IFIH1","RAET1L","ISG20","DHX58","DDX60",
           "BIRC3")
rps <- c("RIF1", "PARPBP", "RAD51", "XRCC5")
degs1 <- c("SERPINE1", "SMARCD3")
degs2 <- c("CYP51A1", "SC5D","ACAT2",  "HMGCR", "HMGCS1", "MVD", "FDPS")
degs3 <- c("SERPINE1", "SMARCD3","CYP51A1", "SC5D","ACAT2",  "HMGCR", "HMGCS1", "MVD", "FDPS")
degs4 <- c("SERPINE1", "SMARCD3","CYP51A1", "SC5D","ACAT2",  "HMGCR", "HMGCS1", "MVD", 
           "FDPS", "KRTCAP3", "SNHG9","SNHG25")
degs5 <- c("YTHDF1", "YTHDF2","YTHDF3")

png("FP_RESISTM.png", width = 1000, height = 1000)
FeaturePlot(HCT116, features = degs3, order = TRUE, pt.size = 1,label.size = 40)
dev.off()

FeaturePlot(bb, features = c("SYT1", "MAPK12","KRTCAP3"), order = FALSE, pt.size = 0.5,label.size = 20)

Idents(bb) <- "Dose"
png("pai1.png", width = 800, height = 600)
VlnPlot(bb, features = degs5, pt.size = 0, log = TRUE, cols = colors) +
  guides(colour = guide_legend(override.aes = list(size=20))) +
  theme(legend.text=element_text(face="bold", size=20))
dev.off()

png("allclusterannotation.png", width = 800, height = 600)
DotPlot(hct116_md, features = degs3, cols = c("lightgrey", "red"), dot.min = 0.1, col.min = 0.5,
        cluster.idents = T) +
  theme(axis.text.x = element_text(angle = 90),panel.grid = element_line(color = "white", size = 0.5, linetype = 2))
dev.off() 

# FINAL Markers for clusters 0-4: "LINC01405","ADM","MGMT","DDIT4","KRTCAP3","STK39",
# "DPEP1","ANXA10","ANGPT2","HAS3","CD33","PPARG","EFNB2","RNF43","NRIP1","ACSL5",
# "CCDC102B","SYTL5","TARID","ANTXR2","GALNT5","LINC01468",
# "SNHG9","SNHG25",
# "CCL5","IFIT2","OASL","IFIT3","IFIT1","RSAD2","RARRES3","HERC5","IFIH1","DDX58"
##################################################################################################
#Merging clusters: MIXING 0,8; MIXING 1,5,9; MIXING 2,3; MIXING 4,7; KEEP 10; KEEP 6: Ran res 0.1
#new obj name bb
bbmarkers<-FindAllMarkers(bb,min.pct = 0.25,logfc.threshold = 0.25)
write.csv(bbmarkers,"./bb.csv",row.names = T)



#Rename idents to seuratclusters
Idents(hct116_md) <- "seurat_clusters"

new.cluster.ids <- c(
  "KRTCAP3 high", #0
  "RESIST-M2 high", #1
  "SMARCD3 high", #2
  "SNHG9/25 high", #3
  "RESIST-M1 high") #4

names(new.cluster.ids) <- levels(hct116_md)
bb <- RenameIdents(hct116_md, new.cluster.ids)
bb$Subtype2 <- Idents(bb)
Idents(bb)
bb$Subtype2

Idents(bb) <- "Subtype2"

png("bbUMAPsubtype2.png", width = 1400, height = 1200)
DimPlot(bb, reduction = "umap", cols = colors,label = TRUE, label.size = 10, pt.size=2) +
  guides(colour = guide_legend(override.aes = list(size=25))) +
  theme(legend.text=element_text(face="bold", size=35))
dev.off()
####################################################################################

#cell proportion
table(bb$Dose)#Count cell number in each group

# PAR HCT116-LD HCT116-MD HCT116-HD 
# 4039      3722      4775      5103

table(bb$Subtype2,bb$Dose)#Count each cell type number in each group

# PAR HCT116-LD HCT116-MD HCT116-HD
# 0     1         9        22      2953
# 1   366       660      1702         2
# 2  1194      1028         0         0
# 3  1158       835         5         0
# 4   293       352       511       644
# 5   251       232       807         3
# 6   400       117       714        47
# 7   215       180       420       357
# 8     0         1         2      1034
# 9   136       301       549         0
# 10   25         7        43        63


# PAR HCT116-LD HCT116-MD HCT116-HD
# RESIST-M2 low  1139      1311      3802        49
# RESIST-M2 high 2374      1873         8         0
# SMARCD3 high      0         8        11      3994
# Unaffected      503       525       912       997
# SERPINE1 high    23         5        42        63

Cellratio <- prop.table(table(bb$Subtype2,bb$Dose),margin = 2)#count cell type ratio in each group
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
#stacked bar chart
png("cellprop_bydose.png", width = 500, height = 500)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()
dev.off()

#facet bar plot
pdf("cellprop_bysubtype2.pdf", width = 8, height = 4)
ggplot(Cellratio, aes(x=Var2, y=Freq))+
  geom_bar(stat='identity', fill="forest green")+  
  facet_wrap(~Var1,ncol = 6, scales = "free_y")+ 
  theme(axis.text.x=element_text(angle=45, hjust=0.8))
dev.off()
# ggsave("./3_Analysis/mergeQC4_proportion.pdf", width = 20, height = 10)
# 
# 
saveRDS(bb, file = "C:/Users/dasmohua/Downloads/pai1scrna/hct116_md.rds")

png("Wnt_FP.png", width = 500, height = 200)
FeaturePlot(HCT116,features = "SERPINE1")

dev.off()
VlnPlot(HCT116,group.by = "Dose2",log = TRUE, features = "SERPINE1")

png("Serpine1_vs_clusters.png", width = 600, height = 300)
VlnPlot(HCT116,group.by = "seurat_clusters",log = TRUE, features = "SERPINE1")
dev.off()


pdf("resistm_dose.pdf", width = 9, height = 6)
DotPlot(hct116_md, features = degs1, cols = c("lightgrey", "red"), dot.min = 0.1, col.min = 0.5,
        cluster.idents = F) + theme(axis.text.x = element_text(angle = 90),panel.grid = element_line(color = "grey", size = 0.1, linetype = 2))
dev.off()

ggplot(hct116_md, aes(x = degs3, y = Idents(hct116_md), 
                        color = `p.adjust`, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("HOLA")


####################################################################################
# The Oncotype DX Colon Cancer Recurrence Score Assay evaluates the expression of 12 genes to predict the risk 
# of recurrence in stage II and III colon cancer. The panel consists of:
# 7 cancer-related genes, grouped into categories:
# Stromal genes: FAP, INHBA, and BGN.
# Cell cycle genes: MKI67 (Ki-67), MYC, and MYBL2. GADD45B, involved in DNA repair and stress response.
# 5 reference genes: These serve to normalize the expression levels of the cancer-related genes, ensuring more 
# accurate interpretation.
# 
# 
# The recombination proficiency score (RPS) is based on the expression levels for four genes involved in DNA repair pathway 
# preference (Rif1, PARI, RAD51, and Ku80), such that high expression of these genes yields a low RPS. 
# Carcinoma cells with low RPS exhibit HR suppression and frequent DNA copy number alterations, which are characteristic of 
# error-prone repair processes that arise in HR-deficient backgrounds. The RPS system was clinically validated in 
# patients with breast or non-small cell lung carcinomas (NSCLCs).

Idents(hct116_md) <- "Dose"

rps <- c("RIF1", "PARPBP", "RAD51", "XRCC5")
rcc <- c("FAP", "INHBA", "BGN",  "MKI67", "MYC", "MYBL2")
rpscc <- c("RIF1", "PARPBP", "RAD51", "XRCC5", "FAP", "INHBA", "BGN",  "MKI67", "MYC", "MYBL2")

pdf("CDX2.pdf", height=3.5, width=5)
VlnPlot(hct116_md, features = "CDX2", pt.size = 0, log = TRUE, cols = colors) +
  guides(colour = guide_legend(override.aes = list(size=20))) +
  theme(legend.text=element_text(face="bold", size=10))
dev.off()


#############
### UCell with rcc and rps ###
#############

signatures <- list("RCC" = c("INHBA", "MKI67", "MYC", "MYBL2"),
                   "RPS" = c("RIF1", "PARPBP", "RAD51", "XRCC5"))

hct116_md <- AddModuleScore_UCell(hct116_md, features = signatures)
hct116_md@meta.data %>% head()
signature.names <- paste0(names(signatures), "_UCell")

pdf("RCC_RPS_UCell.pdf", height=3.5, width=5)
VlnPlot(hct116_md, features = signature.names, group.by = "Dose", pt.size = F)
dev.off()

pdf("R_HCT116_RPSCC_Dotplot.pdf", height=4.5, width=7)
DotPlot(hct116_md, features = rpscc, group.by = "Dose",
        cols = c("lightgrey", "red"),
        dot.min = 0.1, col.min = 0.1) +
  ylab('Dose') +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

sig <- c("SERPINE1", "SMARCD3", "CYP51A1", "SC5D",
         "ACAT2", "HMGCR", "HMGCS1", "MVD", "FDPS")

pdf("R_HCT116_Sig1_Dotplot.pdf", height=4.5, width=7)
DotPlot(hct116_md, features = sig, group.by = "Dose",
        cols = c("lightgrey", "red"),
        dot.min = 0.1, col.min = 0.1) +
  ylab('Dose') +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()
#################