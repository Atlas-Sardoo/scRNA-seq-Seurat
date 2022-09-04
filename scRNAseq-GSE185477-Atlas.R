# Putty Login: ctchpcve007.merck.com port 2022
# go to scratch:
cd /SFS/scratch/mashayek
ls

# R Login:
#Scripts are modified for Rv4
qlogin -l mem_reserve=150G -l mem_free=150G -l h_vmem=150G
module purge
module load R/4.1.2
R

#Ed provided this source for us:
source("/SFS/product/R/site_RLIBS/.Rprofile.R412.Seurat.harmony")

#Required Libraries:
library(Seurat)
library(ggplot2)
library(sctransform)
library(xlsx)
library(magrittr)
library(reticulate)
library(harmony)
library(cowplot)
library(Matrix)
library(devtools)
"%||%" <- function(a, b) if (!is.null(a)) a else b

# User supplied information
# wkdir<-"~/SALAR_Projects/iDART/mESTscRNAseq/Batch2/"
wkdir<-"/SFS/scratch/mashayek/Liver/GSE185477_scRNA/"
setwd(wkdir)

# 1-Load data and create Seurat object
# Note when using Biotouring export need to change features.tsv to genes.tsv
liver.data <- Read10X(data.dir = "/SFS/scratch/mashayek/Liver/GSE185477_snRNA/")

# Initialize the Seurat object with the raw (non-normalized data).
liver  <- CreateSeuratObject(counts = liver.data, project = "GSE185477", min.cells = 3, min.features = 200)

# Make sure to define sample in the metadata
liver@meta.data$sample = liver@meta.data$orig.ident

# 2-Apply sctransform normalization
# store mitochondrial percentage in object meta data. Make sure how MT/mt is mentioned in the metadata.
liver[["percent.mt"]] <- PercentageFeatureSet(object = liver, pattern = "^MT-")
liver <- subset(liver, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10)

# run sctransform
liver <- SCTransform(liver, vars.to.regress = "percent.mt", verbose = FALSE)


# 3-Perform dimensionality reduction by PCA and UMAP embedding and visualization
# These are now standard steps in the Seurat workflow for visualization and clustering

liver <- RunPCA(liver, pc.genes = combined@var.genes, pcs.compute = 30, do.print = FALSE)
liver <- RunHarmony(liver,  "sample", theta = 2, plot_convergence = TRUE, nclust = 50, max.iter.cluster = 20, max.iter.harmony = 5, assay.use="SCT")
liver <- FindNeighbors(liver, dims = 1:30, reduction = "harmony")
liver <- RunUMAP(liver, dims = 1:30, reduction = "harmony") 
liver <- FindClusters(liver, verbose = FALSE)

#make sure to save the reads.
saveRDS(liver,"GSE185477_Harmony.rds")

#Clustering in UMAP space. Remember, Seurat initially generate UMAP with no label and we need to manually label clusters.
DimPlot<-DimPlot(liver, reduction = "umap")
pdf("DimPlot.pdf",width=14 )
print(DimPlot)
dev.off()


p1<-DimPlot(liver, reduction= "harmony", pt.size = 0.1, group.by = "sample")
p2<-DimPlot(liver, reduction="harmony", label = T, pt.size = 0.1)

v<-VlnPlot(liver, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)

pdf("Harmony.pdf",width=14)
print(v)
plot_grid(p1, p2)
dev.off()

# DotPlot: Visualizing how feature (Marker genes) expression changes across different identity classes (clusters)
features <- c( "IGHA1", "DNASE1L3", "CRHBP", "PTPRB", "STAB2", "ARAP3", "PTH1R", 
               "CRISPLD2", "HGF", "ADAMTS2", "ACTA2", "COLEC11", "DCN", "CCBE1", 
               "CALD1", "IGFBP3", "PRSS23", "PECAM1", "CCL5", "IL7R", "RF1", "PTPRC", 
               "C1ORF112", "ALAS1", "TTC7A", "AKR1D1", "CDO1", "ALDH6A1", "NIT2", 
               "TBX3", "BNIP3", "QDPR", "ETNK2", "AKR1C1", "TMEM150C", "PPA1", "HPR", 
               "ADH6", "HDHD3", "GCDH", "BRIP1", "TSKU", "IMPA2", "ECHDC3", "SLC6A12", 
               "GLYATL1", "F10", "DPP4", "ECI2", "CTH", "SLC47A1", "PIPOX", "PBLD", 
               "HBA2", "SLC25A39", "FKBP8", "GYPC", "ANXA4", "ANK3", "PDGFD", "BCL2", 
               "LRRC7", "KRT7", "CD163", "TYR", "BP", "C1QA", "TMSB10", "SAT1", "CD74",
               "C1QB", "CTSS", "AIF1", "S100A9")

DotPlot<-DotPlot(liver, features =features) + RotatedAxis()
pdf("DotPlot.NLabelled.pdf",width=20 )
print(DotPlot)
dev.off()

#Based on expression of marker genes, we can manually label each cluster.
new.cluster.ids <- c("Hepatocyte", "Hepatocyte", "Hepatocyte", "LSEC", "Hepatocyte", 
                     "Hepatocyte", "Hepatocyte", "Stellate", "Macrophage", "LSEC", 
                     "Cholangiocyte", "Hepatocyte", "Hepatocyte", "NK", "Endothelial",
                     "Hepatocyte", "Macrophage", "Hepatocyte", "Hepatocyte", 
                     "Hepatocyte")
names(new.cluster.ids) <- levels(liver)
liver <- RenameIdents(liver, new.cluster.ids)
DimPlot<-DimPlot(liver, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
pdf("DimPlot.Labelled.pdf",width=14 )
print(DimPlot)
dev.off()

#If we generate new Dotplot, it will show assigned cluster names.
DotPlot<-DotPlot(liver, features =features) + RotatedAxis()
pdf("DotPlot.Labelled.pdf",width=20 )
print(DotPlot)
dev.off()


#Define sub-features for hepatocyte marker genes.
Hepatocyte.features <- c( "C1ORF112", "ALAS1", "TTC7A", "AKR1D1", "CDO1", "ALDH6A1", "NIT2", 
                          "TBX3", "BNIP3", "QDPR", "ETNK2", "AKR1C1", "TMEM150C", "PPA1", "HPR", 
                          "ADH6", "HDHD3", "GCDH", "BRIP1", "TSKU", "IMPA2", "ECHDC3", "SLC6A12", 
                          "GLYATL1", "F10", "DPP4", "ECI2", "CTH", "SLC47A1", "PIPOX", "PBLD")


#Generate FeaturePlot for hepatocyte marker genes.
FeaturePlot<-FeaturePlot(liver, features = Hepatocyte.features)
pdf("Hepatocyte.FeaturePlot.pdf",width = 25, height = 45)
print(FeaturePlot)
dev.off()

#Generate RidgePlot for hepatocyte marker genes.
RidgePlot<-RidgePlot(liver, features =Hepatocyte.features, ncol = 2)
pdf("Hepatocyte.RidgePlot.pdf",width = 20, height = 45)
plot(RidgePlot)
dev.off()

#Generate RidgePlot for hepatocyte marker genes.
VlnPlot<-VlnPlot(liver, features = Hepatocyte.features)
pdf("HepatocyteVlnPlot.pdf",width = 25, height = 25 )
print(VlnPlot)
dev.off()

#Generate heatmap for all marker genes.
Heatmap<-DoHeatmap(liver, features = features, size = 3.5, hjust = 0, angle = 90) + NoLegend()
pdf("Heatmap.pdf", width = 10, height = 10)
print(Heatmap)
dev.off()

