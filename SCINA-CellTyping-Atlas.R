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

library(SCINA)
library(MLmetrics)
library(readxl)

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
wkdir<-"/SFS/scratch/mashayek/Liver/GSE185477_snRNA/"
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


# Import Marker Genes
source("Functions-new.R")
marker.long = readxl::read_xlsx("MarkersList.xlsx", col_names = TRUE, sheet = 1)
marker.list = MarkerFormat(dataframe = marker.long, format_type = "long")

# SCINA -------------------------------------------------------------------
markers = marker.list
out_SCINA = Run_SCINA(query = liver[["SCT"]]@data, query_type = "normalized", markers = markers)

# Add annotation and probabilities
liver$SCINA_celltype = out_SCINA$pr
liver$SCINA_probabilities = out_SCINA$pr.score

# Visualize
p1<-DimPlot(liver, group.by = c("SCINA_celltype"), label = TRUE) + FeaturePlot(liver, features = "SCINA_probabilities")
pdf("SCINA.CellType.pdf",width=22, height =15)
print(p1)
plot_grid(p1)
dev.off()


# Thresholding ------------------------------------------------------------

# Set threshold [0 - 1]
prob_threshold = 0.9

# Apply threshold

liver$SCINA_probabilities_thresh = liver$SCINA_probabilities<prob_threshold
liver$SCINA_celltype_thresh = liver$SCINA_celltype
liver$SCINA_celltype_thresh[liver$SCINA_probabilities<prob_threshold] = "unknown"
Plot<-FeaturePlot(liver, features = "SCINA_probabilities")+DimPlot(liver, group.by = "SCINA_probabilities_thresh")+DimPlot(liver, group.by = "SCINA_celltype_thresh", label = TRUE)

pdf("SCINA.Plot.pdf", width=15, height =25)
print(Plot)
dev.off()


