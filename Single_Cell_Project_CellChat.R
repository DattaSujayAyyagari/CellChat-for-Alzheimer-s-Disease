# devtools::install_github("jinworks/CellChat")
# devtools::install_github("jokergoo/ComplexHeatmap")
# .libPaths()
# devtools::install_github("jokergoo/circlize")
# install.packages("digest")
# install.packages('NMF')
# install.packages("presto")
# packages <- c("curl", "glue", "htmltools", "processx", "Rcpp", "promises", "httpuv")
# install.packages(packages)
# install.packages("rlang")
# devtools::install_github("immunogenomics/presto")
# BiocManager::install("Biobase")
# BiocManager::install("BiocNeighbors")
# BiocManager::install("BiocGenerics")
# BiocManager::install('glmGamPoi')
library(CellChat)
library(patchwork)
library(Matrix)
library(Seurat)
library(dplyr)
library(tidyr)

# we are not integrating the data here because we want to look at the cell-cell communication 
# in alzheimer's 

############### FAD (Control) #######################

# Loading the familial Alzheimer's disease dataset that has not been treated with APC
# APC (Activated protein C)
fad.dir <- "C:/Users/sujayayyagari/OneDrive/Desktop/WT/control/GSM7092586_FAD/FAD"
# Load the matrix file using Read10X with gzipped files
fad_cnt_mtx <- Read10X(fad.dir, 
                       gene.column = 2, 
                       cell.column = 1, 
                       unique.features = TRUE, 
                       strip.suffix = FALSE)

# creating a seurat object of fad
fad = CreateSeuratObject(counts = fad_cnt_mtx, 
                         project = "FAD", 
                         min.cells = 3, 
                         min.features = 200)

# creating a new column to store the percentage of mitochondrial counts originating from a set of features
fad[["percent.mt"]] <- PercentageFeatureSet(fad, 
                                            pattern = "^mt-")

# Show QC metrics for the first 5 cells 
head(fad@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(fad, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships.
# fad
plot1 <- FeatureScatter(fad, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot2 <- FeatureScatter(fad, 
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot1 + plot2

# filtering the dataset
fad <- subset(fad, 
              subset = nFeature_RNA > 200 & 
                nFeature_RNA < 4500 & 
                percent.mt < 25)

# Normalizing the data
fad <- NormalizeData(fad, 
                     normalization.method = "LogNormalize", 
                     scale.factor = 10000)
fad[["RNA"]]$data

# finding variable genes
fad <- FindVariableFeatures(fad, 
                            selection.method = "vst", 
                            nfeatures = 2000)

# Identify the 10 most highly variable genes
fad_top10 <- head(VariableFeatures(fad), 10)

# plot variable features with labels
plot1 <- VariableFeaturePlot(fad)
plot2 <- LabelPoints(plot = plot1, 
                     points = fad_top10, 
                     repel = TRUE)
plot2

# Scaling the data
fad.all.genes <- rownames(fad)
fad <- ScaleData(fad, 
                 features = fad.all.genes)

# Perform linear dimensional reduction
fad <- RunPCA(fad, 
              features = VariableFeatures(object = fad))
ElbowPlot(fad)

# finding neighbors by constructing a KNN graph based on the euclidean distance in PCA space
fad <- FindNeighbors(fad, 
                     dims = 1:20)

# grouping cells together
fad <- FindClusters(fad, 
                    resolution = 0.6)
head(Idents(fad), 5)

# Runnning non-linear dimensionality reduction (UMAP)
fad <- RunUMAP(fad, dims = 1:20, 
               # umap.method = 'umap-learn', 
               # metric = 'correlation', 
               verbose = FALSE)
DimPlot(fad, 
        reduction = "umap",
        label =T)

# find markers for every cluster compared to all remaining cells, report only the positive ones
fad_markers <- FindAllMarkers(fad, 
                                  only.pos = TRUE)

fad_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

fad_markers %>%
  dplyr::filter(cluster == 21) %>%
  dplyr::filter(avg_log2FC > 1)

# storing all the markers with 1 log fold change in each cluster
# into a dataframe

# Initialize an empty list to store filtered dataframes
filtered_data <- list()

# Loop through clusters 0 to 21
for (i in 0:21) {
  # Filter fad_markers for the current cluster and avg_log2FC > 1
  filtered <- fad_markers %>%
    filter(cluster == i, avg_log2FC > 1)
  
  # Store the filtered dataframe in the list
  filtered_data[[i+1]] <- filtered
}

filtered_data[[22]]

# Cell Type annotations
fad <- RenameIdents(fad, "0" = "Microglia")
fad <- RenameIdents(fad, "1" = "Microglia")
fad <- RenameIdents(fad, "2" = "Microglia")
fad <- RenameIdents(fad, "3" = "Astrocytes")
fad <- RenameIdents(fad, "4" = "Endothelial cells")
fad <- RenameIdents(fad, "5" = "Oligodendrocytes")
fad <- RenameIdents(fad, "6" = "Neurons")
fad <- RenameIdents(fad, "7" = "Ependymal cells")
fad <- RenameIdents(fad, "8" = "Astrocytes")
fad <- RenameIdents(fad, "9" = "T memory cells")
fad <- RenameIdents(fad, "10" = "Oligodendrocytes")
fad <- RenameIdents(fad, "11" = "Macrophages")
fad <- RenameIdents(fad, "12" = "Unknown")
fad <- RenameIdents(fad, "13" = "Pericytes")
fad <- RenameIdents(fad, "14" = "Unknown")
fad <- RenameIdents(fad, "15" = "Neurons")
fad <- RenameIdents(fad, "16" = "Oligodendrocytes")
fad <- RenameIdents(fad, "17" = "B cells")
fad <- RenameIdents(fad, "18" = "Pericytes")
fad <- RenameIdents(fad, "19" = "Unknown")
fad <- RenameIdents(fad, "20" = "Fibroblasts")
fad <- RenameIdents(fad, "21" = "Neutrophils")

DimPlot(fad, 
        reduction = "umap",
        label = TRUE)



saveRDS(fad, file = "FAD_control")
fad_control <- readRDS(file.choose()) #read in FAD file here





############### FAD APC (Condition) #######################

# Loading the familial Alzheimer's disease dataset that has been treated with APC
fad_apc.dir <- "C:/Users/sujayayyagari/OneDrive/Desktop/WT/condition/GSM7092587_FAD-APC/FAD-APC/"
# Load the matrix file using Read10X with gzipped files
fadapc_cnt_mtx <- Read10X(fad_apc.dir, 
                       gene.column = 2, 
                       cell.column = 1, 
                       unique.features = TRUE, 
                       strip.suffix = FALSE)

# creating a seurat object of fad_apc
fad_apc = CreateSeuratObject(counts = fadapc_cnt_mtx, 
                         project = "FAD_APC", 
                         min.cells = 3, 
                         min.features = 200)

# creating a new column to store the percentage of mitochondrial counts originating from a set of features
fad_apc[["percent.mt"]] <- PercentageFeatureSet(fad_apc, 
                                            pattern = "^mt-")

# Show QC metrics for the first 5 cells 
head(fad_apc@meta.data, 5)

# Visualize QC metrics as a violin plot
VlnPlot(fad_apc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationshi
# fad_apc
plot3 <- FeatureScatter(fad_apc, 
                        feature1 = "nCount_RNA", 
                        feature2 = "percent.mt")
plot4 <- FeatureScatter(fad_apc, 
                        feature1 = "nCount_RNA", 
                        feature2 = "nFeature_RNA")
plot3 + plot4

# filtering the dataset
fad_apc <- subset(fad_apc, 
              subset = nFeature_RNA > 200 & 
                nFeature_RNA < 5000 & 
                percent.mt < 25)

# Normalizing the data
fad_apc <- NormalizeData(fad_apc, 
                    normalization.method = "LogNormalize", 
                    scale.factor = 10000)

fad_apc[["RNA"]]$data

# finding variable genes
fad_apc <- FindVariableFeatures(fad_apc, 
                            selection.method = "vst", 
                            nfeatures = 2000)

# Identify the 10 most highly variable genes
fad_apc_top10 <- head(VariableFeatures(fad_apc), 10)

# plot variable features with labels
plot3 <- VariableFeaturePlot(fad_apc)
plot4 <- LabelPoints(plot = plot3, 
                     points = fad_apc_top10, 
                     repel = TRUE)
plot4

# Scaling the data
fad_apc.all.genes <- rownames(fad_apc)
fad_apc <- ScaleData(fad_apc, 
                features = fad_apc.all.genes)

# wt[["RNA"]]$scale.data (don't print)

# Perform linear dimensional reduction
fad_apc <- RunPCA(fad_apc, 
             features = VariableFeatures(object = fad_apc))
ElbowPlot(fad_apc)

# finding neighbors by constructing a KNN graph based on the euclidean distance in PCA space
fad_apc <- FindNeighbors(fad_apc, 
                     dims = 1:17)

# grouping cells together
fad_apc <- FindClusters(fad_apc, 
                    resolution = 0.3)
head(Idents(fad_apc), 5)

# Runnning non-linear dimensionality reduction (UMAP)
fad_apc <- RunUMAP(fad_apc, dims = 1:20, 
              # umap.method = 'umap-learn', 
              # metric = 'correlation', 
              verbose = FALSE)
DimPlot(fad_apc, 
        reduction = "umap",
        label =T)

# find markers for every cluster compared to all remaining cells, report only the positive ones
fad_apc_markers <- FindAllMarkers(fad_apc, 
                     only.pos = TRUE)

fad_apc_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

fad_apc_markers %>%
  dplyr::filter(cluster == 23) %>%
  dplyr::filter(avg_log2FC > 1)

# Cell Type annotations
fad_apc <- RenameIdents(fad_apc, "0" = "Microglia")
fad_apc <- RenameIdents(fad_apc, "1" = "Oligodendrocytes")
fad_apc <- RenameIdents(fad_apc, "2" = "Astrocytes")
fad_apc <- RenameIdents(fad_apc, "3" = "Endothelial cells")
fad_apc <- RenameIdents(fad_apc, "4" = "Choroid plexus cells")
fad_apc <- RenameIdents(fad_apc, "5" = "Neurons")
fad_apc <- RenameIdents(fad_apc, "6" = "Neurons")
fad_apc <- RenameIdents(fad_apc, "7" = "Ependymal cells")
fad_apc <- RenameIdents(fad_apc, "8" = "Oligodendrocytes")
fad_apc <- RenameIdents(fad_apc, "9" = "Pericytes")
fad_apc <- RenameIdents(fad_apc, "10" = "Macrophages")
fad_apc <- RenameIdents(fad_apc, "11" = "Microglia")
fad_apc <- RenameIdents(fad_apc, "12" = "Fibroblasts")
fad_apc <- RenameIdents(fad_apc, "13" = "Neurons")
fad_apc <- RenameIdents(fad_apc, "14" = "Oligodendrocytes")
fad_apc <- RenameIdents(fad_apc, "15" = "T memory cells")
fad_apc <- RenameIdents(fad_apc, "16" = "Oligodendrocytes")
fad_apc <- RenameIdents(fad_apc, "17" = "Endothelial cells")
fad_apc <- RenameIdents(fad_apc, "18" = "Pericytes")
fad_apc <- RenameIdents(fad_apc, "19" = "Endothelial cells")
fad_apc <- RenameIdents(fad_apc, "20" = "Neutrophils")
fad_apc <- RenameIdents(fad_apc, "21" = "B cells")
fad_apc <- RenameIdents(fad_apc, "22" = "Microglia")
fad_apc <- RenameIdents(fad_apc, "23" = "Unknown")

DimPlot(fad_apc, 
        reduction = "umap",
        label =T)


saveRDS(fad_apc, file = "FAD_condition")
wt <- readRDS(file.choose()) #read in FAD file here

#################################################################################CELLCHAT ANALYSIS######################################################################################################################################
library(dplyr)
library(Seurat)
library(patchwork)
library(CellChat)
path <- "C:\\Users\\sujayayyagari\\OneDrive\\Desktop\\WT"


seurat_object <- readRDS(file.choose())
data.input <- GetAssayData(seurat_object, assay = "RNA", slot = "data")
labels <- Idents(seurat_object)
meta <- data.frame(group=labels, row.names = names(labels))
cellchat <- createCellChat(object = data.input, meta = meta, group.by="group")
CellChatDB <- CellChatDB.mouse 
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "group") 
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
cellchat@DB <- CellChatDB
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, raw.use = FALSE) # use the projected data
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
#For example, showing the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot.
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
pathways.show <- c("ApoE") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
#> Do heatmap based on a single object
netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.ApoE <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.ApoE[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
#> [[1]]
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
#We can also show all the significant interactions (L-R pairs) from some cell groups to other cell groups using netVisual_bubble.
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
# (2) show all the significant interactions (L-R pairs) associated with certain signaling pathways
#Plot the signaling gene expression distribution using violin/dot plot
plotGeneExpression(cellchat, signaling = "ApoE", enriched.only = TRUE, type = "violin")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
