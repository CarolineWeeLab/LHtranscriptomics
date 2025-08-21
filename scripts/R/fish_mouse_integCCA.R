library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(xlsx)
library(enrichR)
library(readxl)
library(ggrepel)
library(ggpubr)
library(ggsignif)
library(tidyverse)
library(sctransform)
library(scCustomize)
library(pheatmap)
library(openxlsx)
library(tidyr)
install.packages('devtools')
devtools::install_github('immunogenomics/presto')
install.packages("Cairo")
install.packages("rgl")

setwd("/Users/vindhyachaganty/Downloads/Use this_250225/scrna")


################################################
#
# Processing mouse dataset
#
################################################

# mouse QC
mouse <- readRDS("Mouse_Mickelsen_MergedSeuratObject.rds") # 7231 cells
mouse[["percent.mt"]] <- PercentageFeatureSet(mouse, pattern = "^mt-")
VlnPlot(mouse, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#mouse <- subset(mouse, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 & percent.mt <= 20)
mouse <- subset(mouse, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 & percent.mt <= 15)
mouse #7231 to 6882 cells for 20%, 6406 for 15% and 4176 for 10%
VlnPlot(mouse, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# mouse clustering
mouse <- NormalizeData(mouse)
mouse <- FindVariableFeatures(mouse, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mouse)
mouse <- ScaleData(mouse, features = all.genes)
mouse <- RunPCA(mouse, features = VariableFeatures(object = mouse))
ElbowPlot(mouse)
mouse <- FindNeighbors(mouse, dims = 1:10)
mouse <- FindClusters(mouse, resolution = 0.5) #17 clusters
mouse <- RunUMAP(mouse, dims = 1:10)
DimPlot(mouse, reduction = "umap", label = TRUE)
DimPlot_scCustom(mouse, aspect_ratio = 1)

mouse.markers <- FindAllMarkers(mouse, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
t <- mouse.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)
mouse <- AddModuleScore(
  object = mouse,
  features = c("Snap25", "Syp", "Tubb3", "Elavl2"),
  name = 'pan_neuronal_markers'
)

FeaturePlot(mouse, features = c("Hes1", "Hes5"))

FeaturePlot(mouse, features = c("pan_neuronal_markers1"))
library(viridis)
pal <- viridis(n = 10, option = "D")
FeaturePlot_scCustom(seurat_object = mouse, features = "pan_neuronal_markers1", aspect_ratio = 1)

new.cluster.ids <- c("Neuron", "Neuron", "Oligodendrocyte", "Neuron", "Neuron",
                     "Endothelial", "Astrocyte", "OPC", "Neuron", 
                     "OPC", "Oligodendrocyte", "Microglia", "Oligodendrocyte", "Endothelial", 
                     "Neuron", "RBC", "Astrocyte")

names(new.cluster.ids) <- levels(mouse)
mouse <- RenameIdents(mouse, new.cluster.ids)
mouse[["celltype"]] <- Idents(object = mouse)
DimPlot_scCustom(mouse, aspect_ratio = 1, label = F)

################################################
#
# Mouse neurons
#
################################################

mouseneurons <- subset(x = mouse, idents = "Neuron")
mouseneurons #4063 cells

mouseneurons <- NormalizeData(mouseneurons)
mouseneurons <- FindVariableFeatures(mouseneurons, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mouseneurons)
mouseneurons <- ScaleData(mouseneurons, features = all.genes)
mouseneurons <- RunPCA(mouseneurons, features = VariableFeatures(object = mouseneurons))
ElbowPlot(mouseneurons)
mouseneurons <- FindNeighbors(mouseneurons, dims = 1:10)
mouseneurons <- FindClusters(mouseneurons, resolution = 0.5)#14 clusters
mouseneurons <- RunUMAP(mouseneurons, dims = 1:10)
DimPlot(mouseneurons, reduction = "umap", label = TRUE)
DimPlot_scCustom(mouseneurons, aspect_ratio = 1)

mouseneurons.markers <- FindAllMarkers(mouseneurons, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
t <- mouseneurons.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.xlsx(t, file = "mouse_seurat_clusters_top20_markers.xlsx")

FeaturePlot_scCustom(seurat_object = mouseneurons, features = c("Slc17a6", "Slc32a1", "Gad1", "Gad2"), aspect_ratio = 1)
Stacked_VlnPlot(seurat_object = mouseneurons, features = c("Slc17a6", "Slc32a1", "Gad1", "Gad2"))

new.cluster.ids <- c("GABA", "GABA", "GABA", "Glut", "Glut",
                     "Glut", "GABA+Glut", "Glut", "Glut", 
                     "GABA+Glut", "GABA+Glut", "GABA", "GABA+Glut", "GABA")

names(new.cluster.ids) <- levels(mouseneurons)
mouseneurons <- RenameIdents(mouseneurons, new.cluster.ids)
mouseneurons[["celltype"]] <- Idents(object = mouseneurons)
DimPlot_scCustom(mouseneurons, aspect_ratio = 1, label = F, colors_use = c("#F6222EFF","#5A5156FF", "lightblue"))
Idents(mouseneurons) <- mouseneurons$celltype
FeaturePlot_scCustom(seurat_object = mouseneurons, features = c("Tac1", "Ghr"), aspect_ratio = 1)
# Tac1+Ghr+ double positive 


################################################
#
# Mouse GABA neurons
#
################################################
mouseGABA <- readRDS("mouseGABA.rds") 

mouseGABA <- subset(x = mouseneurons, idents = "GABA")
mouseGABA #2247 
saveRDS(mouseGABA, file = "mouseGABA.rds")

mouseGABA <- NormalizeData(mouseGABA)
mouseGABA <- FindVariableFeatures(mouseGABA, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mouseGABA)
mouseGABA <- ScaleData(mouseGABA, features = all.genes)
mouseGABA <- RunPCA(mouseGABA, features = VariableFeatures(object = mouseGABA))
ElbowPlot(mouseGABA)
mouseGABA <- FindNeighbors(mouseGABA, dims = 1:20)
mouseGABA <- FindClusters(mouseGABA, resolution = 0.5)#12 clusters
mouseGABA <- RunUMAP(mouseGABA, dims = 1:20)
MouseGABAUMAP <- DimPlot_scCustom(mouseGABA, aspect_ratio = 1, label = F)
# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Mouse/GABAUMAP.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 8, height = 6)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
MouseGABAUMAP

# Turn off the device to finalize the file output
dev.off()

mouseGABA.markers <- FindAllMarkers(mouseGABA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
t.mouse.GABA <- mouseGABA.markers %>% group_by(cluster) %>% top_n(n = 60, wt = avg_log2FC)
write.xlsx(t.mouse.GABA, file = "mouseGABA_seurat_clusters_top40_markers.xlsx")

FeaturePlot_scCustom(seurat_object = mouseGABA, features = c("Tac1", "Ghr"), aspect_ratio = 1)
FeaturePlot_scCustom(seurat_object = mouseGABA, features = c("Pmch"), aspect_ratio = 1)

mouseGABA$Tac1_Ghr_exp <- ifelse(mouseGABA@assays$RNA@data["Tac1", ] > 0 & 
                                   mouseGABA@assays$RNA@data["Ghr", ] > 0 
                                   , "Pos", "Neg")

table(mouseGABA$Tac1_Ghr_exp)

Tac1GhrcoexpUMAP_mouseGABA <- DimPlot_scCustom(mouseGABA, aspect_ratio = 1, pt.size = 2.5, label = F, group.by = "Tac1_Ghr_exp", colors_use = c("grey", "red"))

# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Mouse/Tac1GhrGsx1Adra2Coexp_mouseGABA_UMAP2.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 8, height = 6)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
Tac1GhrcoexpUMAP_mouseGABA
print(Tac1GhrcoexpUMAP_mouseGABA)
# Turn off the device to finalize the file output
dev.off()

split_violin_plot <- function(seurat_obj, features) {
  df <- FetchData(seurat_obj, vars = c(features, "ident"))
  colnames(df) <- c(features, "Cluster")
  
  # Convert to long format for ggplot
  df_long <- reshape2::melt(df, id.vars = "Cluster", variable.name = "Gene", value.name = "Expression")
  
  # Create violin plot
  p <- ggplot(df_long, aes(x = Cluster, y = Expression, fill = Cluster)) +
    geom_violin(scale = "width", alpha = 0.5, size = 0.2) +
    facet_wrap(~Gene, ncol = 1) +  # One column for each gene
    theme_classic() +
    theme(
      legend.position = "none",  # Remove legend
      strip.text = element_blank(), 
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis labels
      axis.text.y = element_text(face = "bold")  # Bold y-axis labels (optional)
    ) +
    xlab("Cluster") +  
    ylab("Expression Level")  
  
  return(p)
}
features_common_GABA_mouse <- c("Zfhx4", "Npas3", "Cabp7", "Fstl5", "Gsx1", "Tac1", "Ghr",
                          "Adra2a", "Lhx6", "Nkx2-2", "Sp9", "Cntn5", "Ctbp2", "Ddc",
                          "Th", "Slc18a2", "Gfra1", "Sst", "Tac2", "Htr1a", "Htr7", "Htr2c", "Penk", "Adcyap1",
                          "Nts", "Gal")

library(dplyr)

# Add seurat cluster info to the expression table
lognorm_expression_mouse <- FetchData(mouseGABA, vars = c("Zfhx4", "Npas3", "Cabp7", "Fstl5", "Gsx1", "Tac1", "Ghr",
  "Adra2a", "Lhx6", "Nkx2-2", "Sp9", "Cntn5", "Ctbp2", "Ddc", "Th", "Slc18a2", "Gfra1", "Sst", "Tac2", "Htr1a", "Htr7", 
  "Htr2c", "Penk", "Adcyap1", "Nts", "Gal"))

# Add cluster info
lognorm_expression$seurat_clusters <- mouseGABA$seurat_clusters

# Group by cluster and take average of each gene
mouse_cluster_avg_expression <- lognorm_expression %>%
  group_by(seurat_clusters) %>%
  summarise(across(everything(), mean), .groups = "drop")  # mean for each gene in cluster

# Save to CSV
write.csv(mouse_cluster_avg_expression, "mouse_average_expression_by_cluster.csv", row.names = FALSE)


split_plot_common_GABA_mouse <- split_violin_plot(mouseGABA, features_common_GABA_mouse)

print(split_plot_common_GABA_mouse)
# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Mouse/commonfeatures_mouseGABA.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 11, height = 15)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
split_plot_common_GABA_mouse

# Turn off the device to finalize the file output
dev.off()

################################################
#
# Mouse Glut neurons
#
################################################

mouseglut <- readRDS("mouseglut.rds") 


mouseglut <- subset(x = mouseneurons, idents = "Glut")
mouseglut #1393 cells
saveRDS(mouseglut, file = "mouseglut.rds")


mouseglut <- NormalizeData(mouseglut)
mouseglut <- FindVariableFeatures(mouseglut, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mouseglut)
mouseglut <- ScaleData(mouseglut, features = all.genes)
mouseglut <- RunPCA(mouseglut, features = VariableFeatures(object = mouseglut))
ElbowPlot(mouseglut)
mouseglut <- FindNeighbors(mouseglut, dims = 1:20)
mouseglut <- FindClusters(mouseglut, resolution = 0.5) #13 clusters
mouseglut <- RunUMAP(mouseglut, dims = 1:20)
MouseGlutUMAP <- DimPlot_scCustom(mouseglut, aspect_ratio = 1, label = F)

# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Mouse/GlutUMAP.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 8, height = 6)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
MouseGlutUMAP

# Turn off the device to finalize the file output
dev.off()

FeaturePlot_scCustom(seurat_object = mouseglut, features = c("Tac1", "Ghr"), aspect_ratio = 1)
FeaturePlot_scCustom(seurat_object = mouseglut, features = c("Pmch"), aspect_ratio = 1)

mouseglut.markers <- FindAllMarkers(mouseglut, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
t.mouse.glut <- mouseglut.markers %>% group_by(cluster) %>% top_n(n = 60, wt = avg_log2FC)
write.xlsx(t.mouse.glut, file = "mouseglut_seurat_clusters_top20_markers.xlsx")

features_common_glut_mouse <- c("Cabp7", "Meis2", "Pbx3", "Nrxn3", "Fezf1", "Cbln1", "Tac1", 
                          "Kctd12", "Bhlhe22", "Neurod2", "Mpped2", "Thsd7b", "Adk" , "Caln1", 
                          "Cck", "Syt2", "Pdyn", "Adcyap1", "Gad1", "Gad2", "Htr1a", "Htr2c", "Htr7", "Nts", "Gal")

# Add seurat cluster info to the expression table
lognorm_expression_mouse_glut <- FetchData(mouseglut, vars = c("Cabp7", "Meis2", "Pbx3", "Nrxn3", "Fezf1", "Cbln1", "Tac1", 
                                                          "Kctd12", "Bhlhe22", "Neurod2", "Mpped2", "Thsd7b", "Adk" , "Caln1", 
                                                          "Cck", "Syt2", "Pdyn", "Adcyap1", "Gad1", "Gad2", "Htr1a", "Htr2c", "Htr7", "Nts", "Gal"))

# Add cluster info
lognorm_expression_mouse_glut$seurat_clusters <- mouseglut$seurat_clusters

# Group by cluster and take average of each gene
mouseglut_cluster_avg_expression <- lognorm_expression_mouse_glut %>%
  group_by(seurat_clusters) %>%
  summarise(across(everything(), mean), .groups = "drop")  # mean for each gene in cluster

# Save to CSV
write.csv(mouseglut_cluster_avg_expression, "mouseglut_average_expression_by_cluster.csv", row.names = FALSE)

split_plot_common_glut_mouse <- split_violin_plot(mouseglut, features_common_glut_mouse)

print(split_plot_common_glut_mouse)
# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Mouse/commonfeatures2_mouseglut.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 11, height = 15)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
split_plot_common_glut_mouse

# Turn off the device to finalize the file output
dev.off()
################################################
#
# Processing zebrafish dataset
#
################################################

#Re-analysing zebrafish data with pre-filtering on Partek, 200<nFeatureRNA<6000, %mt<20, 15 and 10%
samples = c("Voraciously feeding", "Starved");

#2.Create a variable called outdir to specify your output directory:

outdir = "/Users/vindhyachaganty/Downloads/Use this_250225/scrna"
#3.Read in the feature-barcode matrices generated by the cellranger pipeline
LHA = list(); # first declare an empty list in which to hold the feature-barcode matrices
LHA[[1]] <- Read10X_h5( "V3_Fed-SCI7T051-SCI5T051_S13.h5");
LHA[[2]] <- Read10X_h5("V3_Starved-SCI7T039-SCI5T039_S11.h5");

#4.Convert each feature-barcode matrix to a Seurat object
Starvedrefed.list = list(); # First create an empty list to hold the Seurat objects
Starvedrefed.list[[1]] = CreateSeuratObject(counts = LHA[[1]], min.cells=3, min.features=200, project=samples[1]);
Starvedrefed.list[[1]][["DataSet"]] = samples[1];
Starvedrefed.list[[2]] = CreateSeuratObject(counts = LHA[[2]], min.cells=3, min.features=200, project=samples[2]);
Starvedrefed.list[[2]][["DataSet"]] = samples[2];

rm(LHA);

#5.Merge datasets
Starvedrefed <- merge(x=Starvedrefed.list[[1]], y=c(Starvedrefed.list[[2]]), add.cell.ids = c("Voraciously feeding", "Starved"), project="Zf_LHA");
str(Starvedrefed@meta.data) # examine the structure of the Seurat object meta data
saveRDS(Starvedrefed, file = sprintf("%s/MergedSeuratObject.rds", outdir));
Starvedrefed #9049 samples for 20%mt, 8744 cells for 5% mt

Starvedrefed <- Starvedrefed[!is.infinite(rowSums(Starvedrefed)),]
rm(Starvedrefed) # save some memory
fish <- readRDS("MergedSeuratObject.rds") # 9049 cells

fish[["percent.mt"]] <- PercentageFeatureSet(fish, pattern = "^mt-")
VlnPlot(fish, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
fish <- subset(fish, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 & percent.mt <= 15)
fish #9049 cells for 20%mt, 9031 cells for 15%mt, 8974 cells for 10%mt and 8744 cells for 5% mt
table(fish$orig.ident)

# fish clustering

fish <- NormalizeData(fish)
fish <- FindVariableFeatures(fish, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(fish)
fish <- ScaleData(fish, features = all.genes)
fish <- RunPCA(fish, features = VariableFeatures(object = fish))
ElbowPlot(fish)
fish <- FindNeighbors(fish, dims = 1:14) # or 10
fish <- FindClusters(fish, resolution = 0.5) #20 clusters
fish <- RunUMAP(fish, dims = 1:14) # or 10
DimPlot_scCustom(fish, aspect_ratio = 1)
fish <- JoinLayers(fish)
fish.markers <- FindAllMarkers(fish, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

fish <- AddModuleScore(
  object = fish,
  features = list(c("snap25a", "sypa", "tubb5", "elavl4"), 
                  c("snap25a", "sypa", "tubb5", "elavl3", "elavl4", "neurod1", "map2")),
  name = 'pan_neuronal_markers'
)
 
FeaturePlot_scCustom(seurat_object = fish, 
                     features = c("pan_neuronal_markers1", "pan_neuronal_markers2"), aspect_ratio = 1)
FeaturePlot_scCustom(seurat_object = fish_nonLH, 
                     features = c("her5"), aspect_ratio = 1)

t<- fish.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

Idents(fish) <- fish$seurat_clusters
new.cluster.ids <- c("Neurons", "Neurons", "Neurons", "Neurons", "Neurons",
                     "Non-neurons", "Neurons", "Neurons", "Non-neurons", 
                     "Neurons", "Neurons", "Non-neurons", "Neurons", "Non-neurons", 
                     "Non-neurons", "Non-neurons", "Non-neurons", "Non-neurons", 
                     "Non-neurons", "Neurons")

names(new.cluster.ids) <- levels(fish)
fish <- RenameIdents(fish, new.cluster.ids)
fish[["celltype"]] <- Idents(object = fish)
DimPlot_scCustom(fish, aspect_ratio = 1, label = F)
table(Idents(fish))
Plot1 <- FeaturePlot_scCustom(seurat_object = fish, features = "Gal4FFCDS", aspect_ratio = 1)
VlnPlot(fish, features = "Gal4FFCDS")
Plot2 <- FeaturePlot(fish, features = "Gal4FFCDS")

# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Gal4FFCDS_UMAP2.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 8, height = 6)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
Plot2

# Turn off the device to finalize the file output
dev.off()

################################################
#
# Fish neurons
#
################################################

fishLH <- subset(x = fish, idents = c("Neurons"))
fishLH <- subset(x = fishLH, subset = Gal4FFCDS > 0.5)
fishLH #3062 cells

# Get the names of LH cells
LH_cells <- colnames(fishLH)

# Subset all non-LH cells
fish_nonLH <- subset(x = fish, cells = setdiff(colnames(fish), LH_cells))
fish_nonLH #5969 cells out of 9031 cells

fishLH <- NormalizeData(fishLH)
fishLH <- FindVariableFeatures(fishLH, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(fishLH)
fishLH <- ScaleData(fishLH, features = all.genes)
fishLH <- RunPCA(fishLH, features = VariableFeatures(object = fishLH))
ElbowPlot(fishLH)
fishLH <- FindNeighbors(fishLH, dims = 1:20) 
fishLH <- FindClusters(fishLH, resolution = 0.5)#18 clusters
fishLH <- RunUMAP(fishLH, dims = 1:20)
DimPlot_scCustom(fishLH, aspect_ratio = 1)
Stacked_VlnPlot(seurat_object = fishLH, features = c("slc17a6b", "slc32a1", "gad1b", "gad2"))
Stacked_VlnPlot(seurat_object = fishLH, features = c("snap25a", "sypa", "tubb5", "elavl3", "elavl4", "neurod1", "map2", "zgc:136537"))

fishLH.markers <- FindAllMarkers(fishLH, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
t <- fishLH.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
View(t)
write.xlsx(t, file = "fishLH_seurat_clusters_top20_markers.xlsx")


fishLH <- AddModuleScore(
  object = fishLH,
  features = list(c("snap25a", "sypa", "tubb5", "elavl4"), 
                  c("snap25a", "sypa", "tubb5", "elavl3", "elavl4", "neurod1", "map2")),
  name = 'pan_neuronal_markers'
)
FeaturePlot_scCustom(seurat_object = fishLH, 
                     features = c("pan_neuronal_markers1", "pan_neuronal_markers2"), aspect_ratio = 1)
FeaturePlot_scCustom(seurat_object = fishLH, features = c("slc17a6b", "slc32a1", "gad1b", "gad2"), aspect_ratio = 1)

VlnPlot(fishLH, features = c("slc17a6b", "slc32a1", "gad1b", "gad2"), ncol = 2)
VlnPlot(fishLH, features = "her5")

Idents(fishLH) <- fishLH$seurat_clusters
new.cluster.ids <- c("Glut", "Glut", "Unassigned", "GABA", "Glut",
                     "Glut","Glut", "GABA", "Glut", 
                     "GABA","GABA", "GABA", "GABA","Glut", "Glut", "Glut", "GABA", "GABA")

names(new.cluster.ids) <- levels(fishLH)
fishLH <- RenameIdents(fishLH, new.cluster.ids)
fishLH[["celltype"]] <- Idents(object = fishLH)
DimPlot_scCustom(fishLH, aspect_ratio = 1, label = F)
Idents(fishLH) <- fishLH$celltype

#non LH analysis
fish_nonLH <- NormalizeData(fish_nonLH)
fish_nonLH <- FindVariableFeatures(fish_nonLH, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(fish_nonLH)
fish_nonLH <- ScaleData(fish_nonLH, features = all.genes)
fish_nonLH <- RunPCA(fish_nonLH, features = VariableFeatures(object = fish_nonLH))
ElbowPlot(fish_nonLH)
fish_nonLH <- FindNeighbors(fish_nonLH, dims = 1:20) 
fish_nonLH <- FindClusters(fish_nonLH, resolution = 0.5)#24 clusters
fish_nonLH <- RunUMAP(fish_nonLH, dims = 1:20)
DimPlot_scCustom(fish_nonLH, aspect_ratio = 1)
Stacked_VlnPlot(seurat_object = fish_nonLH, features = c("slc17a6b", "slc32a1", "gad1b", "gad2"))
Stacked_VlnPlot(seurat_object = fish_nonLH, features = c("snap25a", "sypa", "tubb5", "elavl3", "elavl4", "neurod1", "map2", "zgc:136537"))

fish_nonLH.markers <- FindAllMarkers(fish_nonLH, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
t.nonLH <- fish_nonLH.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
View(t.nonLH)
write.xlsx(t.nonLH, file = "fishnonLH_seurat_clusters_top20_markers.xlsx")

fish_nonLH <- AddModuleScore(
  object = fish_nonLH,
  features = list(c("snap25a", "sypa", "tubb5", "elavl4"), 
                  c("snap25a", "sypa", "tubb5", "elavl3", "elavl4", "neurod1", "map2")),
  name = 'pan_neuronal_markers'
)
FeaturePlot_scCustom(seurat_object = fish_nonLH, 
                     features = c("pan_neuronal_markers1", "pan_neuronal_markers2"), aspect_ratio = 1, label = TRUE)
FeaturePlot_scCustom(fish_nonLH, features = c("snap25a", "sypa", "tubb5", "elavl3", "elavl4", "neurod1", "map2"), aspect_ratio = 1)
FeaturePlot_scCustom(fish_nonLH, features = c("snap25a", "sypa", "tubb5", "elavl3", "elavl4", "neurod1", "map2"), aspect_ratio = 1)
VlnPlot(fish_nonLH, features = c("pan_neuronal_markers1", "pan_neuronal_markers2"), ncol = 2)
VlnPlot(fish_nonLH, features = c("slc17a6b", "slc32a1", "gad1b", "gad2"), ncol = 2)

new.cluster.ids <- c("Glutamatergic neurons", "GABAergic neurons", "Glutamatergic neurons", "Neuronal progenitors", "Glutamatergic neurons",
                     "GABAergic neurons", "Endothelial cells", "Glutamatergic neurons", "Glutamatergic neurons", "Astrocytes", "OPCs",
                     "GABAergic neurons", "GABAergic neurons", "Retina", "Neuronal progenitors","Microglia", 
                     "Cranial ganglion", "Glutamatergic neurons", "Endothelial cells", "Endothelial cells", "Oligodendrocytes", 
                     "Eye epidermis", "RBCs", "Rostral blood island (myeloid)")


names(new.cluster.ids) <- levels(fish_nonLH)
fish_nonLH <- RenameIdents(fish_nonLH, new.cluster.ids)
fish_nonLH[["celltype"]] <- Idents(object = fish_nonLH)
Plot1 <- DimPlot_scCustom(fish_nonLH, aspect_ratio = 1, label = TRUE, label.size = 3, repel = TRUE)
table(Idents(fish_nonLH))

#Save as SVG
install.packages("svglite")
library(svglite)
# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Zf_Glut/Zfglut_UMAP.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 8, height = 6)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
Plot1

# Turn off the device to finalize the file output
dev.off()


Idents(fish_nonLH) <- fish_nonLH$celltype

table(Idents(fish_nonLH), fish_nonLH$orig.ident)

########Make activity plots for individual clusters

# Load necessary libraries
library(Seurat)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(reshape2)

fish_nonLH[["Average"]] <- rowMeans(FetchData(fish_nonLH, vars = c("fosab", "nr4a1", "npas4a")), na.rm = TRUE)
comparisons <- list("orig.idents")
# Load necessary libraries
library(ggpubr)

# Function to calculate average expression and add it to the Seurat object
fish_nonLH[["Average"]] <- rowMeans(FetchData(fish_nonLH, vars = c("fosab", "nr4a1", "npas4a")), na.rm = TRUE)

# Function to create violin plots for each gene
split_violin_plot <- function(seurat_obj, features) {
  df <- FetchData(seurat_obj, vars = c(features, "ident", "orig.ident"))
  colnames(df) <- c(features, "Cluster", "Condition")
  
  # Convert to long format for ggplot
  df_long <- melt(df, id.vars = c("Cluster", "Condition"), variable.name = "Gene", value.name = "Expression")
  
  # Check if both conditions exist
  print(table(df_long$Condition))  # Debugging step
  
  # Create violin plot
  p <- ggplot(df_long, aes(x = Cluster, y = Expression, fill = Condition)) +
    geom_violin(scale = "width", position = position_dodge(width = 0.8), alpha = 0.5) +
    geom_boxplot(width = 0.1, position = position_dodge(width = 0.8), outlier.shape = NA, color = "black") +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 1, alpha = 0.6) +
    scale_fill_manual(values = c("Starved" = "#FDE725FF", "Voraciously feeding" = "#440154FF")) +  # Fix color mapping
    facet_wrap(~Gene, ncol = 1) +  # One column, four rows
    stat_compare_means(aes(group = Condition), method = "wilcox.test", label = "p.signif") +  # Wilcoxon test
    theme_classic() +
    theme(
      legend.position = "right",
      strip.text = element_text(size = 14, face = "bold"),  # Bold facet labels
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis labels
      axis.text.y = element_text(face = "bold")  # Bold y-axis labels (optional)
    ) +
    xlab("Cluster") +  # Correct x-axis label
    ylab("Expression Level")  # Correct y-axis label
  
  return(p)
}

# Define genes of interest
features <- c("fosab", "nr4a1", "npas4a", "Average")

# Generate the combined plot
Activityplot_ZfnonLH <- split_violin_plot(fish_nonLH, features)

#Save as SVG
install.packages("svglite")
library(svglite)
# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/FishLHandnonLH/Activityplot_ZfnonLH.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 12, height = 6)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
Activityplot_ZfnonLH

# Turn off the device to finalize the file output
dev.off()
################################################
#
#non LH GABA and Glut
#
################################################
nonLH_GABA <- subset(x = fish_nonLH, idents = c("GABAergic neurons")) #1539 neurons
nonLH_glut <- subset(x = fish_nonLH, idents = c("Glutamatergic neurons")) #2806 neurons

nonLH_GABA <- NormalizeData(nonLH_GABA)
nonLH_GABA <- FindVariableFeatures(nonLH_GABA, selection.method = "vst", nfeatures = 2000)
all.genes.GABA <- rownames(nonLH_GABA)
nonLH_GABA <- ScaleData(nonLH_GABA, features = all.genes.GABA)
nonLH_GABA <- RunPCA(nonLH_GABA, features = VariableFeatures(object = nonLH_GABA))
ElbowPlot(nonLH_GABA)
nonLH_GABA <- FindNeighbors(nonLH_GABA, dims = 1:20) 
nonLH_GABA <- FindClusters(nonLH_GABA, resolution = 0.5)#11 clusters
nonLH_GABA <- RunUMAP(nonLH_GABA, dims = 1:20)
DimPlot_scCustom(nonLH_GABA, aspect_ratio = 1)
Stacked_VlnPlot(seurat_object = nonLH_GABA, features = c("slc17a6b", "slc32a1", "gad1b", "gad2", "Gal4FFCDS"))
Stacked_VlnPlot(seurat_object = nonLH_GABA, features = c("snap25a", "sypa", "tubb5", "elavl3", "elavl4", "neurod1", "map2", "zgc:136537"))

nonLH_GABA.markers <- FindAllMarkers(nonLH_GABA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
t.nonLH_GABA<- nonLH_GABA.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
View(t.nonLH_GABA)
write.xlsx(t.nonLH_GABA, file = "fishnonLHGABA_seurat_clusters_top20_markers.xlsx")
new.cluster.ids <- c("Otpa+ population", "Ventral forebrain", "Neurons (telencephalon)", "Neuronal progenitors", "Hypothalamus_1",
                     "Unassigned", "Unassigned","Hindbrain/cranial nerve", "Optic tectum (GABAergic)",
                     "Hypothalamus_2", "Cerebellum/progenitors")
              


names(new.cluster.ids) <- levels(nonLH_GABA)
nonLH_GABA <- RenameIdents(nonLH_GABA, new.cluster.ids)
nonLH_GABA[["celltype"]] <- Idents(object = nonLH_GABA)
Plot1 <- DimPlot_scCustom(nonLH_GABA, aspect_ratio = 1, label = TRUE, label.size = 3, repel = TRUE)
table(Idents(nonLH_GABA))
Plot1
library(svglite)
# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/FishLHandnonLH/UMAP_Zf_NonLH_GABA.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 8, height = 6)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
Plot1

# Turn off the device to finalize the file output
dev.off()

nonLH_glut <- NormalizeData(nonLH_glut)
nonLH_glut <- FindVariableFeatures(nonLH_glut, selection.method = "vst", nfeatures = 2000)
all.genes.glut <- rownames(nonLH_glut)
nonLH_glut <- ScaleData(nonLH_glut, features = all.genes.glut)
nonLH_glut <- RunPCA(nonLH_glut, features = VariableFeatures(object = nonLH_glut))
ElbowPlot(nonLH_glut)
nonLH_glut <- FindNeighbors(nonLH_glut, dims = 1:20) 
nonLH_glut <- FindClusters(nonLH_glut, resolution = 0.5)#15 clusters
nonLH_glut <- RunUMAP(nonLH_glut, dims = 1:20)
DimPlot_scCustom(nonLH_glut, aspect_ratio = 1)
Stacked_VlnPlot(seurat_object = nonLH_glut, features = c("slc17a6b", "slc32a1", "gad1b", "gad2", "Gal4FFCDS"))
Stacked_VlnPlot(seurat_object = nonLH_glut, features = c("snap25a", "sypa", "tubb5", "elavl3", "elavl4", "neurod1", "map2"))

nonLH_glut.markers <- FindAllMarkers(nonLH_glut, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
t.nonLH_glut<- nonLH_glut.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
View(t.nonLH_glut)
write.xlsx(t.nonLH_glut, file = "fishnonLHGlut_seurat_clusters_top20_markers.xlsx")

new.cluster.ids <- c("Neurons (diencephalon, differentiating)", "Dorsal diencephalon (thalamus)", "Neurons", 
                     "Neuronal progenitors", "Midbrain", "Granule cells", "Hypothalamus",
                     "Dorsal telencephalon/dorsal diencephalon", "Diencephalon","Neurons (Glutamatergic)", 
                     "Telencephalon (pallium)", "Dorsal habenula", "Ventral forebrain", 
                     "Cranial ganglion", "Ventral habenula")
names(new.cluster.ids) <- levels(nonLH_glut)
nonLH_glut <- RenameIdents(nonLH_glut, new.cluster.ids)
nonLH_glut[["celltype"]] <- Idents(object = nonLH_glut)
Plot2 <- DimPlot_scCustom(nonLH_glut, aspect_ratio = 1, label = TRUE, label.size = 3, repel = TRUE)
table(Idents(nonLH_glut))
Plot2

library(svglite)
# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/FishLHandnonLH/UMAP_Zf_NonLH_Glut.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 8, height = 6)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
Plot2

# Turn off the device to finalize the file output
dev.off()
#Activity plots
# Load necessary libraries
library(Seurat)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(reshape2)

Idents(nonLH_GABA) <- nonLH_GABA$seurat_clusters
table(Idents(nonLH_GABA), nonLH_GABA$orig.ident)

comparisons <- list("orig.idents")

# Function to calculate average expression and add it to the Seurat object
nonLH_GABA[["Average"]] <- rowMeans(FetchData(nonLH_GABA, vars = c("fosab", "nr4a1", "npas4a")), na.rm = TRUE)

# Function to create violin plots for each gene
split_violin_plot <- function(seurat_obj, features) {
  df <- FetchData(seurat_obj, vars = c(features, "ident", "orig.ident"))
  colnames(df) <- c(features, "Cluster", "Condition")
  
  # Convert to long format for ggplot
  df_long <- melt(df, id.vars = c("Cluster", "Condition"), variable.name = "Gene", value.name = "Expression")
  
  # Check if both conditions exist
  print(table(df_long$Condition))  # Debugging step
  
  # Create violin plot
  p <- ggplot(df_long, aes(x = Cluster, y = Expression, fill = Condition)) +
    geom_violin(scale = "width", position = position_dodge(width = 0.8), alpha = 0.5, size = 0.2) +
    #geom_boxplot(width = 0.1, position = position_dodge(width = 0.8), outlier.shape = NA, color = "grey70") +
    #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 0.3, alpha = 0.2) +
    # Add individual points without jittering
    #geom_point(position = position_dodge(width = 0.8), 
    #size = 1, alpha = 0.3, color = "black") +  # Points with transparency
    scale_fill_manual(values = c("Starved" = "#FDE725FF", "Voraciously feeding" = "#440154FF")) +  # Fix color mapping
    facet_wrap(~Gene, ncol = 1) +  # One column, four rows
    stat_compare_means(aes(group = Condition), method = "wilcox.test", label = "p.signif", label.y = 4.5) +  # Wilcoxon test
    theme_classic() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 14, face = "bold"),  # Increase size of legend title
      legend.text = element_text(size = 12),  # Increase size of legend labels
      strip.text = element_blank(), 
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis labels
      axis.text.y = element_text(face = "bold")  # Bold y-axis labels (optional)
    ) +
    xlab("Cluster") +  # Correct x-axis label
    ylab("Expression Level")  # Correct y-axis label
  
  return(p)
}
# Define genes of interest
features <- c("fosab", "nr4a1", "npas4a", "Average")

# Generate the combined plot
Activityplot_nonLH_GABA <- split_violin_plot(nonLH_GABA, features)
print(Activityplot_nonLH_GABA)
#Save as SVG
#install.packages("svglite")
library(svglite)
# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/FishLHandnonLH/Activityplot_nonLH_GABA.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 12, height = 6)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
Activityplot_nonLH_GABA

# Turn off the device to finalize the file output
dev.off()

################################################
#
# Fish GABA neurons
#
################################################
#fishGABA <- readRDS("fishGABA.RDS")
fishGABA <- subset(x = fishLH, idents = "GABA")
fishGABA #897 cells
saveRDS(fishGABA, file = "fishGABA.rds")
fishGABA <- readRDS("fishGABA.RDS")

fishGABA <- NormalizeData(fishGABA)
fishGABA <- FindVariableFeatures(fishGABA, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(fishGABA)
fishGABA <- ScaleData(fishGABA, features = all.genes)
fishGABA <- RunPCA(fishGABA, features = VariableFeatures(object = fishGABA))
ElbowPlot(fishGABA)
fishGABA <- FindNeighbors(fishGABA, dims = 1:20)
fishGABA <- FindClusters(fishGABA, resolution = 0.5)#11 clusters
fishGABA <- RunUMAP(fishGABA, dims = 1:20)
UMAP_Zf_GABA <- DimPlot_scCustom(fishGABA, aspect_ratio = 1, label = TRUE, label.size = 5, repel = TRUE)

#Save as SVG
install.packages("svglite")
library(svglite)
# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Zf_GABA/UMAP_Zf_GABA.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 8, height = 6)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
UMAP_Zf_GABA

# Turn off the device to finalize the file output
dev.off()

#Find markers
fishGABA <- JoinLayers(fishGABA)
fishGABA.markers <- FindAllMarkers(fishGABA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
filt.fishGABAmarkers <- fishGABA.markers %>% group_by(cluster) %>% top_n(n = 60, wt = avg_log2FC)
write.xlsx(filt.fishGABAmarkers, file = "fishGABA_seurat_clusters_top60_markers.xlsx")
View(filt.fishGABAmarkers)

#Use markers to annotate clusters
FeaturePlot_scCustom(seurat_object = fishGABA, features = c("tac1","ghra"), aspect_ratio = 1)
Idents(fishGABA) <- fishGABA$seurat_clusters
#new.cluster.ids <- c("lhx1a+gata3", "tac1+ghra", "meis2a+cbln2b","adra2da+tac3+man1a1","vip+SYNPR", "penkb+htr1ab+lhx6",
#                   "ghrb+pdyn+vipr1b","th2+slc18a2", "drd3+trhr2", "ccka+bdnf", "tph2+ddc+slc18a2")

#names(new.cluster.ids) <- levels(fishGABA)
#fishGABA <- RenameIdents(fishGABA, new.cluster.ids)
#fishGABA[["celltype"]] <- Idents(object = fishGABA)
#DimPlot_scCustom(fishGABA, aspect_ratio = 1, label = TRUE, label.size = 3, repel = TRUE)
#Idents(fishGABA) <- fishGABA$celltype

#table(Idents(fishGABA), fishGABA$orig.ident)


#Activity plots
# Load necessary libraries
library(Seurat)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(reshape2)

comparisons <- list("orig.idents")

# Function to calculate average expression and add it to the Seurat object
fishGABA[["Average"]] <- rowMeans(FetchData(fishGABA, vars = c("fosab", "nr4a1", "npas4a")), na.rm = TRUE)

# Function to create violin plots for each gene
split_violin_plot <- function(seurat_obj, features) {
  df <- FetchData(seurat_obj, vars = c(features, "ident", "orig.ident"))
  colnames(df) <- c(features, "Cluster", "Condition")
  
  # Convert to long format for ggplot
  df_long <- melt(df, id.vars = c("Cluster", "Condition"), variable.name = "Gene", value.name = "Expression")
  
  # Check if both conditions exist
  print(table(df_long$Condition))  # Debugging step
  
  # Create violin plot
  p <- ggplot(df_long, aes(x = Cluster, y = Expression, fill = Condition)) +
    geom_violin(scale = "width", position = position_dodge(width = 0.8), alpha = 0.5, size = 0.2) +
    #geom_boxplot(width = 0.1, position = position_dodge(width = 0.8), outlier.shape = NA, color = "grey70") +
    #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 0.3, alpha = 0.2) +
    # Add individual points without jittering
    #geom_point(position = position_dodge(width = 0.8), 
               #size = 1, alpha = 0.3, color = "black") +  # Points with transparency
    scale_fill_manual(values = c("Starved" = "#FDE725FF", "Voraciously feeding" = "#440154FF")) +  # Fix color mapping
    facet_wrap(~Gene, ncol = 1) +  # One column, four rows
    stat_compare_means(aes(group = Condition), method = "wilcox.test", label = "p.signif", label.y = 4.5) +  # Wilcoxon test
    theme_classic() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 14, face = "bold"),  # Increase size of legend title
      legend.text = element_text(size = 12),  # Increase size of legend labels
      strip.text = element_blank(), 
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis labels
      axis.text.y = element_text(face = "bold")  # Bold y-axis labels (optional)
    ) +
    xlab("Cluster") +  # Correct x-axis label
    ylab("Expression Level")  # Correct y-axis label
  
  return(p)
}
# Define genes of interest
features <- c("fosab", "nr4a1", "npas4a", "Average")

# Generate the combined plot
Activityplot_ZfGABA <- split_violin_plot(fishGABA, features)

#Save as SVG
#install.packages("svglite")
library(svglite)
# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Zf_GABA/Activityplot_ZfGABA2.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 12, height = 6)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
Activityplot_ZfGABA

# Turn off the device to finalize the file output
dev.off()

######Plot to visualise expression of NMs and receptors expressed in each cluster

# Define genes of interest
features_NMRecs <- c( "adra2a", "adra2b", "ccka", "ghra", "ghrb","htr1ab", "htr1d", "pdyn", "penka", 
                      "penkb", "slc18a2", "sst1.1", "tac1", "th2", "tph1a", "trh", "vip")

# Generate separate plots for Neuromodulators and Receptors
split_plot_NMRecs <- split_violin_plot(fishGABA, features_NMRecs)
#Save as SVG
#install.packages("svglite")
library(svglite)
# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Zf_GABA/NMREC_plot_ZfGABA.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 11, height = 15)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
split_plot_NMRecs

# Turn off the device to finalize the file output
dev.off()

######Plot to visualise top20 markers
# Find all markers
fishGABA.markers <- FindAllMarkers(fishGABA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

filt20.fishGABAmarkers <- fishGABA.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
# Display the top 20 markers (2 per cluster)
print(n =22, filt20.fishGABAmarkers)

split_violin_plot <- function(seurat_obj, features) {
  df <- FetchData(seurat_obj, vars = c(features, "ident"))
  colnames(df) <- c(features, "Cluster")
  
  # Convert to long format for ggplot
  df_long <- reshape2::melt(df, id.vars = "Cluster", variable.name = "Gene", value.name = "Expression")
  
  # Create violin plot
  p <- ggplot(df_long, aes(x = Cluster, y = Expression, fill = Cluster)) +
    geom_violin(scale = "width", alpha = 0.5, size = 0.2) +
    facet_wrap(~Gene, ncol = 1) +  # One column for each gene
    theme_classic() +
    theme(
      legend.position = "none",  # Remove legend
      strip.text = element_blank(), 
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis labels
      axis.text.y = element_text(face = "bold")  # Bold y-axis labels (optional)
    ) +
    xlab("Cluster") +  
    ylab("Expression Level")  
  
  return(p)
}

features_common_GABA <- c("zfhx4", "NPAS3", "CABP7", "fstl5", "gsx1", "tac1", "ghra",
                     "adra2a", "lhx6", "nkx2.2a", "sp9", "cntn5", "ctbp2a", "ddc",
                     "th2", "slc18a2", "gfra1a", "sst1.1", "tac3a", "htr1ab", 
                     "htr7a", "htr2cl1", "adcyap1b", "penkb", "nts", "galn") 

library(dplyr)

# Add seurat cluster info to the expression table
lognorm_expression <- FetchData(fishGABA, vars = c(
  "zfhx4", "NPAS3", "CABP7", "fstl5", "gsx1", "tac1", "ghra",
  "adra2a", "lhx6", "nkx2.2a", "sp9", "cntn5", "ctbp2a", "ddc",
  "th2", "slc18a2", "gfra1a", "sst1.1", "tac3a", "htr1ab", 
  "htr7a", "htr2cl1", "adcyap1b", "penkb", "nts", "galn"
))

# Add cluster info
lognorm_expression$seurat_clusters <- fishGABA$seurat_clusters

# Group by cluster and take average of each gene
fish_cluster_avg_expression <- lognorm_expression %>%
  group_by(seurat_clusters) %>%
  summarise(across(everything(), mean), .groups = "drop")  # mean for each gene in cluster

# Save to CSV
write.csv(fish_cluster_avg_expression, "fishGABA_average_expression_by_cluster.csv", row.names = FALSE)

install.packages("pheatmap")      # For heatmap
install.packages("tidyverse")     # For data wrangling (optional but useful)
# Example of how to input your table manually (truncated for brevity)
expression_data <- read.csv("mouseandfishGABAavgexpforheatmap.csv", row.names = 1, check.names = FALSE)

heatmap_GABA <- pheatmap(expression_data,
         scale = "none",                  # Normalize each gene
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = FALSE,
         color = colorRampPalette(c( "white", "red"))(100),
         fontsize_row = 8,
         fontsize_col = 10,
         main = "Average Log-Normalized Expression per Cluster")

# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/GABAfishmouse_commonfeatures_heatmap.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 11, height = 9)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
heatmap_GABA

# Turn off the device to finalize the file output
dev.off()

split_plot_common_GABA <- split_violin_plot(fishGABA, features_common_GABA)
print(split_plot_common_GABA)
# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Zf_GABA/commonfeatures_ZfGABA.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 11, height = 15)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
split_plot_common_GABA

# Turn off the device to finalize the file output
dev.off()

################################################
#
# Fish Glut neurons
#
################################################
fishglut <- readRDS("fishglut.RDS")
fishglut <- subset(x = fishLH, idents = "Glut")
fishglut #1817 cells
saveRDS(fishglut, file = "fishglut.rds")

fishglut <- NormalizeData(fishglut)
fishglut <- FindVariableFeatures(fishglut, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(fishglut)
fishglut <- ScaleData(fishglut, features = all.genes)
fishglut <- RunPCA(fishglut, features = VariableFeatures(object = fishglut))
ElbowPlot(fishglut)
fishglut <- FindNeighbors(fishglut, dims = 1:20)
fishglut <- FindClusters(fishglut, resolution = 0.5) #10 clusters
fishglut <- RunUMAP(fishglut, dims = 1:20)
UMAP_Zf_glut <- DimPlot_scCustom(fishglut, aspect_ratio = 1, label = F)

# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Zf_Glut/UMAP_Zf_Glut.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 8, height = 6)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
UMAP_Zf_glut

# Turn off the device to finalize the file output
dev.off()

FeaturePlot_scCustom(seurat_object = fishglut, features = c("tac1","ghra"), aspect_ratio = 1)

fishglut <- JoinLayers(fishglut)

fishglut.markers <- FindAllMarkers(fishglut, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
filt.fishglutmarkers <- fishglut.markers %>% group_by(cluster) %>% top_n(n = 60, wt = avg_log2FC)
write.xlsx(filt.fishglutmarkers, file = "fishglut_seurat_clusters_top60_markers.xlsx")
View(filt.fishglutmarkers)
Idents(fishglut) <- fishglut$seurat_clusters

# Load necessary libraries
library(Seurat)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(reshape2)

fishglut[["Average"]] <- rowMeans(FetchData(fishglut, vars = c("fosab", "nr4a1", "npas4a")), na.rm = TRUE)
comparisons <- list("orig.idents")
# Load necessary libraries
library(ggpubr)

# Function to create violin plots for each gene
split_violin_plot <- function(seurat_obj, features) {
  df <- FetchData(seurat_obj, vars = c(features, "ident", "orig.ident"))
  colnames(df) <- c(features, "Cluster", "Condition")
  
  # Convert to long format for ggplot
  df_long <- melt(df, id.vars = c("Cluster", "Condition"), variable.name = "Gene", value.name = "Expression")
  
  # Check if both conditions exist
  print(table(df_long$Condition))  # Debugging step
  
  # Create violin plot
  p <- ggplot(df_long, aes(x = Cluster, y = Expression, fill = Condition)) +
    geom_violin(scale = "width", position = position_dodge(width = 0.8), alpha = 0.5) +
    geom_boxplot(width = 0.1, position = position_dodge(width = 0.8), outlier.shape = NA, color = "black") +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 1, alpha = 0.6) +
    scale_fill_manual(values = c("Starved" = "#FDE725FF", "Voraciously feeding" = "#440154FF")) +  # Fix color mapping
    facet_wrap(~Gene, ncol = 1) +  # One column, four rows
    #stat_compare_means(aes(group = Condition), method = "wilcox.test", label = "p.signif") +  # Wilcoxon test
    theme_classic() +
    theme(
      legend.position = "right",
      strip.text = element_text(size = 14, face = "bold"),  # Bold facet labels
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis labels
      axis.text.y = element_text(face = "bold")  # Bold y-axis labels (optional)
    ) +
    xlab("Cluster") +  # Correct x-axis label
    ylab("Expression Level")  # Correct y-axis label
  
  return(p)
}

# Define genes of interest
features <- c("fosab", "nr4a1", "npas4a", "Average")

# Generate the combined plot
Activityplot_Zf_glut <- split_violin_plot(fishglut, features)

#Save as SVG
install.packages("svglite")
library(svglite)

# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Zf_Glut/Activityplot_ZfGlut.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 12, height = 6)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
Activityplot_Zf_glut

# Turn off the device to finalize the file output
dev.off()

######Plot to visualise expression of NMs and receptors expressed in each cluster

# Define genes of interest
features_NMRecs <- c("adcyap1b","ccka","cckb", "chrm2a", "chrm4a","galr2b", "hrh3","htr1d","htr5ab",
                     "ndnf","pdyn","penka", "penkb","sst1.1","sstr1a", "tac1", "trh","vip")


#Removing ndnf, chrm2a, galr2b, htr5ab, prlhr2a
# Generate separate plots for Neuromodulators and Receptors
split_plot_NMRecs <- split_violin_plot(fishglut, features_NMRecs)

# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Zf_Glut/NMRecs_plot_ZfGlut.svg.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 11, height = 15)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
split_plot_NMRecs

# Turn off the device to finalize the file output
dev.off()
######Plot to visualise top20 markers
# Find all markers
fishglut.markers <- FindAllMarkers(fishglut, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

filt20.fishglutmarkers <- fishglut.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
# Display the top 20 markers (2 per cluster)

features_common_glut <- c("CABP7", "meis2a", "pbx3b", "nrxn3a", "fezf1", "cbln1", "tac1", 
                     "kctd12.1", "bhlhe22", "neurod2", "mpped2", "thsd7ba", "adka" , "caln1", 
                     "cckb","SYT2", "pdyn", "adcyap1b", "gad1b", "gad2","htr1ab", "htr2cl1", "htr7a", "nts", "galn") 

# Add seurat cluster info to the expression table
lognorm_expression_fish_glut <- FetchData(fishglut, vars = c("CABP7", "meis2a", "pbx3b", "nrxn3a", "fezf1", "cbln1", "tac1", 
                                                               "kctd12.1", "bhlhe22", "neurod2", "mpped2", "thsd7ba", "adka" , "caln1", 
                                                               "cckb","SYT2", "pdyn", "adcyap1b", "gad1b", "gad2","htr1ab", "htr2cl1", "htr7a", "nts", "galn"))

# Add cluster info
lognorm_expression_fish_glut$seurat_clusters <- fishglut$seurat_clusters

# Group by cluster and take average of each gene
fishglut_cluster_avg_expression <- lognorm_expression_fish_glut %>%
  group_by(seurat_clusters) %>%
  summarise(across(everything(), mean), .groups = "drop")  # mean for each gene in cluster

# Save to CSV
write.csv(fishglut_cluster_avg_expression, "fishglut_average_expression_by_cluster.csv", row.names = FALSE)

install.packages("pheatmap")      # For heatmap
install.packages("tidyverse")     # For data wrangling (optional but useful)
# Example of how to input your table manually (truncated for brevity)
expression_data <- read.csv("mouseandfishGlutavgexpforheatmap.csv", row.names = 1, check.names = FALSE)

heatmap_glut <- pheatmap(expression_data,
         scale = "none",                  # Normalize each gene
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = FALSE,
         color = colorRampPalette(c("white", "red"))(100),
         fontsize_row = 8,
         fontsize_col = 10,
         main = "Average Log-Normalized Expression per Cluster")

# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Glutfishmouse_commonfeatures_heatmap.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 11, height = 9)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
heatmap_glut

# Turn off the device to finalize the file output
dev.off()

split_plot_common_glut <- split_violin_plot(fishglut, features_common_glut)

print(split_plot_common_glut)
# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Zf_Glut/commonfeatures2.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 11, height = 15)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
split_plot_common_glut

# Turn off the device to finalize the file output
dev.off()


################################################
#
# Mouse and fish integrated (only GABA)
#
################################################

saveRDS(fishGABA, "fishGABA.RDS")
saveRDS(mouseGABA, "mouseGABA.RDS")

mouseGABA <- readRDS("mouseGABA.RDS")
fishGABA <- readRDS("fishGABA.RDS")
fishGABAmapped <- fishGABA

# Mapping_table has two columns: "Zebrafish_genename" and "Mouse_genename"
mapping_table <- read.csv("Mapping_file_BioMart.csv", stringsAsFactors = FALSE)
# Adjust mapping_table such that there is no ghrb, and all ghra maps to Ghr
#mapping_table <- mapping_table[mapping_table$Zebrafish_genename != "ghrb", ]
#mapping_table$Mouse_genename[mapping_table$Zebrafish_genename == "ghra"] <- "Ghr"

#######Vindhya's way of mapping genes to the top %match only
# Check current gene names
current_genes <- rownames(fishGABAmapped[["RNA"]])
head(current_genes)

# Create a named vector for mapping
rename_vector <- setNames(mapping_table$Mouse_genename, mapping_table$Zebrafish_genename)
print(rename_vector['ghrb'])

# Identify rows with empty or missing Zebrafish_genename
mapping_table %>% filter(is.na(Zebrafish_genename) | Zebrafish_genename == "")
# Remove rows with missing or empty Zebrafish_genename
mapping_table <- mapping_table %>% filter(!is.na(Zebrafish_genename) & Zebrafish_genename != "")

# Check the column type
str(mapping_table$`X.id..query.gene.identical.to.target.Mouse.gene`)

# Load the dplyr package
library(dplyr)

# Filter the mapping to get the mouse gene with the highest value in column X
filtered_mapping <- mapping_table %>%
  group_by(Zebrafish_genename) %>%
  slice_max(order_by = `X.id..query.gene.identical.to.target.Mouse.gene`,
            with_ties = FALSE) %>% # Select the row with the highest X
  ungroup()

# Create the named vector for renaming
rename_vector <- setNames(filtered_mapping$Mouse_genename, filtered_mapping$Zebrafish_genename)

print(rename_vector['ghrb'])

rna_data <- GetAssayData(fishGABAmapped, assay = "RNA", layer = "counts")

# Update rownames with the mapping vector
rownames(rna_data) <- ifelse(rownames(rna_data) %in% names(rename_vector), 
                             rename_vector[rownames(rna_data)], 
                             rownames(rna_data))

# Identify rows with missing names
missing_names <- which(rownames(rna_data) == "" | is.na(rownames(rna_data)))
length(missing_names)  # Count of missing names

# Remove rows with missing names
rna_data <- rna_data[!is.na(rownames(rna_data)) & rownames(rna_data) != "", ]

new_assay <- CreateAssayObject(counts = rna_data)
print(new_assay)
fishGABAmapped[["RNA"]] <- new_assay

head(rownames(fishGABAmapped[["RNA"]]))

common_genes <- intersect(rownames(mouseGABA), rownames(fishGABAmapped))
print(length(common_genes))  # 9857

#Kimberle's way of mapping names

VlnPlot(fishGABAmapped, features = c("Tac1", "Ghr"))

VlnPlot(mouseGABA, features = c("Tac1", "Ghr"))
common_genes <- intersect(rownames(mouseGABA), rownames(fishGABAmapped))
length(common_genes) #9857
mouseGABA <- subset(mouseGABA, features = common_genes)
mouseGABA
fishGABAmapped <- subset(fishGABAmapped, features = common_genes)
fishGABAmapped
fishGABAmapped@assays$RNA@meta.features <- data.frame(row.names = rownames(fishGABAmapped@assays$RNA@data))

mouseGABA <- NormalizeData(mouseGABA)
mouseGABA <- FindVariableFeatures(mouseGABA)

fishGABAmapped <- NormalizeData(fishGABAmapped)
fishGABAmapped <- FindVariableFeatures(fishGABAmapped)

features <- SelectIntegrationFeatures(object.list = list(mouseGABA, fishGABAmapped), nfeatures = 10000)
mouseGABA <- ScaleData(mouseGABA, features = features)
mouseGABA <- RunPCA(mouseGABA, features = features)
fishGABAmapped <- ScaleData(fishGABAmapped, features = features)
fishGABAmapped <- RunPCA(fishGABAmapped, features = features)

anchors <- FindIntegrationAnchors(object.list = list(mouseGABA, fishGABAmapped),
                                  anchor.features = features,
                                  reduction = "cca", k.anchor =50) # increase anchor for more integration

merged_obj_GABA <- IntegrateData(anchorset = anchors, k.weight = 40) # decrease k.weight for more integration

merged_obj_GABA <- ScaleData(merged_obj_GABA)
merged_obj_GABA <- RunPCA(merged_obj_GABA)
ElbowPlot(merged_obj_GABA)
merged_obj_GABA <- FindNeighbors(merged_obj_GABA, dims = 1:20) # try different
merged_obj_GABA <- FindClusters(merged_obj_GABA, resolution = 0.5) #12 clusters
merged_obj_GABA <- RunUMAP(merged_obj_GABA, dims = 1:20)
UMAP_GABA_integ <- DimPlot_scCustom(merged_obj_GABA, aspect_ratio = 1, label = T)

#Save as SVG
install.packages("svglite")
library(svglite)
# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Integrated/GABA/Tac1GhrcoexpbarplotZF_GABA_integ.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 8, height = 6)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
Tac1GhrcoexpbarplotZF_GABA_integ

# Turn off the device to finalize the file output
dev.off()

Tac1GhrUMAP_GABA_integ <- FeaturePlot_scCustom(merged_obj_GABA, features = c("Tac1", "Ghr"), aspect_ratio = 1)


VlnPlot(merged_obj_GABA, features = c("Tac1", "Ghr"))
VlnPlot(merged_obj_GABA, features = "Pmch")
merged_obj_GABA$Species <- case_when(
  merged_obj_GABA$DataSet %in% c("A", "B") ~ "Mouse",
  merged_obj_GABA$DataSet %in% c("Starved", "Voraciously feeding") ~ "Zebrafish",
  TRUE ~ "Unknown"  # Optional: Assign "Unknown" to unmatched values
)

# Define custom colors for each species
species_colors <- c("Mouse" = "#00FFFF",  # Cyan for Species1
                    "Zebrafish" = "#FF00FF") # Magenta for Species2
# Create the UMAP plot with the custom colors
SpeciesUMAP_GABA_integ <- DimPlot_scCustom(merged_obj_GABA, 
                                           aspect_ratio = 1, 
                                           label = F, 
                                           group.by = "Species") + 
  scale_color_manual(values = species_colors)

# Print the plot
print(SpeciesUMAP_GABA_integ)

#Save as SVG
install.packages("svglite")
library(svglite)

# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Integrated/GABA/SpeciesUMAP.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 8, height = 6)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
SpeciesUMAP_GABA_integ

# Turn off the device to finalize the file output
dev.off()

fishGABAmapped.markers <- FindAllMarkers(fishGABAmapped, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
t.fishGABAmapped <- fishGABAmapped.markers %>% group_by(cluster) %>% top_n(n = 60, wt = avg_log2FC)

mouseGABA.markers <- FindAllMarkers(mouseGABA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
t.mouseGABA <- mouseGABA.markers %>% group_by(cluster) %>% top_n(n = 60, wt = avg_log2FC)

# Extract the 'gene' column from both tables (assuming 'gene' is a column name)
fishGABA_genes <- t.fishGABAmapped$gene
mouseGABA_genes <- t.mouseGABA$gene

# Find common genes between the two columns
common_genes <- intersect(fishGABA_genes, mouseGABA_genes)

# Display the common genes
print(common_genes)

# If you want to extract rows with common genes from both tables
t.fishGABAmapped_common <- t.fishGABAmapped[t.fishGABAmapped$gene %in% common_genes, ]
t.mouseGABA_common <- t.mouseGABA[t.mouseGABA$gene %in% common_genes, ]

# Optionally, if you want to view the subsetted tables
head(t.fishGABAmapped_common)
head(t.mouseGABA_common)
# Tac1+Ghr+ double positive 

merged_obj_GABA$Tac1_Ghr_exp <- ifelse(merged_obj_GABA@assays$RNA@data["Tac1", ] > 0 & 
                                         merged_obj_GABA@assays$RNA@data["Ghr", ] > 0, "Pos", "Neg")

Tac1GhrcoexpUMAP_GABA_integ <- DimPlot_scCustom(merged_obj_GABA, aspect_ratio = 1, label = F, group.by = "Tac1_Ghr_exp", colors_use = c("grey", "red"))

#Pmch expression
merged_obj_GABA$Pmch_exp <- ifelse(merged_obj_GABA@assays$RNA@data["Pmch", ] > 1
                                   , "Pos", "Neg")

DimPlot_scCustom(merged_obj_GABA, aspect_ratio = 1, label = F, group.by = "Pmch_exp", colors_use = c("grey", "red"))

table(merged_obj_GABA$Tac1_Ghr_exp, merged_obj_GABA$Species, merged_obj_GABA$seurat_clusters)

table(Idents(merged_obj_GABA), merged_obj_GABA$Species)
prop_table <- prop.table(table(Idents(merged_obj_GABA), merged_obj_GABA$Species), margin = 1)
prop_table

# Proportion stacked bar graph
library(ggplot2)
library(dplyr)
library(tidyr)

# Convert table to a data frame while keeping row names as Cluster
prop_df <- as.data.frame.matrix(prop_table)

# Add Cluster names as a new column (assuming row names are cluster labels)
prop_df$Cluster <- rownames(prop_table)

# Convert Cluster column to numeric and remove 0 if needed
prop_df$Cluster <- as.numeric(as.character(prop_df$Cluster))

# Convert to long format (reshape the data)
prop_long <- prop_df %>%
  pivot_longer(cols = c(Mouse, Zebrafish), names_to = "Species", values_to = "Proportion")

# Ensure Cluster is a factor with correct ordering (1 to 11)
prop_long$Cluster <- factor(prop_long$Cluster, levels = sort(unique(prop_long$Cluster)))

# Convert to long format (reshape the data)
prop_long <- prop_df %>%
  pivot_longer(cols = c(Mouse, Zebrafish), names_to = "Species", values_to = "Proportion")

# Ensure Cluster is a factor to retain correct ordering
prop_long$Cluster <- as.factor(prop_long$Cluster)

# Plot stacked proportion bar graph
plot1 <- ggplot(prop_long, aes(x = Cluster, y = Proportion, fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bars
  ylab("Proportion") +
  ggtitle("Proportion Stacked Bar Graph by Cluster and Species") +
  scale_fill_manual(values = c("Mouse" = "#00FFFF", "Zebrafish" = "#FF00FF")) + 
  theme_minimal()
print(plot1)

#Save as SVG
install.packages("svglite")
library(svglite)

# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Integrated/GABA/Species_proportion_plot.svg"

# Save the plot as SVG using svglite
svglite(output_file, width =6, height = 8)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
plot1

# Turn off the device to finalize the file output
dev.off()


# Plot cluster percentages

merged_obj_GABA_mouse <- subset(merged_obj_GABA, subset = Species == "Mouse")
merged_obj_GABA_fish <- subset(merged_obj_GABA, subset = Species == "Zebrafish")

Idents(merged_obj_GABA_mouse) <- "seurat_clusters"
Idents(merged_obj_GABA_fish) <- "seurat_clusters"
all_clusters <- unique(merged_obj_GABA$seurat_clusters)

cluster_pos_percentage_mouse <- merged_obj_GABA_mouse@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(Percent_Pos = mean(Tac1_Ghr_exp == "Pos") * 100) %>%
  right_join(data.frame(seurat_clusters = all_clusters), by = "seurat_clusters") %>%  # Ensure all clusters are included
  mutate(Percent_Pos = ifelse(is.na(Percent_Pos), 0, Percent_Pos))  # Fill missing percentages with 0

Tac1Ghrcoexpbarplotmouse_GABA_integ <- ggplot(cluster_pos_percentage_mouse, aes(x = factor(seurat_clusters), y = Percent_Pos)) +
  geom_bar(stat = "identity", fill = "darkgrey") +
  theme_minimal() +
  labs(x = "Cluster", y = "% Tac1+Ghr+ (mouse)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank() )+
  coord_flip()   # Rotate x-axis labels for readability

cluster_pos_percentage_fish <- merged_obj_GABA_fish@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(Percent_Pos = mean(Tac1_Ghr_exp == "Pos") * 100) %>%
  right_join(data.frame(seurat_clusters = all_clusters), by = "seurat_clusters") %>%  # Ensure all clusters are included
  mutate(Percent_Pos = ifelse(is.na(Percent_Pos), 0, Percent_Pos))  # Fill missing percentages with 0

Tac1GhrcoexpbarplotZF_GABA_integ <- ggplot(cluster_pos_percentage_fish, aes(x = factor(seurat_clusters), y = Percent_Pos)) +
  geom_bar(stat = "identity", fill = "red") +
  theme_minimal() +
  labs(x = "Cluster", y = "% Tac1+Ghr+ (zebrafish)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank() )+
  coord_flip()   # Rotate x-axis labels for readability

merged_obj_GABA.markers <- FindAllMarkers(merged_obj_GABA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
t.GABA.integrated <- merged_obj_GABA.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.xlsx(t.GABA.integrated, file = "GABA_integratedcca_seurat_clusters_2_top20_markers.xlsx")
view(t.GABA.integrated)


# Function to create violin plots for each gene
split_violin_plot <- function(seurat_obj, features) {
  df <- FetchData(seurat_obj, vars = c(features, "ident", "Species"))
  colnames(df) <- c(features, "Cluster", "Species")
  
  # Convert to long format for ggplot
  df_long <- melt(df, id.vars = c("Cluster", "Species"), variable.name = "Gene", value.name = "Expression")
  
  # Check if both conditions exist
  print(table(df_long$Species))  # Debugging step
  
  # Create violin plot
  p <- ggplot(df_long, aes(x = Cluster, y = Expression, fill = Species)) +
    geom_violin(scale = "width", position = position_dodge(width = 0.8), alpha = 0.5, size = 0.2) +
    #geom_boxplot(width = 0.1, position = position_dodge(width = 0.8), outlier.shape = NA, color = "grey70") +
    #geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 0.3, alpha = 0.2) +
    # Add individual points without jittering
    #geom_point(position = position_dodge(width = 0.8), 
    #size = 1, alpha = 0.3, color = "black") +  # Points with transparency
    scale_fill_manual(values = c("Mouse" = "#00FFFF", "Zebrafish" = "#FF00FF")) +  # Fix color mapping
    facet_wrap(~Gene, ncol = 1) +  # One column, four rows
    #stat_compare_means(aes(group = Condition), method = "wilcox.test", label = "p.signif", label.y = 4.5) +  # Wilcoxon test
    theme_classic() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 14, face = "bold"),  # Increase size of legend title
      legend.text = element_text(size = 12),  # Increase size of legend labels
      strip.text = element_blank(), 
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis labels
      axis.text.y = element_text(face = "bold")  # Bold y-axis labels (optional)
    ) +
    xlab("Cluster") +  # Correct x-axis label
    ylab("Expression Level")  # Correct y-axis label
  
  return(p)
}


# Load necessary libraries
library(ggpubr)

# Define genes of interest
GABA_int_NMRecs <- c("Cartpt", "Ghrh", "Penk","Pmch", "Sst", "Slc18a2", "Tac1", "Th",
                     "Adra2a", "Calcr", "Chrm3", "Ghr", "Npy2r", "Sstr1","Sstr3")

GABA_int_NM <- c("Cartpt", "Ghrh", "Nmu", "Oxt", "Penk","Pmch", "Sst", "Slc18a2", "Tac1", "Th")

GABA_int_Recs <- c("Adra2a", "Agtr1a", "Calcr", "Chrm3", "Crhr2", "Drd1", "Ghr", "Lepr", "Mc3r", 
                   "Npy1r",  "Npy2r", "Npy5r", "Nr3c2", "Ntsr1","Oprd1", "Ptger4", "Sstr1", 
                   "Sstr3", "Tacr3", "Thrb")

NMRecsplot_integ_GABA <- split_violin_plot(merged_obj_GABA, GABA_int_NMRecs)
NMplot_integ_GABA <- split_violin_plot(merged_obj_GABA, GABA_int_NM)
Receptorplot_integ_GABA <- split_violin_plot(merged_obj_GABA, GABA_int_recs)

#Save as SVG
install.packages("svglite")
library(svglite)

# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Integrated/GABA/NMRecsplot2_integ_GABA.svg"

# Save the plot as SVG using svglite
svglite(output_file, width =11, height = 15)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
NMRecsplot_integ_GABA

# Turn off the device to finalize the file output
dev.off()


#####Plot to visualise top20 markers
# Find all markers
merged_obj_glut.markers <- FindAllMarkers(merged_obj_glut, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
t.glut.integrated <- merged_obj_glut.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
# Display the top 20 markers (2 per cluster)

features_top20 <- c(
  "Crhbp", "Chat", "Tubb5", "Cdkn1c", "Calb2", "Grb14", "Slc1a6", 
  "Atf3", "Six6", "Vgll2", "Samsn1", "Cenpj", "Nos1", "Pdyn", "Prkcd", 
  "Nupr1", "Cd44", "Hk2", "Ggct", "Syt10", "Cbln3", "Slc17a7", "Gem", 
  "Nuf2", "Npvf", "Stk26", "Gfra1", "Rftn1", "Tcf7", "Spx", "Slc27a5", 
  "Piezo2"
)

split_plot_top20 <- split_violin_plot(merged_obj_glut, features_top20)

# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Integrated/Glut/top20plot_integGlut.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 20, height = 40)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
split_plot_top20

# Turn off the device to finalize the file output
dev.off()

###Sankey plot to compare fish vs integrated
# Set cluster identities
Idents(merged_obj_GABA) <- merged_obj_GABA$seurat_clusters
Idents(fishGABA) <- fishGABA$seurat_clusters

# Find markers for each cluster
merged_obj_GABA.markers <- FindAllMarkers(merged_obj_GABA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
tmerged_obj_GABA <- merged_obj_GABA.markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 60)

fishGABA.markers <- FindAllMarkers(fishGABA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
tfishGABA <- fishGABA.markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 60)

common_genes <- intersect(tfishGABA$gene, tmerged_obj_GABA$gene)
length (common_genes)

# Load necessary libraries
library(dplyr)

# 1. Check the structure of the filtered_mapping data frame to ensure it's correct
head(filtered_mapping)

# 2. Assuming filt.fishGABAmarkers is your data frame and you want to map the gene names
# We will join the filt.fishGABAmarkers with the filtered_mapping file using the zebrafish_gene column
tfishGABA <- tfishGABA %>%
  left_join(filtered_mapping, by = c("gene" = "Zebrafish_genename"))  # Replace "gene" column with mouse gene names

# 3. Now you have a new column "mouse_gene" in filt.fishGABAmarkers, which contains the mouse counterparts
# If you want to replace the original gene column with the mouse_gene column, do this:
tfishGABA$gene <- tfishGABA$Mouse_genename


#Remove rows with NA from the tfishGABA dataframe
tfishGABA <- tfishGABA[complete.cases(tfishGABA), ]

# Print the updated tfishGABA dataframe
print(tfishGABA)
common_genes <- intersect(tfishGABA$gene, tmerged_obj_GABA$gene)
length (common_genes) #103 genes
print(common_genes)

# Create the Stacked_VlnPlot with the filtered gene list
Stacked_VlnPlot(
  seurat_object = merged_obj_GABA, 
  features = sort(common_genes),  # Sorted filtered gene list
  x_lab_rotate = TRUE
)

zebrafish_common_genes <- filtered_mapping %>%
  filter(Mouse_genename %in% common_genes) %>%
  pull(Zebrafish_genename)

print(zebrafish_common_genes)

################################################
 #
# Mouse and fish integrated (only Glut)
#
################################################
fishglut <- readRDS("fishglut.RDS")
mouseglut <- readRDS("mouseglut.RDS")

fishglutmapped <- fishglut

# Mapping_table has two columns: "Zebrafish_genename" and "Mouse_genename"
mapping_table <- read.csv("Mapping_file_BioMart.csv", stringsAsFactors = FALSE)
# Adjust mapping_table such that there is no ghrb, and all ghra maps to Ghr
mapping_table <- mapping_table[mapping_table$Zebrafish_genename != "ghrb", ]
mapping_table$Mouse_genename[mapping_table$Zebrafish_genename == "ghra"] <- "Ghr"

#######Vindhya's way of mapping genes to the top %match only
# Check current gene names
current_genes <- rownames(fishglutmapped[["RNA"]])
head(current_genes)

# Create a named vector for mapping
rename_vector <- setNames(mapping_table$Mouse_genename, mapping_table$Zebrafish_genename)
print(rename_vector['ghra'])

# Identify rows with empty or missing Zebrafish_genename
mapping_table %>% filter(is.na(Zebrafish_genename) | Zebrafish_genename == "")
# Remove rows with missing or empty Zebrafish_genename
mapping_table <- mapping_table %>% filter(!is.na(Zebrafish_genename) & Zebrafish_genename != "")

# Check the column type
str(mapping_table$`X.id..query.gene.identical.to.target.Mouse.gene`)

# Load the dplyr package
library(dplyr)

# Filter the mapping to get the mouse gene with the highest value in column X
filtered_mapping <- mapping_table %>%
  group_by(Zebrafish_genename) %>%
  slice_max(order_by = `X.id..query.gene.identical.to.target.Mouse.gene`,
            with_ties = FALSE) %>% # Select the row with the highest X
  ungroup()

# Create the named vector for renaming
rename_vector <- setNames(filtered_mapping$Mouse_genename, filtered_mapping$Zebrafish_genename)

print(rename_vector['ghra'])

rna_data <- GetAssayData(fishglutmapped, assay = "RNA", layer = "counts")

# Update rownames with the mapping vector
rownames(rna_data) <- ifelse(rownames(rna_data) %in% names(rename_vector), 
                             rename_vector[rownames(rna_data)], 
                             rownames(rna_data))

# Identify rows with missing names
missing_names <- which(rownames(rna_data) == "" | is.na(rownames(rna_data)))
length(missing_names)  # Count of missing names

# Remove rows with missing names
rna_data <- rna_data[!is.na(rownames(rna_data)) & rownames(rna_data) != "", ]

new_assay <- CreateAssayObject(counts = rna_data)
print(new_assay)
fishglutmapped[["RNA"]] <- new_assay

head(rownames(fishglutmapped[["RNA"]]))

common_genes <- intersect(rownames(mouseglut), rownames(fishglutmapped))
print(length(common_genes))  # 9857


VlnPlot(fishglutmapped, features = c("Tac1", "Ghr"))

VlnPlot(mouseglut, features = c("Tac1", "Ghr"))
common_genes <- intersect(rownames(mouseglut), rownames(fishglutmapped))
length(common_genes) #9857
mouseglut <- subset(mouseglut, features = common_genes)
mouseglut
fishglutmapped <- subset(fishglutmapped, features = common_genes)
fishglutmapped
fishglutmapped@assays$RNA@meta.features <- data.frame(row.names = rownames(fishglutmapped@assays$RNA@data))

mouseglut <- NormalizeData(mouseglut)
mouseglut <- FindVariableFeatures(mouseglut)

fishglutmapped <- NormalizeData(fishglutmapped)
fishglutmapped <- FindVariableFeatures(fishglutmapped)

features <- SelectIntegrationFeatures(object.list = list(mouseglut, fishglutmapped), nfeatures = 10000)
mouseglut <- ScaleData(mouseglut, features = features)
mouseglut <- RunPCA(mouseglut, features = features)
fishglutmapped <- ScaleData(fishglutmapped, features = features)
fishglutmapped <- RunPCA(fishglutmapped, features = features)

anchors <- FindIntegrationAnchors(object.list = list(mouseglut, fishglutmapped),
                                  anchor.features = features,
                                  reduction = "rpca", k.anchor = 50) # increase anchor for more integration

merged_obj_glut <- IntegrateData(anchorset = anchors, k.weight = 40) # decrease k.weight for more integration

merged_obj_glut <- ScaleData(merged_obj_glut)
merged_obj_glut <- RunPCA(merged_obj_glut)
ElbowPlot(merged_obj_glut)
merged_obj_glut <- FindNeighbors(merged_obj_glut, dims = 1:20) # try different
merged_obj_glut <- FindClusters(merged_obj_glut, resolution = 0.5)
merged_obj_glut <- RunUMAP(merged_obj_glut, dims = 1:20)
UMAP_Glut_integ <- DimPlot_scCustom(merged_obj_glut, aspect_ratio = 1, label = F)

#Save as SVG
install.packages("svglite")
library(svglite)
# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Integrated/Glut/Tac1GhrcoexpbarplotZF_Glut_integ.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 8, height = 6)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
Tac1GhrcoexpbarplotZF_Glut_integ

# Turn off the device to finalize the file output
dev.off()

Tac1GhrUMAP_Glut_integ <- FeaturePlot_scCustom(merged_obj_glut, features = c("Tac1", "Ghr"), aspect_ratio = 1)

merged_obj_glut$Species <- case_when(
  merged_obj_glut$DataSet %in% c("A", "B") ~ "Mouse",
  merged_obj_glut$DataSet %in% c("Starved", "Voraciously feeding") ~ "Zebrafish",
  TRUE ~ "Unknown"  # Optional: Assign "Unknown" to unmatched values
)

# Define custom colors for each species
species_colors <- c("Mouse" = "#00FFFF",  # Cyan for Species1
                    "Zebrafish" = "#FF00FF") # Magenta for Species2
# Create the UMAP plot with the custom colors
SpeciesUMAP_glut_integ <- DimPlot_scCustom(merged_obj_glut, 
                                           aspect_ratio = 1, 
                                           label = F, 
                                           group.by = "Species") + 
  scale_color_manual(values = species_colors)

# Print the plot
print(SpeciesUMAP_glut_integ)

#Save as SVG
install.packages("svglite")
library(svglite)

# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Integrated/Glut/SpeciesUMAP.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 8, height = 6)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
SpeciesUMAP_glut_integ

# Turn off the device to finalize the file output
dev.off()


# Tac1+Ghr+ double positive 

merged_obj_glut$Tac1_Ghr_exp <- ifelse(merged_obj_glut@assays$RNA@data["Tac1", ] > 0 & 
                                         merged_obj_glut@assays$RNA@data["Ghr", ] > 0, "Pos", "Neg")

Tac1GhrcoexpUMAP_GABA_integ <- DimPlot_scCustom(merged_obj_glut, aspect_ratio = 1, label = F, group.by = "Tac1_Ghr_exp", colors_use = c("grey", "red"))
table(merged_obj_glut$Tac1_Ghr_exp, merged_obj_glut$Species, merged_obj_glut$seurat_clusters)
table(Idents(merged_obj_glut), merged_obj_glut$Species)
prop_table <- prop.table(table(Idents(merged_obj_glut), merged_obj_glut$Species), margin = 1)
prop_table

#Hcrt expression
merged_obj_glut$Hcrt_exp <- ifelse(merged_obj_glut@assays$RNA@data["Hcrt", ] > 1.5
                                   , "Pos", "Neg")

VlnPlot(merged_obj_glut, features = "Hcrt")

DimPlot_scCustom(merged_obj_glut, aspect_ratio = 1, label = F, group.by = "Hcrt_exp", colors_use = c("grey", "red"), split.by = "Species")
table(merged_obj_glut$Hcrt_exp, merged_obj_glut$Species, merged_obj_glut$seurat_clusters)
# Plot cluster percentages

merged_obj_glut_mouse <- subset(merged_obj_glut, subset = Species == "Mouse")
merged_obj_glut_fish <- subset(merged_obj_glut, subset = Species == "Zebrafish")

Idents(merged_obj_glut_mouse) <- "seurat_clusters"
Idents(merged_obj_glut_fish) <- "seurat_clusters"
all_clusters <- unique(merged_obj_glut$seurat_clusters)

cluster_pos_percentage <- merged_obj_glut_mouse@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(Percent_Pos = mean(Tac1_Ghr_exp == "Pos") * 100) %>%
  right_join(data.frame(seurat_clusters = all_clusters), by = "seurat_clusters") %>%  # Ensure all clusters are included
  mutate(Percent_Pos = ifelse(is.na(Percent_Pos), 0, Percent_Pos))  # Fill missing percentages with 0

Tac1Ghrcoexpbarplotmouse_glut_integ <- ggplot(cluster_pos_percentage, aes(x = factor(seurat_clusters), y = Percent_Pos)) +
  geom_bar(stat = "identity", fill = "darkgrey") +
  theme_minimal() +
  labs(x = "Cluster", y = "% Tac1+Ghr+ (mouse)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank() )+
  coord_flip()   # Rotate x-axis labels for readability

cluster_pos_percentage <- merged_obj_glut_fish@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(Percent_Pos = mean(Tac1_Ghr_exp == "Pos") * 100) %>%
  right_join(data.frame(seurat_clusters = all_clusters), by = "seurat_clusters") %>%  # Ensure all clusters are included
  mutate(Percent_Pos = ifelse(is.na(Percent_Pos), 0, Percent_Pos))  # Fill missing percentages with 0

Tac1GhrcoexpbarplotZF_Glut_integ <- ggplot(cluster_pos_percentage, aes(x = factor(seurat_clusters), y = Percent_Pos)) +
  geom_bar(stat = "identity", fill = "red") +
  theme_minimal() +
  labs(x = "Cluster", y = "% Tac1+Ghr+ (zebrafish)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank() )+
  coord_flip()   # Rotate x-axis labels for readability

merged_obj_glut.markers <- FindAllMarkers(merged_obj_glut, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
t.glut.integrated <- merged_obj_glut.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
write.xlsx(t.glut.integrated, file = "Glut_integratedcca_seurat_clusters_2_top60_markers.xlsx")
view(t.glut.integrated)

table(Idents(merged_obj_glut), merged_obj_glut$Species)
prop_table <- prop.table(table(Idents(merged_obj_glut), merged_obj_glut$Species), margin = 1)
prop_table

# Proportion stacked bar graph
library(ggplot2)
library(dplyr)
library(tidyr)

# Convert table to a data frame while keeping row names as Cluster
prop_df <- as.data.frame.matrix(prop_table)

# Add Cluster names as a new column (assuming row names are cluster labels)
prop_df$Cluster <- rownames(prop_table)

# Convert Cluster column to numeric and remove 0 if needed
prop_df$Cluster <- as.numeric(as.character(prop_df$Cluster))

# Convert to long format (reshape the data)
prop_long <- prop_df %>%
  pivot_longer(cols = c(Mouse, Zebrafish), names_to = "Species", values_to = "Proportion")

# Ensure Cluster is a factor with correct ordering (1 to 11)
prop_long$Cluster <- factor(prop_long$Cluster, levels = sort(unique(prop_long$Cluster)))

# Convert to long format (reshape the data)
prop_long <- prop_df %>%
  pivot_longer(cols = c(Mouse, Zebrafish), names_to = "Species", values_to = "Proportion")

# Ensure Cluster is a factor to retain correct ordering
prop_long$Cluster <- as.factor(prop_long$Cluster)

# Plot stacked proportion bar graph
plot1 <- ggplot(prop_long, aes(x = Cluster, y = Proportion, fill = Species)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bars
  ylab("Proportion") +
  ggtitle("Proportion Stacked Bar Graph by Cluster and Species") +
  scale_fill_manual(values = c("Mouse" = "#00FFFF", "Zebrafish" = "#FF00FF")) + 
  theme_minimal()
print(plot1)

#Save as SVG
#install.packages("svglite")
library(svglite)

# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Integrated/Glut/Species_proportion_plot.svg"

# Save the plot as SVG using svglite
svglite(output_file, width =6, height = 8)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
plot1

# Turn off the device to finalize the file output
dev.off()

#Plotting NMs and receptors
comparisons <- list("Species")
# Load necessary libraries
library(ggpubr)
# Function to create violin plots for each gene
split_violin_plot <- function(seurat_obj, features) {
  df <- FetchData(seurat_obj, vars = c(features, "ident", "Species"))
  colnames(df) <- c(features, "Cluster", "Species")
  
  # Convert to long format for ggplot
  df_long <- melt(df, id.vars = c("Cluster", "Species"), variable.name = "Gene", value.name = "Expression")
  
  # Check if both conditions exist
  print(table(df_long$Species))  # Debugging step
  
  # Create violin plot
  p <- ggplot(df_long, aes(x = Cluster, y = Expression, fill = Species)) +
    geom_violin(scale = "width", position = position_dodge(width = 0.8), alpha = 0.5) +
    geom_boxplot(width = 0.1, position = position_dodge(width = 0.8), outlier.shape = NA, color = "black") +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), size = 1, alpha = 0.6) +
    scale_fill_manual(values = c("Mouse" = "#00FFFF", "Zebrafish" = "#FF00FF")) +  # Fix color mapping
    facet_wrap(~Gene, ncol = 1) +  # One column, four rows
    #stat_compare_means(aes(group = Condition), method = "wilcox.test", label = "p.signif") +  # Wilcoxon test
    theme_classic() +
    theme(
      legend.position = "right",
      strip.text = element_text(size = 14, face = "bold"),  # Bold facet labels
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),  # Bold x-axis labels
      axis.text.y = element_text(face = "bold")  # Bold y-axis labels (optional)
    ) +
    xlab("Cluster") +  # Correct x-axis label
    ylab("Expression Level")  # Correct y-axis label
  
  return(p)
}

# Define genes of interest
Glut_int_NMRecs <- c("Cck", "Gad1", "Hcrt", "Npy","Pdyn", "Penk","Pmch",
                     "Slc17a7", "Sst", "Tac1", "Th", "Esrrg", 
                     "Lepr", "Ntsr1",  "Sstr1")
Glut_int_NMs <- c("Adm", "Cck", "Chat", "Crh", "Crhbp", "Gad1", "Gad2", "Hcrt", "Igf2", "Igfbp5", 
                  "Ndnf", "Npb", "Npy", "Oxt","Pdyn", "Penk", "Pomc","Slc17a7", "Sst", "Tac1", "Tac2", "Th" )
Glut_int_recs <- c("Adrb3", "Avpr1a", "Cckbr","Chrna6", "Crhr2", "Drd2", "Esrrg", "Htr1a", "Il13ra1", 
                   "Lepr", "Lrp2", "Mc3r", "Ntsr1", "Pgr", "Rxfp3", "Sstr1", "Tacr3", "Trhr2", "Vipr1")

# Subset fishGABA.markers for only the genes in features_1_nms
subset_Glutrec_markers <- merged_obj_glut.markers[merged_obj_glut.markers$gene %in% Glut_int_recs, c("gene", "avg_log2FC", "p_val_adj")]

# Sort the subset by avg_log2FC in descending order
subset_Glutrec_markers <- subset_Glutrec_markers[order(subset_Glutrec_markers$avg_log2FC, decreasing = TRUE), ]

# Print the sorted table
print(subset_Glutrec_markers)

# Generate the combined plot
NMRecplot_integ_glut <- split_violin_plot(merged_obj_glut, Glut_int_NMRecs)
NMplot_integ_glut <- split_violin_plot(merged_obj_glut, Glut_int_NMs)
Receptorplot_integ_glut <- split_violin_plot(merged_obj_glut, Glut_int_recs)

#Save as SVG
install.packages("svglite")
library(svglite)

# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Integrated/Glut/NMReceptorplot_integ_glut.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 11, height = 15)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
NMRecplot_integ_glut

# Turn off the device to finalize the file output
dev.off()

#####Plot to visualise top20 markers
# Find all markers
merged_obj_GABA.markers <- FindAllMarkers(merged_obj_GABA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
t.GABA.integrated <- merged_obj_GABA.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
# Display the top 20 markers (2 per cluster)

features_top20_GABA <- c(
  "Sp8", "Gsc", "Slc17a6", "Tmem163", "Syt2", "Scn4b", "Rorb", 
  "Ucp2", "Tcf7l2", "Nkd1", "Bhlhe22", "Lhx6", "Cartpt", "Anxa2", 
  "Pthlh", "Crym", "Dlx4", "Ghrh", "Tmem114", "Pcsk5", "Nmu", "Nfe2l2", 
  "Foxa1", "Gadd45b")

split_plot_top20_GABA <- split_violin_plot(merged_obj_GABA, features_top20_GABA)

# Define the output file path
output_file <- "/Users/vindhyachaganty/Downloads/Use this_250225/scrna/Integrated/GABA/top20plot_integGABA.svg"

# Save the plot as SVG using svglite
svglite(output_file, width = 20, height = 40)  # Adjust width and height as needed

# Print the plot (this will send it to the SVG file)
split_plot_top20_GABA

# Turn off the device to finalize the file output
dev.off()


# Set cluster identities
Idents(merged_obj_glut) <- merged_obj_glut$seurat_clusters
Idents(fishglut) <- fishglut$seurat_clusters

# Find markers for each cluster
merged_obj_glut.markers <- FindAllMarkers(merged_obj_glut, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
tmerged_obj_glut <- merged_obj_glut.markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 60)

fishglut.markers <- FindAllMarkers(fishglut, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
tfishglut <- fishglut.markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 60)

common_genes <- intersect(tfishglut$gene, tmerged_obj_glut$gene)
length (common_genes)

# Load necessary libraries
library(dplyr)

# 1. Check the structure of the filtered_mapping data frame to ensure it's correct
head(filtered_mapping)

# 2. Assuming filt.fishGABAmarkers is your data frame and you want to map the gene names

tfishglut <- tfishglut %>%
  left_join(filtered_mapping, by = c("gene" = "Zebrafish_genename"))  # Replace "gene" column with mouse gene names

# 3. Now you have a new column "mouse_gene" in filt.fishGABAmarkers, which contains the mouse counterparts
# If you want to replace the original gene column with the mouse_gene column, do this:
tfishglut$gene <- tfishglut$Mouse_genename

#Remove rows with NA from the tfishGABA dataframe
tfishglut <- tfishglut[complete.cases(tfishglut), ]

# Print the updated tfishGABA dataframe
print(tfishglut)
