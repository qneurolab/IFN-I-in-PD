######################################################################################
## Data process, annotation, and AUCell ##
######################################################################################

#### 01. Data process ####

library(Seurat)
library(data.table)
library(future)
library(ggplot2)
library(tidyverse)

# Set future plan for parallel processing
plan("multisession", workers = 24)
options(future.globals.maxSize = 1000 * 1024^2) # Adjusted maxSize to a realistic value

# Load data
counts <- fread("./data/UMI.tsv")
genes <- fread("./data/genes.tsv")
load("~/ref/annotation.v40.Rdata") # Annotation data

# Create Seurat Object
scRNA <- CreateSeuratObject(counts = counts)
rownames(scRNA) <- genes$gene # Simplified gene renaming

# Define a function for extracting gene names based on prefix
extract_genes <- function(anno, prefix) {
  unique(anno$gene_name[grep(prefix, anno$gene_name, ignore.case = TRUE)])
}

# Extract mito and ribo genes using the function
mito_genes <- extract_genes(mRNA_anno, "^MT-")
ribo_genes <- extract_genes(mRNA_anno, "^Rp[sl]")

# Annotation and cleanup
ids <- bind_rows(mRNA_anno, lnc_anno) %>%
  select(SYMBOL, ENSEMBL) %>%
  distinct() %>%
  mutate(ENSEMBL = str_sub(ENSEMBL, 1, 15)) %>%
  filter(ENSEMBL %in% rownames(scRNA))

# Match and rename based on ENSEMBL IDs
pos <- match(ids$ENSEMBL, rownames(scRNA))
scRNA <- scRNA[pos,]
rownames(scRNA) <- ids$SYMBOL[pos]

# Parse metadata from column names
col_data <- str_split_fixed(colnames(scRNA), "_", n = Inf)
scRNA$orig.ident <- col_data[, 2]
scRNA$Sample <- fread("./data/barcodes.tsv")$patient

# Determine group based on Sample
scRNA$group <- if_else(str_detect(scRNA$Sample, "C"), "CO", "PD")

# Ensure meta.features has row names
scRNA[["RNA"]]@meta.features <- data.frame(row.names = rownames(scRNA[["RNA"]]))

# Save processed data
save(scRNA, file = "./data/scRNA_processed.Rdata")

#### 02. QC ####
library(tidyverse)
library(Seurat)
library(vctrs)
library(viridis)
library(future)

# Set up parallel processing
plan("multisession", workers = 16)
options(future.globals.maxSize = 1000 * 1024^2) # Adjusted to a more realistic value

load("data/scRNA_raw.Rdata")

# Define QC thresholds
nFeature_lower <- 200
nFeature_upper <- 6000
nCount_lower <- 200
nCount_upper <- 20000

# Calculate mitochondrial, hemoglobin, and ribosomal percentages
scRNA <- PercentageFeatureSet(scRNA, pattern = "^MT-", col.name = "pMT")
scRNA <- PercentageFeatureSet(scRNA, pattern = "^HBA|^HBB", col.name = "pHB")
scRNA <- PercentageFeatureSet(scRNA, pattern = "^RPS|^RPL", col.name = "pRP")

# Function to plot violin plots and save them
plot_and_save_vln <- function(scRNA, filename, height = 6, width = 10) {
  group_color <- c("#3979F2", "#DF5233")
  VlnPlot(
    object = scRNA,
    features = c("nFeature_RNA", "nCount_RNA"),
    group.by = "group",
    cols = group_color,
    log = TRUE,
    pt.size = 0
  ) + NoLegend()
  
  dir.create("Figure", showWarnings = FALSE) # Ensure directory exists
  ggsave(filename, height = height, width = width)
}

# Plot before QC
plot_and_save_vln(scRNA, "Figure/boxplot_before_QC.pdf")

# Apply QC filters
scRNA <- subset(
  scRNA,
  subset = nCount_RNA > nCount_lower &
    nCount_RNA < nCount_upper &
    nFeature_RNA > nFeature_lower &
    nFeature_RNA < nFeature_upper
)

# Plot after QC
plot_and_save_vln(scRNA, "Figure/boxplot_after_QC.pdf")

# Perform CellCycleScoring
scRNA <- CellCycleScoring(
  scRNA,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes
)

# Save the filtered data
save(scRNA, file = "data/scRNA_filtered.Rda")

#### 03. Integration ####
library(Seurat)
library(future)

# Set up parallel processing
plan(multisession, workers = 16)
options(future.globals.maxSize = 1000 * 1024^2)  # Adjust to a more realistic value

# Load pre-filtered data
load("data/scRNA.filter.Rda")

# Split the Seurat object by 'Sample'
scRNA <- SplitObject(scRNA, split.by = "Sample")

# Perform SCTransform on each split object
scRNA <- lapply(scRNA, function(x) {
  SCTransform(
    x,
    variable.features.n = 3000,
    vars.to.regress = c("S.Score", "G2M.Score"),
    verbose = FALSE
  )
})

# Integration steps
features <- SelectIntegrationFeatures(object.list = scRNA, nfeatures = 3000)
scRNA <- PrepSCTIntegration(object.list = scRNA, anchor.features = features)

AnchorSet <- FindIntegrationAnchors(object.list = scRNA,
                                    normalization.method = "SCT",
                                    anchor.features = features)

# Save the AnchorSet before integration to avoid losing it in case of errors
save(AnchorSet, file = "data/AnchorSet.Rda")

# Integrate data
scRNA <- IntegrateData(anchorset = AnchorSet, normalization.method = "SCT")

# Save the integrated data
save(scRNA, file = "data/scRNA_SCT.Rda")

#### 04. PCA.FindNeighbors ####
library(Seurat)
library(clustree)
library(ggplot2)

# Load processed data
load(file = "data/scRNA_SCT.Rda")

# Run PCA
scRNA <- RunPCA(scRNA, npcs = 50, verbose = TRUE)
ElbowPlot(scRNA, ndims = 50)

# Dimensionality reduction: t-SNE and UMAP
dims <- 1:20
scRNA <- RunTSNE(scRNA, dims = dims, verbose = FALSE)
scRNA <- RunUMAP(scRNA, dims = dims, verbose = FALSE)

# Visualization function for t-SNE and UMAP
visualize_reduction <- function(scRNA, reduction_type) {
  DimPlot(
    object = scRNA,
    reduction = reduction_type,
    group.by = "Sample",
    dims = c(1, 2),
    shuffle = TRUE,
    label = TRUE,
    label.size = 4,
    label.color = "black",
    label.box = TRUE
  )
}

# Visualize t-SNE and UMAP
visualize_reduction(scRNA, "tsne")
visualize_reduction(scRNA, "umap")

# Neighbors and clustering
scRNA <- FindNeighbors(scRNA, dims = dims)
resolutions <- c(1, 2, 4, 8)

# Loop through resolutions for clustering and visualizing t-SNE plots
for (res in resolutions) {
  scRNA <- FindClusters(scRNA, resolution = res, algorithm = 1, random.seed = 2022)
  DimPlot(scRNA, reduction = "tsne", label = TRUE) +
    labs(title = paste0("Resolution: ", res)) +
    theme_minimal()
}

# clustree visualization
p_cutree <- clustree(scRNA@meta.data, prefix = "integrated_snn_res.")
print(p_cutree)

# Set identity class
Idents(scRNA) <- scRNA$integrated_snn_res.2

# Visualize and save cell clusters
p_dim_1 <- DimPlot(scRNA, reduction = "tsne", group.by = "integrated_snn_res.2", label = TRUE, pt.size = 0.5)
print(p_dim_1)

# Ensure the Figure directory exists
if (!dir.exists("Figure")) dir.create("Figure")
ggsave("Figure/cell_clusters.pdf", plot = p_dim_1, height = 6, width = 8)

# Save final scRNA object
save(scRNA, file = "data/scRNA_cluster.Rda")

#### 05. cell anno ####
library(Seurat)
library(RColorBrewer)
library(SCP)
library(ggplot2)

# Load the clustered data
load(file = "data/scRNA_cluster.Rda")

# Define cell colors and main markers for DotPlot
mycols <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#E31A1C", "#FB9A99",
            "#33A02C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")

mainmarkers <- c("SLC17A6", "GAD2", "GRIK1", "TH", "SLC6A3", "CD74", "CSF3R",
                 "AQP4", "SLC1A2", "VCAN", "OLIG1", "MOBP", "MOG", "FOXJ1",
                 "PDGFRB", "CLDN5", "ABCB1")

# Set default assay
DefaultAssay(scRNA) <- "RNA"

# First DotPlot
DotPlot(scRNA, features = mainmarkers, dot.scale = 5, group.by = "integrated_snn_res.2") +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, face = "bold", hjust = 0.5),
        axis.text.y = element_text(face = "bold"))

# Annotate cell types based on cluster information
table(scRNA$integrated_snn_res.2)
CellType <- data.frame(cluster = 0:c(length(table(scRNA$integrated_snn_res.2)) - 1),
                       CellType = "NA")
CellType[CellType$cluster %in% c(8,13,14), 2] <- "Astrocytes"
CellType[CellType$cluster %in% c(3,12), 2] <- "Microglia"
CellType[CellType$cluster %in% c(1:2,4:7,9,11,15:16,18,20,27:28,31), 2] <- "ODCs"
CellType[CellType$cluster %in% c(0), 2] <- "OPCs"
CellType[CellType$cluster %in% c(22), 2] <- "DaNs"
CellType[CellType$cluster %in% c(19,25), 2] <- "Inhib"
CellType[CellType$cluster %in% c(17,26,29), 2] <- "Excit"
CellType[CellType$cluster %in% c(10), 2] <- "Endothelial"
CellType[CellType$cluster %in% c(23), 2] <- "Ependymal"
CellType[CellType$cluster %in% c(21,24,30), 2] <- "Pericytes"

scRNA$CellType = "NA"
Idents(scRNA) = scRNA$integrated_snn_res.2
for(i in 1:nrow(CellType)){
  scRNA@meta.data[which(scRNA@active.ident == CellType$cluster[i]), "CellType"] <- CellType$CellType[i]}
A <- as.data.frame(table(scRNA$CellType))
scRNA$CellType=factor(scRNA$CellType,
                      levels = c("Excit","Inhib","DaNs","Microglia","Astrocytes",
                                 "OPCs","ODCs", "Ependymal","Pericytes","Endothelial"))

# DimPlot grouped by new CellType
DimPlot(scRNA, reduction = "tsne", group.by = "CellType", label = TRUE, cols = mycols) +
  theme_minimal() +
  NoLegend()

# CellDimPlot using SCP
CellDimPlot(srt = scRNA, group.by = c("CellType"), pt.size = 0.01, split.by = "group",
            reduction = "TSNE", theme_use = "theme_blank", palcolor = mycols)

DotPlot(scRNA, features = mainmarkers, dot.scale = 1, group.by = "CellType") +
  scale_colour_gradient2(low = "gray", high = "black", midpoint = 0.4) +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, face = "bold", hjust = 1),
        axis.text.y = element_text(face = "bold"))

# Save the annotated scRNA-seq object
save(scRNA, file = "data/scRNA_anno.Rda")

# FindAllMarkers
library(Seurat)
library(dplyr)

DefaultAssay(scRNA) <- "RNA"

# Find cluster markers
cluster_markers <- FindAllMarkers(scRNA, min.pct = 0.25, logfc.threshold = 0.25)
save(cluster_markers, file = "data/all_cluster_markers.Rdata")
write.table(cluster_markers, file = "data/all_cluster_markers.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Filter significant markers
degs_sig <- cluster_markers %>%
  filter(pct.1 > 0.5, pct.2 < 0.5, p_val_adj <= 0.01, abs(avg_log2FC) > 2) %>%
  arrange(cluster, desc(avg_log2FC))

# Save significantly differentially expressed genes
write.csv(degs_sig, file = "data/cluster_degs_sig_pct0.5.csv")

# Identify top 10 markers for each cluster and create a heatmap
DefaultAssay(scRNA) <- "integrated"
top10 <- cluster_markers %>%
  group_by(cluster) %>%
  top_n(10, wt = avg_log2FC)

# Subset cells if necessary and plot heatmap
sub_scRNA <- subset(scRNA, downsample = 200)
DoHeatmap(sub_scRNA, features = top10$gene) + NoLegend()

DefaultAssay(scRNA) <- "RNA"
Idents(scRNA) <- scRNA$CellType

# Find cell type markers
cell_markers <- FindAllMarkers(scRNA, min.pct = 0.25, logfc.threshold = 0.25)
save(cell_markers, file = "data/all_cell_markers.Rdata")
write.table(cell_markers, file = "data/all_cell_markers.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# heatmap
library(dplyr)
library(pheatmap)

# Load the annotated single-cell RNA-seq data and cell markers
load("data/scRNA_anno.Rda")
load("data/all_cell_markers.Rdata")

# Filter for genes marked as DEGs in only one cell type to reduce redundancy
cell_markers_filtered <- cell_markers %>%
  filter(is.finite(avg_log2FC)) %>%
  filter(gene %in% names(table(gene)[table(gene) < 2]))

# Select the top 50 DEGs for each cluster based on log2FC
degs_top50 <- cell_markers_filtered %>%
  group_by(cluster) %>%
  top_n(50, wt = avg_log2FC) %>%
  ungroup() %>%
  arrange(cluster, desc(avg_log2FC))

# Prepare data for the heatmap
avgData <- t(sapply(unique(degs_top50$gene), function(gene) {
  rowMeans(scRNA@assays$RNA@data[gene, , drop = FALSE], na.rm = TRUE)
}))

# Normalize and prepare heatmap data
phData <- scale(avgData)  # Auto scales to z-scores
rownames(phData) <- unique(degs_top50$gene)

# Define colors for heatmap
mycols <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#E31A1C", "#FB9A99",
            "#33A02C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")
cluster_colors <- setNames(mycols, levels(scRNA$CellType))

# Create the heatmap
pheatmap(phData, 
         color = colorRampPalette(c("darkblue", "white", "red3"))(99), 
         scale = "row", 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         annotation_row = data.frame(CellType = degs_top50$cluster),
         annotation_colors = list(CellType = cluster_colors))

# Filter significantly differentially expressed genes for further analysis
degs_sig <- cell_markers_filtered %>%
  filter(pct.1 > 0.5, pct.2 < 0.5, p_val_adj <= 0.01, abs(avg_log2FC) > 5) %>%
  arrange(cluster, desc(avg_log2FC))

# Save genes for intersections and significant markers
genes_for_interons <- intersect(rownames(scRNA@assays$integrated), degs_sig$gene)
write.table(genes_for_interons, file = "data/degs_for_interons.txt", row.names = FALSE, col.names = "Gene", quote = FALSE)
save(cell_markers_filtered, degs_sig, genes_for_interons, file = "data/degs_sig_markers.Rdata")

# Save significant DEGs to a CSV file
write.csv(degs_sig, file = "data/output_degs_sig_pct0.5.csv", row.names = FALSE)


# enrichment

# degs_sig for each cell type
for (clustername in levels(scRNA@active.ident)){
  write.table(degs_sig[degs_sig$cluster == clustername, ]$gene, paste0("data/for_enrich/","output_deg_", clustername, ".txt"), row.names = F, col.names = "Gene", quote = F)
}

for (clustername in levels(scRNA@active.ident)){
  write.csv(degs_sig[degs_sig$cluster == clustername, ], paste0("data/for_enrich/","output_deg_pval", clustername, ".csv"), row.names = F, col.names = "Gene", quote = F)
}

# terms used
go_used <- list(
  Excit = c("anterograde trans-synaptic signaling", "cell junction organization", "modulation of chemical synaptic transmission"),
  Inhib = c("modulation of chemical synaptic transmission", "neuron projection morphogenesis", "regulation of membrane potential"),
  DaNs = c("anterograde trans-synaptic signaling", "vesicle-mediated transport in synapse", "neurotransmitter secretion"),
  Microglia = c("leukocyte activation", "regulation of immune system process", "secretion"),
  Astrocytes = c("neuron development", "neuron projection development", "multicellular organismal signaling"),
  OPCs = c("synaptic signaling", "regulation of neuron projection development", "regulation of membrane potential"),
  ODCs = c("axon ensheathment", "ensheathment of neurons", "cytoskeleton organization"),
  Ependymal = c("microtubule cytoskeleton organization", "Golgi vesicle transport", "cilium-dependent cell motility"),
  Pericytes = c("circulatory system development", "blood vessel morphogenesis", "response to wounding"),
  Endothelial = c("vasculature development", "regulation of cell migration", "cell adhesion")
)

degs_n <- table(degs_sig$cluster)

# read enrich results
lapply(levels(scRNA$CellType), function(i){
  tmp <- read.table(paste0("data/for_enrich/output_enrich_", i, ".txt"), sep = "\t", header = T, check.names = F)
  
  cbind(tmp[match(go_used[[i]], tmp$Name),
            c("Name", "p-value", "Hit Count in Query List", "Hit Count in Genome")],
        cluster = i, n_deg = degs_n[[i]])
}) %>% do.call(rbind, .) -> dotData
head(dotData)


# set xmax
(xmax <- round(max(-log10(dotData$`p-value`))))

# dotplot
dotplot <- ggplot(cbind(dotData, Order = nrow(dotData):1)) +
  geom_point(mapping = aes(x = -log10(`p-value`),
                           y = Order,
                           size = `Hit Count in Query List`,
                           fill = `Hit Count in Query List`/n_deg),
             shape = 21) +
  scale_fill_gradientn(colours = c("grey", "gold", "red")) +
  scale_y_continuous(position = "right",
                     breaks = 1:nrow(dotData),
                     labels = Hmisc::capitalize(rev(dotData$Name))) +
  scale_x_continuous(breaks = seq(0, xmax, 10),
                     expand = expansion(mult = c(.15, .1))) +
  
  labs(x = "-log10(P-value)", y = NULL) +
  guides(size = guide_legend(title = "Gene count"),
         fill = guide_colorbar(title = "GeneRatio")) +
  theme_bw() +
  theme(panel.grid =element_blank())

dotplot
ggsave(filename = "./Figure/enrichment.pdf",plot = dotplot,width = 5,height = 6)


# relative proportion
# Average cell number and relative proportion
library(tidyverse)
library(reshape2)
library(patchwork)
library(reshape2)
library(RColorBrewer)
library(scales)
library(ggpubr)
library(ggplotify)

data_plotC <- table(scRNA@meta.data$Sample, scRNA@meta.data$CellType) %>% melt()
table(scRNA@meta.data$orig.ident)
colnames(data_plotC) <- c("Sample", "CellType","Number")

pC2 <- ggplot(data = data_plotC, aes(x = Sample, y = Number, fill = CellType)) +
  geom_bar(stat = "identity", width=0.8,aes(group=CellType),position="fill")+
  scale_fill_manual(values=mycols) +
  theme_bw()+
  theme(panel.grid =element_blank()) +
  labs(x="",y="Cell proportion")+
  theme(axis.text = element_text(size=12, colour = "black"))+
  theme(axis.title.y = element_text(size=12, colour = "black"))+
  theme(panel.border = element_rect(size = 1, linetype = "solid", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 0.8,  vjust = 0.6))

ggsave("./Figure/relative_proportion.pdf", plot = pC2, width = 8, height = 8)


# Data preparation for comparing control and disease groups
data_df <- scRNA@meta.data %>%
  count(CellType, group) %>%
  group_by(group) %>%
  mutate(percent = n / sum(n) * 100)

ggplot(data_df, aes(x = CellType, y = percent, fill = group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8, color = 'black') +
  scale_fill_manual(values = c("CO" = "#3979F2", "PD" = "#DF5233")) +
  theme_bw(base_size = 18) +
  labs(x = NULL, y = "Percent") +
  theme(axis.text.x = element_text(angle = 30, hjust = 0.4, vjust = 0.5))

# Save the cell type and group counts to a file
write.table(table(CellType = scRNA@meta.data$CellType, Group = scRNA@meta.data$group),
            file = "data/celltype_group_counts.csv", sep = "\t", quote = FALSE, col.names = NA)

#### 06. AUCell ####

library(Seurat)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(dplyr)
library(magrittr)
library(AUCell)
library(patchwork)
library(ggplot2)
options(stringsAsFactors = FALSE)

# load data
load("data/scRNA_anno.Rda")
tmp1 <- data.table::fread("data/20220813232518_GeneSearchResults_input.txt") #
tmp <- unique(tmp1$`Gene Name`)
gene_selec <- intersect(tmp,rownames(scRNA@assays$integrated@data)) 
geneSets <- list(ISGs=gene_selec)
scRNA
write.table(gene_selec,file = "ISGs_final.txt",quote = F,col.names = F,row.names = F)

# AUCell_buildRankings
set.seed(1234)
cells_rankings <- AUCell_buildRankings(scRNA@assays$integrated@data, nCores=20)

#AUCell_calcAUC
set.seed(1234)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings,nCores = 20,
                            aucMaxRank = nrow(cells_rankings)*0.2)
getAUC(cells_AUC)[,1:5]

# print
table(scRNA$CellType)
cells_assignment <- AUCell_exploreThresholds(cells_AUC,plotHist=T, assign=TRUE, nCores=20)
cells_assignment$ISGs$aucThr$thresholds

geneSetName <- rownames(cells_AUC)[grep("ISGs", rownames(cells_AUC))]
pdf("Figure/AUCell_plotHist.pdf")
AUCell_plotHist(cells_AUC[geneSetName,], aucThr=0.273)
abline(v=0.273)
dev.off()

geneSet.name <- "ISGs"
AUC_Exp <- as.numeric(getAUC(cells_AUC)[geneSet.name, ])
scRNA$AUC <- AUC_Exp
save(AUC_Exp,file = "data/AUC_score.Rdata")

#load("data/AUC_score.Rdata")
scRNA$AUC <- AUC_Exp
scRNA@meta.data <- scRNA@meta.data %>%
  mutate(ISGscore=if_else(AUC > 0.273,"High_ISGs","Low_ISGs"))
table(scRNA$ISGscore)

DimPlot(scRNA,
        group.by = "ISGscore",
        label = F,pt.size = 0.1,
        label.size = 5,
        cols = c("#CD3242","lightgray"),
        reduction = "tsne")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x =element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =element_blank())

DimPlot(scRNA,
        group.by = "ISGscore",
        split.by = "group" ,
        label = F,pt.size = 0.1,
        label.size = 5,
        cols = c("#CD3242","lightgray"),
        reduction = "tsne")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x =element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y =element_blank())


# data reshape
library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)

# Assuming scRNA, selectedThresholds, and cell_markers are already loaded
selectedThresholds <- 0.273

# Data Preparation for Visualization
plot.df <- scRNA@meta.data %>% 
  mutate(AboveThreshold = if_else(AUC > selectedThresholds, "High_ISGs", "Low_ISGs"))

# High Score Counts
high_counts <- table(plot.df$CellType[plot.df$AUC > selectedThresholds], plot.df$group[plot.df$AUC > selectedThresholds])

# Total Counts
total_counts <- table(plot.df$CellType, plot.df$group)

# Merge and Calculate Percentages
df <- merge(as.data.frame(high_counts), as.data.frame(total_counts), by = c("Var1", "Var2"))
names(df) <- c("celltype", "group", "high_score_num", "total_num")
df <- df %>%
  mutate(
    low_score_num = total_num - high_score_num,
    high_per = round((high_score_num / total_num) * 100, 2),
    low_per = round((low_score_num / total_num) * 100, 2)
  ) %>%
  select(celltype, group, high_per, low_per) %>%
  gather(key = "Score", value = "Percentage", -celltype, -group)

# Visualization
ggplot(df, aes(x = celltype, y = Percentage, fill = Score)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  scale_fill_manual(values = c(high_per = "#D74D35", low_per = "lightgray")) +
  facet_wrap(~ group) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

# Differential Expression Analysis
Idents(scRNA) <- scRNA$ISGscore
ISG_markers <- FindAllMarkers(scRNA, min.pct = 0.25, logfc.threshold = 0.25)
ISG_degs_sig <- ISG_markers %>%
  filter(pct.1 > 0.5, pct.2 < 0.5, p_val_adj <= 0.01, abs(avg_log2FC) > 1) %>%
  arrange(desc(avg_log2FC))

# Save Results
write.csv(ISG_degs_sig, file = "data/ISG_markers.csv")
write.csv(ISG_degs_sig[ISG_degs_sig$cluster == "High_ISGs", ]$gene, file = "data/High_ISG_markers.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


######################################################################################
## Pseudotime and TLR analysis ##
######################################################################################


# 07. pseudotime----
library(future)
library(ggplot2)
library(tidyverse)
library(monocle)

# Prepare for parallel processing
plan(multisession, workers = 20) # Adjust based on your system's capability

# Load the preprocessed scRNAseq dataset
load("data/scRNA_anno.Rda") # Adjust path if needed

# Ensure the output directory exists
if (!dir.exists("pseudotime")) {
  dir.create("pseudotime")
}

# Subset the scRNA-seq data for a specific cell type
scRNAsub <- subset(scRNA, subset = CellType == 'Microglia') # Example for 'Microglia'

# Convert data to Monocle's required format
data <- as(as.matrix(scRNAsub@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNAsub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

# Create a CellDataSet for Monocle
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())

# Pre-process data: size factor estimation, dispersion estimation, and gene detection
set.seed(1234)
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, relative_expr = TRUE)
mycds <- detectGenes(mycds, min_expr = 2)

# Variable features selection for ordering cells
scRNAsub <- FindVariableFeatures(scRNAsub, selection.method = "vst", nfeatures = 3000)
var.genes <- VariableFeatures(scRNAsub)
mycds <- setOrderingFilter(mycds, var.genes)

# Reduce dimension and order cells
mycds <- reduceDimension(mycds, max_components = 2, num_dim = 50, reduction_method = 'DDRTree')
mycds <- orderCells(mycds)

# Visualization of cell trajectories
plot1 <- plot_cell_trajectory(mycds, color_by = "State", cell_size = 0.2)
plot2 <- plot_cell_trajectory(mycds, color_by = "group", show_tree = TRUE, cell_size = 0.2) + 
  scale_color_manual(values = c("#3979F2", "#DF5233"))
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime", cell_size = 0.2)

# Combine plots using patchwork
combined_plots <- plot1 | plot2 | plot3
combined_plots

# Save combined plots and pseudotime data
ggsave("pseudotime/combined_plots.pdf", plot = combined_plots, width = 16, height = 6, units = "in")
write.csv(pData(mycds), "pseudotime/pseudotime.csv")

# BEAM analysis for branch point 1
beam_res <- BEAM(mycds, branch_point = 1)
beam_res <- beam_res[order(beam_res$qval),]
write.csv(beam_res, "pseudotime/beam_results.csv")

# Visualization of top genes along pseudotime
top_genes <- head(rownames(beam_res), 100)
heatmap_plot <- plot_pseudotime_heatmap(mycds[top_genes, ], num_clusters = 4)
ggsave("pseudotime/top_genes_heatmap.pdf", plot = heatmap_plot, width = 8, height = 12)

# enrichment

gene_group=tmp1$annotation_row
gene_group$gene=rownames(gene_group)

library(clusterProfiler)
library(org.Hs.eg.db)
allcluster_go=data.frame()
for (i in unique(gene_group$Cluster)) {
  small_gene_group=filter(gene_group,gene_group$Cluster==i)
  df_name=bitr(small_gene_group$gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  go <- enrichGO(gene         = unique(df_name$ENTREZID),
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.2,
                 readable      = TRUE)
  go_res=go@result
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    allcluster_go=rbind(allcluster_go,go_res)
  }
}
head(allcluster_go[,c("ID","Description","qvalue","cluster")])
go_term2 <- subset(allcluster_go,cluster==2)
go_term1 <- subset(allcluster_go,cluster==1)

write.table(allcluster_go,file = "data/monocle_go_terms_NEW.txt",quote = F,sep = "\t",
            row.names = F,col.names = T)

#Endothelial_pseudotime
load("pseudotime/Endothelial_pseudotime.Rdata")

# disp.genes
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
p3 <- plot_ordering_genes(mycds)
p3

# reduceDimension
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
# orderCells
mycds <- orderCells(mycds)

# color
group_color <- c("#3979F2","#DF5233")
mycols <- c("#A6CEE3","#B2DF8A","#1F78B4","#FB9A99","#E31A1C",
            "#33A02C","#FDBF6F","#FF7F00" , "#CAB2D6" ,"#6A3D9A")
show_col(mycols)

# plot_cell_trajectory (state)
plot1 <- plot_cell_trajectory(mycds,color_by = "State",cell_size=0.2)+
  scale_color_manual(values = mycols)
# plot_cell_trajectory (group)
plot2 <- plot_cell_trajectory(mycds, color_by = "group",
                              show_tree = T,
                              cell_size=0.2,
                              cell_link_size = 0.75,
                              cell_name_size = 8,
                              state_number_size = 8,
                              show_branch_points = F,
                              theta = 0)+ scale_color_manual(values = group_color)
# plot_cell_trajectory(Pseudotime)
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime",cell_size=0.2)
# combine
plotc <- plot1|plot2|plot3
plotc
ggsave("Figure/Endothelial_pseudotime.pdf", plot = plotc, width = 10, height = 3.5)


# Pericytes_pseudotime
load("pseudotime/Pericytes_pseudotime.Rdata")

# disp.genes
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
p3 <- plot_ordering_genes(mycds)
p3

# reduceDimension
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
# orderCells
mycds <- orderCells(mycds)

# color
group_color <- c("#3979F2","#DF5233")
mycols <- c("#A6CEE3","#B2DF8A","#1F78B4","#FB9A99","#E31A1C",
            "#33A02C","#FDBF6F","#FF7F00" , "#CAB2D6" ,"#6A3D9A")
show_col(mycols)

# plot_cell_trajectory (state)
plot1 <- plot_cell_trajectory(mycds,color_by = "State",cell_size=0.2)+
  scale_color_manual(values = mycols)
# plot_cell_trajectory (group)
plot2 <- plot_cell_trajectory(mycds, color_by = "group",
                              show_tree = T,
                              cell_size=0.2,
                              cell_link_size = 0.75,
                              cell_name_size = 8,
                              state_number_size = 8,
                              show_branch_points = F,
                              theta = 0)+ scale_color_manual(values = group_color)
# plot_cell_trajectory(Pseudotime)
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime",cell_size=0.2)
# combine
plotc <- plot1|plot2|plot3
plotc
ggsave("Figure/Pericytes_pseudotime.pdf", plot = plotc, width = 10, height = 3.5)


# 08.TLR analysis----

library(Seurat)
library(magrittr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)

markers <- c("TLR1","TLR2","TLR3","TLR5","TLR7","TLR8","TLR10")
markers %in% rownames(scRNA@assays$RNA@counts)
# FetchData
df <- FetchData(scRNA, vars = c(markers, "group", "CellType")) %>%
  reshape2::melt(id.vars = c("group", "CellType"),
                 variable.name = "gene",
                 value.name = "expression") %>%
  dplyr::mutate(group = factor(group, levels = c("CO","PD")),
                CellType = factor(CellType))

#stat.test
stat.test <- df %>%
  group_by(CellType, gene) %>%
  t_test(expression ~ group) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance()

#color
group_color <- c("#3979F2","#DF5233")
stat.test <- stat.test %>% add_xy_position(x = "group")
class(stat.test)

#plot
v1 <- ggviolin(df, x = "group", y = "expression",
               fill = "group",
               color = "group",
               facet = c("gene","CellType")) +
  scale_fill_manual(values = group_color) +
  scale_color_manual(values = group_color) +
  xlab("") +
  ylab("Expression level (count)") +
  theme(panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0),
        strip.background.x = element_rect(fill = "grey95"),
        strip.background.y = element_blank())
v1

v1 + stat_pvalue_manual(stat.test, hide.ns = TRUE) #
ggsave("./Figure/TLR.pdf", width = 10, height = 10)



######################################################################################
## SCENIC,SEEK and other analysis ##
######################################################################################

# 09.SCENIC ----
library(plyr)
library(permute)
library(data.table,lib.loc = "/usr/local/lib/R/library")
library(SCopeLoomR)
library(SCENIC)
library(arrow)
library(future)

dge <- as.matrix(scRNA@assays$RNA@counts)
dim(dge)

# cellMeta
t# Prepare cell metadata
cell.info <- data.frame(ClusterID = as.integer(scRNA@meta.data$integrated_snn_res.2),
                        Tissue = scRNA@meta.data$Sample,
                        CellType = scRNA@meta.data$CellType,
                        row.names = colnames(scRNA))
cell.info$CellState <- paste(cell.info$Tissue, cell.info$ClusterID, sep = "_")

# Initialize SCENIC options
scenicOptions <- initializeScenic(org="hgnc", dbDir="path/to/cisTarget_databases",
                                  dbs="hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
                                  datasetTitle="SCENIC on PD Brain Cells", nCores=10)

# Gene filtering
genesKept <- geneFiltering(dge, scenicOptions, minCountsPerGene = 1, minSamples = 20)
dgeFiltered <- dge[genesKept, ]

# AvgN 20
avg20.rep1 <- AvgN(dge, cell.info, seed=1)

if(!dir.exists("output")) {
  dir.create("output")
}
saveRDS(cell.info, "output/s1_cell.info.rds")

# saveLoom
saveLoom <- function(exprMat, output){
  ## prepare cell cluster info
  cellInfo <- data.frame(
    row.names = colnames(exprMat),
    cellType = sub("\\.[0-9]+", "", colnames(exprMat))
  )
  ##  save loom
  loom <- build_loom(output, dgem=exprMat)
  loom <- add_cellAnnotation(loom, cellInfo)
  close_loom(loom)
}

saveLoom(avg20.rep1, "output/s1_avg20_rep1.loom")

# add_cellAnnotation
loom <- build_loom("output/s1_exprMat.loom", dgem=dge)
loom <- add_cellAnnotation(loom, cell.info)
close_loom(loom)

#### Linux python work
# s2_runPySCENIC.sh
# bash s2_runPySCENIC.sh avg20_rep1
# nohup bash s2_runPySCENIC.sh avg20_rep1 1>s2_rep1.2022.log 2>&1 &
#
#   # s3_postSCENIC.py
#   f_loom_path_scenic=output/s2_avg20_rep1.pyscenic.loom
# ctx_output=output/s2_avg20_rep1.reg.tsv
# sample_name=output/s3_avg20_rep1
# threads=10
# min_regulon_size=5
# python s3_postSCENIC.py $f_loom_path_scenic $ctx_output $sample_name $threads $min_regulon_size

# load package
library(ggplot2)
library(plyr)
library(ggsci)
library(data.table)
library(pbapply)
library(plyr)
library(philentropy)
library(ggplot2)
library(ggrepel)
library(latex2exp)

# load data
cell.info <- readRDS("output/s1_cell.info.rds")
table(scRNA$integrated_snn_res.2)
df<- data.frame(scRNA@reductions$tsne@cell.embeddings,
                scRNA@reductions$umap@cell.embeddings)
identical(rownames(df),rownames(cell.info))
cell.info <- cbind(cell.info, df[rownames(cell.info), ])

head(cell.info)
saveRDS(cell.info, "output/s4_cell.info.rds")
tmp1 <- cell.info[,c(1,3)]
colnames(tmp1) <- c("Cluster","Cell Type")
write.table(tmp1,file = "../data/CellType.Info.txt",row.names = F,col.names = T,quote = F,sep = "\t")
tmp2 <- cell.info[,c(1,2)]
tmp3  <- tibble::rownames_to_column(tmp2)
colnames(tmp3)[1] <- c("")
write.table(tmp3,file = "../data/Cell.Info.txt",row.names = F,col.names = T,quote = F,sep = "\t")

# rasMat
rasMat <- fread("output/s3_avg20_rep1.AUCell.txt", sep = "\t", header = T, data.table = F) # data.frame
rownames(rasMat) <- rasMat$V1
colnames(rasMat) <- sub("(+)", "", colnames(rasMat), fixed = T)
rasMat <- rasMat[, -1]
saveRDS(rasMat, "output/s5_avg20_rep1.rasMat.rds")

cell.info <- readRDS("output/s4_cell.info.rds")
cell.types <- names(table(cell.info$CellType))
ctMat <- lapply(cell.types, function(i) {
  as.numeric(cell.info$CellType == i)
})
ctMat <- do.call(cbind, ctMat)
colnames(ctMat) <- cell.types
rownames(ctMat) <- rownames(cell.info)

# rssMat
rssMat <- pblapply(colnames(rasMat), function(i) {
  sapply(colnames(ctMat), function(j) {
    1 - JSD(rbind(rasMat[, i], ctMat[, j]), unit = 'log2', est.prob = "empirical")
  })
})
rssMat <- do.call(rbind, rssMat)
rownames(rssMat) <- colnames(rasMat)
colnames(rssMat) <- colnames(ctMat)
saveRDS(rssMat, "output/s5_avg20_rep1.rssMat.rds")

# binMat
rssMat <- readRDS("output/s5_avg20_rep1.rssMat.rds")
binMat <- read.table("output/s3_avg20_rep1.binary_mtx.txt", sep = "\t", header = T, row.names = 1, check.names = FALSE)
colnames(binMat) <- sub("(+)", "", colnames(binMat), fixed = T)

source("utils/plotRegulonRank.R")
colnames(rssMat)

# RegulonRank2
p.list <- lapply(colnames(rssMat)[1:10], function(x)
  PlotRegulonRank(rssMat, cell.type = x, topn = 5))
pdf("../Figure/RegulonRank2.pdf",width = 10,height = 10)
cowplot::plot_grid(plotlist = p.list, ncol = 5)
dev.off()

cell.info <- readRDS("output/s4_cell.info.rds")
cell.info <- cbind(cell.info, binMat[rownames(cell.info), ])

# UMAP for Regulons
source("utils/DimPlot.R")
DimPlot(cell.info, cell.type = "Microglia")
DimPlot(cell.info, regulon = "NFATC2")
DimPlot(cell.info, regulon = "RUNX2")
DimPlot(cell.info, regulon = "IRF5")

fig2Plot <- function(cell.type, regulon) {
  p.list <- list(
    PlotRegulonRank(rssMat, cell.type),
    DimPlot(cell.info, cell.type = cell.type),
    DimPlot(cell.info, regulon = regulon)
  )
  cowplot::plot_grid(plotlist = p.list, ncol = 3,
                     rel_widths = c(3,5,5))
}

p1 <- fig2Plot("Microglia", "NFATC2")
p2 <- fig2Plot("Microglia", "RUNX2")
p3 <- fig2Plot("Microglia", "IRF5")
cowplot::plot_grid(p1,p2,p3, ncol = 1)

# SEEK database
NFATC2.seek <- read.csv("../seek_results/NFATC2_SEEK.txt", sep = "\t", header = T)
RUNX2.seek <- read.csv("../seek_results/RUNX2_SEEK.txt", sep = "\t", header = T)
IRF5.seek <- read.csv("../seek_results/IRF5_SEEK.txt", sep = "\t", header = T)

SeekPlot <- function(seek.res, key.words){
  seek.res$Related <- grepl(key.words, seek.res$Description, ignore.case = T)
  
  data.use <- data.frame(
    Dataset = seek.res$Rank,
    P = -log10(seek.res$Coexpression.PValue),
    Related = seek.res$Related & seek.res$Coexpression.PValue < 0.01
  )
  m = sum(seek.res$Related)
  n = sum(data.use$Related)
  M = nrow(seek.res)
  N = sum(seek.res$Coexpression.PValue < 0.01)
  fisher.res <- fisher.test(matrix(c(n,m-n,N-n,M-N-m+n), ncol = 2), alternative = "greater")
  max.p <- max(data.use$P[is.finite(data.use$P)])
  data.use$P <- ifelse(is.finite(data.use$P), data.use$P, max.p)
  
  ggplot(data.use, aes(Dataset, P)) +
    geom_point(color="#4590CE", size=3) +
    geom_point(inherit.aes = F, data = subset(data.use, Related),
               aes(Dataset, P), color="#E2AE2D", size=3) +
    geom_hline(yintercept=2, color="grey", size=1) +
    annotate("text", x=Inf, y=Inf, hjust=1.3, vjust=8.5, label=paste0("(", n, " out of ", m, ")")) +
    annotate("text", x=Inf, y=Inf, hjust=1.4, vjust=10, label=paste0("p=", signif(fisher.res$p.value, 3))) +
    ylab(TeX("-log_{10}(p-value)")) +
    ggtitle("") +
    theme_bw(base_size = 15) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black"),
          axis.ticks = element_line(color="black"),
          axis.text = element_text(color = "black"),
          plot.title = element_text(hjust = .5)
    )
}

p1.1 <- SeekPlot(NFATC2.seek, "(CNS)|(microglia)|(brain)")
p2.1 <- SeekPlot(RUNX2.seek, "(CNS)|(microglia)|(brain)")
p3.1 <- SeekPlot(IRF5.seek, "(CNS)|(microglia)|(brain)")

p1.2 <- cowplot::plot_grid(p1, p1.1, rel_widths = c(13,3))
p2.2 <- cowplot::plot_grid(p2, p2.1, rel_widths = c(13,3))
p3.2 <- cowplot::plot_grid(p3, p3.1, rel_widths = c(13,3))
cowplot::plot_grid(p1.2, p2.2, p3.2, ncol = 1)

if(!dir.exists("../Figure")) {
  dir.create("../Figure")
}
pdf("../Figure/SEEK.pdf", width = 16, height = 20)
cowplot::plot_grid(p1.2, p2.2,p3.2, ncol = 1)
dev.off()

# Corheatmap
library(grid)
library(pbapply)
library(circlize)
library(ggsci)
library(plyr)
library(ggplot2)
library(ggrepel)
library(dendextend)
library(ComplexHeatmap)

#load data
rasMat <- readRDS("output/s5_avg20_rep1.rasMat.rds");class(rasMat)
regulons <- read.table("output/s3_avg20_rep1.regulons.txt", sep = "\t")
regulon.names <- sapply(strsplit(regulons$V1, split = "\\("), function(x) x[1])
regulon.sizes <- sapply(strsplit(regulons$V3, split = ","), function(x) length(x))
regulon.names <- regulon.names[regulon.sizes>=10]
write.table(regulon.names,"output/s6_regulon_name.txt", sep = "\t",row.names = F,col.names = F,quote = F)
save(regulon.names,file = "regulon.names.Rdata")

# rasMat
rasMat <- rasMat[, regulon.names];class(rasMat)
dim(rasMat)
save(rasMat,file = "rasMat_292.Rdata")
pccMat <- cor(rasMat);class(pccMat)

# CSI
CSI <- function(r1, r2) {
  delta <- pccMat[r1,r2]
  r.others <- setdiff(colnames(pccMat), c(r1,r2))
  N <- sum(pccMat[r1, r.others] < delta) + sum(pccMat[r2, r.others] < delta)
  M <- length(r.others) * 2
  return(N/M)
}
csiMat <- pblapply(rownames(pccMat), function(i) sapply(colnames(pccMat), function(j) CSI(i, j)))
csiMat <- do.call(rbind, csiMat)
rownames(csiMat) <- rownames(pccMat)

round(pccMat[1:10,1:10], 2)
round(csiMat[1:10,1:10], 2)
csiMat.binary <- matrix(as.numeric(csiMat >= 0.7), nrow = nrow(csiMat))
colnames(csiMat.binary) <- colnames(csiMat)
rownames(csiMat.binary) <- rownames(csiMat)
csiMat.binary[1:10,1:10]

saveRDS(csiMat, "output/s6_avg20_rep1.csiMat.rds")
write.table(csiMat.binary, "output/s6_avg20_rep1.csiMat.binary.txt", sep = "\t")

mat = readRDS("output/s6_avg20_rep1.csiMat.rds");class(mat)
h = 6
row_dend = as.dendrogram(hclust(dist(mat), method = "complete"))
clusters <- cutree(row_dend, h = h) # dendextend::cutree()
row_dend = color_branches(row_dend, h = h, col = pal_d3("category20")(20))
plot(row_dend)

col_range = c(0.7, 1)
col_fun <- colorRamp2(col_range, c("#FCF8DE", "#253177"))

ht <- Heatmap(
  matrix = mat,
  col = col_fun,
  name = "ht1",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  show_heatmap_legend = FALSE
)

lgd <- Legend(
  col_fun = col_fun,
  title = "",
  at = col_range,
  labels = c("low", "high"),
  direction = "horizontal",
  legend_width = unit(1, "in"),
  border = FALSE
)

#plot heatmap
pdf("../Figure/corheatmap.pdf", width = 10, height = 10)
draw(ht, heatmap_legend_list = list(lgd), heatmap_legend_side = c("bottom"))
decorate_heatmap_body("ht1", {
  tree = column_dend(ht)
  ind = cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
  first_index = function(l) which(l)[1]
  last_index = function(l) { x = which(l); x[length(x)] }
  clusters <- names(table(ind))
  x1 = sapply(clusters, function(x) first_index(ind == x)) - 1
  x2 = sapply(clusters, function(x) last_index(ind == x))
  x1 = x1/length(ind)
  x2 = x2/length(ind)
  grid.rect(x = x1, width = (x2 - x1), y = 1-x1, height = (x1 - x2),
            hjust = 0, vjust = 0, default.units = "npc",
            gp = gpar(fill=NA, col="#FCB800", lwd=3))
  grid.text(label = paste0("M",clusters),
            x = x2-length(clusters)/length(ind), y = 1-x1-(x2-x1)/2,
            default.units = "npc",
            hjust = 1, vjust = 0.5,
            gp = gpar(fontsize=12, fontface="bold"))
})
decorate_column_dend("ht1", {
  tree = column_dend(ht)
  ind = cutree(as.hclust(tree), h = h)[order.dendrogram(tree)]
  first_index = function(l) which(l)[1]
  last_index = function(l) { x = which(l); x[length(x)] }
  clusters <- names(table(ind))
  x1 = sapply(clusters, function(x) first_index(ind == x)) - 1
  x2 = sapply(clusters, function(x) last_index(ind == x))
  grid.rect(x = x1/length(ind), width = (x2 - x1)/length(ind), just = "left",
            default.units = "npc", gp = gpar(fill = pal_d3("category20")(20), alpha=.5, col = NA))
})
dev.off()

# Heatmap
library(grid)
library(pbapply)
library(circlize)
library(ggsci)
library(plyr)
library(ggplot2)
library(ggrepel)
library(dendextend)
library(ComplexHeatmap)
library(data.table)

#load data
table(scRNA$ISGscore)
rasMat <- fread("output/s3_avg20_rep1.AUCell.txt", sep = "\t", header = T, data.table = F) # data.frame
rownames(rasMat) <- rasMat$V1
colnames(rasMat) <- sub("(+)", "", colnames(rasMat), fixed = T)
rasMat <- rasMat[, -1]
rasMat[1:5,1:5]

# correlation plot
regulons <- read.table("output/s3_avg20_rep1.regulons.txt", sep = "\t")
regulon.names <- sapply(strsplit(regulons$V1, split = "\\("), function(x) x[1])
regulon.sizes <- sapply(strsplit(regulons$V3, split = ","), function(x) length(x))
regulon.names <- regulon.names[regulon.sizes>=10]
rasMat <- rasMat[, regulon.names]
dim(rasMat)

#regulon& correlation
identical(rownames(rasMat),rownames(scRNA@meta.data))
# regulon select
gene_cor <- c("RUNX2","NFATC2","IRF5")
rasMat_cor <- rasMat[,gene_cor]
rasMat_cor$AUCscore <- AUC_Exp
rasMat_cor$group <- scRNA@meta.data$group
rasMat_cor$celltype <- scRNA@meta.data$CellType
table(rasMat_cor$celltype)
rasMat_cor <- rasMat_cor%>%
  dplyr::filter(group=="PD")

# correlation
options (digits = 10)
cor_list <- list()
cor_data_df <- data.frame(colnames(rasMat_cor[1:length(gene_cor)]))

y=rasMat_cor$AUCscore
for (i in 1:length(gene_cor)){
  test <- cor.test(as.numeric(rasMat_cor[,i]),y, method="pearson")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("gene","cor","pvalue")

dir.create("../Figure/")
for (i in 1:length(gene_cor)) {
  cor_gene2 <- cor_data_df %>% dplyr::filter(gene==gene_cor[i])
  labeltext <- paste(paste0("n = 17900"),
                     paste0("r = ",round(cor_gene2[2],3),"(pearson)"),
                     if_else(cor_gene2[3]<0.0001,"p.value<0.0001",
                             paste0("p.value= ",round(cor_gene2[3],3))),
                     sep = ", ")
  ggplot(rasMat_cor,aes(rasMat_cor[,i],AUCscore))+
    geom_point(col="#984ea3",size=0.1)+
    geom_smooth(method=lm, se=T,na.rm=T, fullrange=T,size=2,col="#fdc086")+
    geom_rug(col="#7fc97f")+
    theme_minimal()+
    xlab(as.name(gene_cor[i]))+
    theme(plot.title = element_text(hjust = 0.5),
          plot.margin = margin(1, 1, 1, 1, "cm"))+
    labs(title =labeltext)
  ggsave(paste0("../Figure/S4.",gene_cor[i],"_cor.pdf"))
}


## differentaly expressed regulons in High_ISGs vs Low_ISGs microglia
rasMat_t <- t(rasMat)
rasMat_t[1:5,1:5]
dim(rasMat_t)
scRNA_sub <- subset(scRNA,CellType=="Microglia")
clinical <- as.data.frame(colnames(scRNA_sub@assays$RNA@counts))
colnames(clinical) <- "CellID"
rownames(clinical) <- clinical$CellID
clinical$group <- scRNA_sub$ISGscore #
rasMat_t <- rasMat_t[,rownames(clinical)]
dim(rasMat_t)
clinical=clinical[match(colnames(rasMat_t),clinical$CellID),]
identical(colnames(rasMat_t),rownames(clinical))

group_list <- clinical$group
group_list = factor(group_list,
                    levels = c("Low_ISGs","High_ISGs"))
table(group_list)

#limma
library(limma)
library(dplyr)

design=model.matrix(~group_list)
fit=lmFit(rasMat_t,design)
fit=eBayes(fit)
options(digits = 4)
topTable(fit,coef=2,adjust='BH')
deg=topTable(fit,coef=2,adjust='BH',number = Inf)
table(group_list)
deg <- mutate(deg,symbol=rownames(deg))
head(deg)
deg <- deg[!duplicated(deg$symbol),]
logFC_t=0
P.Value_t = 0.05
k1 = (deg$P.Value < P.Value_t)&(deg$logFC < -logFC_t)
k2 = (deg$P.Value < P.Value_t)&(deg$logFC > logFC_t)
Change = ifelse(k1,"Down",ifelse(k2,"Up","Stable"))
table(Change)
deg$Change <- Change

#regulon cluster
regulons_cluster <- read.table("output/s6_avg20_rep1.regulon_clusters.txt", sep = "\t")
colnames(regulons_cluster) <- regulons_cluster[1,]
regulons_cluster <- regulons_cluster[-1,]
DEGs <- dplyr::left_join(deg,regulons_cluster,by=c("symbol"="regulon"))
deg_mic <- deg
DEGs_mic <- DEGs
save(deg_mic,DEGs_mic,file = "../data/deg_mic_regulon_between_ISGsgroup.Rdata")
write.table(DEGs_mic,file = "../data/deg_mic_regulons_betewwn_ISGsgroup.txt",
            col.names = T,quote = F,row.names = T,sep = "\t")

upregulon_microglia_HighISG <- DEGs_mic[(DEGs_mic$Change =="Up") &
                                          (DEGs_mic$cluster =="M2"|DEGs_mic$cluster =="M3"),]$symbol
upregulon_microglia_HighISG

## differentaly expressed regulons in CO vs PD microglia
clinical$group <- scRNA_sub$group
rasMat_t <- rasMat_t %>% as.data.frame()
rasMat_t <- rasMat_t[,rownames(clinical)]
clinical=clinical[match(colnames(rasMat_t),clinical$CellID),]
identical(colnames(rasMat_t),rownames(clinical))

group_list <- clinical$group
group_list = factor(group_list,
                    levels = c("CO","PD"))
table(group_list)
table(scRNA$group)

#limma
design=model.matrix(~group_list)
fit=lmFit(rasMat_t,design)
fit=eBayes(fit)
options(digits = 4)
topTable(fit,coef=2,adjust='BH')
deg=topTable(fit,coef=2,adjust='BH',number = Inf)
table(group_list)
deg <- mutate(deg,symbol=rownames(deg))
head(deg)
deg <- deg[!duplicated(deg$symbol),]
logFC_t=0.0
P.Value_t = 0.05
k1 = (deg$P.Value < P.Value_t)&(deg$logFC < -logFC_t)
k2 = (deg$P.Value < P.Value_t)&(deg$logFC > logFC_t)
Change = ifelse(k1,"Down",ifelse(k2,"Up","Stable"))
table(Change)
deg$Change <- Change
table(deg$Change)
DEGs <- dplyr::left_join(deg,regulons_cluster,by=c("symbol"="regulon"))

upregulon_microglia_PDCO <- DEGs[(DEGs$Change =="Up") &
                                   (DEGs$cluster =="M2"|DEGs$cluster =="M3"),]$symbol
m2.3 <- DEGs[(DEGs$cluster =="M2"|DEGs$cluster =="M3"),]$symbol

# VennDiagram
library(VennDiagram)
library(tidyverse)
v1 <- "M2/3"
v2 <- "microglia_HighISG"
v3 <- "microglia_PDCO"
x = list(v1=m2.3,
         v2=upregulon_microglia_HighISG,
         v3=upregulon_microglia_PDCO)
names(x) <- c(v1,v2,v3)
mycols <- c("#a6cee3","#b2df8a","#fb9a99") # "#ffffb3",#b2df8a,#fb9a99
venn.plot <-venn.diagram(x,
                         filename = NULL,
                         lty ="dotted",
                         lwd =0.5,
                         col ="black",
                         fill =mycols,
                         alpha =0.60,
                         cat.col =mycols,
                         cat.cex =1,
                         cat.fontface="bold",
                         margin =0.05,
                         cex =1)
pdf("../Figure/S.veen.pdf")
grid.draw(venn.plot)
dev.off()
comregulon <- intersect(upregulon_microglia_HighISG,upregulon_microglia_PDCO)
write.table(comgene,file = "../data/comregulon_HighISGs_PD.txt",quote = F,sep = "\t",
            col.names = T,row.names = F)

# TF gene expression
# TF
markers <- c("RUNX2","NFATC2","IRF5")
markers %in% rownames(scRNA_sub@assays$RNA@counts)

# FetchData
df <- FetchData(scRNA_sub, vars = c(markers, "ISGscore", "CellType")) %>%
  reshape2::melt(id.vars = c("ISGscore", "CellType"),
                 variable.name = "gene",
                 value.name = "expression") %>%
  mutate(group = factor(ISGscore, levels = c("Low_ISGs","High_ISGs")),
         CellType = factor(CellType)) %>%
  dplyr::select(-ISGscore)

ggplot(df, mapping = aes(group, expression, fill = group, color = group)) +
  scale_fill_manual(values = c("#2874C5","#EABF00")) +
  scale_color_manual(values = c("#2874C5","#EABF00")) +
  geom_boxplot(outlier.size = -1,
               show.legend = FALSE,
               color="black", lwd=0.2,
               alpha = 0.7) +
  geom_point(shape = 21, size=.6,
             show.legend = FALSE,
             position = position_jitterdodge(),
             alpha = 0.3)+
  
  facet_grid(gene ~ CellType) +
  theme_bw() +
  xlab("") +
  ylab("Expression level (count)") +
  theme(panel.grid = element_blank(),
        strip.background.x = element_rect(fill = "grey95"),
        strip.text.y = element_text(angle = 0),
        strip.background.y = element_blank())
ggsave("../Figure/TF expression.pdf", width = 6, height = 8)

stat.test <- df %>%
  dplyr::group_by(CellType, gene) %>%
  t_test(expression ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()


