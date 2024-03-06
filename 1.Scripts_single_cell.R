######################################################################################
               ## Data process, annotation, and AUCell ##
######################################################################################

#### 01. Data process ####
# GSE157783
library(Seurat)
library(future)
library(ggplot2)
library(tidyverse)
plan()
plan("multisession", workers = 24)
plan()
options(future.globals.maxSize = 1000 * 1024^24)
counts <- data.table::fread("./data/UMI.tsv")
scRNA <- CreateSeuratObject(counts = counts)

genes <- data.table::fread("./data/genes.tsv")
head(rownames(scRNA))
# RenameGenesSeurat
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = HGNC.updated[[i]]$Suggested.Symbol) {
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA

  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}
scRNA <- RenameGenesSeurat(obj = scRNA,
                           newnames = genes$gene)
head(rownames(scRNA))

# ID annotation
load("~/ref/annotation.v40.Rdata")
length(unique(mRNA_anno$gene_name)) 
mito_genes=mRNA_anno$gene_name[grep("^MT-", mRNA_anno$gene_name)] ;mito_genes
mito_genes=lnc_anno$gene_name[grep("^MT-", lnc_anno$gene_name)] ;mito_genes
ribo_genes=mRNA_anno$gene_name[grep("^Rp[sl]", mRNA_anno$gene_name,ignore.case = T)];ribo_genes
ribo_genes=lnc_anno$gene_name[grep("^Rp[sl]", lnc_anno$gene_name,ignore.case = T)];ribo_genes

ids <- rbind(mRNA_anno,lnc_anno) %>% dplyr::select(-3)
colnames(ids) <- c('SYMBOL','ENSEMBL')
ids$ENSEMBL <- stringr::str_sub(ids$ENSEMBL,1,15)
table(rownames(scRNA) %in% ids$ENSEMBL)
length(unique(ids$SYMBOL))
length(unique(ids$ENSEMBL))

ids=ids[!duplicated(ids$SYMBOL),]
ids=ids[!duplicated(ids$ENSEMBL),]
ids=ids %>% dplyr::filter(ENSEMBL %in% rownames(scRNA))
pos=match(ids$ENSEMBL,rownames(scRNA))
scRNA=scRNA[pos,]
scRNA <- RenameGenesSeurat(obj = scRNA,
                           newnames = ids$SYMBOL)

# metadata
head(stringr::str_split(colnames(scRNA@assays$RNA@data), "_",simplify = T))
scRNA$orig.ident = stringr::str_split(colnames(scRNA@assays$RNA@data), "_",simplify = T)[,2]
table(scRNA$orig.ident)

barcodes <- data.table::fread("./data/barcodes.tsv")
scRNA$Sample = barcodes$patient
table(scRNA$Sample)
sum(table(scRNA$Sample)[1:6]) 
sum(table(scRNA$Sample)[7:11]) 

# group
scRNA$group <- ifelse(stringr::str_detect(scRNA$Sample,"C"),"CO","PD")
table(scRNA$group)
scRNA

# rownames
head(rownames(scRNA[["RNA"]]@meta.features))
scRNA[["RNA"]]@meta.features <- data.frame(row.names = rownames(scRNA[["RNA"]]))
head(rownames(scRNA[["RNA"]]@meta.features))

# save data
save(scRNA, file = "./data/scRNA_raw.Rdata")


#### 02. QC ####
rm(list = ls())

library(tidyverse)
library(Seurat)
library(vctrs)
library(viridis)
library(future)
plan()
plan("multisession", workers = 16)
plan()
options(future.globals.maxSize = 1000 * 1024^24)

load("data/scRNA_raw.Rdata")

nFeature_lower <- 200
nFeature_upper <- 6000
nCount_lower <- 200
nCount_upper <- 20000


scRNA <- PercentageFeatureSet(scRNA,
                              pattern = "^MT-",
                              col.name = "pMT")
scRNA <- PercentageFeatureSet(scRNA,
                              pattern = "^HBA|^HBB",
                              col.name = "pHB")
scRNA <- PercentageFeatureSet(scRNA,
                              pattern = "^RPS|^RPL",
                              col.name = "pRP")
head(scRNA)

# VlnPlot (supplementary Figure 1a)

group_color <- c("#3979F2","#DF5233")

VlnPlot(object = scRNA,
        features = c("nFeature_RNA", "nCount_RNA"),
        group.by  = "group",
        cols =group_color,
        log = T,
        pt.size = 0)
ggsave("Figure/boxplot_before_QC.pdf", height = 6, width = 10)

# filter
scRNA <- subset(scRNA,
                subset = nCount_RNA > nCount_lower &
                  nCount_RNA < nCount_upper  &
                  nFeature_RNA > nFeature_lower &
                  nFeature_RNA < nFeature_upper)

table(scRNA@meta.data$Sample)
scRNA

# VlnPlot
VlnPlot(object = scRNA,
        features = c("nFeature_RNA", "nCount_RNA"),
        group.by  = "group",
        cols =group_color,
        log = T,
        pt.size = 0)
ggsave("Figure/boxplot_after_QC.pdf", height = 6, width = 10)

# CellCycleScoring
scRNA <- CellCycleScoring(
  scRNA,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes)
head(scRNA)
save(scRNA, file = "data/scRNA.filter.Rda")

#### 03. Integration ####
# SCT
rm(list = ls())
load("data/scRNA.filter.Rda")
head(scRNA)
table(scRNA$Sample)
scRNA <- SplitObject(scRNA,
                     split.by = "Sample")
scRNA

plan()
plan(multisession, workers = 16)
plan()
options(future.globals.maxSize = 1000 * 1024^24)

# SCTransform
for(i in 1:length(scRNA)){
  scRNA[[i]] <- SCTransform(
    scRNA[[i]],
    variable.features.n = 3000,
    vars.to.regress = c("S.Score", "G2M.Score"),
    verbose = FALSE)
  print(i)
}

# SelectIntegrationFeatures
features <- SelectIntegrationFeatures(object.list = scRNA,
                                      nfeatures = 3000)
# PrepSCTIntegration
scRNA <- PrepSCTIntegration(object.list = scRNA,
                            anchor.features = features)

# FindIntegrationAnchors
AnchorSet <- FindIntegrationAnchors(object.list = scRNA,
                                    reference = 2,
                                    normalization.method = "SCT",
                                    anchor.features = features)
save(AnchorSet, file = "data/AnchorSet.Rda")

load("data/AnchorSet.Rda")

# IntegrateData
scRNA <- IntegrateData(anchorset = AnchorSet,
                       normalization.method = "SCT")
scRNA


save(scRNA, file = "data/scRNA_SCT.Rda")
memory.size(6*100000000000)


#### 04. PCA.FindNeighbors ####
rm(list = ls())
gc()

load(file = "data/scRNA_SCT.Rda")
scRNA <- scRNA1


scRNA
scRNA <- RunPCA(object = scRNA,
                npcs = 50,
                rev.pca = FALSE,
                weight.by.var = TRUE,
                verbose = TRUE,
                ndims.print = 1:5,
                nfeatures.print = 30,
                reduction.key = "PC_")

ElbowPlot(scRNA,
          ndims = 50)

gene <- as.data.frame(rownames(scRNA@assays$RNA))

dims <- 1:20
scRNA <- RunTSNE(scRNA,
                 dims = dims,
                 verbose = F)
DimPlot(object = scRNA,
        reduction = "tsne",
        group.by = "Sample",
        dims = c(1,2),
        shuffle = TRUE,
        label = TRUE,
        label.size = 4,
        label.color = "black",
        label.box = TRUE,
        sizes.highlight = 1)

# UMAP
scRNA <- RunUMAP(scRNA,
                 dims = dims,
                 verbose = F)
DimPlot(object = scRNA,
        reduction = "umap",
        group.by = "Sample",
        dims = c(1,2),
        shuffle = TRUE,
        label = TRUE,
        label.size = 4,
        label.color = "black",
        label.box = TRUE,
        sizes.highlight = 1)

# KNN
scRNA <- FindNeighbors(scRNA,dims = dims)

# FindClusters
for (i in c( 1, 2, 4, 8)) {
  scRNA <- FindClusters(scRNA,
                        resolution = i,
                        algorithm = 1,
                        random.seed = 2022)
  print(DimPlot(scRNA, reduction = "tsne") +
          labs(title = paste0("resolution: ", i)))
}

DimPlot(scRNA, reduction = "tsne",group.by = "integrated_snn_res.2",
        label = T)

# clustree
library(clustree)
p_cutree = clustree(scRNA@meta.data, prefix = "integrated_snn_res.")
p_cutree

Idents(scRNA)= scRNA$integrated_snn_res.2
levels(Idents(scRNA)) 
table(scRNA@active.ident)

# cell_clusters
p_dim_1 = DimPlot(scRNA,
                  pt.size = 0.5,reduction = "tsne",
                  group.by="integrated_snn_res.2", 
                  label = T)
p_dim_1
ggsave("Figure/cell_clusters.pdf", plot = p_dim_1, height = 6, width = 8)

save(scRNA, file = "data/scRNA_cluster.Rda")


#### 05. cell anno ####
load(file = "data/scRNA_cluster.Rda")

library(RColorBrewer)
# cell colors
mycols <- c("#A6CEE3","#1F78B4","#B2DF8A","#E31A1C","#FB9A99",
            "#33A02C","#FDBF6F","#FF7F00" , "#CAB2D6" ,"#6A3D9A")

# cell anno
mainmarkers <- c("SLC17A6", 
                 "GAD2","GRIK1",
                 "TH","SLC6A3",
                 "CD74","CSF3R",
                 "AQP4","SLC1A2", 
                 "VCAN","OLIG1", 
                 "MOBP","MOG",  
                 "FOXJ1", 
                 "PDGFRB", 
                 "CLDN5", "ABCB1"

)

DefaultAssay(scRNA)  <- "RNA"

# DotPlot
DotPlot(scRNA,
        features = mainmarkers ,
        dot.scale = 5,
        group.by = "integrated_snn_res.2") +
  coord_flip() +
  theme(axis.text.x = element_text(angle = 0, face = "bold", hjust = 0.5),
        axis.text.y = element_text(face = "bold"))


# cell anno
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
table(CellType$CellType)

scRNA$CellType = "NA"
Idents(scRNA) = scRNA$integrated_snn_res.2
for(i in 1:nrow(CellType)){
  scRNA@meta.data[which(scRNA@active.ident == CellType$cluster[i]), "CellType"] <- CellType$CellType[i]}
A <- as.data.frame(table(scRNA$CellType))
scRNA$CellType=factor(scRNA$CellType,
                      levels = c("Excit","Inhib","DaNs","Microglia","Astrocytes",
                                 "OPCs","ODCs", "Ependymal","Pericytes","Endothelial"))

DimPlot(scRNA,group.by = "CellType",
        reduction = "tsne",
        label = T,
        cols =  mycols)

library(SCP)
CellDimPlot(
  srt = scRNA, group.by = c("CellType"),pt.size = 0.01,split.by = "group",
  reduction = "TSNE", theme_use = "theme_blank",palcolor = mycols)

Idents(scRNA) = scRNA$CellType

mainmarkers <- c("CLDN5", "ABCB1",
                 "PDGFRB", 
                 "FOXJ1", 
                 "MOBP","MOG",  
                 "VCAN","OLIG1", 
                 "AQP4","SLC1A2", 
                 "CD74","CSF3R",
                 "TH","SLC6A3",
                 "GAD2", "GRIK1",
                 "SLC17A6"
)

DotPlot(scRNA,
        features = mainmarkers ,
        dot.scale = 1,
        group.by = "CellType") +
  scale_colour_gradient2(low = "gray", high = "black",midpoint = 0.4)+
  coord_flip() +
  theme(axis.text.x = element_text(angle = 60, face = "bold", hjust = 1),
        axis.text.y = element_text(face = "bold"))

save(scRNA, file = "data/scRNA_anno.Rda")

# FindAllMarkers
DefaultAssay(scRNA) <- "RNA"
Idents(scRNA) = scRNA$integrated_snn_res.2
table(scRNA@active.ident)
cluster_markers = FindAllMarkers(scRNA,
                                 min.pct = 0.25,
                                 logfc.threshold = 0.25)
save(cluster_markers,file = "data/all_cluster_markers.Rdata")

write.table(cluster_markers,
            file="data/all_cluster_markers.txt",
            sep="\t",quote = F,row.names = F)

# degs_sig_cluster
library(dplyr)
marker.sig <- cluster_markers %>%
  mutate(Ratio = round(pct.1/pct.2,3)) %>%
  dplyr::filter(p_val_adj <= 0.05)

degs_sig <- cluster_markers %>%
  dplyr::filter(pct.1 > 0.5 & pct.2 < 0.5 & p_val_adj <= 0.01 &
                  abs(avg_log2FC) > 2)  %>%
  arrange(cluster, -avg_log2FC)
table(degs_sig$cluster)
length(unique(degs_sig$gene))

DefaultAssay(scRNA) <- "integrated"
# top10
load("data/all_cluster_markers.Rdata")
cluster_markers <- cell_markers
top10 <- cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
sub_scRNA <- subset(scRNA,downsample = 200)
DoHeatmap(sub_scRNA,features = top10$gene) + NoLegend()
# degs_sig_cluster
write.csv(degs_sig, file = "data/cluster_degs_sig_pct0.5.csv")


# celltype markers
DefaultAssay(scRNA) <- "RNA"
Idents(scRNA) = scRNA$CellType
table(scRNA@active.ident)
cell_markers = FindAllMarkers(scRNA,
                              min.pct = 0.25,
                              logfc.threshold = 0.25)
save(cell_markers,file = "data/all_cell_markers.Rdata")

scRNA$CellType=factor(scRNA$CellType,
                      levels = c("Endothelial","Pericytes","Ependymal","ODCs","OPCs","Astrocytes","Microglia","DaNs","Inhib",
                                 "Excit"))
write.table(cell_markers,
            file="data/all_cell_markers.txt",
            sep="\t",quote = F,row.names = F)


# heatmap
library(dplyr)
library(pheatmap)
load("data/scRNA_anno.Rda")
load("data/all_cell_markers.Rdata")

# select top50 degs for heatmap
cell_markers <- cell_markers %>% dplyr::filter(is.finite(avg_log2FC)==T)
x <- table(cell_markers$gene)
cell_markers <- cell_markers[cell_markers$gene %in% names(x[x < 2]), ]
table(cell_markers$cluster)

degs_top50 <- cell_markers  %>%
  group_by(cluster) %>%
  top_n(50, avg_log2FC) %>%
  top_n(50, avg_log2FC) %>%
  arrange(cluster, -avg_log2FC)

# avgData
avgData <- scRNA@assays$RNA@data[degs_top50$gene,] %>%
  apply(1, function(x){
    tapply(x, scRNA$CellType, mean) # ExpMean
  }) %>% t

phData <- MinMax(scale(avgData), -2, 2) # z-score
rownames(phData) <- 1:nrow(phData)

# color
mycols <- c("#A6CEE3","#1F78B4","#B2DF8A","#E31A1C","#FB9A99",
            "#33A02C","#FDBF6F","#FF7F00" , "#CAB2D6" ,"#6A3D9A")

cluster_colors <- setNames(c("#A6CEE3","#1F78B4","#B2DF8A","#E31A1C","#FB9A99",
                             "#33A02C","#FDBF6F","#FF7F00" , "#CAB2D6" ,"#6A3D9A"),
                           levels(scRNA$CellType))

# heatmap
phres <- pheatmap(
  phData,
  color = colorRampPalette(c("darkblue", "white", "red3"))(99),
  scale = "row",cellwidth = 10,cellheight = 0.5,
  cluster_rows = F,
  cluster_cols = F,
  clustering_method = "complete",
  show_rownames = F,
  annotation_row = data.frame(cluster = degs_top50$cluster),
  annotation_colors = list(cluster = cluster_colors))
phres

degs_sig <- cell_markers %>%
  dplyr::filter(pct.1 > 0.5 & pct.2 < 0.5 & p_val_adj <= 0.01 &
                  abs(avg_log2FC) > 5)  %>%
  arrange(cluster, -avg_log2FC)

genes_for_interons <- intersect(rownames(scRNA@assays$integrated),degs_sig$gene)
write.table(genes_for_interons,file ="data/degs_for_interons.txt" ,row.names = F, col.names = "Gene", quote = F)
save(cell_markers,degs_sig,genes_for_interons,file = "data/degs_sig_markers.Rdata")


# save
write.csv(degs_sig, file = "data/output_degs_sig_pct0.5.csv")

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


df <- scRNA@meta.data %>% dplyr::select(group,CellType)
#control
control <- df %>% dplyr::filter(group=="CO")
control_df <- as.data.frame(table(control$CellType,control$group))
colnames(control_df) <- c("celltype","group","count")
control_df$percent  <- round((control_df$count/nrow(control))*100,2)
#disease
disease <- df %>% dplyr::filter(group=="PD")
disease_df <- as.data.frame(table(disease$CellType,disease$group))
colnames(disease_df) <- c("celltype","group","count")
disease_df$percent  <- round((disease_df$count/nrow(disease))*100,2)

data_df <- rbind(control_df,disease_df)


ggplot(data_df,aes(x=celltype,y=percent,fill=group))+
  geom_bar(stat = 'identity',
           position = 'dodge',
           width = 0.8,
           color='black')+
  labs(x=NULL)+
  theme_bw(base_size = 18)+
  theme(axis.text = element_text(colour = 'black'),
        axis.text.x = element_text(angle=30, hjust=0.4, vjust=.5))+
  scale_fill_manual(values = rep(c("#3979F2","#DF5233"),10))+
  scale_y_continuous(limits = c(0,60),breaks = c(0,20,40,60))

table(scRNA$CellType,scRNA$group)
write.table(table(scRNA$CellType,scRNA$group),file = "data/20.percent.csv",quote = F,sep = "\t")

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
selectedThresholds <- 0.273
plot.df <- scRNA@meta.data
table(plot.df$CellType,plot.df$group)
table(plot.df$CellType[plot.df$AUC>selectedThresholds],plot.df$group[plot.df$AUC>selectedThresholds])
table(plot.df$CellType)
table(plot.df$ISGscore)
df1 <- as.data.frame(table(plot.df$CellType[plot.df$AUC>selectedThresholds],plot.df$group[plot.df$AUC>selectedThresholds]))
colnames(df1) <- c("celltype","group","high_scrore_num")
df2 <- as.data.frame(table(plot.df$CellType,plot.df$group))
colnames(df2) <- c("celltype","group","total_num")
df <- merge(df1,df2)
df$low_score_num <- (df$total_num-df$high_scrore_num)
df$high_per <- round((df$high_scrore_num/df$total_num)*100,2)
df$low_per <- round((df$low_score_num/df$total_num)*100,2)
df <- df[,-(3:5)]
df <-df %>% reshape::melt(id=c("celltype","group"))
colnames(df)[3] <- "Score"
table(scRNA@assays$RNA@var.features %in% cell_markers$gene)

blue <- "#386CB0"
red <- "#E31A1C"
green <- "#1B9E77"

ggplot(df,aes(x=celltype,y=value,
              fill=Score))+
  scale_fill_manual(values = c("#D74D35","lightgray"))+
  geom_bar(stat = 'identity')+
  labs(x=NULL)+
  theme_bw(base_size = 18)+
  theme(axis.text.x = element_text(colour = 'black',angle = 60,hjust = 1))+facet_wrap(vars(group))


# DEGs  High_ISGs vs Low_ISGs
Idents(scRNA) = scRNA$ISGscore
table(scRNA@active.ident)
ISG_markers = FindAllMarkers(scRNA,
                             min.pct = 0.25,
                             logfc.threshold = 0.25)
save(ISG_markers,file = "data/ISGs_markers.Rdata")
ISG_degs_sig <- ISG_markers %>%
  dplyr::filter(pct.1 > 0.5 & pct.2 < 0.5 & p_val_adj <= 0.01 &
                  abs(avg_log2FC) > 1)  %>%
  arrange(cluster, -avg_log2FC)
table(ISG_degs_sig$cluster)
length(unique(ISG_degs_sig$gene))

# save results
write.csv(ISG_degs_sig, file = "data/ISG_markers.csv")
write.csv(ISG_degs_sig[ISG_degs_sig$cluster %in% "High_ISGs",]$gene, file = "data/High_ISG_markers.txt",
          quote = F,row.names = F,col.names = F)
Idents(scRNA) = scRNA$CellType

######################################################################################
                    ## Pseudotime and TLR analysis ##
######################################################################################


# 07. pseudotime----
library(future)
library(ggplot2)
library(tidyverse)
library(monocle)
plan()
plan("multisession", workers = 30)
plan()
options(future.globals.maxSize = 10000 * 1024^20)

# load data
load(".data/scRNA_anno.Rda")
dir.create("pseudotime")

scRNAsub=subset(scRNA, CellType == 'Microglia') #Microglia,Pericytes,Endothelial
scRNAsub
data <- as(as.matrix(scRNAsub@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNAsub@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
set.seed(1234)
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=30, relative_expr = TRUE)
mycds <- detectGenes(mycds, min_expr = 2)

# save
#save(mycds,file = "pseudotime/Microglia_pseudotime_new.Rdata")
#save(mycds,file = "pseudotime/Endothelial_pseudotime.Rdata")
#save(mycds,file = "pseudotime/Pericytes_pseudotime.Rdata")

load("pseudotime/Microglia_pseudotime_new.Rdata")

# FindVariableFeatures
scRNAsub <- FindVariableFeatures(scRNAsub, selection.method = "vst", nfeatures = 3000)
var.genes <- VariableFeatures(scRNAsub)
mycds <- setOrderingFilter(mycds, var.genes)
plot_ordering_genes(mycds)
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
p3 <- plot_ordering_genes(mycds)
p2|p2|p3

#reduceDimension
mycds <- reduceDimension(mycds, max_components = 2,
                         num_dim = 50,
                         reduction_method = 'DDRTree',
                         #residualModelFormulaStr = "~group",
                         verbose = F)
# orderCells
mycds <- orderCells(mycds)

group_color <- c("#3979F2","#DF5233")
mycols <- c("#A6CEE3","#B2DF8A","#1F78B4","#FB9A99","#E31A1C",
            "#33A02C","#FDBF6F","#FF7F00" , "#CAB2D6" ,"#6A3D9A")
scales::show_col(mycols)

plot1 <- plot_cell_trajectory(mycds,color_by = "State",cell_size=0.2)
plot2 <- plot_cell_trajectory(mycds, color_by = "group",
                              show_tree = T,
                              cell_size=0.2,
                              cell_link_size = 0.75,
                              cell_name_size = 8,
                              state_number_size = 8,
                              show_branch_points = F,
                              theta = 0)+ scale_color_manual(values = group_color)

plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime",cell_size=0.2)
plotc <- plot1|plot2|plot3
plotc

# save data
write.csv(pData(mycds), "pseudotime/pseudotime.csv")

table(mycds@phenoData@data$group,mycds@phenoData@data$State)

#BEAM analysis
disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.5&dispersion_empirical >= 1*dispersion_fit)
disp.genes <- as.character(disp.genes$gene_id)
mycds_sub <- mycds[disp.genes,]
plot_cell_trajectory(mycds_sub, color_by = "State")
beam_res <- BEAM(mycds_sub, branch_point = 1, cores = 20)
beam_res <- beam_res[order(beam_res$qval),]
beam_res <- beam_res[,c("gene_short_name", "pval", "qval")]

Time_genes <- top_n(beam_res, n = 100, rev(qval)) %>% pull(gene_short_name) %>% as.character()
p <- plot_pseudotime_heatmap(mycds[Time_genes,], num_clusters=4, show_rownames=T, return_heatmap=T)
ggsave("Figure/Time_heatmapTop100.pdf", p,width = 6, height = 16)

mycds_sub_beam <- mycds_sub[row.names(subset(beam_res, qval < 1e-6)),]
tmp1=plot_genes_branched_heatmap(mycds_sub_beam,
                                 branch_point = 1,
                                 num_clusters = 4,
                                 show_rownames = T,
                                 branch_colors = c("#979797", "#F05662", "#7990C8"),
                                 return_heatmap = T)
tmp1
# save image
pdf("Figure/branched_heatmap.pdf",width = 6,height = 16)
tmp1$ph_res
dev.off()


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
table(scRNA@meta.data$integrated_snn_res.2)
class(scRNA@meta.data$integrated_snn_res.2)
table(as.integer(scRNA@meta.data$integrated_snn_res.2))
table(as.numeric(scRNA@meta.data$integrated_snn_res.2))
cell.info <- as.data.frame(as.integer(scRNA@meta.data$integrated_snn_res.2))
colnames(cell.info) <- "ClusterID"
cell.info$Tissue<- scRNA@meta.data$Sample
cell.info$CellType <- scRNA@meta.data$CellType
rownames(cell.info) <- colnames(scRNA@assays$RNA@data)
cell.info$CellState <- paste(cell.info$Tissue, cell.info$ClusterID, sep = "_")
head(cell.info)

#initializeScenic
scenicOptions <- initializeScenic(
  org="hgnc",
  dbDir="../cisTarget_databases",
  dbs="hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
  datasetTitle="SCENIC on PD Brain Cells",
  nCores=10
)

# filter genes
genesKept <- geneFiltering(
  exprMat=dge,
  scenicOptions=scenicOptions,
  minCountsPerGene = 1,
  minSamples = 20
)
length(genesKept)

dge <- dge[genesKept, ]
dim(dge)

# AvgN 20
source("utils/AvgN.R")
avg20.rep1 <- AvgN(dge, cell.info, seed=1)

if(!dir.exists("output")) {
  dir.create("output")
}
saveRDS(cell.info, "output/s1_cell.info.rds")

# saveLoom
source("utils/add_cellAnnotation.R")
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
