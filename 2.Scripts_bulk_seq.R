######################################################################################
                                  ## Bulk analysis ##
######################################################################################

# Load necessary libraries
library(tidyverse)
library(data.table)
library(GEOquery)
library(limma)
library(hgu133plus2.db)
library(ggsci)
library(cowplot)
library(clusterProfiler)

# Initialization
options(stringsAsFactors = F)
Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)

# Define a function to load GSE data
load_GSE_data <- function(GSE_id) {
  # Construct filename
  filename <- paste0("Raw data/", GSE_id, "_series_matrix.txt.gz")
  # Load data from the GEO database
  GSE_data <- getGEO(GEO = GSE_id, filename = filename, getGPL = F)
  # Extract expression data and sample information
  expr <- exprs(GSE_data)
  pd <- pData(GSE_data)
  # Data processing
  pd_filtered <- pd %>%
    filter(str_detect(title, "Control|BR")) %>%
    mutate(group = if_else(str_detect(title, "Control"), "control", "disease")) %>%
    arrange(group, geo_accession)
  expr_filtered <- expr[, rownames(pd_filtered)] %>%
    as.data.frame()
  # Return list
  list(expr = expr_filtered, pd = pd_filtered)
}

# Load data
GSE49036_data <- load_GSE_data("GSE49036")
GSE49036_expr <- GSE49036_data$expr
GSE49036_pd <- GSE49036_data$pd

# Normalization
GSE49036_norm <- as.data.frame(normalizeBetweenArrays(as.matrix(GSE49036_expr)))

# Plot boxplot after normalization
pdf("1.After_norm_boxplots_49036.PDF", 8, 8)
boxplot(GSE49036_norm, las=2, outline=FALSE, main = "GSE49036")
dev.off()

# Annotation and data merging
gpl_GSE49036 <- GSE49036_sm@annotation
ann_49036 <- toTable(hgu133plus2SYMBOL) %>%
  rename(id = V1, symbol = V2) %>%
  filter(symbol != "", !duplicated(symbol)) %>%
  mutate(median = apply(GSE49036_norm[.$id,], 1, median)) %>%
  arrange(desc(median))

# Perform PCA analysis
perform_PCA <- function(data) {
  library(factoextra)
  pca_result <- PCA(t(data), graph = F)
  fviz_pca_ind(pca_result,
               title = "Principal Component Analysis",
               legend.title = "Groups",
               geom.ind = c("point"),
               pointsize = 1.5,
               labelsize = 4,
               col.ind = GSE49036_pd$group, # Ensure GSE49036_pd includes a 'group' column
               axes.linetype = NA,
               mean.point = F,
               addEllipses = TRUE) + coord_fixed(ratio = 1)
  print(pca_result)
  ggsave("PCA_GSE49036.pdf", width = 7.5, height = 6)
}

# Execute PCA analysis
perform_PCA(GSE49036_ann)


# perform_GSVA
perform_GSVA <- function(expression_data) {
  library(GSVA)
  geneset <- clusterProfiler::read.gmt("h.all.v7.5.1.symbols.gmt")
  gs <- lapply(split(geneset$gene, geneset$term), unique)
  gsva_res <- gsva(as.matrix(expression_data), gs, method = "ssgsea", kcdf = "Poisson")
  write.csv(gsva_res, "GSVA_GSE49036.csv", quote = F)
  return(gsva_res)
}

# Execute GSVA analysis
gsva_res <- perform_GSVA(GSE49036_ann)

# plot_gene_expression
plot_gene_expression <- function(gene, expression_data, metadata) {
  library(ggplot2)
  expression_df <- as.data.frame(t(expression_data[gene, ]))
  names(expression_df) <- c("expression")
  expression_df$group <- metadata$group[match(rownames(expression_df), rownames(metadata))]
  ggplot(expression_df, aes(x = group, y = expression, fill = group)) +
    geom_boxplot() +
    scale_fill_manual(values = c("#3979F2", "#DF5233")) +
    labs(title = paste("Expression of", gene), x = "Group", y = "Expression Level") +
    theme_minimal()
}

# Plot boxplot for a specific gene
plot_gene_expression("RUNX2", GSE49036_ann, GSE49036_pd) # NFATC2,IRF5

# perform_pathway_analysis
perform_pathway_analysis <- function(gsva_scores, group_list) {
  library(limma)
  design <- model.matrix(~ 0 + factor(group_list))
  colnames(design) <- levels(factor(group_list))
  fit <- lmFit(gsva_scores, design)
  contrast.matrix <- makeContrasts(disease-control, levels = design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  top_table <- topTable(fit2, coef = "disease-control", adjust = "BH")
  print(head(top_table))
}

# Execute pathway analysis
perform_pathway_analysis(gsva_res, GSE49036_pd$group)
