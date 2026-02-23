# 1. Load qs file and setup metadata ###########################################
library(qs)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(GSVA)

gse18 <- qread("gse18v2_ana_1227-935104.qs")
DimPlot(gse18, group.by = c("seurat_clusters"), label = TRUE) 

# Update condition metadata based on sample names
gse18@meta.data$condition[gse18$samplename %in% 
                            c("MSC4-Dura", "MSC5-Dura", "S1", "S2")] <- "Dura"
gse18@meta.data$condition[gse18$samplename %in% 
                            c("MSC5-BTI", "MSC6-BTI","BTI935100","BTI936394")] <- "Meningioma" 
gse18@meta.data$condition[gse18$samplename %in% 
                            c("S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10",
                              "Co935100","Co936394")] <- "Meningioma"
table(gse18@meta.data$samplename, gse18@meta.data$condition) 

# Define global color palettes
color1 <- c(  
  "#F5B7B1", 
  "#6b7cc7", 
  "#FF5733", 
  "#FFC300", 
  "#7ab55c", 
  "#00BFFF",
  "#9B59B6"
) 

DimPlot(gse18, group.by = c("seurat_clusters"), label = TRUE, cols = color1) 
DimPlot(gse18, group.by = c("seurat_clusters"), label = FALSE, cols = color1) 
DimPlot(gse18, group.by = c("seurat_clusters"), split.by = "condition", label = FALSE, cols = color1) 

color2 <- c("#E3BBED", "#7EA6D9") 
DimPlot(gse18, group.by = c("condition"), label = FALSE, cols = color2) 

# 2. Dot Plot for Major Cell Types #############################################
marker_genes <- c(
  "SFRP2","CRABP2",  "THSD4","KRT18",    
  "CD14","CD68","FCGR3A","LYZ",          
  "CD3E", "NKG7","CD3D",                 
  "CD79A","MS4A1","CD19",                
  "PLVAP","CLDN5", "VWF",                
  "RGS5","MCAM",  "ACTA2",               
  "PLP1", "MAG", "MBP"                   
) 

celltype_order <- c("1", "0", "2", "12", "4", "7", "11") 

DotPlot(gse18,
        features = marker_genes,
        group.by = "seurat_clusters",
        cols = c("white", "darkred"),
        dot.scale = 8,
        cluster.idents = FALSE) +
  scale_y_discrete(limits = celltype_order) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Marker Genes", y = "Cell Clusters") +
  theme(panel.grid = element_line(color = "grey90")) 

# 3. Calculate Average Counts and Proportions per Cluster by Disease ###########

cluster_colors_manual <- c(
  "0" = "#F5B7B1",  
  "1" = "#6b7cc7",  
  "2" = "#FF5733",  
  "4" = "#FFC300",  
  "7" = "#7ab55c",  
  "11" = "#00BFFF",  
  "12" = "#9B59B6"
) 

# Data preprocessing
plot_data <- gse18@meta.data[, c("condition", "seurat_clusters")] %>% 
  mutate(seurat_clusters = factor(as.numeric(as.character(seurat_clusters)))) %>% 
  filter(!seurat_clusters %in% c("5", "6")) %>%  
  mutate(condition = factor(condition, levels = sort(unique(condition)))) 

# Generate count matrix (only for clusters with data)
count_df <- plot_data %>% 
  count(condition, seurat_clusters, name = "cells") %>% 
  filter(cells > 0) 

write.csv(count_df, "count_df.csv", row.names = TRUE) 

# NOTE: The division by 4 and 14 is hardcoded. Consider making this dynamic based on sample count.
count_df <- count_df %>%
  mutate(
    cell_average = case_when(
      condition == "Dura" ~ cells / 4,
      condition == "Meningioma" ~ cells / 14,
      TRUE ~ NA_real_ 
    )
  ) 

# Stacked Absolute Cell Numbers Bar Plot
ggplot(count_df, aes(x = condition, y = cell_average, fill = seurat_clusters)) +
  geom_col(position = "stack", width = 0.8) +
  scale_fill_manual(
    values = cluster_colors_manual,
    breaks = names(cluster_colors_manual),
    guide = guide_legend(reverse = FALSE)
  ) +
  labs(title = "Stacked Absolute Cell Numbers (Non-zero clusters only)", 
       x = "Condition", y = "Number of Cells", fill = "Cluster") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    axis.text = element_text(color = "black", size = 9)
  ) 

# Relative Proportion Bar Plot
ggplot(count_df, aes(x = condition, y = cells, fill = seurat_clusters)) +
  geom_col(position = "fill", width = 0.8) +
  scale_fill_manual(values = cluster_colors_manual) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Cellular Proportions within Conditions (Non-zero clusters only)", 
       x = "Condition", y = "Proportion", fill = "Cluster") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text = element_text(color = "black", size = 9)) 

# 4. Compare Tumor vs Non-Tumor Cells ##########################################

# Merge vascular and immune clusters into cluster 0; extract others and glia
gse18_sub <- gse18
gse18_sub@meta.data$seurat_clusters[gse18_sub@meta.data$seurat_clusters %in% c(2, 12, 7, 4)] <- 0 
gse18_sub <- subset(gse18_sub, seurat_clusters %in% c(0, 1)) 

DimPlot(gse18_sub, group.by = c("seurat_clusters"),
        split.by = "condition",     
        reduction = "umap",      
        label = FALSE, cols = color1,
        ncol = 3) 

# Calculate expression proportions for gene screening 
SY_genes <- c("GJA1","DNM1L","HOMER3","NCAM1", "CTTN", "ARHGEF7", "CAMKK2","DNM2", 
              "PLCB1", "BAIAP2","ARHGEF9", "AKAP1", "CAMK2G","KALRN", "RASGRF2", 
              "HOMER2","HOMER1", "DLG4", "DNM1", "CNTNAP1", "TIAM1", "CNTNAP2",
              "SRC","NTRK2",  "CAMK2D", "GAP43", "PLCB4",  "NLGN2",  "NRXN2",
              "SYNGAP1", "EPHA4", "NRXN3",  "PLCL2", "GRASP", "NLGN3", "DNM3",
              "CNTN1", "NLGN1", "AKAP2", "GRIA3", "GRID2", "GRIP2", "SHANK2", 
              "CTTNBP2", "GRIA1", "SHANK3", "GRIP1", "CAMK2B", "CAMKK1", "GRIN2D",
              "HCN2", "AKAP5", "GRM7", "NRXN1", "NRCAM", "CNTN2", "TTYH1", "GRIK2", 
              "GRM3", "NOS1","GRIN1", "CAMK2A", "SHANK1", "LRRC7", "GRM1", "GRIA2", 
              "CAMKV", "GRIN2A", "GRM2", "GRIA4","GRM5", "GRIN2B") 

available_genes <- SY_genes[SY_genes %in% rownames(gse18_sub)] 
print(paste("Number of available SY genes:", length(available_genes))) 

cluster_info <- gse18_sub@meta.data$seurat_clusters 
results <- data.frame() 

for(gene in available_genes) {
  expression_data <- FetchData(gse18_sub, vars = gene) 
  for(cluster in unique(cluster_info)) {
    cluster_cells <- which(cluster_info == cluster) 
    # Calculate expression proportion (expression > 0)
    expr_percentage <- mean(expression_data[cluster_cells, 1] > 0) * 100 
    results <- rbind(results, data.frame(Gene = gene, Cluster = cluster, ExpressionPercentage = expr_percentage)) 
  }
}

results_wide <- reshape2::dcast(results, Gene ~ Cluster, value.var = "ExpressionPercentage") 
write.csv(results_wide, "SY_genes_TC-nTC-935104_cluster.csv", row.names = FALSE) 

# 5. Custom Dot Plot for Filtered SY Genes #####################################

gene_set <- c("CTTN", "GJA1", "NTRK2", "ARHGEF7", "CAMK2D", "DNM1L", 
              "DNM2", "NCAM1", "PLCB1", "BAIAP2", "AKAP1", "PLCB4", 
              "ARHGEF9", "CAMK2G", "RASGRF2", "GAP43", "NLGN2", 
              "HOMER3", "KALRN", "NRXN2", "CAMKK2", "HOMER2", 
              "NRXN3", "PLCL2", "EPHA4", "GRASP", "TIAM1", "CTTNBP2") 

existing_genes <- gene_set[gene_set %in% rownames(gse18_sub)] 
expression_data <- FetchData(gse18_sub, vars = c(existing_genes, "seurat_clusters")) 

# Calculate average expression and percentage
avg_exp <- expression_data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(existing_genes), mean)) 

pct_exp <- expression_data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(existing_genes), ~ mean(. > 0) * 100)) 

# Reshape and merge
avg_exp_long <- avg_exp %>% pivot_longer(cols = -seurat_clusters, names_to = "gene", values_to = "avg_exp") 
pct_exp_long <- pct_exp %>% pivot_longer(cols = -seurat_clusters, names_to = "gene", values_to = "pct_exp") 

plot_data <- left_join(avg_exp_long, pct_exp_long, by = c("seurat_clusters", "gene")) %>%
  mutate(seurat_clusters = factor(seurat_clusters), gene = factor(gene, levels = existing_genes)) 

ggplot(plot_data, aes(x = gene, y = seurat_clusters)) +
  geom_point(aes(size = pct_exp, color = avg_exp)) +
  scale_color_gradient(low = "lightgray", high = "darkred", name = "Average Expression") +
  scale_size(range = c(1, 12), name = "Percent Expressed") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_blank(),
        panel.grid.major = element_line(color = "grey90", size = 0.2),
        legend.position = "right") +
  labs(title = "SY Gene Expression by Cell Clusters") +
  guides(size = guide_legend(override.aes = list(color = "black"))) 

# 6. GLU Genes Processing (Repeats similar workflow as SY genes) ###############
GLU_genes <- c("GRIN1", "GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D", "DLGAP1",
               "GRIA1", "GRIA2", "GRIA3", "GRIA4",
               "GRIK1", "GRIK2", "GRIK3", "GRIK4", "GRIK5") 

available_genes <- GLU_genes[GLU_genes %in% rownames(gse18_sub)] 
results <- data.frame() 

for(gene in available_genes) {
  expression_data <- FetchData(gse18_sub, vars = gene) 
  for(cluster in unique(cluster_info)) {
    cluster_cells <- which(cluster_info == cluster) 
    expr_percentage <- mean(expression_data[cluster_cells, 1] > 0) * 100 
    results <- rbind(results, data.frame(Gene = gene, Cluster = cluster, ExpressionPercentage = expr_percentage)) 
  }
}

results_wide <- reshape2::dcast(results, Gene ~ Cluster, value.var = "ExpressionPercentage") 
write.csv(results_wide, "GLU_genes_TC-nTC-935104_cluster.csv", row.names = FALSE) 

gene_set <- c("DLGAP1", "GRIA3", "GRIA1", "GRIK5", "GRIN2D", "GRIK1") 
existing_genes <- gene_set[gene_set %in% rownames(gse18_sub)] 
expression_data <- FetchData(gse18_sub, vars = c(existing_genes, "seurat_clusters")) 

avg_exp <- expression_data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(existing_genes), mean)) 

pct_exp <- expression_data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(existing_genes), ~ mean(. > 0) * 100)) 

avg_exp_long <- avg_exp %>% pivot_longer(cols = -seurat_clusters, names_to = "gene", values_to = "avg_exp") 
pct_exp_long <- pct_exp %>% pivot_longer(cols = -seurat_clusters, names_to = "gene", values_to = "pct_exp") 

plot_data <- left_join(avg_exp_long, pct_exp_long, by = c("seurat_clusters", "gene")) %>%
  mutate(seurat_clusters = factor(seurat_clusters), gene = factor(gene, levels = existing_genes)) 

ggplot(plot_data, aes(x = gene, y = seurat_clusters)) +
  geom_point(aes(size = pct_exp, color = avg_exp)) +
  scale_color_gradient(low = "lightgray", high = "darkred", name = "Average Expression") +
  scale_size(range = c(1, 12), name = "Percent Expressed") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_blank(),
        panel.grid.major = element_line(color = "grey90", size = 0.2),
        legend.position = "right") +
  labs(title = "GLU Gene Expression by Cell Clusters") +
  guides(size = guide_legend(override.aes = list(color = "black"))) 

# 7. Calculate scores uniformly using GSVA #####################################

SY_gene <- c("DLGAP1", "GRIA3", "GRIA1", "GRIK5", "GRIN2D", "GRIK1") 

gene_set_list <- list(SYscore_ssGSEA = SY_gene) 

gse18_TC_sub_human <- subset(gse18_sub, seurat_clusters %in% c(0,1)) 
table(gse18_TC_sub_human$seurat_clusters) 

# Extract expression matrix. (Corrected "RNA3" typo to "RNA")
expr_matrix <- GetAssayData(gse18_TC_sub_human, slot = "data", assay = "RNA") 

# Run GSVA
# Note: kcdf="Gaussian" is for log-normalized data (like Seurat's data slot).
# kcdf="Poisson" is for raw integer counts.
print("Running ssGSEA, this may take a few minutes...") 
ssgsea_scores <- gsva(
  expr = as.matrix(expr_matrix), 
  gset.idx.list = gene_set_list,
  method = "ssgsea",
  kcdf = "Gaussian", 
  verbose = TRUE     
) 

gse18_TC_sub_human <- AddMetaData(
  gse18_TC_sub_human,
  metadata = t(ssgsea_scores),
  col.name = "SYscore_ssGSEA" 
) 

# Violin Plot for GSVA scores
p <- VlnPlot(gse18_TC_sub_human, group.by = "seurat_clusters", features = c("SYscore_ssGSEA"), pt.size = 0) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.size = 0, alpha = 0.6) +
  theme_classic() 

# Adjust labels and calculate statistics
score_range <- range(gse18_TC_sub_human@meta.data$SYscore_ssGSEA, na.rm = TRUE) 
global_label_y <- score_range[2] * 1.15  # Adjust this coefficient to move global labels 

existing_clusters <- as.character(unique(gse18_TC_sub_human@meta.data$seurat_clusters)) 
comparisons <- list() 
if(all(c("0", "1") %in% existing_clusters)) comparisons <- c(comparisons, list(c("0", "1"))) 

if(length(comparisons) > 0) {
  p <- p + stat_compare_means(
    method = "wilcox.test",
    comparisons = comparisons,
    label = "p.format",
    tip.length = 0.02,
    step.increase = 0.12,  
    size = 3,
    vjust = 0.5,
    label.y = global_label_y * 0.9  
  )
} 

p <- p + stat_compare_means(
  method = "kruskal.test",
  label.y = global_label_y,  
  size = 3
) 

p <- p + ylim(score_range[1], global_label_y * 2) 
p + labs(title = "SYscore_ssGSEA across clusters",
         x = "Cluster",
         y = "SYscore_ssGSEA") +
  theme(plot.title = element_text(hjust = 0.5)) 

# Secondary Violin Plot
VlnPlot(gse18_TC_sub_human, features = c("SYscore_ssGSEA"), group.by = "seurat_clusters" ,pt.size = 0, 
        cols = c("#90bff9","#ff8080" )) 

metadata <- gse18_TC_sub_human@meta.data 
score_summary <- metadata %>%
  group_by(seurat_clusters) %>%
  summarise(
    mean_score = mean(SYscore_ssGSEA, na.rm = TRUE),
    median_score = median(SYscore_ssGSEA, na.rm = TRUE),
    sd_score = sd(SYscore_ssGSEA, na.rm = TRUE),
    se_score = sd_score / sqrt(n()),
    n_cells = n(),
    .groups = 'drop'
  ) 

# 8. Focus purely on Tumor Cells ###############################################

gse18_TC_sub <- qread("gse18_TC_sub-935104-3_1227.qs") 
DimPlot(gse18_TC_sub, group.by = c("seurat_clusters"),
        split.by = "condition",     
        reduction = "umap",      
        label = TRUE,
        ncol = 3, cols = color1) 

color3 <- c("#82093B", "#585D5E", "#008A45") 
DimPlot(gse18_TC_sub, group.by = c("seurat_clusters"), label = FALSE, cols = color3) 

color4 <- c("#ee7c7c", "#90b8e2") 
DimPlot(gse18_TC_sub, group.by = c("condition"), label = FALSE, cols = color4) 

# Tumor Subpopulation Dot Plot
marker_genes <- c(
  "MKI67","CENPK","MAD2L1","TUBB2B","CTHRC1", "ADAMTS6",  
  "ATP5F1D","ATP5MC2","SELENOM","SELENOP",                
  "FBLN1","MFAP5","COL8A1","MEG3","FTX", "KCNQ1OT1"       
) 

markers_C1 <- FindMarkers(gse18_TC_sub, ident.1 = 1, max.cells.per.ident = 500) 
write.csv(markers_C1, "markersC1.csv", row.names = TRUE) 

celltype_order <- c("4", "0", "1") 

DotPlot(gse18_TC_sub,
        features = marker_genes,
        group.by = "seurat_clusters",
        cols = c("white", "darkred"),
        dot.scale = 8,
        cluster.idents = FALSE) +
  scale_y_discrete(limits = celltype_order) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Marker Genes", y = "Cell Clusters") +
  theme(panel.grid = element_line(color = "grey90")) 

# 9. Evaluate MKI67 Expression in Subsets ######################################

gse18_TC_sub_BTI <- subset(gse18_TC_sub, condition %in% "MSC_BTI") 
gse18_TC_sub_Core <- subset(gse18_TC_sub, condition %in% "MSC_Core") 

feature_data <- FetchData(gse18_TC_sub, vars = "MKI67") 
new_min <- min(feature_data$MKI67) 
new_max <- 2 # Or adjust dynamically as needed 

# Plot Features with identical color scales for fair comparison
FeaturePlot(gse18_TC_sub_BTI, features = c("MKI67"), pt.size = 0) +  
  scale_color_gradient(low = "darkblue", high = "#FFD700", limits = c(new_min, new_max)) 

FeaturePlot(gse18_TC_sub_Core, features = c("MKI67"), pt.size = 0) +  
  scale_color_gradient(low = "darkblue", high = "#FFD700", limits = c(new_min, new_max))
