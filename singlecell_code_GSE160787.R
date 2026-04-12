#load packages#
library(Seurat)
library(multtest)
library(dplyr)
library(ggplot2)
library(patchwork)
library(SeuratData)
library(tidyverse)
#load data and creat seuratdata
rm(list = ls())
CP <- Read10X(data.dir = "CP/")
CP <- CreateSeuratObject(counts = CP, min.features = 100, project ="CP")
D1 <- Read10X(data.dir = "D1/") 
D1 <- CreateSeuratObject(counts = D1, min.features = 100, project ="D1")
D7 <- Read10X(data.dir = "D7/") 
D7 <- CreateSeuratObject(counts = D7, min.features = 100, project ="D7")
D14 <- Read10X(data.dir = "D14/") 
D14 <- CreateSeuratObject(counts = D14, min.features = 100, project ="D14")
D28 <- Read10X(data.dir = "D28/") 
D28 <- CreateSeuratObject(counts = D28, min.features = 100, project ="D28")
D42 <- Read10X(data.dir = "D42/") 
D42 <- CreateSeuratObject(counts = D42, min.features = 100, project ="D42")
#merge data
merged_seurat <- merge(CP, y = c(D1,D7,D14,D28,D42),
                       add.cell.ids = c('CP',"D1","D7","D14","D28","D42"))
# get percent.mt
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")
#get mata data
metadata <- merged_seurat@meta.data
metadata$cells <- rownames(metadata)
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "CP"))] <- 'CP'
metadata$sample[which(str_detect(metadata$cells, "D1"))] <- 'D1'
metadata$sample[which(str_detect(metadata$cells, "D7"))] <- 'D7'
metadata$sample[which(str_detect(metadata$cells, "D14"))] <- 'D14'
metadata$sample[which(str_detect(metadata$cells, "D28"))] <- 'D28'
metadata$sample[which(str_detect(metadata$cells, "D42"))] <- 'D42'
merged_seurat@meta.data <- metadata
###QC###
filtered_seurat <- subset(merged_seurat,
                        subset = nGene > 200 & 
                          nGene < 7000 & 
                          percent.mt < 5)
#joinlayer
filtered_seurat_QC <- JoinLayers(filtered_seurat)
#RUN and remove genes >3cells
Layers(filtered_seurat_QC)
counts <- LayerData(object = filtered_seurat_QC, layer = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 3
filtered_counts <- counts[keep_genes, ]
# clean mata
matched_meta <- filtered_seurat@meta.data[colnames(filtered_counts), ]
# creat new data
filtered_seurat <- CreateSeuratObject(counts = filtered_counts,meta.data = matched_meta)
###Normalize ###
seurat_phase <- NormalizeData(filtered_seurat)
#add cell cycle information#
data("cc.genes.updated.2019") # load Seurat cell cycle data
seurat_phase <- CellCycleScoring(seurat_phase,
                                  s.features = cc.genes.updated.2019$s.genes,
                                  g2m.features = cc.genes.updated.2019$g2m.genes)
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)
seurat_phase <- ScaleData(seurat_phase)
seurat_phase <- RunPCA(seurat_phase)
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")

##########SCTransform#########
options(future.globals.maxSize = 4000 * 1024^6)
split_seurat <- SplitObject(seurat_phase, split.by = "sample")
split_seurat <- split_seurat[c('CP',"D1","D7","D14","D28","D42")]
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(
    split_seurat[[i]],
    return.only.var.genes = FALSE,
    vars.to.regress = c("percent.mt", "S.Score", "G2M.Score")
  )
}

############################ 6.Integrate ###########################
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

saveRDS(seurat_integrated, file = "TGFD_seurat.rds")
#########analysis#########
#load packages#
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(SingleR)
library(AnnotationHub)
library(ensembldb)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(pheatmap)
library(clustree)
library(remotes)
library(RColorBrewer)
library(nichenetr)
library(tidyverse)
library(CellChat)
library(readxl)
library(patchwork)
#load data#
seurat_integrated <- readRDS("TGFD_seurat.rds")
#check cell number#
metadata <- seurat_integrated@meta.data
cell_counts <- table(metadata$sample)
cell_counts_df <- as.data.frame(cell_counts)
colnames(cell_counts_df) <- c("Sample")
print(cell_counts_df)

###Reducing dimensions ###
seurat_integrated <- RunPCA(seurat_integrated, verbose = FALSE)
DimHeatmap(seurat_integrated, dims = 1:10, cells = 500, balanced = TRUE)
ElbowPlot(object = seurat_integrated, ndims = 40)
pc <- 40
seurat_integrated <- RunTSNE(seurat_integrated, reduction = "pca",dims = 1:pc)
seurat_integrated <- RunUMAP(seurat_integrated, reduction = "pca", dims = 1:pc)
DimPlot(seurat_integrated, reduction = "tsne",group.by = "ident",label=T,pt.size = 0.5)
DimPlot(seurat_integrated,reduction = "umap",label=T,pt.size = 0.3)
###Clustering,resolution chose 0.6###
seurat_integrated <- FindNeighbors(object = seurat_integrated, dims = 1:pc)
seurat_integrated <- FindClusters(object = seurat_integrated,resolution = 0.6)
seurat_integrated <- RunUMAP(seurat_integrated, reduction = "pca", dims = 1:40)
seurat_integrated <- RunTSNE(seurat_integrated,reduction = "pca",dims = 1:40)
DimPlot(seurat_integrated,reduction = "tsne",label = TRUE,label.size = 6)
DimPlot(seurat_integrated,reduction = "umap",label = TRUE,label.size = 6)
####annotation by self###
cell_mark <- c("SOX4", "SOX9",                  # neural crest cells
               "PRRX1",                         # mesenchyme without col
               "SOX2", "OTX2", "OTX1", 'PAX6','NES',    # neuro-
               "MITF", "WNT2B",                        # melano
               "COL9A1", "COL6A1", "COL2A1","COL6A3")  # chondro- 
new.cluster.ids <- c( "mesenchyme_1",    # "PRRX1"
                      "neural_crest_1",  # "SOX4" "SOX9"
                      "chondrocytes_1",  # "COL2A1", "COL6A1"
                      "neurogenic_lineage_1", #"SOX2","OTX2", "OTX1",'PAX6','NES'
                      "neurogenic_lineage_2", #"SOX2",'NES'
                      "chondrocytes_2",  # "COL2A1", "COL9A1"
                      "mesenchyme_2",    # "PRRX1"
                      "melano_1",        # "MITF","WNT2B"
                      "neural_crest_2",  # "SOX4" "SOX9"
                      "neural_crest_3",  # "SOX4" "SOX9"
                      "chondrocytes_3",  # "COL2A1","COL6A1", "COL9A1"
                      "neurogenic_lineage_3",  # "SOX2", "NES"
                      "neural_crest_4",        #"SOX4"
                      "neurogenic_lineage_4", # "SOX2"
                      "neurogenic_lineage_5", #'NES'
                      "neurogenic_lineage_6") #'NES'
names(new.cluster.ids) <- levels(seurat_integrated)
seurat_integrated <- RenameIdents(seurat_integrated, new.cluster.ids)
seurat_integrated <- AddMetaData(object = seurat_integrated, 
                                 metadata = seurat_integrated@active.ident,
                                 col.name = "celltype")
DimPlot(seurat_integrated,reduction = "umap",label = TRUE,label.size = 3)
DimPlot(seurat_integrated,reduction = "tsne",label = TRUE,label.size = 3)
DimPlot(seurat_integrated,reduction = "umap",label = T,label.size =5.8) +theme(legend.position = "none")
#visualize marker#
cell_mark <- c("SOX4", "SOX9",                  # neural crest cells
               "PRRX1",                         # mesenchyme without col
               "SOX2", "OTX2", "OTX1", 'PAX6','NES',    # neuro-
               "MITF", "WNT2B",                        # melano
               "COL9A1", "COL6A1", "COL2A1","COL6A3")  # chondro-  
DotPlot(seurat_integrated, features = cell_mark,assay = "RNA",cluster.idents = T,
        scale.by = "size",scale = TRUE,col.min = -2,col.max = 2) +
  coord_flip() +theme_bw() +  labs(x = NULL , y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1, size=14),
        axis.text.y = element_text(size=15)) +
  scale_color_gradient2(low = "#2166AC", mid = "white",high = "#B2182B")

#visualize chondrocytes marker#
cell_mark <- c("PRRX1","IGFBP5","LECT1", "EPYC","ACAN","COL9A1",
               "COL6A1", "COL2A1","COL6A3","VEGFA")  # chondrocytes 
desired_order <- c("mesenchyme_1",  
                   "mesenchyme_2",
                   "chondrocytes_1",
                   "chondrocytes_2", 
                   "chondrocytes_3", 
                   "neurogenic_lineage_1",
                   "neurogenic_lineage_2", 
                   "neurogenic_lineage_3", 
                   "neurogenic_lineage_4",
                   "neurogenic_lineage_5",
                   "neurogenic_lineage_6",  
                   "neural_crest_1", 
                   "neural_crest_2", 
                   "neural_crest_3",  
                   "neural_crest_4",       
                   "melano_1")
sub_obj <- subset(seurat_integrated, idents = desired_order)
sub_obj$seurat_clusters <- factor(Idents(sub_obj), levels = desired_order)
Idents(sub_obj) <- sub_obj$seurat_clusters
DotPlot(sub_obj , features = cell_mark,assay = "RNA",cluster.idents = F,
        scale.by = "size",scale = TRUE,col.min = -2,col.max = 2) +
  coord_flip() +theme_bw() +  labs(x = NULL , y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1, size=14),
        axis.text.y = element_text(size=15)) +
  scale_color_gradient2(low = "#2166AC", mid = "white",high = "#B2182B")

# visualize CCN family#
cell_mark <- c("CYR61","CTGF","NOV","WISP1","WISP2","WISP3",'ITGB1',
               'ITGA6','ITGAV','ITGB5','PLXNA1','IGF2R','NOTCH1','FGFR2','LRP6','ITGA5','LRP1','SDC4')
rename_genes <- c(CYR61 = "CCN1",CTGF  = "CCN2",NOV   = "CCN3",WISP1 = "CCN4",WISP2 = "CCN5",WISP3 = "CCN6")
DotPlot(sub_obj, features = cell_mark,assay = "RNA",cluster.idents = F,scale.by = "size",
        scale = TRUE,col.min = -2,col.max = 2) +
  coord_flip() +theme_bw() +labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1, size=14),
        axis.text.y = element_text(size=15)) +
  scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B"
  ) +scale_x_discrete(labels = rename_genes)
###
old_genes <- c("CYR61","CTGF","NOV","WISP1","WISP2","WISP3")
new_gene_names <- c(CYR61 = "CCN1",CTGF  = "CCN2",NOV   = "CCN3",WISP1 = "CCN4",
                    WISP2 = "CCN5",WISP3 = "CCN6")
DefaultAssay(seurat_integrated) <- "RNA"
plist <- FeaturePlot(seurat_integrated,features = old_genes,cols = c("lightgrey", "#FF0000"),combine = FALSE)
for (i in seq_along(plist)) {gene <- old_genes[i]
  plist[[i]] <- plist[[i]] + ggtitle(new_gene_names[gene])}
wrap_plots(plist, ncol = 2)

###stacked bar plot###
selected_celltypes <- c("mesenchyme_1",  
                        "mesenchyme_2",
                        "chondrocytes_1",
                        "chondrocytes_2", 
                        "chondrocytes_3", 
                        "neurogenic_lineage_1",
                        "neurogenic_lineage_2", 
                        "neurogenic_lineage_3", 
                        "neurogenic_lineage_4",
                        "neurogenic_lineage_5",
                        "neurogenic_lineage_6",  
                        "neural_crest_1", 
                        "neural_crest_2", 
                        "neural_crest_3",  
                        "neural_crest_4",       
                        "melano_1")  
metadata <-seurat_integrated@meta.data
metadata$celltype <- factor(metadata$celltype, levels = selected_celltypes)
total_per_sample <- metadata %>%dplyr::group_by(sample) %>%dplyr::summarise(total_cells = n())
celltype_percent <- metadata %>%
  dplyr::filter(celltype %in% selected_celltypes) %>% 
  dplyr::group_by(sample, celltype) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  dplyr::left_join(total_per_sample, by = "sample") %>%  
  dplyr::mutate(percent = n / total_cells * 100)  
sample_order <- c("CP", "D1", "D7", "D14", "D28", "D42")
celltype_percent$sample <- factor(celltype_percent$sample, levels = sample_order)
n_colors <- length(selected_celltypes)
celltype_colors <- colorRampPalette(brewer.pal(12, "Set3"))(n_colors)
names(celltype_colors) <- selected_celltypes
ggplot(celltype_percent, aes(x = sample, y = percent, fill = celltype))+ 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),color = "black", width = 0.75) + 
  #geom_text(aes(label = sprintf("%.1f", percent)),position = position_dodge(width = 0.8), vjust = -0.5, size = 3, color = "black") + 
  facet_wrap(~ sample, scales = "free_x") +labs(x = "Sample", y = "Cell Type Proportion (%)",fill = NULL) + 
  theme_classic()+theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),  
                       axis.text.y = element_text(size = 6),
                       axis.title = element_text(size = 15, face = "bold"),
                       legend.position = "top",
                       legend.text = element_text(size = 11.5),         
                       strip.text = element_text(size = 15, face = "bold"),
                       strip.background = element_rect(fill = "gray90", color = NA),
                       panel.grid.major.y = element_line(color = "gray90", size = 0.15))+ 
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
###VlnPlot CCN family###
multi_subset <- subset(seurat_integrated, subset = celltype %in% c("mesenchyme_1", "mesenchyme_2","chondrocytes_1", "chondrocytes_2","chondrocytes_3"))
DefaultAssay(multi_subset) <- "RNA"
celltype_order <- c("mesenchyme_1", "mesenchyme_2","chondrocytes_1", "chondrocytes_2","chondrocytes_3")
multi_subset@meta.data$celltype <- factor(multi_subset@meta.data$celltype, levels = celltype_order)
sample_order <- c("CP", "D1", "D7", "D14", "D28", "D42")
multi_subset@meta.data$sample <- factor(multi_subset@meta.data$sample, levels = sample_order)
sample_colors <- c("#1F77B4",  # CP
                   "#FF7F0E",  # D1
                   "#D62728",  # D7
                   "#2CA02C",  # D14
                   "#9467BD",  # D28
                   "#8C564B"   # D42
)
names(sample_colors) <- sample_order
VlnPlot(multi_subset, features = "NOV",assay = "RNA",
        group.by = "celltype",   
        split.by = "sample",      
        pt.size = 0.1,
        log = FALSE,
        cols = sample_colors     
) + ggtitle("NOV Expression") +
  xlab("Cell Type") + ylab("Expression Level") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.position = "right",legend.title = element_text(face = "bold")) +
  guides(fill = guide_legend(title = "Sample")) 

##########################AUC######################
library(readxl)
library("Matrix")
library("AUCell")
terms<- read_excel("list.xlsx")
print(terms)
geneset<- terms[['CCN_proteins_regulatory_network']]
geneset <- geneset[!is.na(geneset)]
valid_genes <- intersect(geneset, rownames(seurat_integrated))
if (length(valid_genes) < 5) stop
gene_sets <- list("CCN_proteins_regulatory_network" = valid_genes)
seurat_integrated <- JoinLayers(seurat_integrated)
expr_matrix <- LayerData(seurat_integrated, layer = "data")
expr_matrix <- as(expr_matrix, "dgCMatrix")              
cells_rankings <- AUCell_buildRankings(expr_matrix,plotStats = FALSE,  splitByBlocks = TRUE)
cells_AUC <- AUCell_calcAUC(gene_sets,cells_rankings,
                            aucMaxRank = ceiling(0.05 * nrow(cells_rankings)),nCores = 4)
seurat_integrated$CCN_proteins <- as.numeric(getAUC(cells_AUC)["CCN_proteins_regulatory_network", ])
#visualize#
library(ggplot2)
umap_df <- FetchData(seurat_integrated, vars = c("umap_1", "umap_2", "CCN_proteins"))
ggplot(umap_df, aes(x = umap_1, y = umap_2, color = CCN_proteins)) +
  geom_point(size = 0.8, alpha = 0.8) +
  scale_color_viridis_c(option = "viridis",name = "AUC Score" ) +
  labs(x = "UMAP1", y = "UMAP2", 
       #title = "CCN_proteins_regulatory_network (AUCell)" 
  ) +
  theme_classic()
#visualize boxplot#
sample_order <- c("CP", "D1", "D7", "D14", "D28", "D42")
seurat_integrated@meta.data$sample <- factor(
  seurat_integrated@meta.data$sample, 
  levels = sample_order)
umap_df <- FetchData(seurat_integrated, vars = c("umap_1", "umap_2", "CCN_proteins", 'sample')) 
ggplot(umap_df, aes(x = sample, y = CCN_proteins, fill = sample)) +
  geom_boxplot(width = 0.6, outlier.size = 0.5) +
  scale_fill_brewer(palette = "Set3") +  
  labs(x = "sample", y = "AUC Score", 
       #title = "CCN_proteins_regulatory_network "
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 16, angle = 60, hjust = 1),  
        axis.text.y = element_text(size = 15),  
        axis.title = element_text(size = 18),  
        legend.position = "none" )

##visualize_AUC_Bubble_Plot## 
auc_df <- seurat_integrated@meta.data %>%group_by(celltype) %>%
  summarise(auc_mean = mean(CCN_proteins),n_cells = n())
ggplot(auc_df, aes(x = celltype, y = 1,size = n_cells,color = auc_mean)) +
  geom_point(alpha = 0.9) +
  scale_size(range = c(3, 15)) +
  scale_color_gradientn(colors = c("#2166AC", "white", "#B2182B")) +
  theme_bw() +
  theme(axis.title.y = element_blank(),axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),axis.text.x = element_text(angle =90, hjust = 1, size = 15)) +
  labs(x = "Cell Type",color = "AUC Score",size = "Number of Cells",
       title = "AUC Activity Across Cell Types")

###nichenetr###
organism <- "human"
lr_network<-readRDS("lr_network_human_21122021.rds")
ligand_target_matrix<-readRDS("ligand_target_matrix_nsga2r_final.rds")
weighted_networks<-readRDS("weighted_networks_nsga2r_final.rds")
lr_network <- lr_network %>% distinct(from, to)
head(lr_network)
ligand_target_matrix[1:5,1:5] 
head(weighted_networks$lr_sig) 
head(weighted_networks$gr) 
ccn_mapping <- data.frame(
  nichenet_name = c("CCN1", "CCN2", "CCN3", "CCN4", "CCN6"),
  sc_name = c("CYR61", "CTGF", "NOV", "WISP1","WISP3")
)
lr_network_ccn <- lr_network %>%
  dplyr::filter(from %in% ccn_mapping$nichenet_name) %>%
  left_join(ccn_mapping, by = c("from" = "nichenet_name")) %>%
  mutate(from = sc_name) %>%
  dplyr::select(from, to) %>%
  distinct()
seurat_integrated[["RNA"]] <- JoinLayers(seurat_integrated[["RNA"]])
expression_matrix <- LayerData(seurat_integrated, assay = "RNA", layer = "data")
metadata <- seurat_integrated@meta.data
#sender_cells and receiver_cells#
sender_cells = rownames(metadata)[metadata$celltype %in% c("mesenchyme_1", "mesenchyme_2")]  
receiver_cells = rownames(metadata)[metadata$celltype %in% c("chondrocytes_1", "chondrocytes_2","chondrocytes_3")]  
ccn_ligands = c("CYR61", "CTGF", "NOV", "WISP1","WISP3") 
ccn_network <- lr_network_ccn %>%
  dplyr::filter(from %in% ccn_ligands & to %in% weighted_networks$lr_sig$to) %>%
  dplyr::distinct(from, to)
library(igraph)
net <- graph_from_data_frame(ccn_network)
plot(net, edge.arrow.size=0.2, vertex.label.cex=0.7)
potential_receptors <- unique(ccn_network$to)
receptor_expression <- AverageExpression(
  seurat_integrated,
  features = potential_receptors,
  group.by = "celltype",
  assays = "RNA"
)$RNA
expressed_receptors <- rownames(receptor_expression)[
  rowSums(receptor_expression[, c("chondrocytes-1", "chondrocytes-2", "chondrocytes-3")] > 0.05) > 0
]
final_ccn_network <- ccn_network %>%
  dplyr::filter(to %in% expressed_receptors)

Idents(seurat_integrated) <- "celltype"
de_genes <- FindMarkers(seurat_integrated, 
                        ident.1 = "chondrocytes_1", 
                        ident.2 = c("chondrocytes_2", "chondrocytes_3"),  
                        only.pos = TRUE,
                        logfc.threshold = 0.25,min.pct = 0.25) %>% 
  rownames_to_column("gene") %>% 
  dplyr::filter(p_val_adj < 0.05)
#calculate CCN ligand activity#
ccn_ligands_sc <- ccn_mapping$sc_name
ccn_ligand_matrix <- ligand_target_matrix[, ccn_mapping$nichenet_name] 
colnames(ccn_ligand_matrix) <- ccn_ligands_sc
receiver_cell_ids <- rownames(metadata)[metadata$celltype %in% c("chondrocytes_1", "chondrocytes_2", "chondrocytes_3")]
expression_matrix <- LayerData(seurat_integrated, assay = "RNA", layer = "data")  
background_expressed_genes <- rownames(expression_matrix)[
  rowSums(expression_matrix[, receiver_cell_ids] > 0) > 0]
ligand_activities <- nichenetr::predict_ligand_activities(
  geneset = de_genes$gene,                     
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ccn_ligand_matrix,     
  potential_ligands = ccn_ligands              
)
ligand_activities %>% 
  arrange(-aupr_corrected) %>% 
  knitr::kable()
targets <- lapply(ccn_ligands_sc, function(ligand) {
  if (ligand %in% colnames(ccn_ligand_matrix)) {
    ligand_scores <- ccn_ligand_matrix[, ligand]
    top_targets <- names(sort(ligand_scores, decreasing = TRUE))[1:200]
    return(top_targets)
  } else {
    warning(paste("ligand", ligand, "gap"))
    return(character(0))
  }
})
names(targets) <- ccn_ligands_sc 
targets <- targets[sapply(targets, length) > 0]
all_targets <- unique(unlist(targets))
#KEGG#
library(clusterProfiler)
library(org.Hs.eg.db)
entrez_ids <- bitr(all_targets, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
kegg_result <- enrichKEGG(gene = entrez_ids$ENTREZID,organism = 'hsa',pvalueCutoff = 0.05,qvalueCutoff = 0.2)
dotplot(kegg_result, showCategory=15, title="KEGG Pathway Enrichment of CCN Targets")
#visualize strangth sender to receiver#
ligand_expression <- AverageExpression(seurat_integrated, features = ccn_ligands_sc,
                                       group.by = "celltype",assays = "RNA")$RNA
receptor_expression <- AverageExpression(seurat_integrated,features = expressed_receptors,
                                         group.by = "celltype",assays = "RNA")$RNA
#heatmap#
interaction_strength <- matrix(0, nrow = length(ccn_ligands_sc),ncol = length(expressed_receptors),
                               dimnames = list(ccn_ligands_sc, expressed_receptors))
for(lig in ccn_ligands_sc) {
  for(rec in expressed_receptors) {
    if(any(final_ccn_network$from == lig & final_ccn_network$to == rec)) {
      ligand_exp <- mean(ligand_expression[lig, c("mesenchyme-1", "mesenchyme-2")])
      receptor_exp <- mean(receptor_expression[rec, c("chondrocytes-1", "chondrocytes-2")])
      interaction_strength[lig, rec] <- ligand_exp * receptor_exp
    }
  }
}
rename_ligands <- c(CYR61 = "CCN1",CTGF  = "CCN2",NOV   = "CCN3",WISP1 = "CCN4",WISP2 = "CCN5",WISP3 = "CCN6")
nichenetr::make_heatmap_ggplot(
  interaction_strength,
  "CCN ligands", "Receptors",
  color = "purple",
  legend_title = "Interaction\nStrength"
) +
  scale_y_discrete(labels = rename_ligands) +   
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.title.x = element_text(size = 18, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 18),
    legend.position = "bottom"
  )

####Cell chat####
library(CellChat)
library(patchwork)
library(RColorBrewer)
library(ggalluvial)
options(stringsAsFactors = FALSE)
seurat_integrated[["RNA"]] <- JoinLayers(seurat_integrated[["RNA"]])
data.input <- LayerData(seurat_integrated, assay = "RNA", layer = "data")
metadata <- data.frame(cell_type = Idents(seurat_integrated),
                       row.names = colnames(seurat_integrated))
cellchat <- createCellChat(object = data.input)
cellchat <- addMeta(cellchat, meta = metadata)
cellchat <- setIdent(cellchat, ident.use = "cell_type")
cellchat@idents <- factor(cellchat@idents, 
                          levels = levels(seurat_integrated$celltype))
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10) # 常规参数
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
roupSize <- as.numeric(table(cellchat@idents))
par(mar = c(0,0,2,0))
groupSize <- as.numeric(table(cellchat@idents))
colors <- colorRampPalette(brewer.pal(12, "Paired"))(length(groupSize))
netVisual_circle(cellchat@net$weight,
                 color.use = colors,
                 vertex.weight = groupSize,
                 weight.scale = TRUE, 
                 title.name = "Cell-Cell Communication")
######################Findmarkers ###################
DefaultAssay(seurat_integrated) <- "RNA"
seurat_integrated <- JoinLayers(seurat_integrated)
markers <- FindAllMarkers(object = seurat_integrated, only.pos = TRUE,logfc.threshold = 0.25)  
top10_markers <- markers %>% 
  group_by(cluster) %>% 
  arrange(p_val_adj, .by_group = TRUE) %>% 
  slice_head(n = 10)
write.csv(top10_markers, file = "markers_top10.csv", row.names = FALSE)