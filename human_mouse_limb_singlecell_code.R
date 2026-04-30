###human limb code###
library(Seurat)
library(sceasy)
#reticulate::py_install("anndata")
#convert H5AD
#getwd()
#sceasy::convertFormat("skeletalatlas_cellxgene.h5ad",from = "anndata",to = "seurat",outFile = "object.rds")
#read in and ensure filtered
limb1<-readRDS("object.rds")
#
#visualize CCN#
Idents(limb1) <- limb1$fineanno
all_types<- unique(limb1$fineanno)
unique_fine <- unique(limb1$fineanno)
grep("CC", unique_fine, ignore.case = TRUE, value = TRUE)

priority_groups <- c("Proximal Mesenchyme", "Proximal Mesenchyme 2",
                     "Distal Mesenchyme","Distal Mesenchyme 2", 
                     "Distal Mesenchyme 3","ISL1 Mesenchyme","Chondroprogenitor",
                     "Early Resting CC","Resting CC","Proliferating CC",        
                     "Cycling CC","Prehypertrophic CC","Hypertrophic CC","Articular CC",
                     "Pre Osteoblast 2", "Pre Osteoblast 1", "Early Osteoblast", 
                     "Mature Osteoblast", "HEY1 Osteoblast", 
                     "Osteocyte", "Osteoclast Precursor", "Mature Osteoclast"
)
other_groups <- setdiff(all_types, priority_groups)
limb1$fineanno <- factor(
  limb1$fineanno,
  levels = c(priority_groups, other_groups)
)
Idents(limb1) <- limb1$fineanno
#FeaturePlot(limb1, features = "CCN1", reduction = "umap",cols = c("lightgrey", "#B2182B"))
library(ggplot2)
library(dplyr)
cell_mark <- c("CCN1","CCN2","CCN3","CCN4","CCN5","CCN6",'ITGB1') 
DotPlot(limb1, features = cell_mark,assay = "RNA",cluster.idents = F,
        scale.by = "size",scale = TRUE,col.min = -2,col.max = 2) +
  coord_flip() + 
  theme_bw() +  labs(x = NULL , y = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=1, size=10),
        axis.text.y = element_text(size=10)) +
  scale_color_gradient2(low = "#2166AC", mid = "white",high = "#B2182B")
#AUC#
library(readxl)
library(Seurat)
library(UCell)
library(dplyr)
library(ggplot2)
inputlist <- read_excel("HUMANCCNlist.xlsx")
geneset <- inputlist$human
geneset <- geneset[!is.na(geneset)]
valid_genes <- intersect(geneset, rownames(limb1))
if (length(valid_genes) < 5) {
  stop("Valid genes < 5, too few genes for UCell")
}
limb1 <- AddModuleScore_UCell(
  limb1,
  list(CCN = valid_genes),
  name = "CCN"
)
UCell_threshold <- 0.1
auc_df <- limb1@meta.data %>%
  group_by(fineanno) %>%
  summarise(
    auc_mean  = mean(CCNCCN, na.rm = TRUE),
    n_cells   = n(),
    pos_ratio = sum(CCNCCN > UCell_threshold, na.rm = TRUE) / n(),
    .groups   = "drop"
  )
auc_df$fineanno <- factor(
  auc_df$fineanno,
  levels = c(priority_groups, other_groups)
)
ggplot(
  auc_df,
  aes(
    x     = fineanno,
    y     = auc_mean,         
    size  = pos_ratio,        
    color = auc_mean
  )
) +
  geom_point(alpha = 0.9) +
  scale_size_continuous(
    range   = c(1, 6),                     
    name    = "Positive Fraction",
    labels  = scales::percent_format()
  ) +
  scale_color_gradientn(
    colors = c("#2166AC", "white", "#B2182B"),
    name   = "Mean UCell score"
  ) +
  theme_bw() +
  theme(
    axis.text.y  = element_text(size = 12),   
    axis.title.y = element_text(size = 14),
    axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11),
    plot.title   = element_text(size = 16, face = "bold")
  ) +
  labs(
    x     = "Cell type",
    y     = "Mean UCell Score",
    title = "CCN network activity across cell types (human)"
  )
###mouse limb code###
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
E11 <- Read10X(data.dir = "E11/")
E11 <- CreateSeuratObject(counts = E11, min.features = 100, project ="E11")
E13 <- Read10X(data.dir = "E13/") 
E13 <- CreateSeuratObject(counts = E13, min.features = 100, project ="E13")
E15 <- Read10X(data.dir = "E15/") 
E15 <- CreateSeuratObject(counts = E15, min.features = 100, project ="E15")
E18 <- Read10X(data.dir = "E18/") 
E18 <- CreateSeuratObject(counts = E18, min.features = 100, project ="E18")
#merge data
merged_seurat <- merge(E11, y = c(E13,E15,E18),
                       add.cell.ids = c('E11',"E13","E15","E18"))
# get percent.mt
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^mt-")
#get mata data
metadata <- merged_seurat@meta.data
metadata$cells <- rownames(metadata)
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "E11"))] <- 'E11'
metadata$sample[which(str_detect(metadata$cells, "E13"))] <- 'E13'
metadata$sample[which(str_detect(metadata$cells, "E15"))] <- 'E15'
metadata$sample[which(str_detect(metadata$cells, "E18"))] <- 'E18'
merged_seurat@meta.data <- metadata
###QC###
filtered_seurat <- subset(merged_seurat,
                          subset = nGene > 200 & 
                            nGene < 6000 & 
                            percent.mt < 10)
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
### Human → Mouse cell cycle gene orthology ###
library(biomaRt)
data("cc.genes.updated.2019")
human_s <- cc.genes.updated.2019$s.genes
human_g2m <- cc.genes.updated.2019$g2m.genes
mouse_mart <- useMart(
  "ensembl",
  dataset = "mmusculus_gene_ensembl",
  host = "https://dec2021.archive.ensembl.org/"
)
human_mart <- useMart(
  "ensembl",
  dataset = "hsapiens_gene_ensembl",
  host = "https://dec2021.archive.ensembl.org/"
)
s_map <- getLDS(
  attributes  = "external_gene_name",
  filters     = "external_gene_name",
  values      = human_s,
  mart        = human_mart,
  attributesL = "external_gene_name",
  martL       = mouse_mart
)
colnames(s_map) <- c("human", "mouse")
head(s_map)
g2m_map <- getLDS(
  attributes  = "external_gene_name",
  filters     = "external_gene_name",
  values      = human_g2m,
  mart        = human_mart,
  attributesL = "external_gene_name",
  martL       = mouse_mart
)
colnames(g2m_map) <- c("human", "mouse")
head(g2m_map)
library(dplyr)
mouse_s   <- s_map  %>% distinct(human, mouse) %>% pull(mouse) %>% unique()
mouse_g2m <- g2m_map %>% distinct(human, mouse) %>% pull(mouse) %>% unique()
length(mouse_s)
length(mouse_g2m)
### Add mouse cell cycle scores ###
seurat_phase <- CellCycleScoring(seurat_phase,s.features = mouse_s,g2m.features = mouse_g2m,set.ident = TRUE)
#ScaleData#
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
split_seurat <- split_seurat[c('E11',"E13","E15","E18")]
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- SCTransform(
    split_seurat[[i]],
    return.only.var.genes = FALSE,
    vars.to.regress = c("percent.mt", "S.Score", "G2M.Score")
  )
}
#####Integrate ####
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features, 
                                   assay = "SCT")
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
seurat_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")
saveRDS(seurat_integrated, file = "mouselimb_seurat.rds")

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
seurat_integrated <- readRDS("mouselimb_seurat.rds")
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
pc <- 30
seurat_integrated <- RunTSNE(seurat_integrated, reduction = "pca",dims = 1:pc)
seurat_integrated <- RunUMAP(seurat_integrated, reduction = "pca", dims = 1:pc)
Idents(seurat_integrated) <- seurat_integrated$sample
DimPlot(seurat_integrated, reduction = "tsne",group.by = "ident",label=T,pt.size = 0.5)
DimPlot(seurat_integrated,reduction = "umap",label=T,pt.size = 0.3)
###joinLayers combine counts and data###
seurat_integrated [['RNA']] <- JoinLayers(seurat_integrated [['RNA']])
#Dimensionality reduction and clustering using UMAP
seurat_integrated  <- FindNeighbors(seurat_integrated,dims = 1:pc)
# resolution 0.1 - 1 ,and need adjust!!
seurat_integrated  <- FindClusters(seurat_integrated ,resolution = seq(from = 0.1,to = 1.0, by = 0.1))
seurat_integrated  <- RunUMAP(seurat_integrated ,dims = 1:pc)
seurat_integrated  <- RunTSNE(seurat_integrated ,dims = 1:pc)
###to chose resolution 
library(clustree)
clustree(seurat_integrated ) 
DimPlot(seurat_integrated ,reduction = 'umap',group.by = "sample")
# chose 0.5 resolution
seurat_integrated <- FindClusters(object = seurat_integrated,resolution = 0.5)
seurat_integrated <- RunUMAP(seurat_integrated, reduction = "pca", dims = 1:pc)
seurat_integrated <- RunTSNE(seurat_integrated,reduction = "pca",dims = 1:pc)
DimPlot(seurat_integrated,reduction = "tsne",label = TRUE,label.size = 6)
DimPlot(seurat_integrated,reduction = "umap",label = TRUE,label.size = 6)
####annotation###
cell_mark <- c("Hoxd13",              # cartilage_precursors
               "Krt14",               # skin
               "Cdh5",                # vasculature 
               "Lyz2",                # blood
               "Col2a1", "Col1a1",    # mature_cartilage/bone/tendon 
               "Col10a1",             # growth_plate  
               "Myog")                # muscle  
cell_mark <- c("Hoxd13",                 # cartilage_precursors
               "Krt14",                  # skin
               "Cdh5",                   # vasculature 
               "Lyz2",                   # blood
               "Col2a1", "Acan",         # mature_cartilage(Col10a1-)
               "Col10a1",                # growth_plate  
               "Runx2","Sp7", "Bglap",   # mature bone 
               "Scx",                    # tendon_lineage 
               "Col1a1",                 # fibroblast-like (Col1a1⁺)
               "Myog")                   # muscle  
DotPlot(seurat_integrated , features = cell_mark,assay = "RNA",cluster.idents = F,
        scale.by = "size",scale = TRUE,col.min = -2,col.max = 2) +
  coord_flip() + 
  theme_bw() +  labs(x = NULL , y = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=1, size=10),
        axis.text.y = element_text(size=10)) +
  scale_color_gradient2(low = "#2166AC", mid = "white",high = "#B2182B")
#cell.cluster.ids#
new.cluster.ids <- c( "mature_cartilage_1",               # "Col2a1", "Acan"
                      "mature_cartilage_2",               # "Col2a1", "Acan"
                      "cartilage_precursors_1",          # "Hoxd13"
                      "vasculature_1",                    # "Cdh5"
                      "skin_1",                           # "Krt14"
                      "tendon_lineage_1",                 # "Scx","Col1a1", 
                      "Col1a1⁺_mesenchymal_cells_1",                # "Col1a1"
                      "muscle_1",                         # "Myog"
                      "cartilage_precursors_2",           # "Hoxd13"
                      "growth_plate_1",                   # "Col10a1"
                      "cartilage_precursors_3",          # "Hoxd13"
                      "blood_1",                          # "Lyz2"
                      "skin_2",                           # "Krt14"
                      "muscle_2",                         # "Myog"
                      "blood_2",                          # "Lyz2"
                      "Col1a1⁺_mesenchymal_cells_2",                # "Col1a1"
                      "growth_plate_2"   )                # "Col10a1"
names(new.cluster.ids) <- levels(seurat_integrated)
seurat_integrated <- RenameIdents(seurat_integrated, new.cluster.ids)
seurat_integrated <- AddMetaData(object = seurat_integrated, 
                                 metadata = seurat_integrated@active.ident,
                                 col.name = "celltype")
DimPlot(seurat_integrated,reduction = "umap",label = TRUE,label.size = 3)
DimPlot(seurat_integrated,reduction = "tsne",label = TRUE,label.size = 3)
DimPlot(seurat_integrated,reduction = "umap",label = T,label.size =3) +theme(legend.position = "none")


# CCN family DotPlot 
library(Seurat)
library(ggplot2)
library(dplyr)
cell_mark <- c("Cyr61","Ctgf","Nov","Wisp1","Wisp2","Wisp3","Itgb1")
rename_genes <- c(Cyr61 = "Ccn1",Ctgf = "Ccn2",Nov   = "Ccn3",
                  Wisp1 = "Ccn4",Wisp2 = "Ccn5",Wisp3 = "Ccn6",
                  Itgb1 = "Itgb1")

## 
desired_order <- c("Col1a1⁺_mesenchymal_cells_1","Col1a1⁺_mesenchymal_cells_2",
                   "cartilage_precursors_1","cartilage_precursors_2","cartilage_precursors_3",
                   "mature_cartilage_1","mature_cartilage_2",
                   "growth_plate_1","growth_plate_2",
                   "tendon_lineage_1",
                   "muscle_1","muscle_2",
                   "vasculature_1","blood_1","blood_2","skin_1","skin_2")
## 
levels(Idents(seurat_integrated))
Idents(seurat_integrated) <- factor(
  Idents(seurat_integrated),
  levels = desired_order
)
levels(Idents(seurat_integrated))
## 
DotPlot(
  seurat_integrated,
  features = cell_mark,
  assay = "RNA",
  cluster.idents = FALSE,
  scale.by = "size",
  scale = TRUE,
  col.min = -2,
  col.max = 2
) +
  coord_flip() +
  theme_bw() +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 12),
    axis.text.y = element_text(size = 15),
    plot.margin = margin(10, 10, 10, 60)
  ) +
  scale_color_gradient2(
    low  = "#2166AC",
    mid  = "white",
    high = "#B2182B"
  ) +
  scale_x_discrete(labels = rename_genes)

#AUC#
library(readxl)
library("Matrix")
library("AUCell")
inputlist<- read_excel("mouseCCNlist.xlsx")
geneset<- inputlist$mouse
geneset <- geneset[!is.na(geneset)]
geneset_list <- geneset 
valid_genes <- intersect(geneset, rownames(seurat_integrated))
if (length(valid_genes) < 5) {
  stop("Valid genes < 5, too few genes for AUCell")
}
gene_sets <- list(
  CCN_family = valid_genes
)
seurat_integrated[["RNA"]] <- JoinLayers(seurat_integrated[["RNA"]])
expr_matrix <- LayerData(seurat_integrated, layer = "data")
expr_matrix <- as(expr_matrix, "dgCMatrix")
cells_rankings <- AUCell_buildRankings(
  expr_matrix,
  plotStats = FALSE,
  splitByBlocks = TRUE
)
cells_AUC <- AUCell_calcAUC(
  gene_sets,
  cells_rankings,
  aucMaxRank = ceiling(0.05 * nrow(cells_rankings)),
  nCores = 4
)
seurat_integrated$CCN_AUC <- as.numeric(
  getAUC(cells_AUC)["CCN_family", ]
)

library(dplyr)
library(ggplot2)
AUC_threshold <- 0.1
desired_order
auc_df <- seurat_integrated@meta.data %>%
  group_by(celltype) %>%
  summarise(
    auc_mean  = mean(CCN_AUC, na.rm = TRUE),
    n_cells   = n(),
    pos_ratio = sum(CCN_AUC > AUC_threshold, na.rm = TRUE) / n()
  ) %>%
  ungroup()
auc_df$celltype <- factor(auc_df$celltype, levels = desired_order)
ggplot(auc_df, aes(x = celltype, y = auc_mean, size = pos_ratio, color = auc_mean)) +
  geom_point(alpha = 0.9) +
  scale_size_continuous(
    range = c(2, 10),
    name = "Positive Fraction",
    labels = scales::percent_format()
  ) +
  scale_color_gradientn(
    colors = c("#2166AC", "white", "#B2182B"),
    name = "Mean AUC"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    plot.title   = element_text(size = 16, face = "bold"),
    plot.margin  = margin(10, 10, 10, 40)   
  ) +
  labs(
    x = "Cell Type",
    y = "Mean AUCell Score",
    title = "CCN network activity across cell types (mouse)"
  )
