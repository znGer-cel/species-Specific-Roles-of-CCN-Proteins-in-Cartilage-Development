#packages#
library(WGCNA)
library(flashClust)
library(iterators)
library(readr)
library(openxlsx)
#input data#
femData <- read_csv("chicken_data.csv")
femData <- as.data.frame(femData)
duplicates <- duplicated(femData[[1]])
femData <- femData[!duplicates, ]
rownames(femData) <- femData[[1]]
femData <- femData[,-1]
datExpr = as.data.frame(t(femData))
gsg = goodSamplesGenes(datExpr, verbose = 3) 
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}
gsg$allOK 
#outlier sample#
sampleTree = hclust(dist(datExpr), method = "average")
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 400, col = "red") 
clust = cutreeStatic(sampleTree, cutHeight = 400, minSize = 10)
table(clust)   
keepSamples = (clust==1) 
datExpr = datExpr[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
#input Traits#
allTraits = read.csv("Traits.csv", sep=',',header= 1);
dim(allTraits)  
names(allTraits)
femaleSamples = rownames(datExpr);
traitRows = match(femaleSamples, allTraits$sample);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];
collectGarbage()
sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE) 
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
#save#
save(datExpr, datTraits, file = "dataInput.RData")
###WGCNA###
library(WGCNA)
options(stringsAsFactors = FALSE);
lnames = load(file = "dataInput.RData");
lnames
#chose powers#
powers = c(seq(1, 20, by = 1)) 
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
sft$powerEstimate
par(mfrow = c(3,2));
cex1 = 0.7;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red", pos=4);
abline(h=0.7,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#in chicken power we chose 12#
net = blockwiseModules(datExpr, power = 12,
                       TOMType = "signed", minModuleSize = 100,
                       reassignThreshold = 0, mergeCutHeight = 0.3,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "TOM",
                       verbose = 3)
mergedColors = labels2colors(net$colors)
moduleColors <- mergedColors 
table(net$colors)
table(mergedColors)
sizeGrWindow(12, 9)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
#correlation heatmap#
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf(file = "module_trait_relationships.pdf", width = 10, height = 8)
par(mar = c(6, 8, 4, 2) + 0.1) 
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
moduleColors <- labels2colors(net$colors)
geneNames <- colnames(datExpr)
moduleGeneDF <- data.frame(Gene = geneNames, Module = moduleColors)
write.xlsx(moduleGeneDF, file = "WGCNA_module_genes.xlsx", rowNames = FALSE)
# Bootstrap #
set.seed(123)
nIter = 10
jaccard_list = c()
moduleColors_full<-labels2colors(net$colors)
module_full<-moduleColors_full  
for (i in 1:nIter) {
  cat("Bootstrap", i, "\n")
  idx = sample(1:nSamples, size = round(0.8 * nSamples))
  datExpr_boot = datExpr[idx, ]
  net_boot = blockwiseModules(datExpr_boot,
                              power = 12,
                              TOMType = "signed",
                              minModuleSize = 100,
                              mergeCutHeight = 0.3,
                              numericLabels = TRUE)
  
  module_boot = labels2colors(net_boot$colors)
  tab = table(module_full, module_boot)
  jaccard = mean(apply(tab, 1, function(x){ max(x)/sum(x) }))
  
  jaccard_list[i] = jaccard
}
jaccard_list
mean(jaccard_list)

#10% leave-out check#
library(WGCNA)
set.seed(124)
nSamples = nrow(datExpr)
remove_n = round(nSamples * 0.1)
remove_idx = sample(1:nSamples, remove_n)
datExpr_subset = datExpr[-remove_idx, ]
net_subset = blockwiseModules(datExpr_subset,
                              power = 12,
                              TOMType = "signed",
                              minModuleSize = 100,
                              mergeCutHeight = 0.3,
                              numericLabels = TRUE)
moduleColors_full = labels2colors(net$colors)
moduleColors_subset = labels2colors(net_subset$colors)
overlap = table(moduleColors_full, moduleColors_subset)
overlap
#check ccn list#
library(WGCNA)
library(readxl)
CCNlist <- read_excel("CCNlist.xlsx")[[1]]
module_full <- labels2colors(net$colors)
names(module_full) <- colnames(datExpr)
module_subset <- labels2colors(net_subset$colors)
names(module_subset) <- colnames(datExpr_subset)
modules_of_interest <- c("turquoise", "yellow", "greenyellow")
CCN_in_subset <- CCNlist[CCNlist %in% names(module_subset)]
stable_genes <- sum(module_subset[CCN_in_subset] %in% modules_of_interest)
stable_ratio <- stable_genes / length(CCN_in_subset)
cat("percentage of CCN regulatory‑network genes remain preserved within their original module eigengenes (MEs)", r
    ound(stable_ratio * 100, 1), "%\n")