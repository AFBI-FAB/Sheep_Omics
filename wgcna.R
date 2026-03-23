library(WGCNA)
library(flashClust)
library(curl)
library(DESeq2)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(tidyverse)
library(CorLevelPlot)
library(clusterProfiler)
library(pathview)
library(stringr)
library(AnnotationHub)
library(ggridges)
library(enrichplot)

out_dir="WGCNA/"
#load("WGCNA_Norm_counts_workspace.RData")
countData=read.csv("countsdata.csv",header=T,row.names=1,sep=",", check.names = FALSE)
head(countData)
keep_genes <- rowSums( countData > 5 ) >= 7
countData2 <- countData[keep_genes,]
dim(countData)
dim(countData2)

# Read the metadata
sample_metadata = read.csv(file = "metaData.csv", sep=",", row.names = 1)
head(sample_metadata)


#####################################################################################
# Normalization
###################################################################################
# create a deseq2 dataset
data = countData2
data <- data[,unique(rownames(sample_metadata))]
all(colnames(data) == rownames(sample_metadata))
# create dds
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = sample_metadata,
                              design = ~ 1) # not spcifying model
dim(dds)
# perform variance stabilization
dds_norm <- vst(dds)
# get normalized counts
norm.counts <- assay(dds_norm) 
norm.counts.t <- t(norm.counts)

# rerun sample clustering to see any change
sampleTree <- hclust(dist(norm.counts.t), method = "average") #Clustering samples based on distance 
hclust = pdf("0.hclust-samples_on_normalized_counts_with_33_samples.pdf")
par(cex = 0.6);
par(mar = c(0,4,2,0))
# Plotting the cluster dendrogram
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
dev.off()


###########################################################################################
# QC - outlier detection
###########################################################################################

# detect outlier samples - hierarchical clustering - method 1

sampleTree <- hclust(dist(norm.counts.t), method = "average") #Clustering samples based on distance 
sampleConn <- adjacency(t(norm.counts.t), type ="unsigned")
Z.k <- scale(apply(sampleConn, 2, sum))
Z.k


# Outlier sample removal
outlier_samples <- c("7094", "6977", "6764") # based on Z.k less than -2
countData_clean <- norm.counts.t[!(rownames(norm.counts.t) %in% outlier_samples),]
dim(countData_clean)
metaData_clean <- sample_metadata[setdiff(rownames(sample_metadata), outlier_samples), ]
dim(metaData_clean)

# change newly created vars to old count and metadata just for ease of use
data <- countData_clean
sample_metadata <- metaData_clean
dim(data)
dim(sample_metadata)


# detect outlier genes
gsg <- goodSamplesGenes(data)
summary(gsg)
gsg$allOK
table(gsg$goodGenes) # 25030 genes passed
table(gsg$goodSamples)

# # if allOK returen false, remove genes that are detectd as outliers
# data <- data[gsg$goodGenes == TRUE,]
# dim(data)


# rerun sample clustering to see any change
sampleTree <- hclust(dist(data), method = "average") #Clustering samples based on distance 
hclust = pdf("1.hclust-samples_on_normalized_counts.pdf")
par(cex = 0.6);
par(mar = c(0,4,2,0))
# Plotting the cluster dendrogram
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
             cex.axis = 1.5, cex.main = 2)
dev.off()



##################################################################
# PCA
#################################################################

norm.counts <- data
dim(norm.counts)

# detect outlier samples - pca - method 2
pca <- prcomp(norm.counts.t)
pca.dat <- pca$x
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)
pca.dat <- as.data.frame(pca.dat)
pcaplot = file.path(out_dir, "2.PCA_norm_counts.pdf")
plot2=ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))
ggsave(filename = pcaplot, plot = plot2) 



###################################################################################################
# Choosing the soft-thresholding power: analysis of network topology
###################################################################################################
# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts, powerVector = powers, verbose = 5)
# Plot the results:
pdf(paste0(out_dir,"3.soft-power.pdf"),height=7,width=7.5)
par(mfrow = c(1,2))
# Set some parameters
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot3 = plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
             xlab="Soft Threshold (powers)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", 
             main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.85,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot3 = plot(sft$fitIndices[,1], sft$fitIndices[,5], 
             xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
             main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
cex1 = 0.9
print(plot3)
dev.off()

sft

# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)
softPower <- 4 #7 for all samples and 4 for 3 samles removed
# calling adjacency function
adjacency <- adjacency(norm.counts, power = softPower, type="signed")
# TOM
enableWGCNAThreads(nThreads = 16)  # adjust based on your CPU
TOM <- TOMsimilarity(adjacency, TOMType = "signed")#This gives similarity between genes
TOM.dissimilarity <- 1-TOM # get dissimilarity matrix

# Hierarchical Clustering Analysis
#The dissimilarity/distance measures are then clustered using linkage hierarchical clustering and a dendrogram (cluster tree) of genes is constructed.

# creating the dendrogram 
pdf(paste0(out_dir, "4.dendrogram.gene.clustering.TOM.dissimilarity.pdf"))
hclustGeneTree <- hclust(as.dist(TOM.dissimilarity), method = "average")

# Plot the resulting clustering tree (dendogram)

plot4 = plot(hclustGeneTree, xlab = "", sub = "", 
             main = "Gene Clustering on TOM-based disssimilarity", 
             labels = FALSE, hang = 0.04)
print(plot4)
dev.off()

##############################################################################
# identify modules
##############################################################################
# Make the modules larger, so set the minimum higher
minModuleSize <- 50
# Module ID using dynamic tree cut
dynamicMods <- cutreeDynamic(dendro = hclustGeneTree, method="hybrid",cutHeight=0.995,
                             distM = TOM.dissimilarity,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)

table(dynamicMods)#returns a table of the counts of factor levels in an object. In this case how many genes are assigned to each created module.
# Convert numeric labels into colors
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)#returns the counts for each color (aka the number of genes within each module) 
pdf(paste0(out_dir, "5.dendrogram.with.module.colors.pdf"))
# Plot the dendrogram and colors underneath
plot5 = plotDendroAndColors(hclustGeneTree, dynamicColors, "Dynamic Tree Cut", 
                            dendroLabels = FALSE, hang = 0.03, 
                            addGuide = TRUE, guideHang = 0.05, 
                            main = "Gene dendrogram and module colors")
print(plot5)
dev.off()

##############################################################
# Module Eigengenes
#############################################################
# Calculate the module eigengenes
dynamic_MEList <- moduleEigengenes(norm.counts, colors = dynamicColors)
dynamic_MEs <- dynamic_MEList$eigengenes
head(dynamic_MEs)
#To further condense the clusters (branches) into more meaningful modules you can cluster modules based on pairwise eigengene correlations 
#and merge the modules that have similar expression profiles.
# Calculate dissimilarity of module eigengenes
dynamic_MEDiss <- 1-cor(dynamic_MEs) #Calculate eigengene dissimilarity
dynamic_METree <- hclust(as.dist(dynamic_MEDiss), method = "average")#Clustering eigengenes 
pdf(paste0(out_dir, "6.clustering of module eigen genes.pdf"))
# Plot the hclust
plot6 = plot(dynamic_METree, main = "Dynamic Clustering of module eigengenes",
             xlab = "", sub = "",)
abline(h=.25, col = "red") #a height of .25 corresponds to correlation of .75
print(plot6)
dev.off()

########################################################################
# Finding different cutheights
#######################################################################
# Assume you already have datExpr and dynamicColors from WGCNA

# Function to test different cutHeights
test_cutHeight <- function(norm.counts, dynamicColors, cutHeights = c(0.15, 0.20, 0.25, 0.30)) {
  results <- data.frame(cutHeight = numeric(), nModules = integer())
  
  for (h in cutHeights) {
    merge <- mergeCloseModules(norm.counts, dynamicColors,
                               cutHeight = h,
                               verbose = 0)
    mergedColors <- merge$colors
    nModules <- length(unique(mergedColors))
    
    results <- rbind(results, data.frame(cutHeight = h, nModules = nModules))
  }
  
  return(results)
}

# Run the test
results <- test_cutHeight(norm.counts, dynamicColors)
print(results)

# compare module sizes
library(dplyr)

check_module_sizes <- function(norm.counts, dynamicColors, cutHeight) {
  merge <- mergeCloseModules(norm.counts, dynamicColors, cutHeight = cutHeight, verbose = 0)
  mergedColors <- merge$colors
  table(mergedColors) %>% sort(decreasing = TRUE)
}

sizes_015 <- check_module_sizes(norm.counts, dynamicColors, 0.15)
sizes_015
sizes_020 <- check_module_sizes(norm.counts, dynamicColors, 0.20)
sizes_020
sizes_025 <- check_module_sizes(norm.counts, dynamicColors, 0.25)
sizes_025
sizes_030 <- check_module_sizes(norm.counts, dynamicColors, 0.30)
sizes_030

########################################################################
# Merge similar modules
########################################################################
dynamic_MEDissThres <- 0.25
# Call an automatic merging function
merge_dynamic_MEDs <- mergeCloseModules(norm.counts, dynamicColors, cutHeight = dynamic_MEDissThres, verbose = 3)
# The Merged Colors
dynamic_mergedColors <- merge_dynamic_MEDs$colors
# Eigen genes of the new merged modules
mergedMEs <- merge_dynamic_MEDs$newMEs
table(dynamic_mergedColors)
# dendrogram with original and merged modules
pdf(paste0(out_dir, "7.Dendrogram-with-original_and-merged-modules.pdf"))
plot7 = plotDendroAndColors(hclustGeneTree, cbind(dynamicColors, dynamic_mergedColors),
                            c("Dynamic Tree Cut", "Merged dynamic"),
                            dendroLabels = FALSE, hang = 0.03,
                            addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors for original and merged modules")
print(plot7)
dev.off()

# Rename Module Colors 
moduleColors <- dynamic_mergedColors

# Construct numerical labels corresponding to the colors 
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder)-1
MEs <- mergedMEs
#save(MEs, moduleLabels, moduleColors, hclustGeneTree, file = "wgcna-networkConstruction.RData")

###############################################################
# Network visualization
#############################################################
nGenes = ncol(norm.counts)
nSamples = nrow(norm.counts)

## Calculate topological overlap anew: this could be done more efficiently by saving the TOM
## calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(norm.counts, power = softPower);
## Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
## Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
## Call the plot function
sizeGrWindow(9,9)

pdf("Network_heatmap_plot.pdf")
TOMplot(plotTOM, hclustGeneTree, moduleColors, main = "Network heatmap plot, all genes")
dev.off()

nSelect = 400
## For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
## There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
## Open a graphical window
sizeGrWindow(9,9)
## Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
## the color palette; setting the diagonal to NA also improves the clarity of the plot
##png('123.png')
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")


#####################################################################################
# Relating modules to external traits
#####################################################################################
# pull out all continuous traits. Before that convertiing ID column which is a factor to numeric type
sample_metadata_original = sample_metadata
sample_metadata$DMI_acc_BW <- sample_metadata$DMI_performance/sample_metadata$Final_BW
#str(sample_metadata)
allTraits <- sample_metadata[,c("ADG_performance", "DMI_performance", "DMI_acc_BW", "RFI")]
table(rownames(MEs) == rownames(allTraits))

# define numbers of genes and samples
nGenes <- ncol(norm.counts)
nSamples <- nrow(norm.counts)

# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(norm.counts, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
names(MEs) <- substring(names(MEs), 3)
moduleTraitCor <- cor(MEs, allTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
# # Floor p-values to 2 decimal places
# moduleTraitPvalue_floor <- floor(moduleTraitPvalue * 100) / 100
# print(moduleTraitCor)
# print(moduleTraitPvalue)
# print(moduleTraitPvalue_floor)
# 
# # find modules significant for both traits
# sig_modules <- rownames(moduleTraitPvalue)[moduleTraitPvalue[,"DMI_performance"] <= 0.05 & moduleTraitPvalue[,"ADG_performance"] <= 0.05 & moduleTraitPvalue[,"DMI_acc_BW"] <= 0.05]
# sig_modules

# create module-trait heatmap
# Will display correlations and their p-values
pdf(paste0(out_dir, "8.Module-trait-heatmap.pdf"), width = 10, height=11)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap
plot8 = labeledHeatmap(Matrix = moduleTraitCor, 
                       xLabels = names(allTraits),
                       yLabels = names(MEs), 
                       ySymbols = names(MEs), 
                       colorLabels = FALSE, 
                       colors = blueWhiteRed(100),
                       textMatrix = textMatrix, 
                       setStdMargins = FALSE,
                       cex.text = 0.9,
                       zlim = c(-1,1),
                       main = paste("Module-trait Relationships"))
print(plot8)	
dev.off()

table(dynamic_mergedColors)
#Each row corresponds to a module eigengene, and the columns correspond to a trait. 
#Each cell contains a p-value and correlation. Those with strong positive correlations are shaded a darker red while those with stronger negative correlations become more blue. 
# heatmap with significance as stars (*)
heatmap.data <- merge(MEs , allTraits, by = 'row.names')
head(heatmap.data)
heatmap.data <- heatmap.data %>%   column_to_rownames(var = 'Row.names')
colnames(heatmap.data)
pdf(paste0(out_dir, "9.module trait relationships with significance.pdf"), width = 10, height = 7)
colnames(heatmap.data)[18:21] <- c("ADG", "DMI_absolute", "DMI_adjusted", "RFI")
plot9 = CorLevelPlot(heatmap.data,
                     x = names(heatmap.data)[18:21],
                     y = names(heatmap.data)[1:17],
                     col = c("red", "orange", "white", "lightgreen", "green"),
                     signifSymbols = c("***", "**", "*", ""),
                     signifCutpoints = c(0, 0.001, 0.01, 0.05, 1),
                     rotLabX = 30, rotLabY = 10)
print(plot9)
dev.off()
table(dynamic_mergedColors)

###################################
#LIST OF GEnes in each modules
#####################################

topGOgenes <- colnames(norm.counts)[moduleColors =="darkorange"]
write.csv(topGOgenes,"FE_darkorange_genes.txt")

topGOgenes <- colnames(norm.counts)[moduleColors =="darkgrey"]
write.csv(topGOgenes,"FE_darkgrey_genes.txt")

topGOgenes <- colnames(norm.counts)[moduleColors =="black"]
write.csv(topGOgenes,"FE_black_genes.txt")

topGOgenes <- colnames(norm.counts)[moduleColors =="lightcyan"]
write.csv(topGOgenes,"FE_lightcyan_genes.txt")

topGOgenes <- colnames(norm.counts)[moduleColors =="cyan"]
write.csv(topGOgenes,"FE_cyan_genes.txt")

topGOgenes <- colnames(norm.counts)[moduleColors =="purple"]
write.csv(topGOgenes,"FE_purple_genes.txt")


topGOgenes <- colnames(norm.counts)[moduleColors =="lightcyan"]
write.csv(topGOgenes,"FE_lightcyan_genes.txt")

topGOgenes <- colnames(norm.counts)[moduleColors =="cyan"]
write.csv(topGOgenes,"FE_cyan_genes.txt")

topGOgenes <- colnames(norm.counts)[moduleColors =="purple"]
write.csv(topGOgenes,"FE_purple_genes.txt")

topGOgenes <- colnames(norm.counts)[moduleColors =="lightcyan"]
write.csv(topGOgenes,"FE_lightcyan_genes.txt")

topGOgenes <- colnames(norm.counts)[moduleColors =="cyan"]
write.csv(topGOgenes,"FE_cyan_genes.txt")

topGOgenes <- colnames(norm.counts)[moduleColors =="purple"]
write.csv(topGOgenes,"FE_purple_genes.txt")

topGOgenes <- colnames(norm.counts)[moduleColors =="lightgreen"]
write.csv(topGOgenes,"FE_lightgreen_genes.txt")

topGOgenes <- colnames(norm.counts)[moduleColors =="darkgrey"]
write.csv(topGOgenes,"FE_darkgrey_genes.txt")

topGOgenes <- colnames(norm.counts)[moduleColors =="blue"]
write.csv(topGOgenes,"FE_blue_genes.txt")

topGOgenes <- colnames(norm.counts)[moduleColors =="magenta"]
write.csv(topGOgenes,"FE_magenta_genes.txt")

##############################################
#13. correlation values within a heatmap plot #
##############################################
dissTOM = 1-TOMsimilarityFromExpr(norm.counts, power = 4);
plotTOM = dissTOM^7;
diag(plotTOM) = NA;
#Call the plot function
memory.limit(40000)
sizeGrWindow(9,9)


pdf(paste0(out_dir, "11.Network heatmap plot.pdf"), width = 17, height = 10)
TOMplot(plotTOM, hclustGeneTree, dynamicColors, main = "Network heatmap plot")
dev.off()

#####################################################################################################
# Target gene identification
# Gene relationship to trait and important modules
#####################################################################################################
# To Quantify associations of individual gene with trait of interest (methane production)
# For this first define gene significance (GS) as the absolute value of the correlation between the gene and the trait.
# For each module, also define a quantitative measure of module membership (MM) as the correlation of the module eigngene and the gene expression profile. Allows us to quantify the similarity of all genes on the array to every module.
# Define variable methane production 
metpro <- as.data.frame(allTraits)
names(metpro) <- c("ADG_performance", "DMI_performance", "DMI_acc_BW", "RFI") # rename
# Calculate the correlations between modules
geneModuleMembership <- as.data.frame(WGCNA::cor(norm.counts, MEs, use = "p"))

# What's the correlation for the trait?
geneTraitSignificance <- as.data.frame(cor(norm.counts, metpro, use = "p"))
# What are the p-values for each correlation?
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples = nSamples))
names(geneTraitSignificance) <- paste("GS.", names(metpro), sep = "")
names(GSPvalue) <- paste("p.GS.", names(metpro), sep = "")

####################################################
# IIntramodular analysis: identify genes with high GS and MM
####################################################
modNames <- names(geneModuleMembership)
par(mfrow = c(2,3)) 
# Initialize for loop
# NEED: modNames
pdf(paste0(out_dir, "10.Gene_signicance_Module_membership.pdf"))
for (i in names(geneModuleMembership)) {
  # Pull out the module we're working on
  module <- i
  print(module) 
  # Find the index in the column of the dataframe 
  column <- match(module, modNames) #print(column)
  # Pull out the Gene Significance vs module membership of the module
  moduleGenes = moduleColors == module
  genenames = rownames(geneTraitSignificance)
  # Make the plot
  plot10 = verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), 
                              abs(geneTraitSignificance[moduleGenes, 1]),
                              xlab = paste("Module Membership in", module, "module"),
                              ylab = "Gene significance for Total VFA",
                              main = paste("Module membership vs. gene significnace \n"),
                              cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
}    
print(plot10)
dev.off()

# # Print the number of genes in each module
for (i in names(geneModuleMembership)) {
  # Pull out the module we're working on
  module <- i
  # Find the index in the column of the dataframe
  column <- match(module, modNames)
  # Pull out the Gene Significance vs module membership of the module
  moduleGenes = moduleColors == module
  genenames = rownames(geneTraitSignificance)
  print(paste("There are ", length(genenames[moduleGenes]), " genes in the ", module, " module.", sep = ""))

  # NOTE: This makes hidden variables with the gene names
  assign(paste(module, "_genes", sep = ""), genenames[moduleGenes])
}

# merge this important module information with gene annotation information and write out a file that summarizes the results.
# Combine pval, module membership, and gene significance into one dataframe
# Prepare pvalue df
# GSpval <- GSPvalue %>%
#   tibble::rownames_to_column(var = "gene")
# # Prepare module membership df
# gMM_df <- geneModuleMembership %>%
#   tibble::rownames_to_column(var = "gene") %>%
#   gather(key = "moduleColor", value = "moduleMemberCorr", -gene) 
# # Prepare gene significance df
# GS_metprod_df <- geneTraitSignificance %>%
#   data.frame() %>%
#   tibble::rownames_to_column(var = "gene")
# # Put everything together 
# allData_df <- gMM_df %>%
#   left_join(GS_metprod_df, by = "gene") %>%
#   left_join(GSpval, by = "gene") 

# Write a file 
#write.csv(allData_df, file = "allData_df_WGCNA.csv")


#################################################
# Hub gene in each module
#################################################
hub = chooseTopHubInEachModule(norm.counts, moduleColors)
hubgene = file.path(out_dir,"Tophub_genes_in_each_module.csv" )
write.csv(hub, file = hubgene, row.names = T) 

onehub = chooseOneHubInEachModule(norm.counts, moduleColors)
onehubgene = file.path(out_dir,"Onehub_genes_in_each_module.csv" )
write.csv(onehub, file = onehubgene, row.names = T) 


########################################################
# Extracting sig modules and genes
########################################################
# sigModules <- c("darkorange", "black", "lightgreen", "darkgrey", "blue")
# # Gene names (rownames of your expression data)
# allGenes <- colnames(norm.counts)
# # Genes in significant modules
# genesInModules <- allGenes[dynamic_mergedColors %in% sigModules]
# # View first few
# head(genesInModules)
# # Suppose 'moduleTraitCor' is your correlation matrix
# sigModules_ADG <- names(which(moduleTraitCor[,"ADG_performance"] > 0.4))
# sigModules_RFI    <- names(which(moduleTraitCor[,"RFI"] < -0.4))
# 
# genes_ADG <- allGenes[dynamic_mergedColors %in% sigModules_ADG]
# genes_RFI    <- allGenes[dynamic_mergedColors %in% sigModules_RFI]


#########################################################
# Writing all genes into a file
########################################################
# for multiple traits
# Traits you want to analyze
library(ggplot2)

traits <- c("ADG_performance", "DMI_performance", "DMI_acc_BW", "RFI")
geneInfo_list <- list()
# What are the p-values for each correlation?
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", MEs, sep="");
names(MMPvalue) = paste("p.MM", MEs, sep="");

for (trait in traits) {
  
  # Order modules by their significance for this trait
  modOrder <- order(-abs(cor(MEs, metpro[, trait, drop=FALSE], use = "p")))
  
  geneInfo_trait <- geneInfo
  
  for (mod in 1:ncol(geneModuleMembership)) {
    oldNames <- names(geneInfo_trait)
    geneInfo0 <- data.frame(geneInfo_trait,
                            geneModuleMembership[, modOrder[mod]],
                            MMPvalue[, modOrder[mod]])
    
    names(geneInfo0) <- c(oldNames,
                          paste("MM.", modNames[modOrder[mod]], sep=""),
                          paste("p.MM.", modNames[modOrder[mod]], sep=""))
    
    # Order genes by module color and GS for the current trait
    gs_col <- paste0("GS.", trait)
    geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0[[gs_col]]))
    geneInfo_trait <- geneInfo0[geneOrder, ]
  }
  
  # Save result for this trait
  geneInfo_list[[trait]] <- geneInfo_trait
  geneInfowrite <- file.path(out_dir, paste0("geneInfo_", trait, ".csv"))
  write.csv(geneInfo_list[[trait]], file = geneInfowrite, row.names = TRUE)



pdf(paste0(out_dir, "12.Module-Module_Heatmap_correlations.pdf"), width=8, height=10)
plotEigengeneNetworks(MEs,
                      setLabels = "Module-Module Heatmap",
                      plotDendrograms = TRUE,
                      plotHeatmaps = TRUE,
                      excludeGrey = TRUE,
                      signed = TRUE,
                      setMargins = TRUE,
                      # increase margins: marDendro = c(bottom, left, top, right)
                      marDendro = c(8, 4, 2, 2),
                      # marHeatmap = c(bottom, left, top, right)
                      marHeatmap = c(10, 10, 2, 2),
                      cex.adjacency = 0.7)  # or signed = FALSE for unsigned
dev.off()

