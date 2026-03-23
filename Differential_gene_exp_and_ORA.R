
library(dplyr)
library(DESeq2)

# Differential gene expression analyses
metaData <- read.csv("33lambs_all_phenotypes.csv")
countData<-read.csv("33lambs_featureCounts.csv",sep=",", header=T, check.names=F)

orig_names <- names(countData) # keep a back-up copy of the original names
geneID <- countData[,1] # Convert count data to a matrix of appropriate form that DEseq2 can read
countData <- as.matrix(countData[ , -1]) 
sampleIndex <- colnames(countData)
countData <- as.matrix(countData[,sampleIndex])
rownames(countData) <- geneID

rownames(metaData) <- metaData$LambID 
metaData$LambID  <- factor(metaData$LambID)
colnames(countData) == metaData$LambID 
# reorder counts columns based on row order of metadata (coldata)  
countData <- countData[, rownames(metaData)]
all(rownames(metaData) == colnames(countData))

# visualise the distribution of counts using a boxplot and density plot.
rafalib::mypar(1,2,mar=c(6,3,3,2))
boxplot(log2(as.matrix(countData)+1),ylab=expression('Log'[2]~'Read counts'),las=2,main="Raw data")
hist(log2(as.matrix(countData)+1),ylab="",las=2,main="Raw data")
par(mfrow=c(1,1))

# We can check if any samples need to be discarded based on the number of genes detected. 
# We create a barplot of genes detected across samples. I chose 7 because I have 7 replicates in my group with least samples.
barplot(colSums(countData>7),ylab="Number of detected genes",las=2)
abline(h=median(colSums(countData>7))) 

#Across genes
hist(rowSums(countData>7)) 

# Now filter these genes
# remove genes with low counts
keep_genes <- rowSums( countData > 5 ) >= 7
countData2 <- countData[keep_genes,]

# Distribution of the filtered counts. Compare this to the previous boxplot above.
boxplot(log2(as.matrix(countData2)+1),ylab=expression('Log'[2]~'Read counts'),las=2,main="Filtered data")

# In addition, compare the histogram of filtered counts below to the raw data above.

hist(rowSums(countData2>7))

all.equal(colnames(countData2),rownames(metaData))


d1 <- DESeqDataSetFromMatrix(countData=countData2,colData=metaData,design=~ADG)
d1 <- DESeq2::estimateSizeFactors(d1,type="ratio")
cd <- log2( counts(d1,normalized=TRUE) + 1 ) 

boxplot(cd,ylab=expression('Log'[2]~'Read counts'),las=2,main="DESeq2")

deseq2Data <- DESeq(d1)

# extract significant genes based on Log2FC and adj.pvalue thresholds of your wish
pval = 0.05
lfc = 2
results = resultsNames(deseq2Data)
upresultstable = matrix(nrow = length(results), ncol = 1, dimnames = list(results,"upDEGs"))
downresultstable = matrix(nrow = length(results), ncol = 1, dimnames = list(results,"downDEGs"))

for(i in 2:length(results)){
  
  res = results(deseq2Data, 
                name = results[i])# independent filtering occurs in this step to save you from multiple test correction on genes with no power
  resorder <- res[order(res$padj),]
  upDEGs = (length(na.omit(which(res$padj<pval & res$log2FoldChange > lfc))))
  downDEGs = (length(na.omit(which(res$padj<pval & res$log2FoldChange < -lfc))))
  resSig = subset(resorder, padj < pval & log2FoldChange > lfc | padj < pval & log2FoldChange < -lfc)
  #write.csv(resSig , file=paste0(results[i],"_updownDEGs.csv"), row.names = T)
  res_pval = subset(resorder, padj < pval )
  upresultstable[results[i],"upDEGs"] = upDEGs
  downresultstable[results[i],"downDEGs"] = downDEGs 
}

summary(res, alpha=0.05)

# run variance stabilizing transformation on counts
vsd <- vst(deseq2Data)
plotPCA(vsd, intgroup="ADG")

# Find top 500 genes with most variance
# calculate the variance for each gene
rv <- rowVars(assay(vsd))
# Top n genes by variance to keep.
ntop <- 500
# select the ntop genes by variance
select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
# perform a PCA on the data in assay(x) for the selected genes
pca <- prcomp(t(assay(vsd)[select,]))
# Loadings for the first two PCs.
loadings <- pca$rotation[, seq_len(2)]


# Another way to plot PCA
library(PCAtools)
p <- pca(assay(vsd), metadata = metaData, removeVar = 0.1)
screeplot(p, axisLabSize = 18, titleLabSize = 22)
biplot(p, showLoadings = TRUE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 3)


########################### Enrichment Analysis
#extract annotation DB of sheep
library(AnnotationHub)
library(clusterProfiler)
library(pathview)
ah <- AnnotationHub()
AnnotationHub::query(ah, c("Ovis", "aries"))
Oaries <- ah[["AH114633"]] # old one which worked with Biocversion 3.18
columns(Oaries)

res_original = res
resSig_original = resSig
# Taking the resSig variable from the above deseq2 as it has all the genes with corresponding pvalue and log2fc values
resSig$GeneID = rownames(resSig)
resSigEntrez <- AnnotationDbi::select(Oaries, keys =  rownames(resSig),
                                      columns = c('ENTREZID','GENENAME'), keytype = 'SYMBOL') # Oaries is the annotation db created from the geneset_enrichment_analysis.R code
# Now replace the NA values in entrezid column with the values from first column. 
resSigEntrez$ENTREZID <- ifelse(is.na(resSigEntrez$ENTREZID), resSigEntrez$SYMBOL, resSigEntrez$ENTREZID)
colnames(resSigEntrez) = c("GeneID", "ENTREZID")
resSig_with_entrez = merge(as.data.frame(resSig), resSigEntrez, by = "GeneID") #All entrezIDs have been retrieved for the 36 genes with p<0.1 and lfc=0 threshold
# Extract entrezIDs for the res variable
res_allgenes<-res
res_allgenes$GeneID<-rownames(res)
# removing "LOC" info
res_allgenes$GeneID <- stringr::str_remove(res_allgenes$GeneID, "LOC")
rownames(res_allgenes) = stringr::str_remove(rownames(res_allgenes), "LOC")

res_allgenesEntrez <- AnnotationDbi::select(Oaries, keys =  res_allgenes$GeneID,
                                            columns = c('ENTREZID','GENENAME'), keytype = 'SYMBOL')
# Now replace the NA values in entrezid column with the values from first column.
res_allgenesEntrez$ENTREZID <- ifelse(is.na(res_allgenesEntrez$ENTREZID), res_allgenesEntrez$SYMBOL, res_allgenesEntrez$ENTREZID)
colnames(res_allgenesEntrez) = c("GeneID", "ENTREZID")
res_allgenes_with_entrez = merge(as.data.frame(res_allgenes), res_allgenesEntrez, by = "GeneID") #All entrezIDs have been retrieved for the DEGs
dim(res_allgenes_with_entrez)
#we want the log2 fold change 
original_gene_list <- res_allgenes_with_entrez$log2FoldChange
# name the vector
names(original_gene_list) <- res_allgenes_with_entrez$ENTREZID
# omit any NA values 
gene_list<-na.omit(original_gene_list)
# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
# Exctract significant results (padj < 0.05)
sig_genes_df = subset(res_allgenes_with_entrez, padj < 0.05)
# From significant results, we want to filter on log2fold change
genes <- sig_genes_df$log2FoldChange
# Name the vector
names(genes) <- sig_genes_df$ENTREZID
# omit NA values
genes <- na.omit(genes)
# filter on min log2fold change (log2FoldChange > 2)
genes <- names(genes)[abs(genes) > 2]
length(genes)
go_enrich <- enrichGO(gene = genes,
                      universe = names(gene_list),
                      OrgDb = Oaries, 
                      keyType = 'ENTREZID',
                      readable = T,
                      ont = "ALL",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)

as.data.frame(go_enrich)
Ora_go_df = file.path("Raw_ADG_filtering_1_ORA_GO_enrichments.csv")
write.csv(as.data.frame(go_enrich),Ora_go_df, row.names=T)

# Extract significant results from df2
kegg_sig_genes_df = subset(resSig_with_entrez, padj < 0.05)
# From significant results, we want to filter on log2fold change
kegg_genes <- kegg_sig_genes_df$log2FoldChange
# Name the vector with the CONVERTED ID!
names(kegg_genes) <- kegg_sig_genes_df$ENTREZID
# omit NA values
kegg_genes <- na.omit(kegg_genes)
# filter on log2fold change (PARAMETER)
kegg_genes <- names(kegg_genes)[abs(kegg_genes) > 2]
Sys.setenv(R_LIBCURL_SSL_REVOKE_BEST_EFFORT=TRUE)
kegg_organism = "oas"
kk <- enrichKEGG(gene=kegg_genes, universe=names(gene_list),organism=kegg_organism, pvalueCutoff = 0.05)
as.data.frame(kk)
Ora_kegg_df = file.path("Raw_ADG_filtering_1_ORA_KEGG_enrichments.csv")
write.csv(as.data.frame(kk),Ora_kegg_df, row.names=T)

#From	To	Species	David Gene Name
# 100913157	100913157	Ovis aries	integrin linked kinase(ILK)
# 101121460	101121460	Ovis aries	acyl-CoA oxidase 2(ACOX2)
# 101115115	101115115	Ovis aries	platelet glycoprotein 4(LOC101115115)

# Validate DEGs from DESEq2 using Linear Model. If any genes remain significant both in DESEq2 and LM, then those genes have strong assoiations with DMI.

vsd <- vst(deseq2Data, blind = FALSE)
expr_data <- assay(vsd) # extract transformed expression

# Run LM for each gene
lm_results <- apply(expr_data, 1, function(gene_expr){
  model <- lm(gene_expr~ADG, data = metaData)
  summary(model)$coefficients["ADG", c("Estimate", "Pr(>|t|)")]
})

lm_results <- as.data.frame(t(lm_results))
colnames(lm_results) <- c("Estimate", "P.Value")
# Adjust p-values for multiple testing
lm_results$adj.P.Val <- p.adjust(lm_results$P.Value, method ="fdr")

# compare LM with DESeq2
merged_res <- merge(as.data.frame(resSig), lm_results, by ="row.names", all.x = TRUE)



