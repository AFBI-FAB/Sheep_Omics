library(dplyr)
library(plyr)
library(tidyverse)
library(stringr)
library(clusterProfiler)
library(AnnotationDbi)
library(AnnotationHub)
library(enrichplot)
library(KEGGREST)

#extract annotation DB of sheep
ah <- AnnotationHub()
AnnotationHub::query(ah, c("Ovis", "aries"))
Oaries <- ah[["AH114633"]] # old one which worked with Biocversion 3.18
#Oaries <- ah[["AH116933"]]
columns(Oaries)

# ORA
# 1. method 2 to retrieve entrez ids
# Temporarily store the original res into another variable (maybe for later use)
res_original = res
resSig_original = resSig
# check if the number of genes are correct or not
dim(res)
dim(resSig)

# Taking the resSig variable from the above deseq2_MO as it has all the genes with corresponding pvalue and log2fc values
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

# we want the log2 fold change
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
genes <- names(genes)[abs(genes) > 0.584]
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
Ora_go_df = file.path(out_dir, "ORA_GO_enrichments.csv")
write.csv(as.data.frame(go_enrich),Ora_go_df, row.names=T)

# KEGG pathways
# Exctract significant results from df2
kegg_sig_genes_df = subset(resSig_with_entrez, padj < 0.05)
# From significant results, we want to filter on log2fold change
kegg_genes <- kegg_sig_genes_df$log2FoldChange
# Name the vector with the CONVERTED ID!
names(kegg_genes) <- kegg_sig_genes_df$ENTREZID
# omit NA values
kegg_genes <- na.omit(kegg_genes)
# filter on log2fold change (PARAMETER)
kegg_genes <- names(kegg_genes)[abs(kegg_genes) > 0.584]
Sys.setenv(R_LIBCURL_SSL_REVOKE_BEST_EFFORT=TRUE)
kegg_organism = "oas"
kk <- enrichKEGG(gene=kegg_genes, universe=names(gene_list),organism=kegg_organism, pvalueCutoff = 0.05)
as.data.frame(kk)
Ora_kegg_df = file.path(out_dir, "ORA_KEGG_enrichments.csv")
write.csv(as.data.frame(kk),Ora_kegg_df, row.names=T)
