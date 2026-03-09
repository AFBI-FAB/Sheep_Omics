library(dplyr)
library(plyr)
library(tidyverse)
library(ggpubr)
library(EnhancedVolcano)
library(ggplot2)
library(CorLevelPlot)
library(DESeq2)
library(stringr)
library(clusterProfiler)
library(pathview)
library(stringr)
library(AnnotationDbi)
library(AnnotationHub)
library(ggridges)
library(enrichplot)
library(KEGGREST)
library(EnrichmentBrowser)
library(GSVA)
library(pheatmap)
library(car)
library(plotly)
library(rstatix)

countData<-read.csv("20Samples_afterTrimming_6976removed_HISAT2.csv",sep=",", header=T, check.names=F)
orig_names <- names(countData) # keep a back-up copy of the original names
geneID <- countData[,1] # Convert count data to a matrix of appropriate form that DEseq2 can read
countData <- as.matrix(countData[ , -1])
sampleIndex <- colnames(countData)
countData <- as.matrix(countData[,sampleIndex])
rownames(countData) <- geneID
countData2<-countData
# Calculate total number of columns
total_columns <- ncol(countData2)
zero_counts <- rowSums(countData2 == 0)
length(zero_counts)
table(zero_counts)
countData2 <- as.data.frame(countData2) %>% dplyr::filter(zero_counts < 0.7 * total_columns)
mycounts <- countData2

metaData <-read.csv("20Samples_metadata_6976_removed.csv",sep=",",header=T)
metaData %>% arrange(ID)
rownames(metaData) <- metaData$ID
metaData$ID <- factor(metaData$ID)
colnames(countData2) == metaData$ID
# Select the columns from countData2 that are in bno_values
countData2 <- countData2[, colnames(countData2) %in% metaData$ID]

# Function to compute SEM
sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
# Calculate mean and SEM for all numeric variables in 'ad'
summary_stats <- metaData %>%
  summarise(across(where(is.numeric),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        sem  = ~sem(.x))
                   )
            )

summary_stats


#Summary by treatment
sem <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))
# Mean and SEM for all numeric variables grouped by Treatment
summary_by_treatment <- metaData %>%
  group_by(Treatment) %>%
  summarise(across(
    where(is.numeric),
    list(
      mean = ~ mean(.x, na.rm = TRUE),
      sem  = ~ sem(.x)
    ),
    .names = "{.col}_{.fn}"
  ))

summary_by_treatment

# check assumptions
shapiro.test(metaData$CH4production)
shapiro.test(metaData$Ave_MO)

tapply(metaData$CH4production, metaData$Treatment, mean)
tapply(metaData$CH4production, metaData$Treatment, median)
tapply(metaData$CH4production, metaData$Treatment, skewness)
tapply(metaData$CH4production, metaData$Treatment, kurtosis)

# # check variance between groups. If p-value greater than 0.05, vars are same
bartlett.test(metaData$CH4production ~ metaData$Treatment)

# # Run anova
ANO <- aov(metaData$CH4production ~ metaData$Treatment)
summary(ANO)

# normality assumption by groups.
# Computing Shapiro-Wilk test for each group level. If the data is normally distributed, the p-value should be greater than 0.05.
metaData %>%
  group_by(Treatment) %>%
  shapiro_test(CH4production)
ggqqplot(metaData, "CH4production", facet.by = "Treatment")
res.aov <- metaData %>% anova_test(CH4production ~ Treatment)
res.aov

# Pairwise comparisons (exploratory)
pwc <- metaData %>% tukey_hsd(CH4production ~ Treatment)
pwc

# Visualization: box plots with p-values
pwc <- pwc %>% add_xy_position(x = "Treatment")
ggboxplot(metaData, x = "Treatment", y = "CH4production", color="Treatment", palette="jco") +
  stat_pvalue_manual(pwc, hide.ns = TRUE) +
  labs(
    subtitle = get_test_label(res.aov, detailed = TRUE),
    caption = get_pwc_label(pwc)
    )


# Differential analyses
metaData$scaled_CH4production=scale(metaData$CH4production, center=TRUE)
metaData$scaled_Ave_MO = scale(metaData$Ave_MO, center = TRUE)
metaData$MO_acc_for_BW = metaData$Ave_MO/metaData$Final_BW

# For methane production
deseq2Data_meth <- DESeqDataSetFromMatrix(countData=countData2, colData=metaData, design= ~ scaled_CH4production)
deseq2Data_meth <- estimateSizeFactors(deseq2Data_meth)
deseq2Data_meth <- DESeq(deseq2Data_meth)
resultsNames(deseq2Data_meth)
results <- results(deseq2Data_meth, name="scaled_CH4production")

#extract sig genes
### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 0.584
res_tableOE_tb <- results %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
sigOE <- res_tableOE_tb %>%
        filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

# For microalgae oil intake 
deseq2Data_MO <- DESeqDataSetFromMatrix(countData=countData2, colData=metaData, design= ~MO_acc_for_BW)
deseq2Data_MO <- estimateSizeFactors(deseq2Data_MO)
deseq2Data_MO <- DESeq(deseq2Data_MO)
results <- results(deseq2Data_MO, name="MO_acc_for_BW")
#extract sig genes
### Set thresholds
padj.cutoff <- 0.05
lfc.cutoff <- 2
res_tableOE_tb <- results %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
sigOE <- res_tableOE_tb %>%
        filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

# Another for loop method for getting significant genes
pval = 0.05
lfc = 0.584
results = resultsNames(deseq2Data_MO)
upresultstable = matrix(nrow = length(results), ncol = 1, dimnames = list(results,"upDEGs"))
downresultstable = matrix(nrow = length(results), ncol = 1, dimnames = list(results,"downDEGs"))

for(i in 2:length(results)){

  res = results(deseq2Data_MO,
                name = results[i])# independent filtering occurs in this step to save you from multiple test correction on genes with no power
  resorder <- res[order(res$padj),]
  upDEGs = (length(na.omit(which(res$padj<pval & res$log2FoldChange > lfc))))
  downDEGs = (length(na.omit(which(res$padj<pval & res$log2FoldChange < -lfc))))
  resSig = subset(resorder, padj < pval & log2FoldChange > lfc | padj < pval & log2FoldChange < -lfc)
  write.csv(resSig , file=paste0(out_dir,results[i],".0.05P.0.584LFC.updownDEGs.csv"), row.names = T)
  res_pval = subset(resorder, padj < pval )
  upresultstable[results[i],"upDEGs"] = upDEGs
  downresultstable[results[i],"downDEGs"] = downDEGs
}

# Visualisation
# 1. MA plot
# shrinkage of LFC is importnt for visualization
resultsNames(deseq2Data_MO)
resLFC <- lfcShrink(deseq2Data_MO, coef="scaled_CH4production", type="apeglm")
maplot = file.path(out_dir, "MA-Plot.pdf")
p2 = plotMA(resLFC,main = "MA plot", alpha=0.05,colNonSig = "black", colSig = "red")
ggsave(filename = maplot, plot = p2)

# 2. Volcano plot
# change the pvalue and fccutoffs as you wish (lines 124 and 125)
volplot = file.path(out_dir, "Volcano-Plot.pdf")
p3 = EnhancedVolcano(res,
                lab=rownames(res),
                x='log2FoldChange',
                y='padj',
                ylim=c(0,3),
                xlim=c(-4,3),
                title='Volcano plot',
                caption = 'Log2FoldChange cutoff, 0.584; p-value cutoff, 0.05',
                pCutoff=0.05,
                FCcutoff=0.584,
                pointSize = 4.0,
                labSize = 5.0,
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75)
ggsave(filename = volplot, plot = p3)
