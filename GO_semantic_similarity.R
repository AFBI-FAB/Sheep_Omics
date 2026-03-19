library(AnnotationDbi)
library(clusterProfiler)
library(GOSemSim)

hsGOCC <- godata(annoDb = 'org.Hs.eg.db', ont="CC")
hsGOBP <- godata(annoDb = 'org.Hs.eg.db', ont="BP")
hsGOMF <- godata(annoDb = 'org.Hs.eg.db', ont="MF")

# sigOE from previous DGE_analysis.R
hits_entrez <- bitr(sigOE$gene,
                    fromType="SYMBOL", toType="ENTREZID", OrgDb=Oaries)
hits_entrez <- hits_entrez[!is.na(hits_entrez$ENTREZID), ]
hits_entrez <- unique(hits_entrez$ENTREZID)
hits_entrez

# gene semantic similarity
#NME4:4833
#MARCHF3:5365
#PLXNB3: 115123
genes_core <- c("NME4", "MARCHF3", "PLXNB3")
GOSemSim::geneSim("4833","115123", semData=hsGOCC, measure="Wang", combine="BMA")
mgeneSim(genes=c("4833", "5365", "115123"),
         semData=hsGOCC, measure="Wang",verbose=FALSE)

## Heatmap
# Map symbols -> Entrez using org.Hs.eg.db (one-to-one; first hit if multiple)
entrez_map <- AnnotationDbi::mapIds(
  x         = org.Hs.eg.db,
  keys      = genes_core,
  keytype   = "SYMBOL",
  column    = "ENTREZID",
  multiVals = "first"
)

#  Keep only those that mapped (the “mapped ones”)
entrez_vec <- entrez_map[!is.na(entrez_map)] #16 mapped

#  Inspect any symbols that didn’t map
unmapped <- names(entrez_map)[is.na(entrez_map)]

# Print the final named vector you can drop straight into GOSemSim / enrichment
entrez_vec

sim_mat <- mgeneSim(
  unname(entrez_vec),
  semData = hsGOCC,
  measure = "Wang",        # best for BP
  combine = "BMA"          # Best-Match Average
)


keep <- which(rowSums(!is.na(sim_mat)) > 0 & colSums(!is.na(sim_mat)) > 0)
sim_mat_filt <- sim_mat[keep, keep]
sim_mat_filt

# Convert similarity to distance
#dist_mat <- as.dist(1 - sim_mat_filt)
# Build label vectors aligned to the matrix order
symbol_map <- AnnotationDbi::mapIds(
  x         = org.Hs.eg.db,
  keys      = rownames(sim_mat_filt),  # Entrez IDs
  keytype   = "ENTREZID",
  column    = "SYMBOL",
  multiVals = "first"
)

labels <- ifelse(is.na(symbol_map), rownames(sim_mat_filt), symbol_map)
labels_unique <- make.unique(labels, sep = "_")

pheatmap(
  sim_mat_filt,
  clustering_distance_rows = as.dist(1 - sim_mat_filt),
  clustering_distance_cols = as.dist(1 - sim_mat_filt),
  clustering_method = "average",
  color = colorRampPalette(c("#313695","#74add1","#ffffbf",
                             "#f46d43","#a50026"))(100),
  main = "GO Semantic Similarity (P.adjusted < 0.05, absolute LFC > 0.584)",
  fontsize = 10,
  border_color = NA,
  labels_row = labels_unique,
  labels_col = labels_unique
)
