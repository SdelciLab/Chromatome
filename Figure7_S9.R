#title: "Differential expression analysis IMPDH2"
#author: "Alisa Schmidt"
#date: "22/09/2025"

#libraries
library(DESeq2)
library(ggplot2)
library(dbplyr)
library(tidyverse)
library(apeglm)
library(EnsDb.Hsapiens.v86)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(fgsea)
library(msigdbr)
library(GSEABase)
library(Biobase)
library(GSVA)
library(rrvgo)

#load dds file -------------------------------------------
#setup: MDA-MB-231 cells in 6 conditions
#conditions: WT, K0WT, K0N0G, K0SIG, K0NES, K0NLS
load("/Users/aschmidt/Desktop/Chromatome_paper/DeSeq R_old/DDS_file.RData")
resultsNames(dds)

##Differential expression ----
#change reference to KOWT
dds$Condition <- relevel(dds$Condition, ref = "K0WT")
# drop samples belonging to WT
dds <- dds[, !(dds$Condition %in% c("WT"))]
# drop unused factor levels
dds$Condition <- droplevels(dds$Condition)

# re-run DESeq to fit the model without those conditions
dds <- DESeq(dds)
resultsNames(dds)

##Plot PCA plot ---------
  #Supplementary Figure 9A ----
rld <- rlog(dds, blind = TRUE)
rld_counts <- assay(rld)
pca <- prcomp(t(rld_counts)) #transpose the data
eigs <- pca$sdev^2
percentVar <-  round(100 * eigs / sum(eigs))
condtions <- as.data.frame(dds@colData@listData[["Condition"]],)
colnames(condtions) <- "condition"
pcaData <- cbind(condtions, pca$x)
pcaData$sample <- rownames(pcaData)

ggplot(pcaData, aes(x=PC1, y=PC2, color=condition, label=sample)) + geom_point(size =4) + 
  geom_text_repel(data=pcaData,aes(label=sample),size=4, max.overlaps =100) +
  scale_color_manual(values = c("#990099", "#303030", "#006000", "#ffd320", "blue")) + xlab(paste0("Principle Component 1: ", percentVar[1], "%")) +
  ylab(paste0("Principle Component 2: ", percentVar[2], "%")) + theme_bw()

##plot sample distances ------------
  #Supplementary Figure 9B
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rld$sample
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


#GO term analysis with PCs ------------------------------------------
#FGSEA 
#Ranked stats with Ensembl IDs
stats_pc1 <- pca$rotation[, "PC1"]
names(stats_pc1) <- rownames(pca$rotation)   
stats_pc1 <- sort(stats_pc1, decreasing = TRUE)

stats_pc2 <- pca$rotation[, "PC2"]
names(stats_pc2) <- rownames(pca$rotation)
stats_pc2 <- sort(stats_pc2, decreasing = TRUE)

#Get GO:BP gene sets 
msig_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")

# Make list: pathway -> vector of Ensembl IDs
go_bp_ens <- msig_bp %>%
  split(.$gs_name) %>%
  lapply(function(d) unique(d$ensembl_gene))

##Run fgsea -----------
set.seed(1)
BiocParallel::register(BiocParallel::SnowParam(workers = 4, type = "SOCK"))
fg_pc1 <- fgsea(pathways = go_bp_ens,
                stats    = stats_pc1,
                minSize  = 10,
                maxSize  = 500) %>%
  arrange(padj, desc(NES))

fg_pc2 <- fgsea(pathways = go_bp_ens,
                stats    = stats_pc2,
                minSize  = 10,
                maxSize  = 500) %>%
  arrange(padj, desc(NES))

##Plot GOs as Dotplot -----
  #Supplementary Figure 9D ----
plot_fgsea_pc1_bubble <- function(fg, padj_cutoff = 0.01, k = 15,
                                  use_leading_edge = TRUE, title = "fgsea (GO:BP)") {
  df <- fg %>%
    dplyr::filter(padj < padj_cutoff) %>%
    mutate(
      nLE = lengths(leadingEdge),
      gene_count = if (use_leading_edge) nLE else size,
      Direction = ifelse(NES >= 0, "Positive NES", "Negative NES")
    ) %>%
    arrange(desc(abs(NES))) %>%
    slice_head(n = k) %>%
    mutate(pathway = factor(pathway, levels = rev(pathway)))
  
  ggplot(df, aes(x = NES, y = pathway)) +
    geom_point(aes(size = gene_count, color = padj)) +
    scale_size_area(max_size = 10,
                    name = "Contributing genes")+
    scale_color_viridis_c(name = "padj") +
    labs(x = "Normalized Enrichment Score (NES)",
         y = NULL,
         title = title,
         subtitle = paste0("padj < ", padj_cutoff, ", top ", k, " by |NES|")) +
    theme_bw() +
    theme(legend.position = "right")
}


plot_fgsea_pc1_bubble(fg_pc1, padj_cutoff = 0.01, k = 16, title = "PC1 – GO:BP enrichment", use_leading_edge = TRUE)
plot_fgsea_pc1_bubble(fg_pc2, padj_cutoff = 0.01, k = 16, use_leading_edge = TRUE,  title = "PC2 – GO:BP enrichment")

fg_pc2_plot <- fg_pc2[1:8,]
fg_pc2_plot$leadingEdge <- sapply(fg_pc2_plot$leadingEdge, function(x) paste(x, collapse = ", "))
#write.csv(fg_pc2_plot, "SFig9D.csv")


#rrvgo ----------------
#Cluster similar terms for PC1
  #Supplementary Figure 9C ----

id_map <- msig_bp %>%
  distinct(gs_name,
           go_id = dplyr::coalesce(gs_exact_source,
                                   sub(".*/(GO:\\d+).*", "\\1", gs_url))) %>%
  dplyr::filter(grepl("^GO:\\d{7}$", go_id))

# 2) Keep only significant PC1 terms and get GO IDs
fg_pc1_sif <- fg_pc1[fg_pc1$padj < 0.01,]
sig <- fg_pc1_sif %>%
  transmute(gs_name = pathway, padj) %>%
  inner_join(id_map, by = "gs_name") %>%
  arrange(padj) %>%
  distinct(go_id, .keep_all = TRUE)

stopifnot(nrow(sig) > 0)
terms  <- sig$go_id
scores <- setNames(-log10(sig$padj), terms)

# 3) Semantic similarity (BP) and matrix (NA→0)
semData <- tryCatch(
  GOSemSim::godata('org.Hs.eg.db', ont='BP'),
  warning = function(w) GOSemSim::godata(annoDb='org.Hs.eg.db', ont='BP')
)
sim_mat <- GOSemSim::mgoSim(terms, terms, semData=semData, measure="Wang", combine=NULL)
sim_mat[is.na(sim_mat)] <- 0
diag(sim_mat) <- 1


reducedTerms <- reduceSimMatrix(sim_mat,
                                scores,
                                threshold=0.6,
                                orgdb="org.Hs.eg.db")



# rank clusters by size, then total score, keep top 10
rank_tbl <- reducedTerms %>%
  dplyr::filter(!is.na(parentTerm)) %>%
  group_by(parentTerm) %>%
  summarise(size = n(), total_score = sum(score), .groups = "drop") %>%
  # arrange(desc(size), desc(total_score))
  arrange(desc(total_score))
keep_parents <- head(rank_tbl$parentTerm, 9)
reduced_top10 <- reducedTerms %>% dplyr::filter(parentTerm %in% keep_parents)

#write.csv(reduced_top10, "SFig9C.csv")

# Make sure the GO IDs match the matrix row/col names
terms_keep <- intersect(reduced_top10$go, rownames(sim_mat))

# Subset the similarity matrix symmetrically
sim_mat_top10 <- sim_mat[terms_keep, terms_keep, drop = FALSE]


scatterPlot(sim_mat_top10, reduced_top10) + theme_bw()
p <- scatterPlot(sim_mat_top10, reduced_top10) + theme_bw()
p
p + theme(legend.position = "right") +
  guides(size = guide_legend(title = "Semantic similarity"))

#Plot Differential expression volcanos -------------------
#K0NESvsK0WT ---- 
  #Figure 9C
res <- results(dds, contrast=c("Condition","K0NES","K0WT"))
coefName <- "Condition_K0NES_vs_K0WT"
resLFC <- lfcShrink(dds, coef=coefName, type="apeglm", res=res)

# Convert Ensembl IDs to gene symbols
ensembl_ids <- rownames(resLFC)
hgnc_genenames <- ensembldb::select(EnsDb.Hsapiens.v86, keys=ensembl_ids, keytype="GENEID", columns=c("SYMBOL", "GENEID"))
resLFC$symbol <- hgnc_genenames[match(ensembl_ids, hgnc_genenames$GENEID), "SYMBOL"]

highlight_genes <- c("DCAF8L1", "CCL2", "OLFML1", "IL13RA2", "IMPDH2")

volcano_df_KONES <- as.data.frame(resLFC) %>%
  dplyr::filter(!is.na(padj)) %>%
  mutate(neglog10 = -log10(padj),
         sig = case_when(
           padj < 0.01 & log2FoldChange > 1  ~ "Up",
           padj < 0.01 & log2FoldChange < -1 ~ "Down",
           TRUE                              ~ "NotSig"
         ))

xmax <- max(abs(volcano_df_KONES$log2FoldChange), na.rm = TRUE)

ggplot(volcano_df_KONES, aes(
  x = log2FoldChange,
  y = neglog10,
  color = sig
)) +
  geom_point(aes(size = neglog10 * abs(log2FoldChange), alpha=0.5)) +
  scale_color_manual(values = c(Up = "red", Down = "blue", NotSig = "grey")) +
  scale_size_continuous(range = c(1, 17)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  labs(
    x = "log2 fold change",
    y = expression(-log[10]("adj p-value")),
    title = "K0NES vs K0WT: Volcano plot",
    legend = F
  ) + 
  coord_cartesian(xlim = c(-xmax, xmax)) +
  theme_bw() +
  ggrepel::geom_text_repel(
    data = subset(volcano_df_KONES, symbol %in% highlight_genes),
    aes(label = symbol),
    size = 6,              
    max.overlaps = Inf,
    box.padding = 1,
    point.padding = 0.6,
    force = 5,
    color = "black"
  ) +theme(legend.position = "none")

#ggsave("volcanoKONES_KOWT_genes_big.pdf",width = 9, height = 9, , limitsize = FALSE)

#K0NLS vs K0WT ----
  #Figure 7B ----
res <- results(dds, contrast=c("Condition","K0NLS","K0WT"))
coefName <- "Condition_K0NLS_vs_K0WT"
resLFC <- lfcShrink(dds, coef=coefName, type="apeglm", res=res)

# Convert Ensembl IDs to gene symbols
ensembl_ids <- rownames(resLFC)
hgnc_genenames <- ensembldb::select(EnsDb.Hsapiens.v86, keys=ensembl_ids, keytype="GENEID", columns=c("SYMBOL", "GENEID"))
resLFC$symbol <- hgnc_genenames[match(ensembl_ids, hgnc_genenames$GENEID), "SYMBOL"]


volcano_df_KONLS <- as.data.frame(resLFC) %>%
  dplyr::filter(!is.na(padj)) %>%
  mutate(neglog10 = -log10(padj),
         sig = case_when(
           padj < 0.01 & log2FoldChange > 1  ~ "Up",
           padj < 0.01 & log2FoldChange < -1 ~ "Down",
           TRUE                              ~ "NotSig"
         ))
xmax <- max(abs(volcano_df_KONLS$log2FoldChange), na.rm = TRUE)
# Volcano plot
ggplot(volcano_df_KONLS, aes(
  x = log2FoldChange,
  y = neglog10,
  color = sig
)) +
  geom_point(aes(size = neglog10 * abs(log2FoldChange), alpha=0.5)) +
  scale_color_manual(values = c(Up = "red", Down = "blue", NotSig = "grey")) +
  scale_size_continuous(range = c(1, 17)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  labs(
    x = "log2 fold change",
    y = expression(-log[10]("adj p-value")),
    title = "K0NLS vs K0WT",
    legend = F
  ) + 
  coord_cartesian(xlim = c(-xmax, xmax)) +
  theme_bw() +
  ggrepel::geom_text_repel(
    data = subset(volcano_df_KONLS, symbol %in% highlight_genes),
    aes(label = symbol),
    size = 6,              
    max.overlaps = Inf,
    box.padding = 1,
    point.padding = 0.6,
    force = 5,
    color = "black"
  ) +theme(legend.position = "none")

#ggsave("volcanoKONLS_KOWT_genes_big.pdf",width = 9, height = 9, limitsize = FALSE)

#K0SIGvsK0WT ----
  #Supplementary Figure 9E ----
res <- results(dds, contrast=c("Condition","K0SIG","K0WT"))
coefName <- "Condition_K0SIG_vs_K0WT"
resLFC <- lfcShrink(dds, coef=coefName, type="apeglm", res=res)

# Convert Ensembl IDs to gene symbols
ensembl_ids <- rownames(resLFC)
hgnc_genenames <- ensembldb::select(EnsDb.Hsapiens.v86, keys=ensembl_ids, keytype="GENEID", columns=c("SYMBOL", "GENEID"))
resLFC$symbol <- hgnc_genenames[match(ensembl_ids, hgnc_genenames$GENEID), "SYMBOL"]


volcano_df_K0SIG <- as.data.frame(resLFC) %>%
  dplyr::filter(!is.na(padj)) %>%
  mutate(neglog10 = -log10(padj),
         sig = case_when(
           padj < 0.01 & log2FoldChange > 1  ~ "Up",
           padj < 0.01 & log2FoldChange < -1 ~ "Down",
           TRUE                              ~ "NotSig"
         ))
xmax <- max(abs(volcano_df_K0SIG$log2FoldChange), na.rm = TRUE)

# Volcano plot
ggplot(volcano_df_K0SIG, aes(
  x = log2FoldChange,
  y = neglog10,
  color = sig
)) +
  geom_point(aes(size = neglog10 * abs(log2FoldChange), alpha=0.5)) +
  scale_color_manual(values = c(Up = "red", Down = "blue", NotSig = "grey")) +
  scale_size_continuous(range = c(1, 17)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  labs(
    x = "log2 fold change",
    y = expression(-log[10]("adj p-value")),
    title = "KOSIG vs K0WT",
    legend = F
  ) + 
  coord_cartesian(xlim = c(-xmax, xmax)) +
  theme_bw() +
  ggrepel::geom_text_repel(
    data = subset(volcano_df_K0SIG, symbol %in% highlight_genes),
    aes(label = symbol),
    size = 6,              
    max.overlaps = Inf,
    box.padding = 1,
    point.padding = 0.6,
    force = 5,
    color = "black"
  ) +theme(legend.position = "none")

#ggsave("volcanoK0SIG_KOWT_genes_big.pdf",width = 9, height = 9, limitsize = FALSE)

#K0NOGvsK0WT ----
  #Supplementary Figure 9F ----
res <- results(dds, contrast=c("Condition","K0N0G","K0WT"))
coefName <- "Condition_K0N0G_vs_K0WT"
resLFC <- lfcShrink(dds, coef=coefName, type="apeglm", res=res)

# Convert Ensembl IDs to gene symbols
ensembl_ids <- rownames(resLFC)
hgnc_genenames <- ensembldb::select(EnsDb.Hsapiens.v86, keys=ensembl_ids, keytype="GENEID", columns=c("SYMBOL", "GENEID"))
resLFC$symbol <- hgnc_genenames[match(ensembl_ids, hgnc_genenames$GENEID), "SYMBOL"]


volcano_df_K0N0G <- as.data.frame(resLFC) %>%
  dplyr::filter(!is.na(padj)) %>%
  mutate(neglog10 = -log10(padj),
         sig = case_when(
           padj < 0.01 & log2FoldChange > 1  ~ "Up",
           padj < 0.01 & log2FoldChange < -1 ~ "Down",
           TRUE                              ~ "NotSig"
         ))
xmax <- max(abs(volcano_df_K0N0G$log2FoldChange), na.rm = TRUE)

#volcano
ggplot(volcano_df_K0N0G, aes(
  x = log2FoldChange,
  y = neglog10,
  color = sig
)) +
  geom_point(aes(size = neglog10 * abs(log2FoldChange), alpha=0.5)) +
  scale_color_manual(values = c(Up = "red", Down = "blue", NotSig = "grey")) +
  scale_size_continuous(range = c(1, 17)) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +
  labs(
    x = "log2 fold change",
    y = expression(-log[10]("adj p-value")),
    title = "KONOG vs K0WT",
    legend = F
  ) + 
  coord_cartesian(xlim = c(-xmax, xmax)) +
  theme_bw() +
  ggrepel::geom_text_repel(
    data = subset(volcano_df_K0N0G, symbol %in% highlight_genes),
    aes(label = symbol),
    size = 6,              
    max.overlaps = Inf,
    box.padding = 1,
    point.padding = 0.6,
    force = 5,
    color = "black"
  ) +theme(legend.position = "none")

#ggsave("volcanoK)NOG_KOWT_genes_big.pdf",width = 9, height = 9, limitsize = FALSE)

#check guanosine dependence ---------------------
#Supplementary Figure 9G ----
# 1) Keep symbol + sig, one row per gene per condition
df_sig  <- volcano_df_K0SIG %>% dplyr::select(symbol, sig) %>% distinct(symbol, .keep_all = TRUE)
df_n0g  <- volcano_df_K0N0G %>% dplyr::select(symbol, sig) %>% distinct(symbol, .keep_all = TRUE)

# 2) Merge and classify
KO_df <- full_join(df_sig, df_n0g, by = "symbol", suffix = c("_K0SIG", "_K0N0G")) %>%
  mutate(
    sig_K0SIG = coalesce(sig_K0SIG, "NotSig"),
    sig_K0N0G = coalesce(sig_K0N0G, "NotSig"),
    category = case_when(
      sig_K0SIG == "Up"   & sig_K0N0G == "Up"    ~ "Same Up",
      sig_K0SIG == "Down" & sig_K0N0G == "Down"  ~ "Same Down",
      sig_K0SIG == "Up"   & sig_K0N0G == "NotSig"~ "K0SIG Up only",
      sig_K0SIG == "Down" & sig_K0N0G == "NotSig"~ "K0SIG Down only",
      sig_K0N0G == "Up"   & sig_K0SIG == "NotSig"~ "K0N0G Up only",
      sig_K0N0G == "Down" & sig_K0SIG == "NotSig"~ "K0N0G Down only",
      sig_K0N0G == "Up"   & sig_K0SIG == "Down"~ "opposite",
      sig_K0N0G == "Down" & sig_K0SIG == "Up"~ "opposite",
      TRUE ~ "Other"  # (e.g., Opposite or Both NotSig)
    )
  )
#pie chart
counts <- KO_df %>%
  dplyr::filter(category != "Other") %>%
  count(category) %>%
  arrange(desc(n))

ggplot(counts, aes(x = "", y = n, fill = category)) +
  geom_col(width = 1) +
  coord_polar(theta = "y") +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 4) +
  theme_void() +
  labs(title = "Overlap of differential genes: K0SIG vs K0N0G") +
  guides(fill = guide_legend(title = NULL))




##Check SIGNATURES--------------------------------------
  #Figure 7H ----
#Extract Signatures of interest with ENSEMBLE IDs

# BRCA/HR-focused sets from  GO 
go_hr <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>%
  dplyr::filter(grepl("HOMOLOGOUS_RECOMBINATION|DOUBLE[-_ ]STRAND_BREAK|BRCA|REPAIR", gs_name, ignore.case = TRUE))

#combine
brca_like_sets2 <- bind_rows(go_hr) %>%
  distinct(gs_name, ensembl_gene)

brca_like_sets_filtered <- brca_like_sets2 %>%
  # remove protein/plasma membrane/mitochondrial terms
  dplyr::filter(!grepl("PROTEIN_REPAIR|PLASMA_MEMBRANE_REPAIR|MITOCHONDRIAL", gs_name, ignore.case = TRUE)) %>%
  # remove regulation (positive/negative/unspecified) terms
  dplyr::filter(!grepl("REGULATION", gs_name, ignore.case = TRUE)) %>%
  # remove disease terms
  dplyr::filter(!grepl("DISEASE", gs_name, ignore.case = TRUE)) %>%
  # keep only DDR-related terms
  dplyr::filter(grepl("DNA|REPAIR|HOMOLOGOUS|NONHOMOLOGOUS|FANCONI|NER|MMR|EXCISION", gs_name, ignore.case = TRUE))

# Named list: pathway -> vector of unique symbols
genesets_list <- brca_like_sets_filtered %>%
  split(.$gs_name) %>%
  lapply(function(df) unique(df$ensembl_gene))

# check redundancy
# Compute Jaccard similarity matrix
sets <- genesets_list
set_names <- names(sets)

jaccard_mat <- matrix(NA, nrow = length(sets), ncol = length(sets),
                      dimnames = list(set_names, set_names))

for (i in seq_along(sets)) {
  for (j in seq_along(sets)) {
    inter <- length(intersect(sets[[i]], sets[[j]]))
    uni   <- length(union(sets[[i]], sets[[j]]))
    jaccard_mat[i, j] <- inter / uni
  }
}

# VST
vst_mat <- assay(vst(dds, blind = TRUE))
gsvaPar <- gsvaParam(vst_mat,genesets_list )
ss <- gsva(gsvaPar)
# Average across replicates to get one column per condition
ph <- as.data.frame(colData(dds))
stopifnot("Condition" %in% names(ph))
ph <- ph[colnames(ss), , drop=FALSE]

conds <- levels(factor(ph$Condition))
mean3_all <- sapply(conds, function(k){
  idx <- rownames(ph)[ph$Condition == k]
  rowMeans(ss[, idx, drop=FALSE], na.rm=TRUE)
})
colnames(mean3_all) <- conds   # pathways x conditions

desired_order <- c("K0N0G","K0WT", "K0NES", 
                   "K0NLS", "K0SIG")
ss_means_z <- t(scale(t(mean3_all)))
ss_means_z <- ss_means_z[, desired_order]  


#order based on similarity of GO terms
common <- intersect(rownames(ss_means_z), rownames(jaccard_mat))
d      <- as.dist(1 - jaccard_mat[common, common])
hc     <- hclust(d, method = "average")

row_ord <- c(hc$labels[hc$order], setdiff(rownames(ss_means_z), common))  # append any unmatched

ss_means_z_sorted <- ss_means_z[row_ord, , drop = FALSE]

pheatmap::pheatmap(
  ss_means_z_sorted,
  cluster_rows = FALSE, cluster_cols = FALSE,
  show_rownames = TRUE, show_colnames = TRUE,
  color = colorRampPalette(c("blue","white","red"))(100),
  main = "Pathway activity (ordered by Jaccard similarity)"
)

#save df as csv
go_df <- as.data.frame(ss_means_z_sorted)

# add a column with rownames
go_df$GO_term <- rownames(go_df)

# add the ensembl IDs as a new column using the GO term names
go_df$ensembl_ids <- sapply(go_df$GO_term, function(term) {
  ids <- genesets_list[[term]]
  if (is.null(ids)) return(NA)             # if term not found in list
  paste(ids, collapse = ", ")              # collapse into a single string
})

# move GO_term to the first column for readability
go_df <- go_df[, c("GO_term", setdiff(names(go_df), "GO_term"))]

write.csv(go_df, "data45_Fig9H.csv")

#check BRCA1 counts ----------------
# Supplementary Figure 9H

#ENSG00000012048
rld <- rlog(dds, blind = F) #normalize the counts first
rld_counts <- assay(rld)
gene_symbols <- c("ENSG00000012048")
filtered_counts <- as.data.frame(rld_counts[rownames(rld_counts) %in% gene_symbols, ])
colnames(filtered_counts) <- "counts"

filtered_counts$condition <- rownames(filtered_counts)

filtered_counts <- filtered_counts %>%
  mutate(condition = sub("[0-9]+$", "", condition))
filtered_counts <- filtered_counts[!filtered_counts$condition %in% c("K0SIG", "K0N0G"),]

filtered_counts_mean <- filtered_counts %>%
  group_by(condition) %>%
  summarise(mean(counts), sd(counts))

colnames(filtered_counts_mean) <- c("condition", "Mean", "Std")



order <- c("K0WT", "K0NES", "K0NLS")
filtered_counts_mean$condition <- factor(filtered_counts_mean$condition, levels = order)
#plot counts
ggplot(filtered_counts_mean, aes(x = condition, y = Mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = Mean - Std, ymax = Mean + Std), width = 0.2) +
  theme_minimal() +
  labs(title = "Mean BRCA1 Counts with Standard Deviation",
       y = "Mean BRCA1 Counts")  + ylim(0,12)

#ggsave("BRCA1_counts.pdf",width = 8, height = 8)

#check ERK signatures ----------------------------------------------------------
  #Figure 7I ----
#extract signatures
#Extract Signatures of interest with ENSEMBLE IDs
# 1) ERK GO
go_erk <- msigdbr(species="Homo sapiens", category="C5", subcategory="GO:BP") %>%
  dplyr::filter(grepl("ERK1_AND_ERK2_CASCADE", gs_name))
# Reactome ERK-related
reactome_erk <- msigdbr(species="Homo sapiens", category="C2", subcategory="CP:REACTOME") %>%
  dplyr::filter(grepl("ERK", gs_name, ignore.case=TRUE))

#combine
ERK_sets <- bind_rows(go_erk, reactome_erk) %>%
  distinct(gs_name, ensembl_gene)

# Named list: pathway -> vector of unique symbols
genesets_list2 <- ERK_sets %>%
  split(.$gs_name) %>%
  lapply(function(df) unique(df$ensembl_gene))


#merge redundant gene signatures
gene_col <- "ensembl_gene"

merge_map <- list(
  ERK1_2_CASCADE_AND_REGULATION = c(
    "GOBP_ERK1_AND_ERK2_CASCADE",
    "GOBP_NEGATIVE_REGULATION_OF_ERK1_AND_ERK2_CASCADE",
    "GOBP_POSITIVE_REGULATION_OF_ERK1_AND_ERK2_CASCADE"
  ),
  ERK1_2_ACTIVATION = c(
    "REACTOME_MAPK1_ERK2_ACTIVATION",
    "REACTOME_MAPK3_ERK1_ACTIVATION"
  )
)

remove_sets <- c("REACTOME_ERKS_ARE_INACTIVATED")

# Build merged sets (union of genes; no duplicates)
merged_df <- lapply(names(merge_map), function(new_name) {
  old_names <- merge_map[[new_name]]
  ERK_sets %>%
    dplyr::filter(gs_name %in% old_names) %>%
    distinct(!!sym(gene_col)) %>%        # ensure each ID only once
    mutate(gs_name = new_name)
}) %>% bind_rows()

# Keep everything else except the sets being merged and the one to remove
kept_df <- ERK_sets %>%
  dplyr::filter(!gs_name %in% c(unlist(merge_map), remove_sets))

# Final merged table
ERK_sets_merged <- bind_rows(kept_df, merged_df) %>%
  dplyr::select(gs_name, all_of(gene_col)) %>%
  arrange(gs_name, !!sym(gene_col))
ERK_sets_merged_grouped <- ERK_sets_merged %>%
  group_by(gs_name) %>%
  summarise(ensembl_gene = paste(unique(ensembl_gene), collapse = ","), .groups = "drop")


# Named list: pathway -> vector of unique symbols
genesets_list2 <- ERK_sets_merged %>%
  split(.$gs_name) %>%
  lapply(function(df) unique(df$ensembl_gene))

#take all replicates and average later after scoring 
gsvaPar <- gsvaParam(vst_mat, genesets_list2) 
ss <- gsva(gsvaPar)                             
# 3) Average across replicates to get one column per condition
ph <- as.data.frame(colData(dds))
stopifnot("Condition" %in% names(ph))
ph <- ph[colnames(ss), , drop=FALSE]

conds <- levels(factor(ph$Condition))
mean3_all <- sapply(conds, function(k){
  idx <- rownames(ph)[ph$Condition == k]
  rowMeans(ss[, idx, drop=FALSE], na.rm=TRUE)
})
colnames(mean3_all) <- conds   # pathways x conditions



# 4) Z-score per pathway across conditions and plot
ss_means_z <- t(scale(t(mean3_all)))
desired_order <- c("K0WT", "K0N0G", 
                   "K0NES", "K0NLS", "K0SIG")

ss_means_z <- ss_means_z[, desired_order]  
pheatmap(ss_means_z,
         cluster_rows=FALSE, cluster_cols=FALSE,
         show_rownames=TRUE, show_colnames=TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main="ERK pathway activity") 

#save df as csv
go_df <- as.data.frame(ss_means_z)

# add a column with rownames
go_df$GO_term <- rownames(go_df)

# add the ensembl IDs as a new column using the GO term names
go_df$ensembl_ids <- sapply(go_df$GO_term, function(term) {
  ids <- genesets_list2[[term]]
  if (is.null(ids)) return(NA)             # if term not found in list
  paste(ids, collapse = ", ")              # collapse into a single string
})

# move GO_term to the first column for readability
go_df <- go_df[, c("GO_term", setdiff(names(go_df), "GO_term"))]

#write.csv(go_df, "data46_Fig9I.csv")

#Heatmap of ERK1/2 --------------------------------
  #Figure 7J ----
#Heatmap of ERK1/2 --------------------------------
# get pathway genes and z-score expression
sig_name <- "ERK1_2_ACTIVATION"
g <- intersect(genesets_list2[[sig_name]], rownames(vst_mat))
Ez_sig <- t(scale(t(vst_mat[g, , drop = FALSE])))

# map Ensembl -> SYMBOL
gm <- ensembldb::select(
  EnsDb.Hsapiens.v86,
  keys = rownames(Ez_sig),
  keytype = "GENEID",
  columns = c("GENEID","SYMBOL")
)
gm <- gm[!duplicated(gm$GENEID), ]
sym_lookup <- setNames(gm$SYMBOL, gm$GENEID)

# replace rownames with SYMBOLs (fallback: keep Ensembl if SYMBOL missing)
rownames(Ez_sig) <- ifelse(
  is.na(sym_lookup[rownames(Ez_sig)]) | sym_lookup[rownames(Ez_sig)] == "",
  rownames(Ez_sig),
  sym_lookup[rownames(Ez_sig)]
)

# column annotation
ann_col <- data.frame(
  Condition = ph[colnames(Ez_sig), "Condition", drop = TRUE],
  row.names = colnames(Ez_sig)
)


# group vector per sample
cond_vec <- ph[colnames(Ez_sig), "Condition", drop = TRUE]

# average replicates (columns) within each condition
Ez_cond <- sapply(split(colnames(Ez_sig), cond_vec),
                  function(s) rowMeans(Ez_sig[, s, drop = FALSE], na.rm = TRUE))


row_order <- order(Ez_cond[, "K0NLS"], decreasing = TRUE)
Ez_cond_sorted <- Ez_cond[row_order, ]

pheatmap::pheatmap(
  Ez_cond_sorted,
  cluster_rows = F, cluster_cols = FALSE,
  show_rownames = TRUE, show_colnames = TRUE,
  main = paste0(sig_name, " — condition means"),
  color = colorRampPalette(c("blue", "white", "red"))(100),
  border_color = NA
)
