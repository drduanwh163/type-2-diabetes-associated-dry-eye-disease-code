# DEGs----
library(DESeq2)
library(FactoMineR)
library(factoextra)
library(limma)
library(ggrepel)
library(showtext)
library(pheatmap)

DATA_DIR <- "path/to/data"
OUTPUT_DIR <- "path/to/output"

COUNT_FILE <- file.path(DATA_DIR, "gene_count.xls")
SAMPLE_FILE <- file.path(DATA_DIR, "sample.group.xls")
TPM_FILE <- file.path(DATA_DIR, "gene_tpm.xls")
GO_FILE <- file.path(DATA_DIR, "Unigene.GO.xls")
KEGG_FILE <- file.path(DATA_DIR, "Unigene.KEGG.xls")

GROUP_LEVELS <- c("GROUP_1", "GROUP_2", "GROUP_3", "CONTROL")
GROUP_1 <- "GROUP_1"
GROUP_2 <- "CONTROL"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

count <- read.table(COUNT_FILE, sep = "\t", row.names = 1, check.names = FALSE,
                    stringsAsFactors = FALSE, header = TRUE)
sample <- read.table(SAMPLE_FILE, sep = "\t", check.names = FALSE,
                     stringsAsFactors = FALSE, header = TRUE)
tpm <- read.table(TPM_FILE, sep = "\t", row.names = 1, check.names = FALSE,
                  stringsAsFactors = FALSE, header = TRUE)

exp <- log2(tpm[2:length(colnames(tpm))] + 1)
expr1 <- count

group_list <- factor(sample$Group, levels = GROUP_LEVELS)

data <- as.data.frame(t(expr1))
data.pca <- PCA(data, graph = FALSE)

pca_plot <- fviz_pca_ind(
  data.pca,
  geom.ind = "point",
  col.ind = group_list,
  pointsize = 1.3,
  alpha.ind = 0.6,
  addEllipses = TRUE,
  alpha.ellipse = 0.3,
  legend.title = "group",
  repel = TRUE
) +
  theme_bw() +
  ggtitle(paste(GROUP_LEVELS, collapse = "-")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13)
  )

ggsave(
  pca_plot,
  filename = file.path(OUTPUT_DIR, "PCA.pdf"),
  width = 5.5,
  height = 4
)

sample_sub <- sample[sample$Group %in% c(GROUP_1, GROUP_2), , drop = FALSE]
expr_sub <- expr1[, rownames(sample_sub), drop = FALSE]

sample_sub$Group <- factor(
  sample_sub$Group,
  levels = c(GROUP_1, GROUP_2)
)

dds <- DESeqDataSetFromMatrix(
  countData = expr1,
  colData = data.frame(condition = factor(sample$Group)),
  design = ~ condition
)

dds <- DESeq(dds)
res <- results(dds)

logFC_threshold <- log2(2)
pval_threshold <- 0.05

DEG_deseq2 <- as.data.frame(res)
DEG_deseq2$gene_id <- rownames(DEG_deseq2)
DEG_deseq2$Regulation <- "ns"

DEG_deseq2$Regulation[
  (DEG_deseq2$pvalue < pval_threshold) &
    (DEG_deseq2$log2FoldChange < -logFC_threshold)
] <- "Down"

DEG_deseq2$Regulation[
  (DEG_deseq2$pvalue < pval_threshold) &
    (DEG_deseq2$log2FoldChange > logFC_threshold)
] <- "Up"

DEG_deseq2$FoldChange <- 2^DEG_deseq2$log2FoldChange

comparison_name <- paste(GROUP_1, "vs", GROUP_2, sep = "-")
comparison_dir <- file.path(OUTPUT_DIR, comparison_name)
dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)

write.table(
  DEG_deseq2,
  file = file.path(
    comparison_dir,
    paste0(comparison_name, "-diff-p-0.05-FC-2.gene.xls")
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

comparison_files <- c(
  file.path(OUTPUT_DIR, "GROUP_1-vs-CONTROL/1.GROUP_1-vs-CONTROL-diff-p-0.05-FC-2.gene.xls"),
  file.path(OUTPUT_DIR, "GROUP_2-vs-CONTROL/1.GROUP_2-vs-CONTROL-diff-p-0.05-FC-2.gene.xls"),
  file.path(OUTPUT_DIR, "GROUP_3-vs-CONTROL/1.GROUP_3-vs-CONTROL-diff-p-0.05-FC-2.gene.xls"),
  file.path(OUTPUT_DIR, "COMBINATION-vs-GROUP_3/1.COMBINATION-vs-GROUP_3-diff-p-0.05-FC-2.gene.xls"),
  file.path(OUTPUT_DIR, "COMBINATION-vs-GROUP_2/1.COMBINATION-vs-GROUP_2-diff-p-0.05-FC-2.gene.xls")
)

DbXe <- read.table(comparison_files[1], sep = "\t", header = TRUE, check.names = FALSE)
Xe <- read.table(comparison_files[2], sep = "\t", header = TRUE, check.names = FALSE)
Db <- read.table(comparison_files[3], sep = "\t", header = TRUE, check.names = FALSE)
DbXe.Db <- read.table(comparison_files[4], sep = "\t", header = TRUE, check.names = FALSE)
DbXe.Xe <- read.table(comparison_files[5], sep = "\t", header = TRUE, check.names = FALSE)

DbXe$FoldChange <- 2^DbXe$log2FoldChange
Xe$FoldChange <- 2^Xe$log2FoldChange
Db$FoldChange <- 2^Db$log2FoldChange

for (x in list(DbXe, Xe, Db)) {
  x$Regulation <- "ns"
  x$Regulation[
    (x$`p-value` < pval_threshold) &
      (x$log2FoldChange < -logFC_threshold)
  ] <- "Down"
  x$Regulation[
    (x$`p-value` < pval_threshold) &
      (x$log2FoldChange > logFC_threshold)
  ] <- "Up"
}

DbXe$Regulation <- "ns"
DbXe$Regulation[
  (DbXe$`p-value` < pval_threshold) &
    (DbXe$log2FoldChange < -logFC_threshold)
] <- "Down"
DbXe$Regulation[
  (DbXe$`p-value` < pval_threshold) &
    (DbXe$log2FoldChange > logFC_threshold)
] <- "Up"

Xe$Regulation <- "ns"
Xe$Regulation[
  (Xe$`p-value` < pval_threshold) &
    (Xe$log2FoldChange < -logFC_threshold)
] <- "Down"
Xe$Regulation[
  (Xe$`p-value` < pval_threshold) &
    (Xe$log2FoldChange > logFC_threshold)
] <- "Up"

Db$Regulation <- "ns"
Db$Regulation[
  (Db$`p-value` < pval_threshold) &
    (Db$log2FoldChange < -logFC_threshold)
] <- "Down"
Db$Regulation[
  (Db$`p-value` < pval_threshold) &
    (Db$log2FoldChange > logFC_threshold)
] <- "Up"

write.table(DbXe, comparison_files[1], sep = "\t", row.names = FALSE)
write.table(Xe, comparison_files[2], sep = "\t", row.names = FALSE)
write.table(Db, comparison_files[3], sep = "\t", row.names = FALSE)

DEG_deseq2 <- DbXe.Xe

top_up <- DEG_deseq2[DEG_deseq2$Regulation == "Up", ]
top_up <- top_up[order(-top_up$log2FoldChange), ][1:5, ]

top_down <- DEG_deseq2[DEG_deseq2$Regulation == "Down", ]
top_down <- top_down[order(top_down$log2FoldChange), ][1:5, ]

top_genes <- rbind(top_up, top_down)

p2 <- ggplot(
  DEG_deseq2,
  aes(x = log2FoldChange, y = -log10(`p-value`), colour = Regulation)
) +
  geom_point(alpha = 0.8, size = 3) +
  geom_vline(
    xintercept = c(-log2(2), log2(2)),
    lty = 4,
    col = "black",
    lwd = 0.8
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    lty = 4,
    col = "black",
    lwd = 0.8
  ) +
  labs(
    title = comparison_name,
    x = "log2FoldChange",
    y = "-log10(P.Value)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.title = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 20)
  ) +
  geom_text_repel(
    data = top_genes,
    aes(label = gene_name),
    size = 6,
    color = "black",
    max.overlaps = 100,
    box.padding = 1,
    point.padding = 0.3,
    segment.size = 0.6,
    segment.alpha = 0.5,
    min.segment.length = 0
  )

ggsave(
  p2,
  filename = file.path(OUTPUT_DIR, paste0(comparison_name, ".Volcano.pdf")),
  width = 8,
  height = 6
)

DbXe <- read.table(comparison_files[1], sep = "\t", header = TRUE, check.names = FALSE)
Xe <- read.table(comparison_files[2], sep = "\t", header = TRUE, check.names = FALSE)
Db <- read.table(comparison_files[3], sep = "\t", header = TRUE, check.names = FALSE)
DbXe.Db <- read.table(comparison_files[4], sep = "\t", header = TRUE, check.names = FALSE)
DbXe.Xe <- read.table(comparison_files[5], sep = "\t", header = TRUE, check.names = FALSE)

DEG_deseq2 <- DbXe.Xe

up_genes <- DEG_deseq2[DEG_deseq2$Regulation == "Up", ]
down_genes <- DEG_deseq2[DEG_deseq2$Regulation == "Down", ]

top10_up <- up_genes[order(-up_genes$log2FoldChange), ][1:10, "gene_id"]
top10_down <- down_genes[order(down_genes$log2FoldChange), ][1:10, "gene_id"]

selected_genes <- c(top10_up, top10_down)

TOP10.sig.DEGs <- DEG_deseq2[DEG_deseq2$gene_id %in% selected_genes, ]

write.table(
  TOP10.sig.DEGs,
  file = file.path(
    OUTPUT_DIR,
    paste0(comparison_name, ".TOP10.sig.DEGs.FC2.xls")
  ),
  sep = "\t"
)

keep <- sample$Group %in% c(GROUP_1, GROUP_2)
sample_sub <- sample[keep, ]

heatmap_data <- exp[, sample_sub$Sample][selected_genes, , drop = FALSE]

id2name <- setNames(
  TOP10.sig.DEGs$gene_name,
  TOP10.sig.DEGs$gene_id
)

rownames(heatmap_data) <- id2name[rownames(heatmap_data)]

group_list <- factor(
  sample_sub$Group,
  levels = c(GROUP_1, GROUP_2)
)

annotation_col <- data.frame(Group = group_list)
rownames(annotation_col) <- colnames(heatmap_data)

pdf(
  file = file.path(
    OUTPUT_DIR,
    paste0(comparison_name, ".heatmap.FC2.pdf")
  ),
  width = 8,
  height = 6
)

pheatmap(
  mat = heatmap_data,
  scale = "row",
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  angle_col = 45,
  border_color = "black",
  fontsize = 12,
  fontsize_col = 8,
  fontsize_row = 10
)

dev.off()

top_up <- DEG_deseq2[DEG_deseq2$Regulation == "Up", ]
top_down <- DEG_deseq2[DEG_deseq2$Regulation == "Down", ]
top_genes <- rbind(top_up, top_down)

exp.sig <- exp[rownames(top_genes), ]
write.table(
  exp.sig,
  file.path(OUTPUT_DIR, "significant_gene_expression.xls"),
  sep = "\t"
)

go <- read.table(
  GO_FILE,
  sep = "\t",
  header = TRUE
)

kegg <- read.table(
  KEGG_FILE,
  sep = "\t",
  header = TRUE
)

DEG_deseq2 <- read.table(
  file.path(OUTPUT_DIR, "comparison-diff-p-0.05-FC-2.gene.xls"),
  sep = "\t",
  header = TRUE
)

rownames(DEG_deseq2) <- DEG_deseq2$ID

top_up <- DEG_deseq2[DEG_deseq2$Regulation == "up", ]
top_up <- top_up[order(-top_up$log2FoldChange), ][1:5, ]

top_down <- DEG_deseq2[DEG_deseq2$Regulation == "down", ]
top_down <- top_down[order(top_down$log2FoldChange), ][1:5, ]

top_genes <- rbind(top_up, top_down)

p2 <- ggplot(
  DEG_deseq2,
  aes(x = log2FoldChange, y = -log10(pvalue), colour = Regulation)
) +
  geom_point(alpha = 0.8, size = 3) +
  geom_vline(
    xintercept = c(-log2(2), log2(2)),
    lty = 4,
    col = "black",
    lwd = 0.8
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    lty = 4,
    col = "black",
    lwd = 0.8
  ) +
  labs(
    title = "Comparison",
    x = "log2FoldChange",
    y = "-log10(P.Value)"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    legend.title = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 20)
  ) +
  geom_text_repel(
    data = top_genes,
    aes(label = rownames(top_genes)),
    size = 6,
    color = "black",
    max.overlaps = 100,
    box.padding = 1,
    point.padding = 0.3,
    segment.size = 0.6,
    segment.alpha = 0.5,
    min.segment.length = 0
  )

ggsave(
  p2,
  filename = file.path(OUTPUT_DIR, "Comparison.Volcano.pdf"),
  width = 10.8,
  height = 8.8
)

up_genes <- DEG_deseq2[DEG_deseq2$Regulation == "up", ]
down_genes <- DEG_deseq2[DEG_deseq2$Regulation == "down", ]

top10_up <- up_genes[order(-up_genes$log2FoldChange), ][1:10, "ID"]
top10_down <- down_genes[order(down_genes$log2FoldChange), ][1:10, "ID"]

selected_genes <- c(top10_up, top10_down)

TOP10.sig.DEGs <- DEG_deseq2[selected_genes, ]

write.table(
  TOP10.sig.DEGs,
  file.path(OUTPUT_DIR, "TOP10.sig.DEGs.xls"),
  sep = "\t"
)

heatmap_data <- fpkm[selected_genes, ]

annotation_col <- data.frame(Group = group_list)
rownames(annotation_col) <- colnames(heatmap_data)

pdf(
  file = file.path(OUTPUT_DIR, "TOP10.heatmap.pdf"),
  width = 8,
  height = 6
)

pheatmap(
  mat = heatmap_data,
  scale = "row",
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  angle_col = 45,
  border_color = "black",
  fontsize = 12,
  fontsize_col = 8,
  fontsize_row = 10
)

dev.off()

pdf(
  file = file.path(OUTPUT_DIR, "ALL.heatmap.pdf"),
  width = 8,
  height = 15
)

pheatmap(
  mat = heatmap_data,
  scale = "row",
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  angle_col = 45,
  border_color = "black",
  fontsize = 12,
  fontsize_col = 8,
  fontsize_row = 10
)

dev.off()

#GSEA----
library(clusterProfiler)
library(org.Hs.eg.db)
library(fgsea)
library(enrichplot)
library(ggplot2)
library(dplyr)

setwd("PATH_TO_WORKING_DIRECTORY")

deg_file <- "PATH_TO_DEG_FILE"
gmt_file <- "PATH_TO_GMT_FILE"
output_file <- "GSEA_RESULTS_FILE"

deg <- read.table(
  deg_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

geneList <- deg$log2FoldChange
names(geneList) <- deg$gene_name

geneList <- geneList[!is.na(names(geneList))]
geneList <- geneList[order(abs(geneList), decreasing = TRUE)]
geneList <- geneList[!duplicated(names(geneList))]
geneList <- sort(geneList, decreasing = TRUE)

pathways <- gmtPathways(gmt_file)

fgseaRes <- fgsea(
  pathways = pathways,
  stats = geneList,
  minSize = 15,
  maxSize = 500,
  nperm = 10000
)

fgseaRes <- fgseaRes %>%
  arrange(pval)

fgsea_out <- as.data.frame(fgseaRes)

fgsea_out$leadingEdge <- sapply(
  fgsea_out$leadingEdge,
  paste,
  collapse = ","
)

write.table(
  fgsea_out,
  output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

top_n <- 10
plot_dir <- "GSEA_plots"

dir.create(
  plot_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

top10 <- fgseaRes %>%
  filter(pval < 0.05) %>%
  arrange(desc(abs(NES))) %>%
  slice_head(n = top_n)

for (i in seq_len(nrow(top10))) {
  
  pathway_name <- gsub(
    "[/\\\\:*?\"<>|]",
    "_",
    top10$pathway[i]
  )
  
  pdf(
    file.path(plot_dir, paste0(pathway_name, ".pdf")),
    width = 2,
    height = 2
  )
  
  print(
    plotEnrichment(
      pathways[[top10$pathway[i]]],
      geneList
    ) +
      labs(title = pathway_name) +
      theme_bw(base_size = 9)
  )
  
  dev.off()
}

#GSVA----
library(GSVA)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(GSEABase)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

setwd("PATH_TO_WORKING_DIRECTORY")

expression_file <- "PATH_TO_EXPRESSION_FILE"
group_file <- "PATH_TO_GROUP_FILE"
gmt_file <- "PATH_TO_GMT_FILE"

expression_data <- read.delim(
  expression_file,
  header = TRUE,
  stringsAsFactors = FALSE
)

group_data <- read.table(
  group_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

expression_data <- expression_data %>%
  group_by(gene_name) %>%
  summarise(
    across(where(is.numeric), mean),
    .groups = "drop"
  )

expression_data[, -1] <- log2(expression_data[, -1] + 1)

gmt <- getGmt(gmt_file)

pathway_names <- sapply(gmt, setName)

target_pathways <- c(
  "PATHWAY_1",
  "PATHWAY_2",
  "PATHWAY_3",
  "PATHWAY_4",
  "PATHWAY_5"
)

gene_sets <- gmt[pathway_names %in% target_pathways]

expr_matrix <- as.data.frame(expression_data)
rownames(expr_matrix) <- expr_matrix$gene_name
expr_matrix <- expr_matrix[, -1]
expr_matrix <- as.matrix(expr_matrix)

param <- gsvaParam(
  exprData = expr_matrix,
  geneSets = gene_sets
)

gsva_scores <- gsva(param)

write.table(
  gsva_scores,
  "GSVA_scores.xls",
  sep = "\t",
  quote = FALSE
)

heatmap_data <- t(gsva_scores)
heatmap_data <- scale(heatmap_data)

color_fun <- colorRamp2(
  c(-1, 0, 1),
  c("COLOR_LOW", "white", "COLOR_HIGH")
)

row_annotation <- rowAnnotation(
  Group = group_data$Group,
  col = list(
    Group = c(
      GROUP_1 = "COLOR_1",
      GROUP_2 = "COLOR_2",
      GROUP_3 = "COLOR_3",
      GROUP_4 = "COLOR_4"
    )
  )
)

pdf(
  "GSVA_heatmap.pdf",
  width = 5,
  height = 13
)

Heatmap(
  heatmap_data,
  na_col = "white",
  col = color_fun,
  name = "GSVA score",
  show_column_names = TRUE,
  show_row_names = FALSE,
  row_split = group_data$Group,
  row_title = NULL,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  left_annotation = row_annotation,
  column_names_gp = gpar(fontsize = 15),
  row_names_gp = gpar(fontsize = 10),
  row_title_gp = gpar(fontsize = 15),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 15),
    labels_gp = gpar(fontsize = 13)
  ),
  rect_gp = gpar(col = "black", lwd = 1),
  height = unit(nrow(heatmap_data) * 5, "mm"),
  width = unit(ncol(heatmap_data) * 10, "mm"),
  row_gap = unit(0, "mm"),
  column_gap = unit(1.5, "mm")
)

dev.off()

plot_data <- as.data.frame(t(gsva_scores))

plot_data$Group <- group_data$Group
plot_data$Sample <- rownames(plot_data)

plot_data <- melt(
  plot_data,
  id.vars = c("Sample", "Group")
)

colnames(plot_data) <- c(
  "Sample",
  "Group",
  "Pathway",
  "Score"
)

ggplot(
  plot_data,
  aes(
    x = Pathway,
    y = Score,
    fill = Group
  )
) +
  geom_boxplot(width = 0.6) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

gene_ids <- expression_data$gene_id

expression_values <- log2(
  expression_data[, -1] + 1
)

expression_data_ssgsea <- data.frame(
  gene_id = gene_ids,
  expression_values,
  check.names = FALSE
)

expr_matrix_ssgsea <- expression_data_ssgsea
rownames(expr_matrix_ssgsea) <- toupper(
  expr_matrix_ssgsea$gene_id
)

expr_matrix_ssgsea <- expr_matrix_ssgsea[, -1]

GENE_SET_1 <- toupper(
  as.character(GENE_SET_1[[1]])
)

GENE_SET_2 <- toupper(
  as.character(GENE_SET_2[[1]])
)

GENE_SET_1 <- GENE_SET_1[
  GENE_SET_1 != "" &
    !is.na(GENE_SET_1)
]

GENE_SET_2 <- GENE_SET_2[
  GENE_SET_2 != "" &
    !is.na(GENE_SET_2)
]

gene_sets_custom <- GeneSetCollection(
  list(
    GeneSet(
      GENE_SET_1,
      setName = "GENE_SET_1",
      geneIdType = SymbolIdentifier()
    ),
    GeneSet(
      GENE_SET_2,
      setName = "GENE_SET_2",
      geneIdType = SymbolIdentifier()
    )
  )
)

param_ssgsea <- gsvaParam(
  exprData = as.matrix(expr_matrix_ssgsea),
  geneSets = gene_sets_custom
)

ssgsea_scores <- gsva(param_ssgsea)

write.table(
  ssgsea_scores,
  "ssGSEA_scores.xls",
  sep = "\t",
  quote = FALSE
)

ssgsea_heatmap_data <- as.data.frame(
  t(ssgsea_scores)
)

ssgsea_heatmap_data <- ssgsea_heatmap_data %>%
  mutate_all(as.numeric)

ssgsea_heatmap_data <- scale(
  ssgsea_heatmap_data
)

ssgsea_heatmap_data[
  ssgsea_heatmap_data > 2
] <- 2

ssgsea_heatmap_data[
  ssgsea_heatmap_data < -2
] <- -2

color_fun <- colorRamp2(
  c(-2, 0, 2),
  c("COLOR_LOW", "white", "COLOR_HIGH")
)

group_colors <- c(
  GROUP_1 = "COLOR_1",
  GROUP_2 = "COLOR_2"
)

row_annotation <- rowAnnotation(
  Group = group_data$Group,
  col = list(
    Group = group_colors
  )
)

pdf(
  "GSVA_heatmap.pdf",
  width = 8,
  height = 10
)

Heatmap(
  ssgsea_heatmap_data,
  na_col = "white",
  col = color_fun,
  name = "ssGSEA score",
  show_column_names = TRUE,
  show_row_names = FALSE,
  row_split = group_data$Group,
  row_title = NULL,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  left_annotation = row_annotation,
  column_names_gp = gpar(fontsize = 15),
  row_names_gp = gpar(fontsize = 10),
  row_title_gp = gpar(fontsize = 15),
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 15),
    labels_gp = gpar(fontsize = 13)
  ),
  rect_gp = gpar(col = "black", lwd = 1),
  height = unit(
    nrow(ssgsea_heatmap_data) * 4,
    "mm"
  ),
  width = unit(
    ncol(ssgsea_heatmap_data) * 50,
    "mm"
  ),
  row_gap = unit(0, "mm"),
  column_gap = unit(1.5, "mm")
)

dev.off()

# SVM-RFE----
library(e1071)
library(openxlsx)
library(caret)
library(dplyr)

source("PATH_TO_SVM_RFE_SCRIPT")

nfold <- 5

set.seed(1999)

expression_file <- "PATH_TO_EXPRESSION_MATRIX"
label_file <- "PATH_TO_GROUP_INFORMATION"

expression_data <- read.table(
  expression_file,
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

labels <- read.table(
  label_file,
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

expression_data <- expression_data[, rownames(labels)]

expression_data <- as.data.frame(
  t(log2(expression_data + 1))
)

expression_data$group <- labels$Group

expression_data <- expression_data %>%
  select(group, everything())

colnames(expression_data) <- make.names(
  colnames(expression_data)
)

input <- expression_data
input$group <- as.factor(input$group)

nrows <- nrow(input)

folds <- rep(
  1:nfold,
  len = nrows
)[sample(nrows)]

folds <- lapply(
  1:nfold,
  function(x) which(folds == x)
)

results <- lapply(
  folds,
  svmRFE.wrap,
  input,
  k = 1,
  halve.above = 100
)

top.features <- WriteFeatures(
  results,
  input,
  save = FALSE
)

featsweep <- lapply(
  1:32,
  FeatSweep.wrap,
  results,
  input
)

write.table(
  top.features,
  "SVM_RFE_feature_ranking.xls",
  sep = "\t"
)

no.info <- min(
  prop.table(table(input[, 1]))
)

errors <- sapply(
  featsweep,
  function(x) ifelse(
    is.null(x),
    NA,
    x$error
  )
)

pdf(
  "SVM_RFE_error_rate.pdf",
  height = 3.5,
  width = 4
)

PlotErrors(
  errors,
  no.info = no.info
)

dev.off()

min.error.idx <- which.min(errors)

best.features <- top.features[
  1:min.error.idx,
]

write.table(
  best.features,
  "SVM_RFE_selected_features.xls",
  sep = "\t",
  col.names = NA,
  quote = FALSE
)

# RandomForest----

library(randomForest)
library(fastshap)
library(shapviz)
library(ggplot2)

setwd("PATH_TO_WORKING_DIRECTORY")

expression_file <- "PATH_TO_EXPRESSION_MATRIX"
label_file <- "PATH_TO_GROUP_INFORMATION"

expression_data <- read.table(
  expression_file,
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

labels <- read.table(
  label_file,
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

expression_data <- expression_data[, rownames(labels)]

expression_data <- t(
  log2(expression_data + 1)
)

set.seed(1888)

rf_model <- randomForest(
  x = expression_data,
  y = as.factor(labels$Group),
  ntree = 100,
  importance = TRUE
)

pdf(
  "RandomForest_OOB_error.pdf",
  width = 5,
  height = 5
)

plot(
  rf_model,
  main = "Random Forest",
  lwd = 2
)

dev.off()

optimal_trees <- which.min(
  rf_model$err.rate[, 1]
)

optimal_trees

rf_model_final <- randomForest(
  x = expression_data,
  y = as.factor(labels$Group),
  ntree = optimal_trees,
  importance = TRUE
)

gene_importance <- importance(
  rf_model_final
)

gene_importance_df <- data.frame(
  Gene = colnames(expression_data),
  Importance = gene_importance[, "MeanDecreaseGini"]
)

gene_importance_df <- gene_importance_df[
  order(
    gene_importance_df$Importance,
    decreasing = TRUE
  ),
]

write.table(
  gene_importance_df,
  "gene_importance.xls",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

head(
  gene_importance_df,
  30
)

pdf(
  "geneImportance.pdf",
  width = 5,
  height = 6
)

varImpPlot(
  rf_model_final,
  n.var = 32,
  main = ""
)

dev.off()

selected_genes <- gene_importance_df$Gene[
  gene_importance_df$Importance > 0
]

X <- expression_data[, selected_genes, drop = FALSE]

rf_model_shap <- randomForest(
  x = X,
  y = as.factor(labels$Group)
)

pred_wrapper <- function(object, newdata) {
  predict(
    object,
    newdata = newdata,
    type = "prob"
  )[, "TARGET_GROUP"]
}

shap_values <- fastshap::explain(
  rf_model_shap,
  X = X,
  pred_wrapper = pred_wrapper,
  nsim = 100
)

shap_object <- shapviz(
  as.matrix(shap_values),
  X = X
)

pdf(
  "SHAP_barplot.pdf",
  width = 6,
  height = 5
)

print(
  sv_importance(
    shap_object,
    kind = "bar",
    max_display = 30
  ) +
    ggtitle(
      "SHAP mean(|value|) - Random Forest"
    )
)

dev.off()

pdf(
  "SHAP_beeswarm.pdf",
  width = 6,
  height = 5
)

print(
  sv_importance(
    shap_object,
    kind = "both",
    max_display = 30
  ) +
    ggtitle(
      "SHAP beeswarm - Random Forest"
    )
)

dev.off()

# wilcox.test----
library(ggplot2)
library(ggpubr)
library(reshape2)
library(dplyr)
library(tidyverse)
library(gghalves)
library(ggsignif)
library(rstatix)
library(limma)

setwd("PATH_TO_WORKING_DIRECTORY")

expression_file <- "PATH_TO_EXPRESSION_MATRIX"
group_file <- "PATH_TO_GROUP_INFORMATION"

gene_list <- c(
  "KRT86",
  "MUC1",
  "RGL3",
  "LRRC32",
  "MSC",
  "PF4V1",
  "CRYBB1",
  "OLFM1"
)

expression_data <- read.delim(
  expression_file,
  row.names = 1
)

expression_data <- expression_data[
  rownames(expression_data) %in% gene_list,
]

expression_data <- log2(
  expression_data + 1
)

expression_data <- as.data.frame(
  t(expression_data)
)

expression_data$Sample <- rownames(
  expression_data
)

group_data <- read.table(
  group_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

merged_data <- merge(
  expression_data,
  group_data,
  by.x = "Sample",
  by.y = "Sample"
)

merged_data$Group <- factor(
  merged_data$Group,
  levels = c(
    "GROUP_1",
    "GROUP_2",
    "GROUP_3",
    "GROUP_4"
  )
)

long_data <- melt(
  merged_data,
  id.vars = c("Sample", "Group"),
  variable.name = "Gene",
  value.name = "Expression"
)

long_data$Gene <- factor(
  long_data$Gene,
  levels = gene_list
)

plot_groups <- c(
  "GROUP_1",
  "GROUP_4"
)

plot_data <- long_data[
  long_data$Group %in% plot_groups,
]

plot_data$Group <- factor(
  plot_data$Group,
  levels = plot_groups
)

pdf(
  "Key_gene_expression.pdf",
  width = 8,
  height = 8
)

ggviolin(
  plot_data,
  x = "Gene",
  y = "Expression",
  fill = "Group",
  add = "boxplot",
  palette = c(
    "GROUP_1" = "COLOR_1",
    "GROUP_4" = "COLOR_2"
  ),
  position = position_dodge(width = 1),
  scale = "width"
) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 16
    ),
    axis.text.y = element_text(
      size = 16
    ),
    axis.title = element_text(
      face = "bold",
      size = 18
    ),
    legend.text = element_text(
      size = 15
    ),
    legend.title = element_text(
      size = 18
    ),
    plot.title = element_text(
      face = "bold"
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(
      color = "black"
    ),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 1
    )
  ) +
  stat_compare_means(
    aes(group = Group),
    method = "wilcox.test",
    label = "p.signif",
    size = 7,
    label.y = 6.5
  )

dev.off()

# gene.roc----
library(pROC)

expression_file_1 <- "PATH_TO_EXPRESSION_MATRIX_1"
expression_file_2 <- "PATH_TO_EXPRESSION_MATRIX_2"
group_file <- "PATH_TO_GROUP_INFORMATION"

genes <- c(
  "KRT86",
  "MUC1",
  "RGL3",
  "LRRC32",
  "MSC",
  "PF4V1",
  "CRYBB1",
  "OLFM1"
)

group_levels_1 <- c(
  "GROUP_1",
  "GROUP_2"
)

group_levels_2 <- c(
  "GROUP_3",
  "GROUP_2"
)

plot_colors <- c(
  "#CA2C2C",
  "#2775AC",
  "#2A9D8F",
  "#F2692F",
  "#8E44AD",
  "#E67E22",
  "#16A085",
  "#34495E"
)

prepare_data <- function(
    expression_file,
    group_file,
    group_levels
) {
  
  expression_data <- read.delim(
    expression_file,
    row.names = 1
  )
  
  expression_data <- log2(
    expression_data + 1
  )
  
  expression_data <- as.data.frame(
    t(expression_data)
  )
  
  expression_data$Sample <- rownames(
    expression_data
  )
  
  group_data <- read.table(
    group_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
  )
  
  group_data <- group_data[
    group_data$Group %in% group_levels,
  ]
  
  group_data$Group <- factor(
    group_data$Group,
    levels = group_levels
  )
  
  merged_data <- merge(
    expression_data,
    group_data,
    by = "Sample"
  )
  
  rownames(merged_data) <- merged_data$Sample
  merged_data$Sample <- NULL
  
  return(merged_data)
}

plot_roc <- function(
    data,
    genes,
    output_file,
    colors
) {
  
  pdf(
    output_file,
    width = 3.5,
    height = 3.5
  )
  
  roc_model <- roc(
    data[, "Group"],
    data[[genes[1]]],
    smooth = TRUE
  )
  
  auc_text <- paste0(
    genes[1],
    ": ",
    sprintf(
      "%0.3f",
      auc(roc_model)
    )
  )
  
  plot(
    roc_model,
    col = colors[1],
    main = "",
    lwd = 2
  )
  
  if (length(genes) > 1) {
    
    for (i in 2:length(genes)) {
      
      roc_model <- roc(
        data[, "Group"],
        data[[genes[i]]],
        smooth = TRUE
      )
      
      lines(
        roc_model,
        col = colors[i],
        lwd = 2
      )
      
      auc_text <- c(
        auc_text,
        paste0(
          genes[i],
          ": ",
          sprintf(
            "%0.3f",
            auc(roc_model)
          )
        )
      )
    }
  }
  
  legend(
    "bottomright",
    auc_text,
    lwd = 2,
    bty = "n",
    col = colors[seq_along(genes)],
    cex = 0.56
  )
  
  dev.off()
}

data_1 <- prepare_data(
  expression_file = expression_file_1,
  group_file = group_file,
  group_levels = group_levels_1
)

plot_roc(
  data = data_1,
  genes = genes,
  output_file = "ROC_comparison_1.pdf",
  colors = plot_colors
)

data_2 <- prepare_data(
  expression_file = expression_file_2,
  group_file = group_file,
  group_levels = group_levels_2
)

plot_roc(
  data = data_2,
  genes = genes,
  output_file = "ROC_comparison_2.pdf",
  colors = plot_colors
)

# bulk.ssGSEA----
library(GSVA)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(GSEABase)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(ggsci)

setwd("PATH_TO_WORKING_DIRECTORY")

expression_file <- "PATH_TO_EXPRESSION_MATRIX"
immune_gene_set_file <- "PATH_TO_IMMUNE_GENE_SET"
group_file <- "PATH_TO_GROUP_INFORMATION"

ssgsea_score_file <- "ssGSEA_scores.xls"
immune_gmt_file <- "immune_cell_type.gmt"
heatmap_file <- "ssGSEA_heatmap.pdf"
boxplot_file <- "ssGSEA_boxplot.pdf"

expression_data <- read.delim(
  expression_file,
  row.names = 1,
  check.names = FALSE
)

immune_gene_sets <- read.table(
  immune_gene_set_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)

expression_data <- log2(expression_data + 1)

gene_set_list <- immune_gene_sets %>%
  group_by(Cell.type) %>%
  summarise(genes = list(unique(Metagene)), .groups = "drop") %>%
  rowwise() %>%
  mutate(
    gene_set = list(
      GeneSet(
        setName = Cell.type,
        geneIds = genes,
        geneIdType = SymbolIdentifier(),
        collectionType = NullCollection()
      )
    )
  ) %>%
  pull(gene_set)

gene_set_collection <- GeneSetCollection(gene_set_list)

toGmt(
  gene_set_collection,
  con = immune_gmt_file
)

ssgsea_param <- ssgseaParam(
  exprData = as.matrix(expression_data),
  geneSets = gene_set_collection
)

ssgsea_scores <- gsva(ssgsea_param)

write.table(
  ssgsea_scores,
  ssgsea_score_file,
  sep = "\t"
)

group_data <- read.table(
  group_file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  stringsAsFactors = FALSE
)

heatmap_data <- as.data.frame(t(ssgsea_scores))

stopifnot(identical(rownames(heatmap_data), rownames(group_data)))

heatmap_data <- heatmap_data %>%
  mutate(across(everything(), as.numeric))

heatmap_data <- t(heatmap_data)

stopifnot(identical(colnames(heatmap_data), rownames(group_data)))

color_fun <- colorRamp2(
  c(0, 0.5, 1),
  c("#194E7A", "white", "#CA2C2C")
)

group_levels <- c(
  "GROUP_1",
  "GROUP_2",
  "GROUP_3",
  "GROUP_4"
)

group_colors <- c(
  "GROUP_1" = "#EC6356",
  "GROUP_2" = "#B45CAF",
  "GROUP_3" = "#65A2DE",
  "GROUP_4" = "#40A6A0"
)

group_data$Group <- factor(
  group_data$Group,
  levels = group_levels
)

top_annotation <- HeatmapAnnotation(
  Group = group_data$Group,
  col = list(Group = group_colors)
)

column_order <- colnames(heatmap_data)[
  order(
    factor(
      group_data$Group,
      levels = group_levels
    )
  )
]

pdf(
  heatmap_file,
  height = 5,
  width = 10
)

Heatmap(
  heatmap_data,
  na_col = "white",
  col = color_fun,
  show_row_names = TRUE,
  show_column_names = FALSE,
  column_order = column_order,
  column_split = group_data$Group,
  gap = unit(0, "mm"),
  border = FALSE,
  column_title = NULL,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  top_annotation = top_annotation,
  row_names_gp = gpar(fontsize = 12),
  column_names_gp = gpar(fontsize = 0),
  name = "ssGSEA"
)

dev.off()

plot_data <- as.data.frame(t(ssgsea_scores))

stopifnot(identical(rownames(plot_data), rownames(group_data)))

plot_data$Group <- group_data$Group
plot_data$Sample <- rownames(plot_data)

plot_data <- melt(
  plot_data,
  id.vars = c("Sample", "Group")
)

colnames(plot_data) <- c(
  "Sample",
  "Group",
  "imm.type",
  "ssGSEA_scores"
)

plot_data <- plot_data[
  plot_data$Group %in% c("GROUP_1", "GROUP_4"),
]

plot_data$Group <- factor(
  plot_data$Group,
  levels = c("GROUP_1", "GROUP_4")
)

boxplot_ssgsea <- ggplot(
  plot_data,
  aes(
    x = imm.type,
    y = ssGSEA_scores
  )
) +
  labs(
    y = "imm ssGSEA Scores",
    x = ""
  ) +
  geom_boxplot(
    aes(fill = Group),
    position = position_dodge(0.5),
    width = 0.5
  ) +
  scale_fill_manual(
    values = c(
      "GROUP_1" = "#CA2C2C",
      "GROUP_4" = "#194E7A"
    )
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(
      size = 15,
      color = "black"
    ),
    axis.text = element_text(
      size = 13,
      color = "black"
    ),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(
      hjust = 1,
      angle = 45
    ),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 15)
  ) +
  stat_compare_means(
    aes(group = Group),
    label = "p.signif",
    method = "wilcox.test",
    hide.ns = TRUE,
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05, 1),
      symbols = c("***", "**", "*", "ns")
    ),
    size = 9
  )

ggsave(
  boxplot_file,
  boxplot_ssgsea,
  width = 15,
  height = 8
)

# Cell.Gene.corr----
library(limma)
library(psych)
library(reshape2)
library(pheatmap)
library(dplyr)

expression_file <- "PATH_TO_EXPRESSION_MATRIX"
ssgsea_file <- "PATH_TO_SSGSEA_SCORE"
group_file <- "PATH_TO_GROUP_INFORMATION"

gene_correlation_file <- "gene_immune_correlation.xls"
gene_correlation_heatmap <- "gene_immune_correlation_heatmap.pdf"
cell_correlation_file <- "immune_cell_correlation_within_group.xls"

genes <- c(
  "MUC1",
  "PF4V1",
  "CRYBB1",
  "OLFM1"
)

group_levels <- c(
  "GROUP_1",
  "GROUP_2",
  "GROUP_3",
  "GROUP_4"
)

expr1 <- read.delim(
  expression_file,
  row.names = 1
)

expr1 <- t(as.matrix(log2(expr1 + 1)))

expr1 <- expr1[, genes, drop = FALSE]

expr2 <- read.table(
  ssgsea_file,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

expr2 <- t(expr2)

stopifnot(identical(rownames(expr1), rownames(expr2)))

cor_res <- corr.test(
  expr1,
  expr2,
  method = "spearman",
  adjust = "fdr"
)

cor_matrix <- cor_res$r
p_matrix <- cor_res$p

cor_data <- melt(
  cor_matrix,
  value.name = "corr"
)

cor_data$pvalue <- as.vector(p_matrix)

colnames(cor_data) <- c(
  "Gene",
  "Cell",
  "corr",
  "pvalue"
)

write.table(
  cor_data,
  gene_correlation_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

sig_levels <- function(p) {
  if (is.na(p)) {
    return("")
  }
  if (p < 0.001) {
    return("***")
  }
  if (p < 0.01) {
    return("**")
  }
  if (p < 0.05) {
    return("*")
  }
  return("")
}

sig_matrix <- matrix(
  "",
  nrow = nrow(p_matrix),
  ncol = ncol(p_matrix),
  dimnames = dimnames(p_matrix)
)

for (i in seq_len(nrow(p_matrix))) {
  for (j in seq_len(ncol(p_matrix))) {
    sig_matrix[i, j] <- sig_levels(p_matrix[i, j])
  }
}

pdf(
  gene_correlation_heatmap,
  width = 8,
  height = 8
)

pheatmap(
  cor_matrix,
  display_numbers = sig_matrix,
  number_color = "black",
  color = colorRampPalette(
    c("#194E7A", "#F2F2F2", "#CA2C2C")
  )(20),
  fontsize_number = 16,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "Correlation Heatmap",
  angle_col = 45,
  cellwidth = 32,
  cellheight = 32,
  border_color = "black"
)

dev.off()

ssgsea_data <- read.table(
  ssgsea_file,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

ssgsea_data <- t(ssgsea_data)

clinical_data <- read.table(
  group_file,
  sep = "\t",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

stopifnot(all(rownames(ssgsea_data) %in% rownames(clinical_data)))

clinical_data <- clinical_data[
  rownames(ssgsea_data),
  ,
  drop = FALSE
]

clinical_data$Group <- factor(
  clinical_data$Group,
  levels = group_levels
)

expr2_df <- as.data.frame(ssgsea_data)
expr2_df$Group <- clinical_data$Group

group_list <- unique(
  as.character(expr2_df$Group)
)

group_list <- group_list[!is.na(group_list)]

cor_list <- list()

for (g in group_list) {
  tmp <- expr2_df %>%
    filter(Group == g) %>%
    select(-Group)
  
  tmp <- tmp[
    complete.cases(tmp),
    ,
    drop = FALSE
  ]
  
  cor_res <- corr.test(
    tmp,
    method = "spearman",
    adjust = "none"
  )
  
  cor_matrix <- cor_res$r
  p_matrix <- cor_res$p
  
  q_matrix <- matrix(
    p.adjust(
      as.vector(p_matrix),
      method = "fdr"
    ),
    nrow = nrow(p_matrix),
    ncol = ncol(p_matrix),
    dimnames = dimnames(p_matrix)
  )
  
  cor_data <- melt(
    cor_matrix,
    value.name = "corr"
  )
  
  cor_data$pvalue <- as.vector(p_matrix)
  cor_data$qvalue <- as.vector(q_matrix)
  cor_data$Group <- g
  
  colnames(cor_data) <- c(
    "Cell1",
    "Cell2",
    "corr",
    "pvalue",
    "qvalue",
    "Group"
  )
  
  cor_list[[g]] <- cor_data
}

cor_all <- bind_rows(cor_list)

write.table(
  cor_all,
  cell_correlation_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# SigGene.GSEA----
library(clusterProfiler)
library(enrichplot)
library(future.apply)
library(ggplot2)
library(GseaVis)
library(jjAnno)
library(dplyr)
library(stringr)

setwd("PATH_TO_WORKING_DIRECTORY")

expression_file <- "PATH_TO_EXPRESSION_MATRIX"
gmt_file <- "PATH_TO_GMT_FILE"
output_directory <- "PATH_TO_OUTPUT_DIRECTORY"

key_genes <- c(
  "MUC1",
  "PF4V1",
  "CRYBB1",
  "OLFM1"
)

expression_data <- read.delim(
  expression_file,
  row.names = 1
)

expression_data <- log2(expression_data + 1)

hallmark_gene_sets <- read.gmt(gmt_file)

batch_cor <- function(gene) {
  y <- as.numeric(expression_data[gene, ])
  gene_names <- rownames(expression_data)
  
  do.call(
    rbind,
    future_lapply(
      gene_names,
      function(x) {
        test_result <- cor.test(
          as.numeric(expression_data[x, ]),
          y,
          method = "spearman"
        )
        
        data.frame(
          gene = gene,
          mRNAs = x,
          cor = test_result$estimate,
          p.value = test_result$p.value
        )
      }
    )
  )
}

run_keygene_GSEA <- function(
    gene_symbol,
    gene_sets,
    outdir = output_directory
) {
  cor_results <- batch_cor(gene_symbol)
  
  gene_df <- data.frame(
    cor = cor_results$cor,
    SYMBOL = cor_results$mRNAs
  )
  
  gene_list <- gene_df$cor
  names(gene_list) <- gene_df$SYMBOL
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  gsea_result <- GSEA(
    gene_list,
    TERM2GENE = gene_sets,
    pvalueCutoff = 1,
    pAdjustMethod = "BH"
  )
  
  significant_results <- subset(
    gsea_result@result,
    pvalue < 0.05 & qvalue < 0.25
  )
  
  write.table(
    significant_results,
    paste0(
      outdir,
      gene_symbol,
      "_GSEA_results.txt"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  plot_data <- as.data.frame(gsea_result@result) %>%
    filter(pvalue < 0.05) %>%
    arrange(desc(abs(NES))) %>%
    slice_head(n = 5)
  
  plot_data$Description <- str_replace(
    plot_data$Description,
    "^HALLMARK_",
    ""
  )
  
  wrap_fun <- label_wrap_gen(width = 25)
  
  plot_data$Description_wrapped <- plot_data$Description %>%
    str_replace_all("_", " ") %>%
    wrap_fun() %>%
    as.character()
  
  xmin <- min(
    plot_data$enrichmentScore,
    na.rm = TRUE
  ) * 1.2
  
  xmax <- max(
    plot_data$enrichmentScore,
    na.rm = TRUE
  ) * 1.2
  
  pdf(
    paste0(
      outdir,
      gene_symbol,
      "_GSEA_bubble_custom.pdf"
    ),
    width = 6,
    height = 4
  )
  
  print(
    ggplot(
      plot_data,
      aes(
        x = enrichmentScore,
        y = reorder(
          Description_wrapped,
          setSize
        )
      )
    ) +
      geom_point(
        aes(
          size = setSize,
          color = p.adjust
        )
      ) +
      scale_color_gradient(
        low = "#2775AC",
        high = "#EC6356"
      ) +
      scale_size(
        range = c(8, 15)
      ) +
      coord_cartesian(
        xlim = c(xmin, xmax)
      ) +
      labs(
        title = paste0(
          "GSEA of ",
          gene_symbol
        ),
        x = "enrichmentScore",
        y = NULL,
        color = "p.adjust",
        size = "setSize"
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(
          size = 20,
          face = "bold",
          hjust = 0.5
        ),
        axis.text = element_text(
          size = 13
        ),
        axis.title = element_text(
          size = 16
        ),
        legend.title = element_text(
          size = 18
        ),
        legend.text = element_text(
          size = 16
        )
      )
  )
  
  dev.off()
  
  return(gsea_result)
}

plot_GSEA_gseaNb <- function(
    gsea_result,
    gene_symbol,
    outdir = output_directory,
    topN = 5,
    legend.position,
    curveCol
) {
  result_data <- gsea_result@result
  
  gene_set_id <- result_data %>%
    filter(pvalue < 0.05) %>%
    arrange(desc(abs(NES))) %>%
    slice_head(n = topN) %>%
    pull(ID)
  
  plot <- gseaNb(
    object = gsea_result,
    geneSetID = gene_set_id,
    subPlot = 3,
    rmHt = TRUE,
    legend.position = legend.position,
    curveCol = curveCol
  )
  
  pdf(
    paste0(
      outdir,
      gene_symbol,
      "_GseaVis.pdf"
    ),
    width = 5,
    height = 5
  )
  
  print(plot)
  
  dev.off()
}

legend_positions <- list(
  MUC1 = c(0.30, 0.35),
  PF4V1 = c(0.36, 0.43),
  CRYBB1 = c(0.35, 0.35),
  OLFM1 = c(0.32, 0.35)
)

curve_colors <- list(
  MUC1 = c(
    "#CA2C2C",
    "#F2692F",
    "#00B4D7",
    "#EEA63A",
    "#AA3183"
  ),
  PF4V1 = c(
    "#EC6353",
    "#41A990",
    "#00B4D7",
    "#EEA63A",
    "#AA3183"
  ),
  CRYBB1 = c(
    "#EC6353",
    "#41A990",
    "#00B4D7",
    "#EEA63A",
    "#AA3183"
  ),
  OLFM1 = c(
    "#EC6353",
    "#41A990",
    "#00B4D7",
    "#EEA63A",
    "#AA3183"
  )
)

for (gene in key_genes) {
  gsea_result <- run_keygene_GSEA(
    gene_symbol = gene,
    gene_sets = hallmark_gene_sets
  )
  
  plot_GSEA_gseaNb(
    gsea_result = gsea_result,
    gene_symbol = gene,
    legend.position = legend_positions[[gene]],
    curveCol = curve_colors[[gene]]
  )
}

# Clinical.Logistic----
library(broom)
library(showtext)
library(forestploter)
library(grid)


clinical_file <- "PATH_TO_CLINICAL_DATA"
raw_result_file <- "Clinical_Univariate_Logistic_Raw.xls"
formatted_result_file <- "Clinical_Univariate_Logistic.xls"
forest_plot_file <- "Clinical_Logistic_Forest.pdf"

data <- read.delim(
  clinical_file,
  check.names = FALSE,
  row.names = 1
)

train <- data[, 1:11, drop = FALSE]

categorical_vars <- c(
  "DbDE",
  "Sex",
  "Symtom_dryness",
  "Symtom_foreign_body_sensation",
  "Symtom_burning",
  "Symtom_redness",
  "Symtom_ocular_fatigue"
)

continuous_vars <- c(
  "Age",
  "Diabetes_Duration",
  "OSDI",
  "BUT_mean"
)

for (v in categorical_vars) {
  train[[v]] <- as.factor(train[[v]])
}

train$DbDE <- as.factor(train$DbDE)

results <- data.frame()

vars <- setdiff(
  colnames(train),
  "DbDE"
)

for (v in vars) {
  formula <- as.formula(
    paste("DbDE ~ `", v, "`", sep = "")
  )
  
  fit <- tryCatch(
    glm(
      formula,
      data = train,
      family = binomial
    ),
    error = function(e) NULL
  )
  
  if (!is.null(fit) && fit$converged) {
    coef_table <- coef(summary(fit))
    
    if (nrow(coef_table) >= 2) {
      beta <- coef_table[2, 1]
      se <- coef_table[2, 2]
      p_value <- coef_table[2, 4]
      
      results <- rbind(
        results,
        data.frame(
          Variable = v,
          OR = exp(beta),
          CI_Low = exp(beta - 1.96 * se),
          CI_High = exp(beta + 1.96 * se),
          P.Value = p_value
        )
      )
    }
  }
}

results <- results[
  order(-results$OR),
]

write.table(
  results,
  raw_result_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

results$OR <- round(
  results$OR,
  2
)

results$CI_Low <- round(
  results$CI_Low,
  2
)

results$CI_High <- round(
  results$CI_High,
  2
)

results$`OR (95% CI)` <- paste0(
  sprintf("%.2f", results$OR),
  " (",
  sprintf("%.2f", results$CI_Low),
  " - ",
  sprintf("%.2f", results$CI_High),
  ")"
)

results$P.Value.Display <- ifelse(
  results$P.Value < 0.001,
  "<0.001",
  sprintf("%.3f", results$P.Value)
)

write.table(
  results,
  formatted_result_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

showtext_auto()

results$` ` <- paste(
  rep(" ", 25),
  collapse = " "
)

forest_theme <- forest_theme(
  ci_pch = 15,
  ci_lty = 1,
  ci_lwd = 2,
  ci_Theight = 0.3,
  refline_lwd = 1,
  refline_lty = "dashed",
  base_size = 13,
  vertical_spacing = 1.5
)

forest_data <- results[
  ,
  c(
    "Variable",
    " ",
    "OR (95% CI)",
    "P.Value.Display"
  )
]

forest_plot <- forest(
  forest_data,
  est = results$OR,
  lower = results$CI_Low,
  upper = results$CI_High,
  sizes = 0.8,
  ci_column = 3,
  ref_line = 1,
  xlim = c(0, 5),
  ticks_at = c(1),
  theme = forest_theme
)

forest_plot <- insert_text(
  forest_plot,
  text = "Clinical Logistic Regression",
  col = 1:4,
  part = "header",
  gp = gpar(fontface = "bold")
)

forest_plot <- add_border(
  forest_plot,
  part = "header",
  where = "bottom"
)

forest_plot <- edit_plot(
  forest_plot,
  col = 1:4,
  which = "text",
  hjust = unit(0.5, "npc"),
  x = unit(0.5, "npc")
)

forest_plot <- edit_plot(
  forest_plot,
  row = which(results$P.Value < 0.05),
  which = "background",
  gp = gpar(fill = "#FDEFDD")
)

forest_plot <- edit_plot(
  forest_plot,
  row = which(results$P.Value >= 0.05),
  which = "background",
  gp = gpar(fill = "white")
)

pdf(
  forest_plot_file,
  height = 5,
  width = 8
)

print(forest_plot)

dev.off()

# Nomogram----
library(rms)
library(regplot)
library(pROC)

clinical_file <- "PATH_TO_CLINICAL_DATA"
expression_file <- "PATH_TO_EXPRESSION_MATRIX"

nomogram_file <- "Nomogram.pdf"
regplot_file <- "Nomogram_regplot.pdf"
roc_file <- "Nomogram_ROC.pdf"
calibration_file <- "Nomogram_Calibration.pdf"

gene_list <- c(
  "MUC1",
  "PF4V1",
  "CRYBB1",
  "OLFM1"
)

sig_vars <- c(
  "Sex",
  "Diabetes_Duration",
  "BUT_mean",
  "MUC1",
  "PF4V1",
  "CRYBB1",
  "OLFM1"
)

data <- read.delim(
  clinical_file,
  check.names = FALSE,
  row.names = 1
)

train <- data[, 1:11, drop = FALSE]

categorical_vars <- c(
  "DbDE",
  "Sex",
  "Symtom_dryness",
  "Symtom_foreign_body_sensation",
  "Symtom_burning",
  "Symtom_redness",
  "Symtom_ocular_fatigue"
)

for (v in categorical_vars) {
  train[[v]] <- as.factor(train[[v]])
}

train <- train[
  ,
  c(
    "DbDE",
    "Sex",
    "Diabetes_Duration",
    "BUT_mean"
  ),
  drop = FALSE
]

expression_data <- read.delim(
  expression_file,
  row.names = 1
)

expression_data <- expression_data[
  rownames(expression_data) %in% gene_list,
  ,
  drop = FALSE
]

expression_data <- t(
  log2(expression_data + 1)
)

train <- merge(
  train,
  expression_data,
  by = "row.names",
  all = FALSE
)

rownames(train) <- train$Row.names
train$Row.names <- NULL

train$DbDE <- as.factor(train$DbDE)

dd <- datadist(train)
options(datadist = "dd")

formula_multi <- as.formula(
  paste(
    "DbDE ~",
    paste(
      sprintf("`%s`", sig_vars),
      collapse = " + "
    )
  )
)

fit_rms <- lrm(
  formula_multi,
  data = train,
  x = TRUE,
  y = TRUE
)

nomogram_model <- nomogram(
  fit_rms,
  fun = plogis,
  fun.at = c(0.01, 0.5, 0.99),
  funlabel = "Probability of DbDE",
  lp = TRUE
)

pdf(
  nomogram_file,
  height = 10,
  width = 15
)

plot(
  nomogram_model,
  lplabel = "Linear Predictor",
  xfrac = 0.2,
  tcl = -0.2,
  lmgp = 0.2,
  points.label = "Points",
  total.points.label = "Total Points",
  cap.labels = FALSE,
  cex.var = 1.5,
  cex.axis = 1.3,
  col.grid = gray(c(0.8, 0.95))
)

dev.off()

colnames(train) <- make.names(
  colnames(train)
)

sig_vars2 <- make.names(sig_vars)

formula_regplot <- as.formula(
  paste(
    "DbDE ~",
    paste(
      sig_vars2,
      collapse = " + "
    )
  )
)

fit_regplot <- lrm(
  formula_regplot,
  data = train,
  x = TRUE,
  y = TRUE
)

pdf(
  regplot_file,
  width = 8,
  height = 8
)

regplot(
  fit_regplot,
  dencol = "#F2692F",
  boxcol = "#D670B5",
  plots = c("violin", "boxes"),
  observation = train[8, ],
  center = TRUE,
  subticks = TRUE,
  droplines = TRUE,
  title = "Nomogram",
  points = TRUE,
  odds = FALSE,
  showP = TRUE,
  rank = "sd",
  interval = "confidence",
  clickable = FALSE
)

dev.off()

train$pred_prob <- predict(
  fit_rms,
  newdata = train,
  type = "fitted"
)

roc_train <- roc(
  train$DbDE,
  train$pred_prob
)

auc_value <- auc(roc_train)

ci_auc <- ci.auc(
  roc_train
)

pdf(
  roc_file,
  width = 3.5,
  height = 3.5
)

plot(
  roc_train,
  col = "#CA2C2C",
  lwd = 2,
  main = paste(
    "Nomogram ROC (AUC =",
    round(auc_value, 3),
    ")"
  )
)

legend(
  "bottomright",
  cex = 0.8,
  legend = paste(
    "AUC =",
    round(auc_value, 3),
    "\n95% CI: [",
    round(ci_auc[1], 3),
    ", ",
    round(ci_auc[3], 3),
    "]"
  ),
  bty = "n",
  col = "#CA2C2C",
  lwd = 2
)

dev.off()

cal_train <- calibrate(
  fit_rms,
  method = "boot",
  B = 1000
)

pdf(
  calibration_file,
  width = 5,
  height = 5
)

plot(
  cal_train,
  xlab = "Predicted probability",
  ylab = "Observed probability",
  main = "Calibration Curve - Nomogram"
)

abline(
  0,
  1,
  lty = 2,
  col = "gray"
)

dev.off()

# Clinical.Gene_Corr----
library(psych)
library(reshape2)
library(pheatmap)

expression_file <- "PATH_TO_EXPRESSION_MATRIX"
clinical_file <- "PATH_TO_CLINICAL_DATA"

correlation_file <- "correlation_results.xls"
heatmap_file <- "correlation_heatmap.pdf"

gene_list <- c(
  "MUC1",
  "PF4V1",
  "CRYBB1",
  "OLFM1"
)

clinical_vars <- c(
  "Diabetes_Duration",
  "BUT_mean",
  "Disease_Severity"
)

expr1 <- read.delim(
  expression_file,
  row.names = 1
)

expr1 <- t(
  as.matrix(
    log2(expr1 + 1)
  )
)

expr1 <- expr1[, gene_list, drop = FALSE]

expr2 <- read.delim(
  clinical_file,
  check.names = FALSE,
  row.names = 1
)

expr2 <- expr2[, clinical_vars, drop = FALSE]
expr2 <- expr2[rownames(expr1), , drop = FALSE]

stopifnot(identical(rownames(expr1), rownames(expr2)))

cor_result <- corr.test(
  expr1,
  expr2,
  method = "spearman",
  adjust = "fdr"
)

cor_matrix <- cor_result$r
p_matrix <- cor_result$p

cor_df <- melt(
  cor_matrix,
  value.name = "corr"
)

cor_df$pvalue <- as.vector(p_matrix)

colnames(cor_df) <- c(
  "Gene",
  "Clinical",
  "corr",
  "pvalue"
)

write.table(
  cor_df,
  correlation_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

sig_levels <- function(p) {
  if (is.na(p)) {
    return("")
  }
  
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  }
  
  return("")
}

sig_matrix <- matrix(
  "",
  nrow = nrow(p_matrix),
  ncol = ncol(p_matrix),
  dimnames = dimnames(p_matrix)
)

for (i in seq_len(nrow(p_matrix))) {
  for (j in seq_len(ncol(p_matrix))) {
    sig_matrix[i, j] <- sig_levels(
      p_matrix[i, j]
    )
  }
}

pdf(
  heatmap_file,
  width = 6,
  height = 5
)

pheatmap(
  cor_matrix,
  display_numbers = sig_matrix,
  number_color = "black",
  color = colorRampPalette(
    c(
      "#194E7A",
      "#F2F2F2",
      "#CA2C2C"
    )
  )(20),
  fontsize_number = 16,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "Correlation Heatmap",
  angle_col = 45,
  cellwidth = 32,
  cellheight = 32,
  border_color = "black"
)

dev.off()

