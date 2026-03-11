########################################################################
# Bioinformatics Course Keio 2025 day3 前半
# Bulk RNA-seq DEG解析（RStudio用）
########################################################################

# ============================================================
# 作業ディレクトリとデータパスの設定
# ============================================================
# RStudioで開いたら、まず作業ディレクトリを設定してください。
# 以下のパスをご自身の環境に合わせて変更してください。
# （例：デスクトップに bioinformatics-course-keio-2025 フォルダを作成した場合）

# --- macOS の場合 ---
# project_dir <- "~/Downloads/bioinformatics-course-keio-2025"

# --- Windows の場合 ---
# project_dir <- "C:/Users/あなたのユーザー名/Downloads/bioinformatics-course-keio-2025"

project_dir <- "~/Downloads/bioinformatics-course-keio-2025"  # ← 自分の環境に合わせて変更
setwd(project_dir)

data_dir     <- file.path(project_dir, "data", "bulk")
meta_file    <- file.path(project_dir, "data", "bulk", "sample_meta.txt")
ifn_genes_file     <- file.path(project_dir, "data", "bulk", "IFNgenes100.txt")
img_dir      <- file.path(project_dir, "img")

# ============================================================
# ライブラリの読み込み
# ============================================================
library(ggplot2)
library(dplyr)
library(edgeR)
library(ggrepel)
library(gplots)
library(grid)
library(png)
library(variancePartition)

# ============================================================
# 使用するデータについて
# ============================================================
# 本実習では、末梢血から25種類の免疫細胞をFACSで分離しBulk RNA-seqを
# 行った公開データを使用します。
# 健常者（HC）とSLE（全身性エリテマトーデス）患者のTh1細胞について、
# DEG（発現変動遺伝子）解析を行います。

# ============================================================
# メタデータの読み込み
# ============================================================
meta <- read.table(meta_file, header = T)
table(meta$disease)

########################################################################
# PCAとは
########################################################################
# PCA（主成分分析）は、高次元データの次元を削減し、データの構造を
# 可視化するための手法です。

# --- PCA概念図 ---
if (file.exists(file.path(img_dir, "pca_concept.png"))) {
  img <- readPNG(file.path(img_dir, "pca_concept.png"))
  grid.newpage()
  grid.raster(img)
}

########################################################################
# DEG解析（SLE vs HC、Th1細胞）
########################################################################

# ============================================================
# Th1のカウントデータ読み込み
# ============================================================
file <- file.path(data_dir, "Th1_count.txt")
data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
dim(data)

HC_SLE_samples <- meta %>% filter(disease %in% c("HC", "SLE")) %>% pull(id)

# Th1が欠損していないサンプルに絞る
HC_SLE_samples <- intersect(colnames(data), HC_SLE_samples)
# HC,SLEサンプルだけに絞る
HC_SLE_cols <- intersect(colnames(data), c("Gene_id", "Gene_name", HC_SLE_samples))
data <- data %>% select(all_of(HC_SLE_cols))
dim(data)

head(data)

# ============================================================
# TMM正規化
# ============================================================
# TMM (Trimmed Mean of M-values) はサンプル間のライブラリ組成の違いを
# 補正する正規化手法。
# Reference: Robinson MD, Oshlack A.
#   "A scaling normalization method for differential expression analysis
#    of RNA-seq data."
#   Genome Biology. 2010;11:R25. doi:10.1186/gb-2010-11-3-r25
gene_info <- data %>% select(Gene_id, Gene_name)
count_matrix <- data %>% select(-Gene_id, -Gene_name)
rownames(count_matrix) <- data$Gene_id

dge <- DGEList(counts = as.matrix(count_matrix))

# フィルタリング（低発現遺伝子の除去）
sample_percentage <- 0.1
min_samples <- ceiling(ncol(dge) * sample_percentage)
keep <- rowSums(cpm(dge) > 1) >= min_samples
dge_filtered <- dge[keep, ]
print(paste0("フィルタリング後の遺伝子数: ", nrow(dge_filtered), " / 元の遺伝子数: ", nrow(dge)))

dge_filtered <- calcNormFactors(dge_filtered, method = "TMM")
log_cpm <- cpm(dge_filtered, log = TRUE)

# ============================================================
# PCA（HC vs SLE）
# ============================================================
gene_variance <- apply(log_cpm, 1, var)
top_genes_indices <- order(gene_variance, decreasing = TRUE)[1:min(5000, length(gene_variance))]
expression_top <- log_cpm[top_genes_indices, ]

expression_top_t <- t(expression_top)
pca_result <- prcomp(expression_top_t, center = TRUE, scale = TRUE)

var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
pc1_var <- round(var_explained[1] * 100, 1)
pc2_var <- round(var_explained[2] * 100, 1)

pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$id <- rownames(pca_df)
pca_df$disease <- left_join(pca_df, meta, by = c("id" = "id")) %>% pull(disease)
pca_df$disease <- factor(pca_df$disease)

# --- PCAプロット（HC vs SLE）---
pca_plot_disease <- ggplot(pca_df, aes(x = PC1, y = PC2, color = disease)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text_repel(
    aes(label = id),
    size = 3,
    box.padding = 0.35,
    point.padding = 0.5,
    segment.color = "grey50",
    max.overlaps = 2
  ) +
  labs(
    title = "PCA Analysis (Th1: HC vs SLE)",
    subtitle = paste0("Top ", length(top_genes_indices), " Most Variable Genes (TMM normalized & log2CPM data)"),
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)"),
    color = "Disease"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12)
  ) +
  coord_cartesian(expand = TRUE)

print(pca_plot_disease)

# ============================================================
# variancePartition：疾患が遺伝子発現の分散をどの程度説明するか
# ============================================================
# variancePartitionは各遺伝子の発現分散を、指定した要因（ここではdisease）と
# 残差（個人差など）に分解する。
# PCAでは全体的な傾向を見るが、variancePartitionは遺伝子ごとに
# 「diseaseがどれだけ発現変動を説明するか」を定量化できる。
#
# Reference: Hoffman GE, Schadt EE.
#   "variancePartition: interpreting drivers of variation in complex
#    gene expression studies."
#   BMC Bioinformatics. 2016;17(1):483. doi:10.1186/s12859-016-1323-z

# メタデータをlog_cpmのサンプル順に合わせる
vp_meta <- data.frame(
  id = colnames(log_cpm),
  disease = factor(meta$disease[match(colnames(log_cpm), meta$id)])
)

# モデル式：diseaseの効果 + 残差
# diseaseはカテゴリ変数なので (1|disease) としてランダム効果で指定
vp_form <- ~ (1|disease)

# variancePartition の実行
vp_result <- fitExtractVarPartModel(log_cpm, vp_form, vp_meta)

# 結果の要約
print("各遺伝子の分散に対するdiseaseの寄与率（中央値）:")
print(summary(vp_result))

# Violin plotで可視化
vp_plot <- plotVarPart(vp_result)
print(vp_plot)

# diseaseの寄与が大きい上位遺伝子を確認
vp_sorted <- sortCols(vp_result)
head(vp_sorted, 20)

# ============================================================
# edgeRによるDEG解析
# ============================================================
#
# edgeRについて:
#   RNA-seqカウントデータのDEG解析のための代表的なBioconductorパッケージ。
#   負の二項分布モデルを用いて、サンプル間の遺伝子発現の差を統計的に検定する。
#   Reference: Robinson MD, McCarthy DJ, Smyth GK.
#     "edgeR: a Bioconductor package for differential expression analysis
#      of digital gene expression data."
#     Bioinformatics. 2010;26(1):139-140. doi:10.1093/bioinformatics/btp616
#
#   公式ユーザーガイド（式の書き方の詳細やさまざまな解析デザインの例が豊富）:
#   https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf
#
# 補足: DESeq2について
#   DEG解析のもう一つの代表的なツールとしてDESeq2がある。
#   edgeRと同様に負の二項分布モデルを使うが、分散の縮小推定の方法が異なる。
#   Reference: Love MI, Huber W, Anders S.
#     "Moderated estimation of fold change and dispersion for RNA-seq data
#      with DESeq2."
#     Genome Biology. 2014;15:550. doi:10.1186/s13059-014-0550-8
#
# 今回の解析で使用するアプローチ: Quasi-likelihood (QL) F-test
#   edgeRにはLRT（尤度比検定）とQL F-testの2つの検定方法がある。
#   QL F-testは分散推定の不確実性を考慮するため、特にサンプル数が少ない場合に
#   偽陽性の制御が優れており、edgeR公式でも推奨されている。
#   各ステップの意味:
#     estimateDisp()  : 負の二項分布の分散パラメータを推定
#     glmQLFit()      : 準尤度（QL）負の二項GLMをフィットし、
#                        経験ベイズ法で遺伝子ごとのQL分散を縮小推定
#     glmQLFTest()    : QL F検定により、指定した係数についてDEGを検出

print("疾患グループ：")
print(levels(pca_df$disease))

dge_filtered$samples$group <- factor(meta$disease[match(colnames(dge_filtered), meta$id)])

print("各疾患グループのサンプル数:")
print(table(dge_filtered$samples$group))

disease_groups <- levels(dge_filtered$samples$group)
reference_group <- disease_groups[1]
treatment_group <- disease_groups[2]
print(paste0("比較: ", treatment_group, " vs ", reference_group))

# デザインマトリックスの設定
# ~group は「切片 + グループ差」のモデル。SLE_vs_HC の係数がSLEとHCの差を表す。
design <- model.matrix(~group, data = dge_filtered$samples)
colnames(design) <- c("Intercept", paste0(treatment_group, "_vs_", reference_group))
head(design)

# Step 1: 負の二項分布の分散を推定
dge_filtered <- estimateDisp(dge_filtered, design)

# Step 2: 準尤度GLMをフィット（経験ベイズによる分散の縮小推定を含む）
fit <- glmQLFit(dge_filtered, design)

# Step 3: QL F検定でDEGを検出
qlf <- glmQLFTest(fit, coef = paste0(treatment_group, "_vs_", reference_group))

# 結果の取得
results <- topTags(qlf, n = Inf, sort.by = "PValue")
results_df <- as.data.frame(results)

# 有意なDEGの数をカウント（FDR < 0.05）
de_genes_up <- sum(results_df$logFC > 0 & results_df$FDR < 0.05)
de_genes_down <- sum(results_df$logFC < 0 & results_df$FDR < 0.05)

print(paste0(treatment_group, "で発現上昇している遺伝子数（FDR < 0.05）: ", de_genes_up))
print(paste0(treatment_group, "で発現低下している遺伝子数（FDR < 0.05）: ", de_genes_down))

# 結果にGene_nameを追加
results_df$Gene_id <- rownames(results_df)
results_df <- left_join(results_df, gene_info, by = "Gene_id")

print("DEG解析結果の上位:")
head(results_df, 10)

# ============================================================
# Volcano Plotの作成
# ============================================================
volcano_data <- results_df
volcano_data$significance <- "Not Significant"
volcano_data$significance[volcano_data$logFC >= 0 & volcano_data$FDR < 0.05] <- "Up-regulated"
volcano_data$significance[volcano_data$logFC < 0 & volcano_data$FDR < 0.05] <- "Down-regulated"
volcano_data$significance <- factor(volcano_data$significance,
                                    levels = c("Up-regulated", "Down-regulated", "Not Significant"))

top_genes_label <- volcano_data %>%
  filter(FDR < 0.05) %>%
  arrange(FDR) %>%
  head(10)

volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = -log10(FDR), color = significance)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgray") +
  scale_color_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue", "Not Significant" = "gray")) +
  geom_text_repel(
    data = top_genes_label,
    aes(label = Gene_name),
    size = 3,
    max.overlaps = 20,
    box.padding = 0.5
  ) +
  labs(
    title = paste0("Volcano Plot: ", treatment_group, " vs ", reference_group),
    x = "log2 Fold Change",
    y = "-log10 FDR",
    color = "Regulation"
  ) +
  theme_bw() +
  xlim(-6, 6) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

print(volcano_plot)

# ============================================================
# DEG解析結果のヒートマップ（上位50DEG）
# ============================================================
top_degs <- results_df %>%
  filter(FDR < 0.05) %>%
  arrange(FDR) %>%
  head(50)

deg_expression <- log_cpm[top_degs$Gene_id, HC_SLE_samples]
rownames(deg_expression) <- top_degs$Gene_name

sample_order <- meta %>%
  filter(id %in% HC_SLE_samples) %>%
  arrange(disease) %>%
  pull(id)

deg_expression <- deg_expression[, sample_order]

my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 100)

disease_colors <- c("darkgreen", "darkred")
names(disease_colors) <- c("HC", "SLE")

sample_disease <- meta$disease[match(colnames(deg_expression), meta$id)]
col_colors <- disease_colors[sample_disease]

heatmap.2(as.matrix(deg_expression),
          col = my_palette,
          scale = "row",
          trace = "none",
          density.info = "none",
          key = TRUE,
          keysize = 1.5,
          ColSideColors = col_colors,
          margins = c(8, 10),
          cexRow = 0.8,
          cexCol = 0.8,
          labCol = FALSE,
          main = "Top 50 DEGs: SLE vs HC",
          dendrogram = "both")

legend("topright",
       legend = c("HC", "SLE"),
       fill = disease_colors,
       border = FALSE,
       bty = "n",
       cex = 0.8)

# ============================================================
# Pathway解析の例（IFN関連遺伝子のエンリッチメント）
# ============================================================

# IFN関連遺伝子グループの読み込み
IFNgenes <- read.table(ifn_genes_file, header = F)[, 1]
head(IFNgenes)
length(IFNgenes)

# 発現データに含まれるIFN遺伝子に絞る
expressed_genes_ids <- rownames(dge_filtered)
expressed_genes_names <- gene_info$Gene_name[match(expressed_genes_ids, gene_info$Gene_id)]
expressed_genes_names <- expressed_genes_names[!is.na(expressed_genes_names) & expressed_genes_names != ""]

filtered_IFNgenes <- intersect(IFNgenes, expressed_genes_names)
print(paste0("発現データに含まれるIFN遺伝子数: ", length(filtered_IFNgenes)))

# 発現上昇/低下DEGを取得
up_DEGs <- results_df %>%
  filter(FDR < 0.05 & logFC > 0) %>%
  pull(Gene_name)

down_DEGs <- results_df %>%
  filter(FDR < 0.05 & logFC < 0) %>%
  pull(Gene_name)

print(paste0("発現上昇DEG数: ", length(up_DEGs)))
print(paste0("発現低下DEG数: ", length(down_DEGs)))

# IFN遺伝子と発現上昇/低下DEGの重複を確認
ifn_in_up_DEGs <- intersect(up_DEGs, filtered_IFNgenes)
print(paste0("発現上昇したIFN関連遺伝子数: ", length(ifn_in_up_DEGs)))

ifn_in_down_DEGs <- intersect(down_DEGs, filtered_IFNgenes)
print(paste0("発現低下したIFN関連遺伝子数: ", length(ifn_in_down_DEGs)))

# upDEGでのFisherの正確確率検定によるエンリッチメント解析
total_genes <- length(expressed_genes_names)

contingency_up <- matrix(c(
  length(ifn_in_up_DEGs),
  length(up_DEGs) - length(ifn_in_up_DEGs),
  length(filtered_IFNgenes) - length(ifn_in_up_DEGs),
  total_genes - length(up_DEGs) - (length(filtered_IFNgenes) - length(ifn_in_up_DEGs))
), nrow = 2)

rownames(contingency_up) <- c("IFN", "non-IFN")
colnames(contingency_up) <- c("upDEG", "non-upDEG")

print("発現上昇DEGに対するIFN遺伝子の2x2Table:")
contingency_up

# Fisherの正確確率検定
fisher_test_up <- fisher.test(contingency_up, alternative = "greater")
print("発現上昇DEGに対するIFN遺伝子のエンリッチメント分析結果:")
print(fisher_test_up)
