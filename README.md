# Bulk RNA-seq DEG解析 ワークショップ

Bioinformatics Course Keio 2025 Day3 の教材です。

## 概要

末梢血免疫細胞のBulk RNA-seqデータを用いて、以下を実習します：

- **PCA**：Th1細胞のHC vs SLEの主成分分析による可視化
- **variancePartition**：疾患（SLE vs HC）が各遺伝子の発現分散をどの程度説明するかの定量化
- **DEG解析**：Th1細胞における健常者（HC）vs SLE患者の発現変動遺伝子の同定（edgeR）
- **Volcano Plot / Heatmap**：DEG結果の可視化
- **Pathway解析**：IFN関連遺伝子のエンリッチメント解析（Fisher正確確率検定）

## ファイル構成

| ファイル | 内容 |
|---------|------|
| [setup_guide_bulk_jp.md](setup_guide_bulk_jp.md) | 環境構築ガイド（Mac / Windows対応） |
| [day3_bulk_DEG.R](day3_bulk_DEG.R) | 解析スクリプト（RStudio用） |
| [day3_bulk_DEG.html](day3_bulk_DEG.html) | 解析レポート（HTML版、コード・結果・図を含む） |
| img/pca_concept.png | PCA概念説明図 |

## セットアップ

詳細は [setup_guide_bulk_jp.md](setup_guide_bulk_jp.md) を参照してください。

```bash
conda create -n bulkworkshop -c conda-forge r-base=4.4.3 -y
conda activate bulkworkshop
conda install -c bioconda -c conda-forge r-ggplot2 r-dplyr r-ggrepel r-gplots bioconductor-edger bioconductor-variancepartition r-rmarkdown r-knitr -y
```

## データ

解析に必要なデータは [Google Drive](https://drive.google.com/drive/folders/1uO1c8wSZVccdQdvWx8liTaJ-e_l1cddL?usp=drive_link) からダウンロードしてください（慶應アドレスでのログインが必要です）。

ダウンロード後、`~/Downloads/bioinformatics-course-keio-2025-day3/` に展開し、以下の構成になるようにしてください：

```
bioinformatics-course-keio-2025-day3/
├── day3_bulk_DEG.R
├── img/
│   └── pca_concept.png
└── data/
    └── bulk/
        ├── sample_meta.txt
        ├── IFNgenes100.txt
        └── Th1_count.txt
```

## 使用パッケージ

| パッケージ | 用途 |
|-----------|------|
| ggplot2 | 描画 |
| dplyr | データ操作 |
| edgeR | DEG解析 |
| ggrepel | ラベル描画 |
| gplots | ヒートマップ |
| variancePartition | 分散分解解析 |
| rmarkdown | HTMLレポート生成 |
| knitr | レポート生成エンジン |
