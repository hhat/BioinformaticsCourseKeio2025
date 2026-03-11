# Bulk RNA-seq DEG解析 ワークショップ

慶應義塾大学 微生物・免疫学実習 Day3 前半の教材です。

## 概要

末梢血免疫細胞のBulk RNA-seqデータを用いて、以下を実習します：

- **PCA**：CD4+ T細胞サブセットの主成分分析による可視化
- **DEG解析**：Th1細胞における健常者（HC）vs SLE患者の発現変動遺伝子の同定（edgeR）
- **Volcano Plot / Heatmap**：DEG結果の可視化
- **Pathway解析**：IFN関連遺伝子のエンリッチメント解析（Fisher正確確率検定）

## ファイル構成

| ファイル | 内容 |
|---------|------|
| [setup_guide_bulk_jp.md](setup_guide_bulk_jp.md) | 環境構築ガイド（Mac / Windows対応） |
| [day3_bulk_DEG.R](day3_bulk_DEG.R) | 解析スクリプト（RStudio用） |
| img/pca_concept.png | PCA概念説明図 |

## セットアップ

[setup_guide_bulk_jp.md](setup_guide_bulk_jp.md) を参照してください。

```bash
conda create -n bulkworkshop -c conda-forge r-base=4.3.3 -y
conda activate bulkworkshop
conda install -c bioconda -c conda-forge r-ggplot2 r-dplyr r-ggrepel r-ggsci r-gplots bioconductor-edger -y
```

## データ

解析に必要なデータは講義時に別途共有します。
