# セットアップガイド：Bulk RNA-seq DEG解析環境の構築

このガイドでは、Bulk RNA-seq DEG（発現変動遺伝子）解析の実習に必要な環境を構築する手順を説明します。**macOS** と **Windows** の両方に対応しています。

> **⚠️ 重要：RStudioの起動方法について**
> 既にRStudioをお持ちの方も、本ガイドの手順に沿って **conda環境を有効化したターミナルから** RStudioを起動してください。
> 普段の起動方法（Applicationsやスタートメニューから開く等）では、必要なパッケージが使えない場合があります。

---

## 0. 必要なPCスペック

| 項目 | 推奨スペック |
|------|------------|
| OS | macOS 13以降 / Windows 10 (64-bit) 以降 |
| メモリ (RAM) | **8 GB以上** |
| ディスク空き容量 | **5 GB以上** |
| CPU | Intel / Apple Silicon (M1/M2/M3/M4) |

---

## 1. コマンドライン（ターミナル）の使い方

### macOS の場合

1. **Spotlight検索**（`Command ⌘` + `Space`）で「**ターミナル**」または「**Terminal**」と入力
2. 表示された **ターミナル.app** を開きます

> **💡 ヒント**：Launchpad → その他 → ターミナル からも開けます。

### Windows の場合

本ガイドでは **Miniforge Prompt** を使用します（セクション2でインストールします）。

Miniforgeインストール前に動作確認したい場合は、**Windows キー** を押して「**PowerShell**」と入力し、開くことができます。

> **⚠️ 注意（Windows）**：セクション2以降では、必ず **Miniforge Prompt** を使用してください。

### コマンドの実行方法

コマンドラインでは、プロンプト（`$` や `>` の後）にコマンドを入力し、**Enterキー** を押して実行します。
本ガイドのコードブロック内のテキストを**1行ずつ**コピー＆ペーストして実行してください。

---

## 2. Conda（Miniforge）のインストール

> **💡 既にcondaをインストール済みの方**（稲毛先生の実習等でインストール済みの場合）：この手順はスキップしてセクション3に進んでください。
> ターミナル / Miniforge Prompt で `conda --version` を実行し、バージョンが表示されればインストール済みです。

### macOS の場合

ターミナルを開き、以下を実行します。

#### Apple Silicon (M1/M2/M3/M4) Mac の場合

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh"
bash Miniforge3-MacOSX-arm64.sh
```

#### Intel Mac の場合

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh"
bash Miniforge3-MacOSX-x86_64.sh
```

インストール後、**ターミナルを閉じて再度開いてください**。
プロンプトの先頭に `(base)` と表示されていればインストール成功です。

### Windows の場合

1. [Miniforgeダウンロードページ](https://github.com/conda-forge/miniforge/releases/latest) にアクセスします
2. **`Miniforge3-Windows-x86_64.exe`** をダウンロードして実行します
3. インストール途中の設定画面で：
   - **「Add Miniforge3 to my PATH environment variable」** に ☑️ チェックを入れてください
   - **「Register Miniforge3 as my default Python」** にもチェックを入れてください
4. インストール完了後、**スタートメニュー** で「**Miniforge Prompt**」と検索して開きます

---

## 3. Conda環境の作成とパッケージのインストール

以下のコマンドを **1行ずつ** ターミナル (macOS) または Miniforge Prompt (Windows) に入力して実行してください。

### Step A. 環境の作成（R のみ）

```bash
conda create -n bulkworkshop -c conda-forge r-base=4.4.3 -y
```

> **⏱ 注意：初回のインストールには10〜20分程度かかることがあります。**

> **⚠️ バージョン指定でエラーが出る場合**：`=4.4.3` を削除して再実行してください。
> ```bash
> conda create -n bulkworkshop -c conda-forge r-base -y
> ```

### Step B. 環境の有効化

```bash
conda activate bulkworkshop
```

プロンプトの先頭が `(bulkworkshop)` に変わったことを確認してください。

> **⚠️ 重要**：ターミナル / Miniforge Prompt を開くたびに `conda activate bulkworkshop` を実行する必要があります。

### Step C. パッケージのインストール

```bash
conda install -c bioconda -c conda-forge \
  r-ggplot2 \
  r-dplyr \
  r-ggrepel \
  r-gplots \
  bioconductor-edger \
  bioconductor-variancepartition \
  -y
```

> **💡 Windows での注意**：`\`（行の継続記号）が動作しない場合は、すべてを **1行** にまとめて実行してください：
> ```bash
> conda install -c bioconda -c conda-forge r-ggplot2 r-dplyr r-ggrepel r-gplots bioconductor-edger bioconductor-variancepartition -y
> ```

> **💡 解説**：edgeR は Bioconductor のパッケージで、DEG解析の中核を担います。依存パッケージの limma も自動的にインストールされます。

---

## 4. RStudio のインストールと起動

### RStudio のインストール

RStudio が未インストールの場合は、以下からインストールしてください（既にインストール済みの方はスキップ）。

- [posit.co/download/rstudio-desktop/](https://posit.co/download/rstudio-desktop/) から無料版をダウンロード・インストール

> **💡 注意**：RStudio は conda ではなく、通常のアプリケーションとしてインストールします。

### RStudio の起動方法

**どのような方法でRStudioをインストールした場合でも**、以下の手順でconda環境を有効化してからRStudioを起動してください。

#### macOS の場合

ターミナルを開き、以下を実行します：

```bash
conda activate bulkworkshop
export RSTUDIO_WHICH_R=$(which R)
open -a RStudio
```

> **💡 `open -a RStudio` でエラーが出る場合**：RStudio のアプリ名が異なる可能性があります。
> 以下を試してください：
> ```bash
> # アプリを探す
> ls /Applications/ | grep -i rstudio
> # 見つかったアプリ名で起動（例）
> open -a "RStudio"
> ```
> それでもダメな場合は、ターミナルから直接パスを指定して起動できます：
> ```bash
> /Applications/RStudio.app/Contents/MacOS/RStudio &
> ```

#### Windows の場合

Miniforge Prompt を開き、**毎回以下を実行**してください：

```bash
conda activate bulkworkshop
set RSTUDIO_WHICH_R=%CONDA_PREFIX%\Scripts\R.exe
"C:\Program Files\RStudio\rstudio.exe"
```

> **⚠️ パスでエラーが出る場合**：以下も試してください：
> ```bash
> "C:\Program Files\RStudio\bin\rstudio.exe"
> ```
> それでも見つからない場合は、スタートメニューで「RStudio」を右クリック →「ファイルの場所を開く」でパスを確認してください。

#### 起動後の確認

RStudio が起動したら、コンソール（左下パネル）で以下を実行してください：

```r
R.home()
```

出力に `miniforge3` や `bulkworkshop` を含むパスが表示されていれば正しく設定されています。

- OK の例：`/Users/xxx/miniforge3/envs/bulkworkshop/lib/R`
- NG の例：`/Library/Frameworks/R.framework/...` や `C:\Program Files\R\...`

> **⚠️ NG の場合**：RStudio を一度閉じて、上記の手順（conda activate → RSTUDIO_WHICH_R設定 → 起動）をやり直してください。
> デスクトップやスタートメニューのショートカットから直接RStudioを開くと、conda環境のRが使われません。

#### （代替案）RStudioなしでRを直接使う

どうしても RStudio がうまくいかない場合は、ターミナル / Miniforge Prompt から R を直接起動できます：

```bash
conda activate bulkworkshop
R
```

---

## 5. データの準備

### A. プロジェクトフォルダの作成

#### macOS の場合（ターミナル）
```bash
cd ~/Downloads
mkdir bioinformatics-course-keio-2025-day3
cd bioinformatics-course-keio-2025-day3
mkdir -p data/bulk
mkdir img
```

#### Windows の場合（Miniforge Prompt）
```bash
cd %USERPROFILE%\Downloads
mkdir bioinformatics-course-keio-2025-day3
cd bioinformatics-course-keio-2025-day3
mkdir data\bulk
mkdir img
```

#### フォルダ構成

```
bioinformatics-course-keio-2025-day3/
├── day3_bulk_DEG.R                    (解析スクリプト)
├── img/
│   └── pca_concept.png                (PCA概念図)
└── data/
    └── bulk/
        ├── sample_meta.txt            (サンプルメタデータ)
        ├── IFNgenes100.txt            (IFN関連遺伝子リスト)
        └── Th1_count.txt              (Th1カウントデータ)
```

### B. データのダウンロード

本実習では、公開データベースから取得した末梢血免疫細胞のBulk RNA-seqデータ（健常者およびSLE患者）を使用します。

以下のリンクからデータをダウンロードし、上記のフォルダ構成に従って配置してください：

- **データ一式**：[Google Drive](https://drive.google.com/drive/folders/1uO1c8wSZVccdQdvWx8liTaJ-e_l1cddL?usp=drive_link)（慶應アドレスでのログインが必要です）

> **💡 解説**：本実習で使用するデータは、健常者（HC）および全身性エリテマトーデス（SLE）患者の末梢血から25種類の免疫細胞をFACSで分離し、RNA-seqを行ったものです。DEG解析ではTh1細胞のHC vs SLEの比較を行います。

---

## 6. インストールの確認

RStudioを起動（セクション4参照）し、コンソールで以下のスクリプトを実行してください。

```r
cat("=== R バージョン ===\n")
cat(R.version.string, "\n\n")

cat("=== パッケージバージョン確認 ===\n")

packages <- c("ggplot2", "dplyr", "edgeR", "ggrepel", "gplots", "limma", "variancePartition")

for (pkg in packages) {
  tryCatch({
    library(pkg, character.only = TRUE)
    cat(sprintf("  OK  %-12s : %s\n", pkg, packageVersion(pkg)))
  }, error = function(e) {
    cat(sprintf("  NG  %-12s : インストールされていません！\n", pkg))
  })
}

cat("\n=== 結果 ===\n")
cat("上記のパッケージがすべて OK であれば、準備完了です！\n")
```

### 期待される出力（参考）

| パッケージ | 説明 |
|-----------|------|
| R         | 4.4.x |
| ggplot2   | 描画 |
| dplyr     | データ操作 |
| edgeR     | DEG解析 |
| ggrepel   | ラベル描画 |
| gplots    | ヒートマップ |
| limma     | edgeRの依存パッケージ |
| variancePartition | 分散分解解析 |

---

## 7. トラブルシューティング

### Q1. `conda` コマンドが見つからない（"command not found"）

**macOS の場合：**
ターミナルを閉じて再度開いてください。それでも解決しない場合：
```bash
~/miniforge3/bin/conda init
```
その後、ターミナルを再起動してください。

**Windows の場合：**
通常のコマンドプロンプトではなく、**Miniforge Prompt** を使用してください。

### Q2. `conda activate bulkworkshop` が動作しない

```bash
conda init
```
を実行後、ターミナル / Miniforge Prompt を**閉じて再度開いて**ください。

### Q3. RStudioで「conda環境のR」が使われていない

RStudioのコンソールで `R.home()` を実行し、パスを確認してください。

- `miniforge3` や `bulkworkshop` を含む → OK
- `/Library/Frameworks/R.framework/...` や `C:\Program Files\R\...` → conda環境のRではない

**対処法**：RStudioを閉じて、ターミナル / Miniforge Prompt から以下の手順でやり直してください：
1. `conda activate bulkworkshop`
2. `export RSTUDIO_WHICH_R=$(which R)` (macOS) / `set RSTUDIO_WHICH_R=%CONDA_PREFIX%\Scripts\R.exe` (Windows)
3. RStudio を起動

### Q4. パッケージのインストールでエラーが出る

conda での解決が難しい場合は、R内から直接インストールを試みてください：

```bash
conda activate bulkworkshop
R
```

R起動後：

```r
# CRANパッケージの場合
install.packages(c("ggplot2", "dplyr", "ggrepel", "gplots"))

# edgeR（Bioconductor）の場合
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("edgeR", "variancePartition"))
```

**macOS でコンパイルエラーが出る場合：**
```bash
xcode-select --install
```

**Windows でコンパイルエラーが出る場合：**
```bash
conda install -c conda-forge m2w64-toolchain -y
```

### Q5. メモリ不足エラーが発生する

- 不要なアプリケーションを閉じる
- RStudioの Environment タブで不要なオブジェクトを削除（`rm(object_name); gc()`）

---

## 8. 代替手段：Condaを使わない環境構築

どうしてもCondaでの環境構築がうまくいかない場合の手順です。

### A. Rのインストール
1. [CRAN](https://cran.r-project.org/) からお使いのOSに合わせてRをダウンロード・インストールします

### B. RStudio のインストール
1. [posit.co/download/rstudio-desktop/](https://posit.co/download/rstudio-desktop/) から無料版をダウンロード・インストールします

### C. パッケージの手動インストール

RStudioを起動し、コンソールで以下を実行します：

```r
# CRANパッケージ
install.packages(c("ggplot2", "dplyr", "ggrepel", "gplots"))

# Bioconductorパッケージ（edgeR）
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(c("edgeR", "variancePartition"))
```

インストール後、セクション6の確認スクリプトを実行してください。
