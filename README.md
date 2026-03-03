# Single-Cell-Learning
Learning scRNA-seq analysis with Seurat
# Single-Cell RNA-seq Analysis Learning Log
**Day 1: Data Loading and Quality Control (QC) using Seurat**

## 1. Background & Objective
小児神経科医として、将来的にGlut1欠損症（Glut1-DS）をはじめとする神経疾患の病態解明（iPS細胞由来脳オルガノイドやBBBモデルの解析）を行うことを見据え、シングルセルRNA-seq（scRNA-seq）のドライ解析スキルの習得を開始した。

初日の本日は、R環境における標準的ツールである `Seurat` パッケージの基礎的な使い方と、シングルセル解析において最も重要となる「細胞の品質評価（Quality Control: QC）」の可視化手法について学習した。

## 2. Dataset
* **Data:** `pbmc_small` (Seuratパッケージに内蔵されている軽量の末梢血単核球テストデータ)
* **Goal:** 巨大な公共データを扱う前の基礎固めとして、まずは数MBクラスのデータでSeuratオブジェクトの構造と可視化の基本を理解する。

## 3. Code & Explanation
今回実行したRスクリプトの解説は以下の通り。

```R
# ① パッケージの呼び出し
# シングルセル解析における事実上の標準（デファクトスタンダード）ツールであるSeuratを読み込む。
library(Seurat)

# ② データの読み込み
# Seuratに最初から組み込まれている練習用データ「pbmc_small」を呼び出す。
data("pbmc_small")

# ③ 品質評価（QC）グラフの出力
# VlnPlot関数を用いて、細胞の品質を評価するためのバイオリンプロットを描画する。
# - nFeature_RNA: 各細胞が発現している「遺伝子の種類数」
# - nCount_RNA: 各細胞が持っている「RNAの総量」
VlnPlot(pbmc_small, features = c("nFeature_RNA", "nCount_RNA"))
## 4. Biological Interpretation (生物学的な解釈)
出力されたバイオリンプロットの各ドットは「1つの細胞」を示している。実際の疾患モデル（脳オルガノイド等）の解析においては、このグラフを確認し、以下の基準で低品質な細胞（ゴミデータ）をフィルタリングして取り除く必要がある。

1. **極端に数値が高い細胞（上振れ）:** 2つ以上の細胞が誤って同じ液滴に入ってしまった「ダブレット（Doublet）」の可能性が高いため除外する。
2. **極端に数値が低い細胞（下振れ）:** 細胞膜が破綻し、中のRNAが漏れ出してしまった「死にかけている細胞」、あるいはRNAを含まない「空の液滴」である可能性が高いため除外する。

## 5. Next Steps
* クラウド環境（Google Colab等）を用いた大容量メモリの確保
* 実際の公開データ（GEO）を用いた、脳オルガノイドデータの読み込みとミトコンドリア遺伝子割合（percent.mt）の算出



**Day 2: Full scRNA-seq Pipeline (From Raw Data to UMAP & Annotation)**

## 1. Background & Objective
前回の基礎的な可視化に続き、今回はより実践的な解析フローを習得した。
無料クラウド環境（Google Colab: 12GB RAM）における巨大な結合データ（.rds）の展開によるクラッシュを経験したことから、完成品ではなく「シーケンス施設から納品される標準的な生データ形式（10x Genomicsのカウントマトリックス等）」から自力でSeuratオブジェクトを構築し、軽量かつ効率的に一連の標準パイプラインを走らせる手法へと切り替えた。
将来の疾患モデリング（患者由来iPS細胞を用いた脳オルガノイド等）を想定し、データの読み込みから細胞種のアノテーション（意味づけ）までを単独で完走することを目的とした。

## 2. Dataset
* **Data:** 3k PBMCs from a Healthy Donor (10x Genomics公式データ)
* **Input format:** `.mtx`, `barcodes.tsv`, `features.tsv` （シーケンス外注時の標準納品フォーマット）

## 3. Workflow & Code
生のカウントデータから細胞地図の作成まで、以下の標準パイプラインを実行した。

    # 1. 生データからのSeuratオブジェクト構築
    counts <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
    seurat_obj <- CreateSeuratObject(counts = counts, project = "PBMC3k", min.cells = 3, min.features = 200)

    # 2. 品質評価 (QC) とフィルタリング
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
    seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

    # 3. データの正規化と高変動遺伝子の抽出
    seurat_obj <- NormalizeData(seurat_obj)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

    # 4. 次元圧縮 (PCA) と次元数の見極め
    seurat_obj <- ScaleData(seurat_obj)
    seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

    # 5. クラスタリングとUMAP可視化
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

## 4. Biological Interpretation & Annotation
次元圧縮とクラスタリングの結果、2638個の細胞は9つの独立したクラスター（0〜8）に分類された。
各クラスターの生物学的な同一性を確認するため、既知のマーカー遺伝子の発現を地図上にマッピング（FeaturePlot）した。

* **T cells:** クラスター0, 2, 4において、T細胞マーカーである *CD3D* の特異的かつ強い発現を確認した。
* **B cells:** クラスター3において、B細胞マーカーである *MS4A1* (CD20) の特異的な発現を確認した。

将来的に神経疾患モデルを解析する際も、この手法を応用して *SOX2* (神経幹細胞) や *GFAP* (アストロサイト) 等のマーカーを評価し、細胞集団ごとの病態や代謝異常の差異を解明していく。
