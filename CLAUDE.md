# mdspace パッケージ開発ルール

## 概要

MDS (多次元尺度構成法) の各種アルゴリズム実装 + Shiny GUI パッケージ。
科研費プロジェクト 24K06481 の成果物。CRAN公開を目指す。

## コーディング規約

### 命名規則

- 関数名: snake_case
- MDS手法の関数: {手法名}\_mds()（例: classical_mds(), kruskal_mds()）
- S3クラス名: “mdspace\_{手法}”（例: “mdspace_classical”）
- 内部関数: .で始める（例: .majorize_step()）

### コードスタイル

- R/ 以下（パッケージ本体）: base R のみ。tidyverse 依存禁止
- develop/ 以下（開発用スクリプト）: tidyverse 使用可
- Rcpp は使わない（純R実装）

### 返り値の共通構造

全MDS関数は以下のスロットを持つS3オブジェクト（list）を返す: - \$conf:
布置行列（n x ndim） - \$stress: ストレス値 - \$call: 関数呼び出し -
\$niter: 反復回数（反復法の場合）

### S3メソッド

全クラスに print(), plot(), summary() を実装する。

### ドキュメント

- roxygen2 使用
- 各関数に @references を必ず記載

### テスト

- testthat 3rd edition
- smacof の結果との照合テストを含める

### Shiny

- golem 不使用。素の Shiny + モジュール化
- inst/shiny/ に配置。run_app() で起動
- 外部パッケージ連携（smacof, exametrika 等）は requireNamespace()
  で判定

## 実装予定モデル

### 基礎モデル（テキスト基礎編対応）

| 関数名          | モデル                            | 参考文献                           | 状態   |
|-----------------|-----------------------------------|------------------------------------|--------|
| classical_mds() | 古典的MDS（Torgerson法）          | Torgerson (1952)                   | 未実装 |
| kruskal_mds()   | 非計量MDS（Kruskal法）            | Kruskal (1964)                     | 未実装 |
| indscal_mds()   | 個人差MDS（INDSCAL）              | Carroll & Chang (1970)             | 未実装 |
| prefmap_mds()   | 選好マッピング（PREFMAP）         | Carroll (1972)                     | 未実装 |
| dual_scaling()  | 双対尺度法 / コレスポンデンス分析 | Nishisato (1980) / Benzécri (1973) | 未実装 |
| mca()           | 多重対応分析                      | Greenacre (1984)                   | 未実装 |
| som()           | 自己組織化写像（SOM）             | Kohonen (1982)                     | 未実装 |

### 発展モデル（テキスト発展編対応）

| 関数名         | モデル                  | 参考文献                      | 状態   |
|----------------|-------------------------|-------------------------------|--------|
| hfm_mds()      | 階層的因子モデル（HFM） | 千野 (Chino)                  | 未実装 |
| dynascal_mds() | DYNASCAL                | 千野 (Chino)                  | 未実装 |
| ellipse1_mds() | 楕円モデル I            | 岡太・今泉 (Okada & Imaizumi) | 未実装 |
| ellipse2_mds() | 楕円モデル II           | 岡太・今泉 (Okada & Imaizumi) | 未実装 |
| ellipse3_mds() | 楕円モデル III          | 岡太・今泉 (Okada & Imaizumi) | 未実装 |
| vonmises_mds() | フォンミーゼスモデル    | 荘島 (Shojima)                | 未実装 |

### Shiny連携（外部パッケージ経由）

| パッケージ | 用途                                                 | Shiny内で利用 |
|------------|------------------------------------------------------|---------------|
| smacof     | スライドベクトルモデル等、mdspaceで自作しないMDS手法 | オプション    |
| exametrika | Biclustering                                         | オプション    |
| igraph     | ネットワーク分析・可視化                             | オプション    |
| kohonen    | SOM（参照実装との比較用）                            | オプション    |

※ 今後追加の可能性あり ※
双対尺度法とコレスポンデンス分析は数学的に同等（表現形が異なる）

## ディレクトリ構成

- R/: パッケージ関数（base R のみ）
- inst/shiny/: Shiny アプリ本体
- inst/shiny/modules/: Shiny モジュール（mod\_\*.R）
- tests/testthat/: テスト
- develop/: 開発用スクリプト（tidyverse 可）
- man/: roxygen2 生成ドキュメント
