# mdspace パッケージ開発ルール

## 概要

MDS (多次元尺度構成法) の各種アルゴリズム実装 + Shiny GUI パッケージ。
科研費プロジェクト 24K06481 の成果物。CRAN公開を目指す。

### 現在の状態（2026-04-14更新）

- SMACOF アルゴリズム（ratio MDS）の開発用コードが
  `develop/2026MMDSchap8.R` で動作確認済み
- 次: `classical_mds()` を R/ に実装 → `kruskal_mds()`
  のパッケージ関数化

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

| 関数名 | モデル | 参考文献 | 状態 |
|----|----|----|----|
| classical_mds() | 古典的MDS（Torgerson法） | Torgerson (1952) | 未実装 |
| kruskal_mds() | 非計量MDS（Kruskal法） | Kruskal (1964) | 未実装 |
| indscal_mds() | 個人差MDS（INDSCAL） | Carroll & Chang (1970) | 未実装 |
| prefmap_mds() | 選好マッピング（PREFMAP） | Carroll (1972) | 未実装 |
| dual_scaling() | 双対尺度法 / コレスポンデンス分析 | Nishisato (1980) / Benzécri (1973) | 未実装 |
| mca() | 多重対応分析 | Greenacre (1984) | 未実装 |
| som() | 自己組織化写像（SOM） | Kohonen (1982) | 未実装 |

### 発展モデル（テキスト発展編対応）

| 関数名 | モデル | 参考文献 | 状態 |
|----|----|----|----|
| hfm_mds() | 階層的因子モデル（HFM） | 千野 (Chino) | 未実装 |
| dynascal_mds() | DYNASCAL | 千野 (Chino) | 未実装 |
| ellipse1_mds() | 楕円モデル I | 岡太・今泉 (Okada & Imaizumi) | 未実装 |
| ellipse2_mds() | 楕円モデル II | 岡太・今泉 (Okada & Imaizumi) | 未実装 |
| ellipse3_mds() | 楕円モデル III | 岡太・今泉 (Okada & Imaizumi) | 未実装 |
| vonmises_mds() | フォンミーゼスモデル | 荘島 (Shojima) | 未実装 |
| contour_mds() | 等高線マップモデル | Abelson (1954-55); 小杉・藤原 (2004) | 未実装 |

### Shiny連携（外部パッケージ経由）

| パッケージ | 用途 | Shiny内で利用 |
|----|----|----|
| smacof | スライドベクトルモデル等、mdspaceで自作しないMDS手法 | オプション |
| exametrika | Biclustering | オプション |
| igraph | ネットワーク分析・可視化 | オプション |
| kohonen | SOM（参照実装との比較用） | オプション |

※ 今後追加の可能性あり ※
双対尺度法とコレスポンデンス分析は数学的に同等（表現形が異なる）

## 参考文献・リソースメモ

### 千野の非対称MDSコード（develop/に格納済み、2026-03-08）

`~/Dropbox/Labo/Archive：研究会資料/CognitiveMiser/KAKEN24730510/R/`
からコピー:

| ファイル | 内容 |
|----|----|
| hfm.R | HFM（Hermitian Form Model）— エルミート形式による非対称性分解。固有値分解ベース |
| hfm.plot.R | HFMの布置プロット |
| gstatis.R | GSTATIS（一般化STATIS） |
| Abelson.map.R | Abelsonのマップ |
| rotConf.R | 布置の回転 |
| power.R | べき乗関連ユーティリティ |

### 岡太・今泉の円モデル・楕円モデル

#### 主要文献

- **円モデル（radius-distance model）初出**: Okada, A. & Imaizumi, T.
  (1987). Nonmetric multidimensional scaling of asymmetric proximities.
  *Behaviormetrika*, 14(21), 81-96. DOI: 10.2333/bhmk.14.21_81 —
  **PDF所持**
- **楕円モデル（ellipse-distance model）初出**: Okada, A. (1988).
  Asymmetric multidimensional scaling of car switching data. In Gaul &
  Schader (Eds.), *Data, Expert Knowledge and Decisions*, pp.279-290.
  Springer. — Springer有料
- **楕円モデル一般化**: Okada, A. (1990). A generalization of asymmetric
  multidimensional scaling. In Schader & Gaul (Eds.), *Knowledge, Data
  and Computer-Assisted Decisions*, NATO ASI Series vol.61, Springer. —
  Springer有料
- **書籍（包括的）**: Okada, A. & Imaizumi, T. (2024). *Applied
  Multidimensional Scaling of Asymmetric Relationships*. Springer,
  Behaviormetrics vol.19, 193pp. — **所持**

#### モデルの数式（判明分）

**円モデル（Okada & Imaizumi, 1987）**: - 各対象 i に座標
x_i（p次元）と半径 r_i \>= 0 を割り当てる - 非対称距離: g_jk = d_jk -
r_j + r_k（d_jk = \|\|x_j - x_k\|\| ユークリッド距離） - g_jk ≠ g_kj
なので非対称を表現 - 目的関数: STRESS最小化（非計量: 単調変換つき） -
MAXSCAL（Takane, 1981）ベースのアルゴリズム（ASYMMAXSCAL-OI）

**楕円モデル（Okada, 1988 / 1990）**: - 円モデルのスカラー半径 r_i
を次元ごとの半軸 a_ik に一般化 -
各対象を「点＋楕円（超楕円体）」で表現 -
半軸の長さが次元ごとの非対称性（支配性・魅力度）を表す - ellipse1/2/3
の区別の詳細は書籍（Okada & Imaizumi, 2024）を参照のこと

### 等高線マップモデル（Contour Map Model）

#### 文献

- **原著**: Abelson, R.P. (1954-55). A technique and a model for
  multi-dimensional scaling. *Public Opinion Quarterly, Winter*,
  405-418.
- **応用論文**: 小杉考司・藤原武弘 (2004).
  等高線マッピングによる態度布置モデル. *行動計量学*, 31(1), 17-24. —
  **PDF所持**（`~/Dropbox/Labo/MyLibrary/MyArticles/2004等高線.pdf`）

#### モデルの数式

MDS布置上の任意の座標 P における誘因価（valence）:

V(P) = Σ\_{j=1}^{n} V(j) / (1 + d²\_{Pj})

- V(j): 対象 j に対する評価値（好意度等）
- d\_{Pj}: 座標 P と対象 j の MDS 布置上のユークリッド距離
- 電位のメタファー: 力の伝播は距離に反比例

#### 特徴

- MDS布置（類似性の認知的成分）に好意度（感情的成分）を重ねて等高線地図として表現
- PREFMAPと異なり、好意度を距離ではなく力場の強度として表現 →
  当てはまりが良い
- 等高線の勾配が心理的緊張の指標、複数の丘は態度の非統合を示す
- 3次元曲面（誘因価曲面）としても可視化可能
- 入力: (1) MDS布置（類似性データから）、(2) 各個人の各対象への評定値

#### その他の関連Behaviormetrika論文（応用論文）

- Okada & Imaizumi (2012). Asymmetric MDS of brand switching among
  margarine brands. *Behaviormetrika*, 39(1), 111-. — SVDベースの応用
- Okada (2012). A brief survey of asymmetric MDS and some open problems.
  *Behaviormetrika*, 39(1), 127-. — サーベイ、円・楕円含む概観
- Okada & Imaizumi (2016). Asymmetric MDS of N-mode M-way categorical
  data using a log-linear model. *Behaviormetrika*, 43(1), 103-.

#### J-STAGE公開の大会抄録（develop/にPDF保存済み）

- 岡太 (2018) 対称関係から非対称関係へ（招待講演）
- 岡太・今泉 (2018) 単相2元非対称MDS比較（距離モデルの比較）
- 今泉・岡太 (2018) Imaizumi slide vector model
- 岡太・今泉 (2006) 非対称MDS外部分析（行動計量学 33(2)）
- 岡太・今泉 (2021) 非対称性の評価 I / II
- 岡太・今泉 (2023) Two-mode three-way I / II

## ディレクトリ構成

- R/: パッケージ関数（base R のみ）
- inst/shiny/: Shiny アプリ本体
- inst/shiny/modules/: Shiny モジュール（mod\_\*.R）
- tests/testthat/: テスト
- develop/: 開発用スクリプト（tidyverse 可）
- man/: roxygen2 生成ドキュメント
