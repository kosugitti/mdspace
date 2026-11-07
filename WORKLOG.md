# WORKLOG

## 2026-04-14

### やったこと

- Okada & Imaizumi (1987) “Nonmetric MDS of Asymmetric Proximities”
  (bhmk.14.21_81.pdf) を日本語翻訳し、lualatex
  でコンパイル（~/Desktop/bhmk_14_21_81_ja.tex / .pdf）
- `develop/2026MMDSchap8.R`: SMACOF
  アルゴリズムの実装（テキスト第8章に準拠）
  - B(Z) 行列の計算関数 `BZ()` を実装
  - V 行列の構成（式 8.17-8.18 の A 行列による累積）
  - ストレス関数 `sigma_function()` を実装（η²_δ + η²(X) - 2ρ の分解）
  - Guttman 変換による反復ループ（while + 収束判定 + 最大反復数）
  - テキスト Table 8.4 の結果（35反復で σ_r = 0.0174
    に収束）と一致を確認
- 実装中に発見・修正したバグ:
  - V 行列: ループ範囲 `j in 2:4` → `j in (i+1):n`、`V[i,j]` だけでなく
    `V <- V + W[i,j] * A` で行列全体を加算
  - ストレス計算: `sum(Delta^2)` →
    `sum(Delta[upper.tri(Delta)]^2)`（上三角のみ）
  - ρ の計算: `trace(B)` → `sum(diag(t(Z) %*% B %*% Z))`
  - ループ内で `sigma_old` の更新漏れ
  - `sigma_function` 内のグローバル変数参照 → 引数 `D` に統一

### 次回への引き継ぎ

- `classical_mds()` の実装（ダブルセンタリング + 固有値分解）→ SMACOF
  の初期値生成に必要
- `develop/2026MMDSchap8.R` をベースに `R/kruskal_mds.R`
  としてパッケージ関数化
  - S3 クラス `mdspace_kruskal` の設計（\$conf, \$stress, \$call,
    \$niter）
  - print/plot/summary メソッド
  - smacof パッケージとの照合テスト
- 実装順序: classical_mds() → kruskal_mds() → indscal_mds() → …

## 2026-11-07(土・午前)

- 設計: `develop/設計メモ_20261107.md`。既存 R(smacof・asymmetry 2.0.7
  など)でできることと自作するものの切り分け，工程(週末1コマずつ・非対称系から)，nASYMS(岡太・今泉の公式
  Shiny)のマニュアルの読み取り
- 決定(先生): 古典的MDS・Kruskal
  は先生が自分で書く／楕円は1関数で引数切り替え／CRAN
  は非対称系がそろってから検討／GSTATIS
  を工程に入れる／スライドベクトル・ASYMSCAL は asymmetry を包む
- 実装:
  [`hfm_mds()`](https://kosugitti.github.io/mdspace/reference/hfm_mds.md)(回転は
  std・varimax・none)・共通の骨組み・[`as_proximity()`](https://kosugitti.github.io/mdspace/reference/as_proximity.md)・[`asym_decomp()`](https://kosugitti.github.io/mdspace/reference/asym_decomp.md)。テスト
  24 件が通過。コミット 3f1e659(push はしていない)
- 2001年の VB 版
  GSTATIS(`Labo/Study：バックアップ/Chikage/非対称ＭＤＳスペシャル/GSTATIS.bas`)を発見。標準型の回転の向きが逆の疑い(設計メモ
  6)
- 確認用スクリプト `develop/20261107.R`(4x4
  の例・2019年版との一致・回転の比較・studentmigration・socioKG)
- DESCRIPTION の ORCID を
  0000-0001-5816-0099(先生の正しい番号)に訂正してコミット
- 次: 先生が `develop/20261107.R` で確かめる→
  GSTATIS([`gstatis_mds()`](https://kosugitti.github.io/mdspace/reference/gstatis_mds.md))か円モデル

## 2026-11-07(土・午前の続き)

- HFM の適合指標を Stress から Unexplained SS(1 −
  累積GOF)に変更(先生の指示)
- [`gstatis_mds()`](https://kosugitti.github.io/mdspace/reference/gstatis_mds.md)
  を実装し，千野先生の FORTRAN 版の出力(1999-12-30
  のメール・`Labo/Study：バックアップ/修士論文/Gstatis/GSTATISのoutput.txt`)と照合。出力は縦長の
  CSV 8本にした(`develop/references/chino1999_gstatis/`・git
  管理外)。入力 19 行列は socioKG と完全一致
- 原論文 Chino, Grorud & Yoshino「A revised version of GSTATIS」(BSJ
  抄録 pp.409-412)を発見し，第3段を改訂版(行ベクトルの射影
  Ψ^{-1/2}O\*h・負の固有値は i 倍)に合わせて直した。重みは GRV
  の固有ベクトルそのもの
- 照合の結果は全段一致(RV・妥協の固有値・軌跡 582 値・Lp ノルム 2,340
  値・距離)。`develop/20261107_chino_check.R`。未特定は欠測の埋め方だけ(設計メモ
  6)
- コミット 9c3d113〜3f47502(push は 9c3d113 以降まだ)
- 次: 欠測の扱い(千野版の再現か，mdspace
  独自の補完か)を先生と決める→円モデル
