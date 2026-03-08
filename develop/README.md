# develop/ — 開発用スクリプト・参考資料

mdspaceパッケージの実装にあたり、過去のリポジトリから収集した種コード・参考資料。
ここのコードはそのままパッケージに入るものではなく、アルゴリズムの参考・検証用。

## ファイル一覧と出自

### PMDS（~/Dropbox/Git/PMDS/）由来

「パソコン多次元尺度構成法」（岡田・今泉, 1994, 共立出版）のR移植を目指したパッケージ。
mdspaceに吸収されたため、元リポジトリは削除済み。

| ファイル | 元パス | 内容 | 作成時期 |
|---------|--------|------|---------|
| smacof20240209.R | develop/smacof20240209.R | SMACOFアルゴリズム実装（B行列・ストレス関数・ダブルセンタリング）。classical_mds()の種 | 2024-02 |
| kruskal_dev.R | develop/kruskal_dev.R | Kruskal NMDSの勾配計算実装。kruskal_mds()の種 | 2023-12 |
| PCKRUS_dev.R | develop/PCKRUS_dev.R | kruskal_dev.Rの改良版。コメント追加、複数の勾配計算方法の比較 | 2023-12 |
| BF02289565.pdf | develop/BF02289565.pdf | 参考論文PDF | 2023-11 |
| BF02289694.pdf | develop/BF02289694.pdf | 参考論文PDF | 2023-11 |

### BMDS（~/Dropbox/Git/BMDS/）由来

Bayesian MDSパッケージの試験的実装（Stan使用）。2016年で開発停止。
mdspaceへの参考コード抽出後、元リポジトリは削除済み。

| ファイル | 元パス | 内容 | 作成時期 |
|---------|--------|------|---------|
| bmds_example_code.R | develop/example_code.R | Bayesian古典的MDSの実装例。Procrustes回転・反射補正ロジック（54-93行）が特に有用 | 2016-03 |
| bmds_vonmises.R | R/20190126.R | von Mises（方向統計）モデルの実装試験。circularパッケージ使用。vonmises_mds()の種 | 2019-01 |

### JSCP_MDS_2019（~/Dropbox/Git/JSCP_MDS_2019/）由来

日本認知心理学会・研究法研究部会 第1回研究会（2019-03-03）「MDS、使ってる？」の教材。
元リポジトリは教材サイトとして独立価値があるため残存。

| ファイル | 元パス | 内容 | 作成時期 |
|---------|--------|------|---------|
| jscp_lesson2.Rmd | lesson2.Rmd | INDSCAL（smacof使用）、Prefmap（手動実装）、Abelsonマップ（自作関数）、HFM（複素固有値分解による非対称MDS自作実装）の教材。hfm_mds(), prefmap_mds()の種 | 2019-03 |

### imKmeans（~/Dropbox/Git/imKmeans/）由来

マハラノビス距離を用いた改良K-means法。通常のK-means（ユークリッド距離）の後、
マハラノビス距離ベースの尤度で再割り当てを反復し、AIC/BICでクラスタ数を評価する。
2016年作成。元リポジトリは削除済み。

| ファイル | 元パス | 内容 | 作成時期 |
|---------|--------|------|---------|
| imKmeans.R | R/imKmeans.R | マハラノビス距離K-means。mvtnormで混合正規尤度を計算しAIC/BIC算出。SOM・クラスタリング関連の参考 | 2016-02 |

### OrdinalScale（~/Dropbox/Git/OrdinalScale/）由来

順序尺度データを適切に扱うためのパッケージ試作。双対尺度法（Dual Scaling）の
固有値分解法・相互平均法・順序データ版の3実装を含む。Nishisato (1980) に基づく。
2019年作成。元リポジトリは削除済み。

| ファイル | 元パス | 内容 | 作成時期 |
|---------|--------|------|---------|
| os_DualScaling.R | R/DualScaling.R | 双対尺度法（固有値分解法）。クロス表の行・列スコア、特異値、寄与率を算出。dual_scaling()の主要な種 | 2019-07 |
| os_DualScalingMA.R | R/DualScalingMA.R | 双対尺度法（相互平均法, Mutual Averages）。反復アルゴリズムによる解法 | 2019-07 |
| os_DualScalingOrdinal.R | R/DualScalingOrdinal.R | 順序データ用双対尺度法。支配行列（dominance matrix）経由で順序尺度を処理 | 2019-07 |
| os_make_dominance.R | R/make_dominance.R | 支配行列の生成関数。DualScalingOrdinalの前処理 | 2019-07 |
| os_make_dummy.R | R/make_dummy.R | ダミー変数行列の生成関数 | 2019-07 |
| os_developing.R | R/developing.R | 開発用スクリプト（テスト・検証コード） | 2019-07 |
| os_CE.csv | R/CE.csv | テスト用データ | 2019-07 |

---

## 注意事項

- R/（パッケージ本体）はbase Rのみだが、develop/ではtidyverse使用可
- 古いコードにはtidyverse依存が混在しているので、パッケージ化時に除去すること
- BMDSのStan依存コードはアルゴリズムの参考のみ（mdspaceはStan不使用）
- smacofパッケージの結果との照合テストを必ず実施すること
