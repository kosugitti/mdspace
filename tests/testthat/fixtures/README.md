# テスト用の照合データ(円モデル・circle_mds())

- `okada_imaizumi_2024_table1_19.csv`: Okada & Imaizumi (2024) 表1.19，10職業間の類似度(行から列へ・対角は NA)。nASYMS マニュアル(`develop/references/nASYMS/nASYMS_user_manual.pdf`)7ページの Table 1(書籍表1.19をスプレッドシート形式にしたもの)から 2026-11-17 に読み取った。nASYMS に読み込ませた Excel は `develop/references/nASYMS/data_occupation_NA.xlsx`。
- `nasyms_radius_occupation_all.csv`: 上のデータを nASYMS(https://nasyms.shinyapps.io/nASYMS-radius/)の円モデルで 2026-11-17 に解析した出力。1〜5次元の布置と半径(小数3桁・画面の表から取得)。設定はマニュアルと同じ(最大5次元・反復200回・初期布置 Torgerson(rational)・Data)。
- nASYMS の最小化された S: 1次元 0.499，2次元 0.372，3次元 0.276，4次元 0.202，5次元 0.183(すべて Converged)。マニュアルの Figure 7 は5次元が 0.182 で，今回の実行と小数3桁目が1違う。2次元の Primary の D2 もマニュアルは -0.147，今回は -0.146。
- 書籍の数値は複数の初期値から最良のものを選ぶ手順(書籍 4.1.1)なので，nASYMS の1回の解析とは一致しないことがある。
