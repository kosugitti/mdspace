# テスト用の照合データ(円モデル・circle_mds())

- `okada_imaizumi_2024_table1_19.csv`: Okada & Imaizumi (2024) 表1.19，10職業間の類似度(行から列へ・対角は NA)。nASYMS マニュアル(`develop/references/nASYMS/nASYMS_user_manual.pdf`)7ページの Table 1(書籍表1.19をスプレッドシート形式にしたもの)から 2026-11-17 に読み取った。nASYMS に読み込ませた Excel は `develop/references/nASYMS/data_occupation_NA.xlsx`。
- `nasyms_radius_occupation_all.csv`: 上のデータを nASYMS(https://nasyms.shinyapps.io/nASYMS-radius/)の円モデルで 2026-11-17 に解析した出力。1〜5次元の布置と半径(小数3桁・画面の表から取得)。設定はマニュアルと同じ(最大5次元・反復200回・初期布置 Torgerson(rational)・Data)。
- nASYMS の最小化された S: 1次元 0.499，2次元 0.372，3次元 0.276，4次元 0.202，5次元 0.183(すべて Converged)。マニュアルの Figure 7 は5次元が 0.182 で，今回の実行と小数3桁目が1違う。2次元の Primary の D2 もマニュアルは -0.147，今回は -0.146。
- 書籍の数値は複数の初期値から最良のものを選ぶ手順(書籍 4.1.1)なので，nASYMS の1回の解析とは一致しないことがある。

# 二相三元の円モデル(未実装)

- `okada_imaizumi_2024_table1_14.csv`: Okada & Imaizumi (2024) 表1.14，アジア10か国間の類似度(行から列へ)を2つのソースについて縦に積んだもの(列 `source`・`from`)。nASYMS マニュアル45ページの Table 7 から 2026-11-17 に読み取った。マニュアルの例(セル E13＝ソース2の Indonesia→Malaysia が 3.16)と一致を確認
- nASYMS(https://nasyms.shinyapps.io/nASYMS-twomode/)で 2026-11-17 に解析(最大5次元・反復200回・Torgerson(rational)・Data＝マニュアルと同じ設定)。Excel は `develop/references/nASYMS/data_nations_2sources_NA.xlsx`
  - `nasyms_twomode_nations_weights.csv`: 1〜5次元の S(全体と各ソース)，非対称重み u_it(Asymmetry weight)と対称重み w_i(Symmetry weight)
  - `nasyms_twomode_nations_conf.csv`: 共通布置と半径。2〜5次元。**1次元は nASYMS の画面が "Error: An error has occurred" になり表示されない**(アプリ側の不具合)
  - 同じ設定で2回実行して S・2次元の布置・重みが完全に一致(再現性あり)
  - 式の確認: 2次元の解から m_jki = w_i d_jk − v_jki (r_j − r_k)，v_jki = w_i d_jk / √Σ_t((x_jt − x_kt)/u_it)² で計算すると，nASYMS の表示する m(11セルで照合)と丸めの差 0.001 以内で一致
  - 画面の表は Shiny の描画に時間がかかり，次元を切り替えた直後は前の次元の表が残る。取り直すときは1次元ずつ十分に待つこと
