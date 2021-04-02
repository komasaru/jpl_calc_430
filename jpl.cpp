#include "jpl.hpp"

namespace jpl_calc_430 {

// 定数
static constexpr char         kFBin[]    = "JPLEPH";        // バイナリファイル名
static constexpr unsigned int kKsize     = 2036;            // KSIZE
static constexpr unsigned int kRecl      =    4;            // 1レコード = KSIZE * 4
                                                            // ヘッダは2レコードで構成
static constexpr unsigned int kPosTtl    =    0;            // 位置: TTL
static constexpr unsigned int kPosCnam   =  252;            // 位置: CNAM
static constexpr unsigned int kPosSs     = 2652;            // 位置: SS
static constexpr unsigned int kPosNcon   = 2676;            // 位置: NCON
static constexpr unsigned int kPosAu     = 2680;            // 位置: AU
static constexpr unsigned int kPosEmrat  = 2688;            // 位置: EMRAT
static constexpr unsigned int kPosIpt    = 2696;            // 位置: IPT(IPT13以外)
static constexpr unsigned int kPosNumde  = 2840;            // 位置: NUMDE
static constexpr unsigned int kPosIpt2   = 2844;            // 位置: IPT13
static constexpr unsigned int kPosCnam2  = 2856;            // 位置: CNAM2(CNAMの続き)
static constexpr unsigned int kPosCval   = kKsize * kRecl;  // 位置: CVAL
static constexpr unsigned int kReclTtl   =   84;            // レコード長: TTL
static constexpr unsigned int kReclCnam  =    6;            // レコード長: CNAM
static constexpr unsigned int kReclSs    =    8;            // レコード長: SS
static constexpr unsigned int kReclNcon  =    4;            // レコード長: NCON
static constexpr unsigned int kReclAu    =    8;            // レコード長: AU
static constexpr unsigned int kReclEmrat =    8;            // レコード長: EMRAT
static constexpr unsigned int kReclIpt   =    4;            // レコード長: IPT
static constexpr unsigned int kReclNumde =    4;            // レコード長: IPT
static constexpr unsigned int kReclCval  =    8;            // レコード長: CVAL
static constexpr unsigned int kCntTtl    =    3;            // 件数: TTL
static constexpr unsigned int kCntCnam   =  400;            // 件数: CNAM
static constexpr unsigned int kCntSs     =    3;            // 件数: SS
static constexpr unsigned int kCntIpt    =   12;            // 件数: IPT（IPT13以外）
static constexpr unsigned int kKind      =    1;            // 計算区分
                                                            //   0: 計算しない
                                                            //   1: 位置・速度を計算
static constexpr double       kSecDay    = 86400.0;         // Seconds in a day

/*
 * @brief      コンストラクタ
 *
 * @param[in]  ユリウス日 (double)
 * @param[in]  単位フラグ (bool; optional)
 *             (true: km, km/sec, false: AU, AU/day)
 * @param[in]  基準フラグ (bool; optional)
 *             (true: 太陽系重心が基準, false: 太陽が基準)
 */
Jpl::Jpl(double jd, const bool is_km, const bool is_bary) {
  this->jd      = jd;
  this->is_km   = is_km;
  this->is_bary = is_bary;
  unsigned int i;
  for (i = 0; i < 3; ++i) {
    pos[i] = 0.0;
    vel[i] = 0.0;
  }
  ifs.open(kFBin, std::ios::binary);
  if (!ifs) {
    std::cout << "[ERROR] " << kFBin
              << " could not be found in this directory!" << std::endl;
    exit(EXIT_FAILURE);
  }
  // vector 用メモリ確保
  ttls.reserve(3);
  cnams.reserve(800);
  sss.reserve(3);
  ipts.reserve(13);
  cvals.reserve(572);  // DE430 で NCON = 572 であることを前提に
  coeffs.reserve(13);
}

/*
 * @brief  デストラクタ
 *
 * @param  <none>
 */
Jpl::~Jpl() { ifs.close(); }

/*
 * @brief   バイナリファイル読み込み
 *
 * @param   <none>
 * @return  <none>
 */
void Jpl::read_bin() {

  try {
    // ヘッダ（1レコード目）
    get_ttl(ttls);      // TTL  (タイトル)
    get_cnam(cnams);    // CNAM (定数名)(400+400件)
    get_ss(sss);        // SS   (ユリウス日(開始,終了),分割日数)
    get_ncon(ncon);     // NCON (定数の数)
    get_au(au);         // AU   (天文単位)
    get_emrat(emrat);   // EMRAT(地球と月の質量比)
    get_ipt(ipts);      // IPT  (オフセット,係数の数,サブ区間数)(水星〜月の章動,月の秤動)
    get_numde(numde);   // NUMDE(DEバージョン番号)
    // ヘッダ（2レコード目）
    get_cval(cvals);    // CVAL (定数値)
    // レコードインデックス取得
    idx = static_cast<int>(jd - sss[0]) / sss[2];
    // 係数取得（対象のインデックス分を取得）
    get_coeff(coeffs);  // COEFF（係数）
  } catch (...) {
    throw;
  }
}

/*
 * @brief      位置・速度(Positions(Radian), Velocities(Radian/Day)) 計算
 *
 * @param[in]  天体番号: 対象 (unsigned int)
 * @param[in]  天体番号: 基準 (unsigned int)
 * @return     <none>
 */
void Jpl::calc_pv(unsigned int t, unsigned int c) {
  unsigned int i;
  unsigned int j;

  try {
    // 計算結果初期化
    for (i = 0; i < 3; ++i) {
      p_sun[i] = 0.0;
      v_sun[i] = 0.0;
      p_nut[i] = 0.0;
      v_nut[i] = 0.0;
    }
    for (i = 0; i < 11; ++i) {
      for (j = 0; j < 3; ++j) {
        ps[i][j] = 0.0;
        vs[i][j] = 0.0;
      }
    }
    for (i = 0; i < 13; ++i) {
      for (j = 0; j < 3; ++j) {
        ps_2[i][j] = 0.0;
        vs_2[i][j] = 0.0;
      }
    }

    // 計算対象フラグ一覧取得
    get_list(t, c, list);

    // 補間（11: 太陽）
    interpolate(11, p_sun, v_sun);

    // 補間（1:水星〜10:月）
    for (i = 0; i < 10; ++i) {
      if (list[i] == 0) { continue; }
      interpolate(i + 1, ps[i], vs[i]);
      if (i > 8) { continue; }
      if (is_bary) { continue; }
      for (j = 0; j < 3; ++j) {
        ps[i][j] = ps[i][j] - p_sun[j];
        vs[i][j] = vs[i][j] - v_sun[j];
      }
    }

    // 補間（14:地球の章動）
    if (list[10] > 0 && ipts[11][1] > 0) {
      interpolate(14, p_nut, v_nut);
    }

    // 補間（15:月の秤動）
    if (list[11] > 0 && ipts[12][1] > 0) {
      interpolate(15, ps[10], vs[10]);
    }

    // 対象天体と基準天体の差
    if (t == 14) {
      if (ipts[11][1] > 0) {
        for (i = 0; i < 3; ++i) {
          pos[i] = p_nut[i];
          vel[i] = v_nut[i];
        }
      }
    } else if (t == 15) {
      if (ipts[12][1] > 0) {
        for (i = 0; i < 3; ++i) {
          pos[i] = ps[10][i];
          vel[i] = vs[10][i];
        }
      }
    } else {
      for (i = 0; i < 10; ++i) {
        for (j = 0; j < 3; ++j) {
          ps_2[i][j] = ps[i][j];
          vs_2[i][j] = vs[i][j];
        }
      }
      if (t == 11 || c == 11) {
        for (i = 0; i < 3; ++i) {
          ps_2[10][i] = p_sun[i];
          vs_2[10][i] = v_sun[i];
        }
      }
      if (t == 12 || c == 12) {
        for (i = 0; i < 3; ++i) {
          ps_2[11][i] = 0.0;
          vs_2[11][i] = 0.0;
        }
      }
      if (t == 13 || c == 13) {
        for (i = 0; i < 3; ++i) {
          ps_2[12][i] = ps[2][i];
          vs_2[12][i] = vs[2][i];
        }
      }
      if (t * c == 30 || t + c == 13) {
        for (i = 0; i < 3; ++i) {
          ps_2[2][i] = 0.0;
          vs_2[2][i] = 0.0;
        }
      } else {
        if (list[2] != 0) {
          for (i = 0; i < 3; ++i) {
            ps_2[2][i] = ps[2][i] - ps[9][i] / (1.0 + emrat);
            vs_2[2][i] = vs[2][i] - vs[9][i] / (1.0 + emrat);
          }
        }
        if (list[9] != 0) {
          for (i = 0; i < 3; ++i) {
            ps_2[9][i] = ps_2[2][i] + ps[9][i];
            vs_2[9][i] = vs_2[2][i] + vs[9][i];
          }
        }
      }
      for (i = 0; i < 3; ++i) {
        pos[i] = ps_2[t - 1][i] - ps_2[c - 1][i];
        vel[i] = vs_2[t - 1][i] - vs_2[c - 1][i];
      }
    }
  } catch (...) {
    throw;
  }
}

/*
 * @brief       取得: TTL
 *              - 84 byte * 3
 *
 * @param[ref]  TTL 一覧 (vector<string>)
 */
void Jpl::get_ttl(std::vector<std::string>& ttls) {
  try {
    get_str_list(kPosTtl, kReclTtl, kCntTtl, ttls);
  } catch (...) {
    throw;
  }
}

/*
 * @brief       取得: CNAM
 *              - 6 byte * (400 + 400)
 *
 * @param[ref]  CNAM 一覧 (vector<string>)
 */
void Jpl::get_cnam(std::vector<std::string>& cnams) {
  std::vector<std::string> cnam2s;   // CNAM2 (6 byte * 400)

  try {
    get_str_list(kPosCnam,  kReclCnam, kCntCnam, cnams );  // 最初の400件
    get_str_list(kPosCnam2, kReclCnam, kCntCnam, cnam2s);  // 後部の400件
    std::copy(cnam2s.begin(), cnam2s.end(), std::back_inserter(cnams));
  } catch (...) {
    throw;
  }
}

/*
 * @brief       取得: SS
 *              - 8 byte * 3
 *
 * @param[ref]  SS 一覧 (vector<double>)
 */
void Jpl::get_ss(std::vector<double>& sss) {
  try {
    get_dbl_list(kPosSs, kReclSs, kCntSs, sss);
  } catch (...) {
    throw;
  }
}

/*
 * @brief       取得: NCON
 *              - 4 byte * 1
 *
 * @param[ref]  NCON (unsigned int)
 */
void Jpl::get_ncon(unsigned int& ncon) {
  try {
    get_val<unsigned int>(kPosNcon, kReclNcon, ncon);
  } catch (...) {
    throw;
  }
}

/*
 * @brief       取得: AU
 *              - 8 byte * 1
 *
 * @param[ref]  AU (double)
 */
void Jpl::get_au(double& au) {
  try {
    get_val<double>(kPosAu, kReclAu, au);
  } catch (...) {
    throw;
  }
}

/*
 * @brief       取得: EMRAT
 *              - 8 byte * 1
 *
 * @param[ref]  EMRAT (double)
 */
void Jpl::get_emrat(double& emrat) {
  try {
    get_val<double>(kPosEmrat, kReclEmrat, emrat);
  } catch (...) {
    throw;
  }
}

/*
 * @brief       取得: NUMDE
 *              - 4 byte * 1
 *
 * @param[ref]  NUMDE (unsigned int)
 */
void Jpl::get_numde(unsigned int& numde) {
  try {
    get_val<unsigned int>(kPosNumde, kReclNumde, numde);
  } catch (...) {
    throw;
  }
}

/*
 * @brief       取得: IPT
 *
 * @param[ref]  値一覧 (vector<vector<unsigned int>>)
 */
void Jpl::get_ipt(std::vector<std::vector<unsigned int>>& vals) {
  unsigned int i;
  unsigned int j;
  unsigned int pos  = kPosIpt;
  unsigned int recl = kReclIpt;
  std::vector<unsigned int> ary;
  unsigned int buf;

  try {
    // IPT13 以外
    for (i = 0; i < kCntIpt; ++i) {
      ary.clear();
      for (j = 0; j < 3; ++j) {
        ifs.seekg(pos);
        ifs.read((char*)&buf, recl);
        ary.push_back(buf);
        pos += recl;
      }
      vals.push_back(ary);
    }
    // IPT13
    ary.clear();
    pos = kPosIpt2;
    for (j = 0; j < 3; ++j) {
      ifs.seekg(pos);
      ifs.read((char*)&buf, recl);
      ary.push_back(buf);
      pos += recl;
    }
    vals.push_back(ary);
  } catch (...) {
    throw;
  }
}

/*
 * @brief       取得: CVAL
 *              - 8 byte * NCON
 *
 * @param[ref]  CVAL 一覧 (vector<double>)
 */
void Jpl::get_cval(std::vector<double>& cvals) {
  try {
    get_dbl_list(kPosCval, kReclCval, ncon, cvals);
  } catch (...) {
    throw;
  }
}

/*
 * @brief       取得: COEFF
 *              - 8 byte * ?
 *              - 地球・章動のみ要素数が 2 で、その他の要素数は 3
 *
 * @param[ref]  値一覧 (vector<vector<vector<vector<double>>>>)
 */
void Jpl::get_coeff(
    std::vector<std::vector<std::vector<std::vector<double>>>>& vals) {
  unsigned int i;
  unsigned int j;
  unsigned int k;
  unsigned int l;
  double buf;
  std::vector<double> ary_a;                 // 該当インデックス分全て
  std::vector<double> ary;                   // 該当惑星分のみ
  std::vector<double> ary_w_1;               // 作業用
  std::vector<std::vector<double>> ary_w_2;  // 作業用
  std::vector<std::vector<std::vector<double>>> ary_p;  // 作業用
  unsigned int offset;
  unsigned int cnt_coeff;
  unsigned int cnt_sub;
  unsigned int n;
  unsigned int pos  = kKsize * kRecl * (2 + idx);
  unsigned int recl = 8;

  try {
    // 該当インデックス分全て取得
    for (i = 0; i < kKsize / 2; ++i) {
      ifs.seekg(pos);
      ifs.read((char*)&buf, recl);
      ary_a.push_back(buf);
      pos += recl;
    }

    // Julian Day (start, end)
    for (i = 0; i < 2; ++i) { jds[i] = ary_a[i]; }

    // 全惑星分
    for (i = 0; i < 13; ++i) {
      // 該当惑星分のみ抽出
      offset    = ipts[i][0];
      cnt_coeff = ipts[i][1];
      cnt_sub   = ipts[i][2];
      n = 3;
      if ((i + 1) == 12) { n = 2; }
      ary.clear();
      for (j = offset - 1; j < offset - 1 + cnt_coeff * n * cnt_sub; ++j) {
        ary.push_back(ary_a[j]);
      }

      // 3次元配列化
      // [サブ区間数, 要素数(3 or 2), 係数の数]
      ary_p.clear();
      for (j = 0; j < ary.size() / (cnt_coeff * n); ++j) {
        ary_w_2.clear();
        for (k = 0; k < n; ++k) {
          ary_w_1.clear();
          for (l = 0; l < cnt_coeff; ++l) {
            ary_w_1.push_back(ary[j * cnt_coeff * n + k * cnt_coeff + l]);
          }
          ary_w_2.push_back(ary_w_1);
        }
        ary_p.push_back(ary_w_2);
      }
      vals.push_back(ary_p);
    }
  } catch (...) {
    throw;
  }
}

/*
 * @brief       取得: (unsigned int|double) 型 1件 template
 *
 * @param[in]   レコード位置 (unsigned int)
 * @param[in]   レコード長 (unsigned int)
 * @param[ref]  値 (class T)
 */
template <class T>
void Jpl::get_val(unsigned int pos, unsigned int recl, T& val) {
  try {
    ifs.seekg(pos);
    ifs.read((char*)&val, recl);
  } catch (...) {
    throw;
  }
}

/*
 * @brief       取得: vector<double> 型
 *
 * @param[in]   レコード位置 (unsigned int)
 * @param[in]   レコード長 (unsigned int)
 * @param[in]   件数 (unsigned int)
 * @param[ref]  値一覧 (vector<double>)
 */
void Jpl::get_dbl_list(
    unsigned int pos, unsigned int recl, unsigned int cnt,
    std::vector<double>& vals) {
  unsigned int i;
  double       buf;

  try {
    for (i = 0; i < cnt; ++i) {
      ifs.seekg(pos);
      ifs.read((char*)&buf, recl);
      vals.push_back(buf);
      pos += recl;
    }
  } catch (...) {
    throw;
  }
}

/*
 * @brief       取得: vector<string> 型
 *              (後のスペースは trim)
 *
 * @param[in]   レコード位置 (unsigned int)
 * @param[in]   レコード長 (unsigned int)
 * @param[in]   件数 (unsigned int)
 * @param[ref]  値一覧 (vector<string>)
 */
void Jpl::get_str_list(
    unsigned int pos, unsigned int recl, unsigned int cnt,
    std::vector<std::string>& vals) {
  unsigned int i;
  char*        buf;
  std::string  str;

  try {
    for (i = 0; i < cnt; ++i) {
      ifs.seekg(pos);
      buf = new char[recl];
      ifs.read(buf, recl);
      str = buf;
      vals.push_back(str.erase(str.find_last_not_of(" ") + 1));
      pos += recl;
      delete buf;
    }
  } catch (...) {
    throw;
  }
}

/*
 * @brief       計算対象フラグ一覧（係数データの並びに対応）取得
 *              （計算区分 0: 計算しない、1: 位置・速度を計算）
 *
 * @param[in]   天体番号: 対象 (unsigned int)
 * @param[in]   天体番号: 基準 (unsigned int)
 * @param[ref]  計算対象フラグ一覧 (unsigned int[12])
 */
void Jpl::get_list(unsigned int t, unsigned int c, unsigned int(&list)[12]) {
  unsigned int i;

  try {
    for (i = 0; i < 12; ++i) { list[i] = 0; }  // 0 で初期化
    if (t == 14) {
      if (ipts[11][1] > 0) { list[10] = kKind; }
      return;
    }
    if (t == 15) {
      if (ipts[12][1] > 0) { list[11] = kKind; }
      return;
    }
    if (t <= 10) { list[t - 1] = kKind; }
    if (t == 10) { list[2]     = kKind; }
    if (t ==  3) { list[9]     = kKind; }
    if (t == 13) { list[2]     = kKind; }
    if (c <= 10) { list[c - 1] = kKind; }
    if (c == 10) { list[2]     = kKind; }
    if (c ==  3) { list[9]     = kKind; }
    if (c == 13) { list[2]     = kKind; }
  } catch (...) {
    throw;
  }
}

/*
 * @brief       補間
 *              * 使用するチェビシェフ多項式の係数は、
 *              * 天体番号が 1 〜 13 の場合は、 x, y, z の位置・速度（6要素）、
 *                天体番号が 14 の場合は、 Δψ, Δε の角位置・角速度（4要素）、
 *                天体番号が 15 の場合は、 φ, θ, ψ の角位置・角速度（6要素）。
 *              * 天体番号が 12 の場合は、 x, y, z の位置・速度の値は全て 0.0 とする。
 *
 * @param[in]   天体番号 (unsigned int)
 * @param[ref]  位置(x, y, z) (double[3])
 * @param[ref]  速度(x, y, z) (double[3])
 *              * 14（地球の章動）の場合、
 *                  位置 = [Δψ の角位置, Δε の角位置, 0.0]
 *                  速度 = [Δψ の角速度, Δε の角速度, 0.0]
 *              * 15（月の秤動）の場合、
 *                  位置 = [ φ の角位置, θ の角位置, ψ の角位置]
 *                  速度 = [ φ の角速度, θ の角速度, ψ の角速度]
 */
void Jpl::interpolate(
    unsigned int astr, double(&pos)[3], double(&vel)[3]) {
  unsigned int n_item = 3;         // 要素数
  unsigned int i_ipt  = astr - 1;  // インデックス（ipts 用）
  unsigned int i_coef = astr - 1;  // インデックス（coeffs 用）
  std::vector<double> wk_pos;      // 作業用 vector （位置）
  std::vector<double> wk_vel;      // 作業用 vector （速度）
  unsigned int        s;           // 作業用 vector サイズ
  double              v;           // 作業用
  unsigned int        i;           // ループインデックス
  unsigned int        j;           // ループインデックス

  try {
    // チェビシェフ時間、サブインデックス等
    norm_time(astr, tc, idx_s);
    if (astr == 14) { n_item = 2; }
    if (astr > 13) {
      i_ipt  = astr - 3;
      i_coef = astr - 3;
    }

    // 位置
    wk_pos.push_back(1.0);
    wk_pos.push_back(tc);
    for (i = 2; i < ipts[i_ipt][1]; ++i) {
      s = wk_pos.size();
      wk_pos.push_back(2.0 * tc * wk_pos[s - 1] - wk_pos[s - 2]);
    }
    for (i = 0; i < n_item; ++i) {
      v = 0;
      for (j = 0; j < ipts[i_ipt][1]; ++j) {
        v += coeffs[i_coef][idx_s][i][j] * wk_pos[j];
      }
      if (!is_km && astr < 14) { v /= au; }
      pos[i] = v;
    }

    // 速度
    wk_vel.push_back(0.0);
    wk_vel.push_back(1.0);
    wk_vel.push_back(2.0 * 2.0 * tc);
    for (i = 3; i < ipts[i_ipt][1]; ++i) {
      s = wk_vel.size();
      wk_vel.push_back(
          2.0 * tc * wk_vel[s - 1] + 2.0 * wk_pos[i - 1] - wk_vel[s - 2]);
    }
    for (i = 0; i < n_item; ++i) {
      v = 0;
      for (j = 0; j < ipts[i_ipt][1]; ++j) {
        v += coeffs[i_coef][idx_s][i][j] * wk_vel[j] * 2.0 * ipts[i_ipt][2]
           / sss[2];
      }
      if (astr < 14) {
        if (is_km) { v /= kSecDay; } else { v /= au; }
      }
      vel[i] = v;
    }
  } catch (...) {
    throw;
  }
}

/*
 * @brief       チェビシェフ多項式用に時刻を正規化、サブ区間のインデックス算出
 *
 * @param[in]   天体番号 (unsigned int)
 * @param[ref]  チェビシェフ時間 (double)
 * @param[ref]  サブ区間のインデックス (unsigned int)
 */
void Jpl::norm_time(unsigned int astr, double& tc, unsigned int& idx_s) {
  double tmp;

  try {
    idx_s = astr;
    if (astr > 13) { idx_s = astr - 2; }
    tc = (jd - jds[0]) / sss[2];
    tmp = tc * ipts[idx_s - 1][2];
    idx_s = static_cast<int>(tmp - static_cast<int>(tc));
    tc = (fmod(tmp, 1.0) + static_cast<int>(tc)) * 2 - 1;
  } catch (...) {
    throw;
  }
}

}  // namespace jpl_calc_430

