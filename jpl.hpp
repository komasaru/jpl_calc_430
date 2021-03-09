#ifndef JPL_READ_430_JPL_HPP_
#define JPL_READ_430_JPL_HPP_

#include <cmath>
#include <cstdlib>   // for EXIT_XXXX
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace jpl_calc_430 {

class Jpl {
  unsigned int  astr_t;       // 天体番号: 対象
  unsigned int  astr_c;       // 天体番号: 基準
  double        jd;           // ユリウス日
  bool          is_km;        // 単位フラグ
                              // (true: km, km/sec, false: AU, AU/day)
  bool          is_bary;      // 基準フラグ
                              // (true: 太陽系重心が基準, false: 太陽が基準)
  std::ifstream ifs;          // バイナリファイル
  double        p_sun[3];     // 位置: 11:太陽
  double        v_sun[3];     // 速度: 11:太陽
  double        ps[11][3];    // 位置: 1:水星〜10:月, 15:月の秤動
  double        vs[11][3];    // 速度: 1:水星〜10:月, 15:月の秤動
  double        p_nut[3];     // 位置: 14:地球の章動
  double        v_nut[3];     // 速度: 14:地球の章動（実質、無使用）
  double        ps_2[13][3];  // 位置: 対象天体と基準天体の差算出用
  double        vs_2[13][3];  // 速度: 対象天体と基準天体の差算出用
  unsigned int  list[12];     // 計算対象フラグ一覧
  double        tc;           // チェビシェフ時間
  unsigned int  idx_s;        // サブ区間のインデックス

  void get_ttl(std::vector<std::string>&);       // 取得: TTL
  void get_cnam(std::vector<std::string>&);      // 取得: CNAM
  void get_ss(std::vector<double>&);             // 取得: SS
  void get_ncon(unsigned int&);                  // 取得: NCON
  void get_au(double&);                          // 取得: AU
  void get_emrat(double&);                       // 取得: EMRAT
  void get_numde(unsigned int&);                 // 取得: NUMDE
  void get_ipt(std::vector<std::vector<unsigned int>>&);   // 取得: IPT
  void get_cval(std::vector<double>&);           // 取得: CVAL
  void get_coeff(std::vector<std::vector<std::vector<std::vector<double>>>>&);
                                                 // 取得: COEFF
  template <class T>
  void get_val(unsigned int, unsigned int, T&);  // 取得: 値1件(template)
  void get_dbl_list(
      unsigned int, unsigned int, unsigned int,
      std::vector<double>&);                     // 取得: vector<double> 型
  void get_str_list(
      unsigned int, unsigned int, unsigned int,
      std::vector<std::string>&);                // 取得: vector<string> 型
  void get_list(unsigned int, unsigned int, unsigned int(&)[12]);
                                                 // 計算対象フラグ一覧取得
  void interpolate(unsigned int, double(&)[3], double(&)[3]);  // 補間
  void norm_time(unsigned int, double&, unsigned int&);
                                                 // チェビシェフ多項式用に時刻を正規化、
                                                 // サブ区間のインデックス算出

public:
  std::vector<std::string>               ttls;    // TTL   (84 byte *   3)
  std::vector<std::string>               cnams;   // CNAM  ( 6 byte * 800)
  std::vector<double>                    sss;     // SS    ( 8 byte *   3)
  unsigned int                           ncon;    // NCON  ( 4 byte *   1)
  double                                 au;      // AU    ( 8 byte *   1)
  double                                 emrat;   // EMRAT ( 8 byte *   1)
  unsigned int                           numde;   // NUMDE ( 4 byte *   1)
  std::vector<std::vector<unsigned int>> ipts;    // IPT   ( 4 byte * 13 * 3)
  std::vector<double>                    cvals;   // SS    ( 8 byte * NCON)
  unsigned int                           idx;     // レコードインデックス
  std::vector<std::vector<std::vector<std::vector<double>>>> coeffs;
                                         // COEFF ( 8 byte *   ?)
                                         // 全惑星分の cnt_sub * (2 or 3) * cnt_coeff
                                         // （4次元配列）
  double                                 jds[2];  // JD (開始、終了)
  double                                 pos[3];  // 計算結果: 位置
  double                                 vel[3];  // 計算結果: 速度

  Jpl(double, const bool = false, const bool = true);  // コンストラクタ
                       // (引数: ユリウス日, [単位フラグ, [基準フラグ]])
  ~Jpl();                                              // デストラクタ
  void read_bin();                                     // バイナリファイル読み込み
  void calc_pv(unsigned int, unsigned int);            // 位置・速度計算
};

}  // namespace jpl_calc_430

#endif

