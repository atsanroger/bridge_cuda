/*!
     @file    vsimd_double-inc.h
     @brief
     @author  Hideo Matsufuru (matufuru)
     $LastChangedBy: matufuru $
     @date    $LastChangedDate:: 2023-12-24 19:34:22 #$
     @version $LastChangedRevision: 2562 $
*/

#ifndef QXS_VSIMD_INCLUDED
#define QXS_VSIMD_INCLUDED

typedef struct{
  double v[VLEND];
} Vsimd_t;

#ifdef __ARM_FEATURE_SVE
#include <arm_sve.h>

typedef svfloat64_t svreal_t;

typedef svint64_t svint_t;

typedef svuint64_t svuint_t;

typedef int64_t int_t;

typedef uint64_t uint_t;

namespace {

  inline svbool_t set_predicate(){
    return svptrue_b64();
  }

  inline svbool_t set_predicate_false(){
    return svpfalse();
  }

  inline svbool_t set_predicate_whilelt(int range){
    return svwhilelt_b64(0, range);
  }

  inline void set_predicate_xp(svbool_t& pg1, svbool_t& pg2){
    svbool_t pg0 = svptrue_b64();
    pg2 = svpfalse();
    for(int iy = VLENYD; iy > 0; --iy){
      pg1 = svwhilelt_b64(0, iy * VLENXD);
      pg2 = sveor_z(pg0, pg2, pg1);
      pg1 = svwhilelt_b64(0, iy * VLENXD - 1);
      pg2 = sveor_z(pg0, pg2, pg1);
    }
    pg1 = svnot_z(pg0, pg2);
  }

  inline void set_predicate_xm(svbool_t& pg1, svbool_t& pg2){
    svbool_t pg0 = svptrue_b64();
    pg2 = svwhilelt_b64(0, VLEND - VLENXD + 1);
    for(int iy = VLENYD-1; iy > 0; --iy){
      pg1 = svwhilelt_b64(0, iy * VLENXD);
      pg2 = sveor_z(pg0, pg2, pg1);
      pg1 = svwhilelt_b64(0, (iy-1) * VLENXD + 1);
      pg2 = sveor_z(pg0, pg2, pg1);
    }
    pg1 = svnot_z(pg0, pg2);
  }

  inline void set_predicate_xp_eo(svbool_t& pg1, svbool_t& pg2,
                                  svbool_t& pg3, int ieo) {
    svbool_t pg0 = svptrue_b64();
    pg2 = svpfalse();
    pg3 = svpfalse();
    for(int iy = VLENYD; iy > 0; --iy){
      if(iy % 2 == ieo){
        pg1 = svwhilelt_b64(0, iy * VLENXD);
        pg2 = sveor_z(pg0, pg2, pg1);
        pg1 = svwhilelt_b64(0, iy * VLENXD - 1);
        pg2 = sveor_z(pg0, pg2, pg1);
      }else{
        pg1 = svwhilelt_b64(0, iy * VLENXD);
        pg3 = sveor_z(pg0, pg3, pg1);
        pg1 = svwhilelt_b64(0, (iy-1) * VLENXD);
        pg3 = sveor_z(pg0, pg3, pg1);
      }
    }
    pg1 = sveor_z(pg0, pg2, pg3);
    pg1 = svnot_z(pg0, pg1);
  }

  inline void set_predicate_xm_eo(svbool_t& pg1, svbool_t& pg2,
                                  svbool_t& pg3, int ieo) {
    svbool_t pg0 = svptrue_b64();
    pg1 = svpfalse();
    pg3 = svpfalse();
    for(int iy = VLENYD; iy > 0; --iy){
      if(iy % 2 == ieo){
        pg2 = svwhilelt_b64(0, iy * VLENXD);
        pg3 = sveor_z(pg0, pg3, pg2);
        pg2 = svwhilelt_b64(0, (iy-1) * VLENXD);
        pg3 = sveor_z(pg0, pg3, pg2);
      }else{
        pg2 = svwhilelt_b64(0, iy * VLENXD);
        pg1 = sveor_z(pg0, pg1, pg2);
        pg2 = svwhilelt_b64(0, (iy-1) * VLENXD + 1);
        pg1 = sveor_z(pg0, pg1, pg2);
      }
    }
    pg2 = sveor_z(pg0, pg1, pg3);
    pg2 = svnot_z(pg0, pg2);
  }

  inline void set_predicate_yp(svbool_t& pg1, svbool_t& pg2){
    svbool_t pg0 = svptrue_b64();
    pg1 = svwhilelt_b64(0, VLENXD * (VLENYD-1));
    pg2 = svnot_z(pg0, pg1);
  }

  inline void set_predicate_ym(svbool_t& pg1, svbool_t& pg2){
    svbool_t pg0 = svptrue_b64();
    pg2 = svwhilelt_b64(0, VLENXD);
    pg1 = svnot_z(pg0, pg2);
  }


  // Nitadori-san's implementation
  inline void set1_at(const int i, svbool_t &pg){
    svbool_t pg0 = svptrue_b64();
    svbool_t pg1 = svwhilelt_b64(0, i);
    svbool_t pg2 = svwhilelt_b64(0, i+1);
    pg = sveor_z(pg0, pg, pg1);
    pg = sveor_z(pg0, pg, pg2);
  }

  inline void rot1_R(uint64_t *u, const int len = VLENXD)
  {
    uint64_t tmp = u[len-1]; // tail
    for(int i=len-1; i>=1; i--){
      u[i] = u[i-1];
    }
    u[0] = tmp;
  }

  inline void rot1_L(uint64_t *u, const int len = VLENXD)
  {
    uint64_t tmp = u[0]; // head
    for(int i = 0; i < len-1; ++i){
      u[i] = u[i+1];
    }
    u[len-1] = tmp;
  }

  // Nitadori-san's implementation
  inline void set_idx_predicate_xp_eo(svbool_t  &pg,
                                      svuint64_t &idx,
                                      const int ieo)
  {
    uint64_t u[VLEN];
    for(int i = 0; i < VLEN; ++i) u[i] = i;

    pg = svpfalse();
    if(0 == ieo){
      // L-shift odd rows
      for(int i = VLENXD; i < VLEN; i += 2*VLENXD){
        set1_at(i, pg);
        rot1_L(u+i, VLENXD);
      }
    }
    if(1 == ieo){
      // L-shift env rows
      for(int i = 0; i < VLEN; i += 2*VLENXD){
	set1_at(i, pg);
	rot1_L(u+i, VLENXD);
      }
    }
    idx = svld1_u64(svptrue_b64(), u);
  }

  inline void set_idx_predicate_xm_eo(svbool_t  &pg,
                                      svuint64_t &idx,
                                      const int ieo)
  {
    uint64_t u[VLEN];
    for(int i = 0; i<VLEN; i++) u[i] = i;

    pg = svpfalse();
    if(0 == ieo){
      // R-shift evn rows
      for(int i = 0; i < VLEN; i += 2*VLENXD){
	set1_at(i + VLENXD-1, pg); // 3, 11
	rot1_R(u+i, VLENXD);
      }
    }
    if(1 == ieo){
      // R-shift odd rows
      for(int i = VLENXD; i < VLEN; i += 2*VLENXD){
	set1_at(i + VLENXD-1, pg); // 7, 15
	rot1_R(u+i, VLENXD);
      }
    }
    idx = svld1_u64(svptrue_b64(), u);
  }

  inline void set_idx_predicate_yp(svbool_t& pg1, svuint64_t &idx)
  {
    pg1 = svwhilelt_b64(0, VLENXD);
    uint64_t u[VLEND];
    for(int i = 0; i < VLENXD*(VLENYD-1); ++i){
      u[i] = i + VLENXD;
    }
    for(int i = 0; i < VLENXD; ++i){
      u[i+VLENXD*(VLENYD-1)] = i;
    }
    idx = svld1_u64(svptrue_b64(), u);
  }

  inline void set_idx_predicate_ym(svbool_t& pg1, svuint64_t &idx)
  {
    svbool_t pg2 = svwhilelt_b64(0, VLENXD*(VLENYD-1));
    pg1 = svnot_z(svptrue_b64(), pg2);
    uint64_t u[VLEND];
    for(int i = 0; i < VLENXD; ++i){
      u[i] = i + VLENXD*(VLENYD-1);
    }
    for(int i = VLENXD; i < VLEND; ++i){
      u[i] = i - VLENXD;
    }
    idx = svld1_u64(svptrue_b64(), u);
  }

  inline void load_vec(svbool_t pg, svfloat64_t& v, const float64_t* vp){
    v = svld1_f64(pg, vp);
  }

  inline void load_add(svbool_t pg, svfloat64_t& v, float64_t* vp){
    svfloat64_t v2;
    v2 = svld1_f64(pg, vp);
    v  = svsel_f64(pg, v2, v);
    //    v  = svmul_m(pg, v, 0.0);
    //    v  = svadd_m(pg, v, v2);
  }

  inline void set_vec(svbool_t pg, svfloat64_t& v, float64_t a){
    v  = svdup_f64_m(v, pg, a);
  }

  inline void clear_vec(svbool_t pg, svfloat64_t& v){
    v  = svdup_f64_m(v, pg, 0.0);
  }

  inline void load_svint(svbool_t pg, svint64_t& v, int64_t* vp){
    v = svld1_s64(pg, vp);
  }

  inline void load_svuint(svbool_t pg, svuint64_t& v, uint64_t* vp){
    v = svld1_u64(pg, vp);
  }

  inline void save_vec(svbool_t pg, double* v, svfloat64_t vt){
    svst1(pg, v, vt);
  }

  inline void load_vec_gather(svbool_t pg, svfloat64_t& v, float64_t* vp,
                              svint64_t index){
    v = svld1_gather_s64index_f64(pg, vp, index);
  }

  inline void load_add_gather(svbool_t pg, svfloat64_t& v, float64_t* vp,
                              svint64_t index){
    svfloat64_t v2;
    v2 = svld1_gather_s64index_f64(pg, vp, index);
    v  = svsel_f64(pg, v2, v);
    //    v  = svmul_m(pg, v, 0.0);
    //    v  = svadd_m(pg, v, v2);
  }

  inline void flip_sign(svbool_t pg, svfloat64_t& v){
    v = svneg_f64_m(v, pg, v);
  }

  inline void flip_sign(svbool_t pg, svfloat64_t& v1, svfloat64_t &v2){
    v1 = svneg_f64_m(v2, pg, v2);
  }

}

#endif  // __ARM_FEATURE_SVE

#endif
