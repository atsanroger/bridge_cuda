/*!
        @file    vsimd_half-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
	$LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-12-24 19:34:22 #$
        @version $LastChangedRevision: 2562 $
*/

#ifndef QXS_VSIMD_INCLUDED
#define QXS_VSIMD_INCLUDED

#ifdef __ARM_FEATURE_SVE
#include <arm_sve.h>

#include <stdlib.h>

typedef struct{
  _Float16  v[VLENH];
  //float16_t v[VLENH];
} Vsimd_t;

//typedef float16_t half;

typedef svfloat16_t svreal_t;

typedef svint16_t svint_t;

typedef svuint16_t svuint_t;

typedef int16_t int_t;

typedef uint16_t uint_t;

namespace {

  inline svbool_t set_predicate(){
    return svptrue_b16();
  }

  inline svbool_t set_predicate_false(){
    return svpfalse();
  }

  inline svbool_t set_predicate_whilelt(int range){
    return svwhilelt_b16(0, range);
  }

  inline void set_predicate_xp(svbool_t& pg1, svbool_t& pg2){
    svbool_t pg0 = svptrue_b16();
    pg2 = svpfalse();
    for(int iy = VLENYH; iy > 0; --iy){
      pg1 = svwhilelt_b16(0, iy * VLENXH);
      pg2 = sveor_z(pg0, pg2, pg1);
      pg1 = svwhilelt_b16(0, iy * VLENXH - 1);
      pg2 = sveor_z(pg0, pg2, pg1);
    }
    pg1 = svnot_z(pg0, pg2);
  }

  inline void set_predicate_xm(svbool_t& pg1, svbool_t& pg2){
    svbool_t pg0 = svptrue_b16();
    pg2 = svwhilelt_b16(0, VLENH - VLENXH + 1);
    for(int iy = VLENYH-1; iy > 0; --iy){
      pg1 = svwhilelt_b16(0, iy * VLENXH);
      pg2 = sveor_z(pg0, pg2, pg1);
      pg1 = svwhilelt_b16(0, (iy-1) * VLENXH + 1);
      pg2 = sveor_z(pg0, pg2, pg1);
    }
    pg1 = svnot_z(pg0, pg2);
  }

  inline void set_predicate_xp_eo(svbool_t& pg1, svbool_t& pg2,
                                  svbool_t& pg3, int ieo) {
    svbool_t pg0 = svptrue_b16();
    pg2 = svpfalse();
    pg3 = svpfalse();
    for(int iy = VLENYH; iy > 0; --iy){
      if(iy % 2 == ieo){
        pg1 = svwhilelt_b16(0, iy * VLENXH);
        pg2 = sveor_z(pg0, pg2, pg1);
        pg1 = svwhilelt_b16(0, iy * VLENXH - 1);
        pg2 = sveor_z(pg0, pg2, pg1);
      }else{
        pg1 = svwhilelt_b16(0, iy * VLENXH);
        pg3 = sveor_z(pg0, pg3, pg1);
        pg1 = svwhilelt_b16(0, (iy-1) * VLENXH);
        pg3 = sveor_z(pg0, pg3, pg1);
      }
    }
    pg1 = sveor_z(pg0, pg2, pg3);
    pg1 = svnot_z(pg0, pg1);
  }

  inline void set_predicate_xm_eo(svbool_t& pg1, svbool_t& pg2,
                                  svbool_t& pg3, int ieo) {
    svbool_t pg0 = svptrue_b16();
    pg1 = svpfalse();
    pg3 = svpfalse();
    for(int iy = VLENYH; iy > 0; --iy){
      if(iy % 2 == ieo){
        pg2 = svwhilelt_b16(0, iy * VLENXH);
        pg3 = sveor_z(pg0, pg3, pg2);
        pg2 = svwhilelt_b16(0, (iy-1) * VLENXH);
        pg3 = sveor_z(pg0, pg3, pg2);
      }else{
        pg2 = svwhilelt_b16(0, iy * VLENXH);
        pg1 = sveor_z(pg0, pg1, pg2);
        pg2 = svwhilelt_b16(0, (iy-1) * VLENXH + 1);
        pg1 = sveor_z(pg0, pg1, pg2);
      }
    }
    pg2 = sveor_z(pg0, pg1, pg3);
    pg2 = svnot_z(pg0, pg2);
  }

  inline void set_predicate_yp(svbool_t& pg1, svbool_t& pg2){
    svbool_t pg0 = svptrue_b16();
    pg1 = svwhilelt_b16(0, VLENXH * (VLENYH-1));
    pg2 = svnot_z(pg0, pg1);
  }

  inline void set_predicate_ym(svbool_t& pg1, svbool_t& pg2){
    svbool_t pg0 = svptrue_b16();
    pg2 = svwhilelt_b16(0, VLENXH);
    pg1 = svnot_z(pg0, pg2);
  }

  // Nitadori-san's implementation
  inline void set1_at(const int i, svbool_t &pg){
    svbool_t pg0 = svptrue_b16();
    svbool_t pg1 = svwhilelt_b16(0, i);
    svbool_t pg2 = svwhilelt_b16(0, i+1);
    pg = sveor_z(pg0, pg, pg1);
    pg = sveor_z(pg0, pg, pg2);
  }

  inline void rot1_R(uint16_t *u, const int len = VLENXH)
  {
    uint16_t tmp = u[len-1]; // tail
    for(int i = len-1; i >= 1; --i){
      u[i] = u[i-1];
    }
    u[0] = tmp;
  }

  inline void rot1_L(uint16_t *u, const int len = VLENXH)
  {
    uint16_t tmp = u[0]; // head
    for(int i = 0; i < len-1; ++i){
      u[i] = u[i+1];
    }
    u[len-1] = tmp;
  }

  inline void set_idx_predicate_xp_eo(svbool_t  &pg,
                                      svuint16_t &idx,
                                      const int ieo)
  {
    uint16_t u[VLEN];
    for(int i=0; i<VLEN; i++) u[i] = i;

    pg = svpfalse();
    if(0 == ieo){
      // L-shift odd rows
      for(int i = VLENXH; i < VLEN; i += 2*VLENXH){
	set1_at(i, pg); // 3, 11
	rot1_L(u+i, VLENXH);
      }
    }
    if(1 == ieo){
      // L-shift env rows
      for(int i = 0; i < VLEN; i += 2*VLENXH){
	set1_at(i, pg); // 7, 15
	rot1_L(u+i, VLENXH);
      }
    }
    idx = svld1_u16(svptrue_b16(), u);
  }

  inline void set_idx_predicate_xm_eo(svbool_t  &pg,
                                      svuint16_t &idx,
                                      const int ieo)
  {
    uint16_t u[VLEN];
    for(int i = 0; i < VLEN; ++i) u[i] = i;

    pg = svpfalse();
    if(0 == ieo){
      // R-shift evn rows
      for(int i=0; i<VLEN; i+=2*VLENXH){
	set1_at(i + VLENXH-1, pg); // 3, 11
	rot1_R(u+i, VLENXH);
      }
    }
    if(1 == ieo){
      // R-shift odd rows
      for(int i = VLENXH; i < VLEN; i += 2*VLENXH){
	set1_at(i + VLENXH-1, pg); // 7, 15
	rot1_R(u+i, VLENXH);
      }
    }
    idx = svld1_u16(svptrue_b16(), u);
  }

  inline void set_idx_predicate_yp(svbool_t& pg1, svuint16_t &idx)
  {
    pg1 = svwhilelt_b16(0, VLENXH);
    uint16_t u[VLENH];
    for(int i = 0; i < VLENXH*(VLENYH-1); ++i){
      u[i] = i + VLENXH;
    }
    for(int i = 0; i < VLENXH; ++i){
      u[i+VLENXH*(VLENYH-1)] = i;
    }
    idx = svld1_u16(svptrue_b16(), u);
  }

  inline void set_idx_predicate_ym(svbool_t& pg1, svuint16_t &idx)
  {
    svbool_t pg2 = svwhilelt_b16(0, VLENXH*(VLENYH-1));
    pg1 = svnot_z(svptrue_b16(), pg2);
    uint16_t u[VLENH];
    for(int i = 0; i < VLENXH; ++i){
      u[i] = i + VLENXH*(VLENYH-1);
    }
    for(int i = VLENXH; i < VLENH; ++i){
      u[i] = i - VLENXH;
    }
    idx = svld1_u16(svptrue_b16(), u);
  }


  inline void set_vec(svbool_t pg, svfloat16_t& v, float16_t a){
    v  = svdup_f16_m(v, pg, a);
  }

  inline void clear_vec(svbool_t pg, svfloat16_t& v){
    v  = svdup_f16_m(v, pg, 0.0);
  }

  inline void load_vec(svbool_t pg, svfloat16_t& v, const float16_t* vp){
    v = svld1_f16(pg, vp);
  }

  inline void load_vec(svbool_t pg, svfloat16_t& v, const _Float16* vp){
    v = svld1_f16(pg, (float16_t*)vp);
  }

  inline void load_add(svbool_t pg, svfloat16_t& v, float16_t* vp){
    svfloat16_t v2;
    v2 = svld1_f16(pg, vp);
    v = svsel_f16(pg, v2, v);
  }

  inline void load_add(svbool_t pg, svfloat16_t& v, _Float16* vp){
    svfloat16_t v2;
    v2 = svld1_f16(pg, (float16_t*)vp);
    v = svsel_f16(pg, v2, v);
  }

  inline void load_svint(svbool_t pg, svint16_t& v, int16_t* vp){
    v = svld1_s16(pg, vp);
  }

  inline void load_svuint(svbool_t pg, svuint16_t& v, uint16_t* vp){
    v = svld1_u16(pg, vp);
  }

  inline void save_vec(svbool_t pg, _Float16* v, svfloat16_t vt){
    svst1(pg, (float16_t*)v, vt);
  }

  inline void save_vec(svbool_t pg, float16_t* v, svfloat16_t vt){
    svst1(pg, v, vt);
  }

  // Note that gather-load is not defined for FP16
  // inline void load_vec_gather(svbool_t pg, svfloat16_t& v, float16_t* vp,
  //                             svint16_t index) {
  //   v = svld1_gather_s16index_f16(pg, vp, index);
  // }

  // inline void load_add_gather(svbool_t pg, svfloat16_t& v, float16_t* vp,
  //                             svint16_t index){
  //   svfloat16_t v2;
  //   v2 = svld1_gather_s16index_f16(pg, vp, index);
  //   v = svsel_f16(pg, v2, v);
  // }

  // Note that svcompact is not defined for FP16
  inline svreal_t svcompact(svbool_t pg, svreal_t& yt)
  { abort(); }  // dummy to succeed compilation

  inline void flip_sign(svbool_t pg, svfloat16_t& v){
    v = svneg_f16_m(v, pg, v);
  }

  inline void flip_sign(svbool_t pg, svfloat16_t& v1, svfloat16_t &v2){
    v1 = svneg_f16_m(v2, pg, v2);
  }


}

#endif // __ARM_FEATURE_SVE

#endif
