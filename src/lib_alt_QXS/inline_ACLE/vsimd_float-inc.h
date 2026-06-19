/*!
        @file    vsimd_float-inc.h
        @brief
        @author  Hideo Matsufuru (matufuru)
	$LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-12-24 19:34:22 #$
        @version $LastChangedRevision: 2562 $
*/

#ifndef QXS_VSIMD_INCLUDED
#define QXS_VSIMD_INCLUDED

typedef struct{
  float v[VLENS];
} Vsimd_t;

#ifdef __ARM_FEATURE_SVE
#include <arm_sve.h>

typedef svfloat32_t svreal_t;

typedef svint32_t svint_t;

typedef svuint32_t svuint_t;

typedef int32_t int_t;

typedef uint32_t uint_t;

namespace {

  inline svbool_t set_predicate(){
    return svptrue_b32();
  }

  inline svbool_t set_predicate_false(){
    return svpfalse();
  }

  inline svbool_t set_predicate_whilelt(int range){
    return svwhilelt_b32(0, range);
  }

  inline void set_predicate_xp(svbool_t& pg1, svbool_t& pg2){
    svbool_t pg0 = svptrue_b32();
    pg2 = svpfalse();
    for(int iy = VLENYS; iy > 0; --iy){
      pg1 = svwhilelt_b32(0, iy * VLENXS);
      pg2 = sveor_z(pg0, pg2, pg1);
      pg1 = svwhilelt_b32(0, iy * VLENXS - 1);
      pg2 = sveor_z(pg0, pg2, pg1);
    }
    pg1 = svnot_z(pg0, pg2);
  }

  inline void set_predicate_xm(svbool_t& pg1, svbool_t& pg2){
    svbool_t pg0 = svptrue_b32();
    pg2 = svwhilelt_b32(0, VLENS - VLENXS + 1);
    for(int iy = VLENYS-1; iy > 0; --iy){
      pg1 = svwhilelt_b32(0, iy * VLENXS);
      pg2 = sveor_z(pg0, pg2, pg1);
      pg1 = svwhilelt_b32(0, (iy-1) * VLENXS + 1);
      pg2 = sveor_z(pg0, pg2, pg1);
    }
    pg1 = svnot_z(pg0, pg2);
  }

  inline void set_predicate_xp_eo(svbool_t& pg1, svbool_t& pg2,
                                  svbool_t& pg3, int ieo) {
    svbool_t pg0 = svptrue_b32();
    pg2 = svpfalse();
    pg3 = svpfalse();
    for(int iy = VLENYS; iy > 0; --iy){
      if(iy % 2 == ieo){
        pg1 = svwhilelt_b32(0, iy * VLENXS);
        pg2 = sveor_z(pg0, pg2, pg1);
        pg1 = svwhilelt_b32(0, iy * VLENXS - 1);
        pg2 = sveor_z(pg0, pg2, pg1);
      }else{
        pg1 = svwhilelt_b32(0, iy * VLENXS);
        pg3 = sveor_z(pg0, pg3, pg1);
        pg1 = svwhilelt_b32(0, (iy-1) * VLENXS);
        pg3 = sveor_z(pg0, pg3, pg1);
      }
    }
    pg1 = sveor_z(pg0, pg2, pg3);
    pg1 = svnot_z(pg0, pg1);
  }

  inline void set_predicate_xm_eo(svbool_t& pg1, svbool_t& pg2,
                                  svbool_t& pg3, int ieo) {
    svbool_t pg0 = svptrue_b32();
    pg1 = svpfalse();
    pg3 = svpfalse();
    for(int iy = VLENYS; iy > 0; --iy){
      if(iy % 2 == ieo){
        pg2 = svwhilelt_b32(0, iy * VLENXS);
        pg3 = sveor_z(pg0, pg3, pg2);
        pg2 = svwhilelt_b32(0, (iy-1) * VLENXS);
        pg3 = sveor_z(pg0, pg3, pg2);
      }else{
        pg2 = svwhilelt_b32(0, iy * VLENXS);
        pg1 = sveor_z(pg0, pg1, pg2);
        pg2 = svwhilelt_b32(0, (iy-1) * VLENXS + 1);
        pg1 = sveor_z(pg0, pg1, pg2);
      }
    }
    pg2 = sveor_z(pg0, pg1, pg3);
    pg2 = svnot_z(pg0, pg2);
  }

  inline void set_predicate_yp(svbool_t& pg1, svbool_t& pg2){
    svbool_t pg0 = svptrue_b32();
    pg1 = svwhilelt_b32(0, VLENXS * (VLENYS-1));
    pg2 = svnot_z(pg0, pg1);
  }

  inline void set_predicate_ym(svbool_t& pg1, svbool_t& pg2){
    svbool_t pg0 = svptrue_b32();
    pg2 = svwhilelt_b32(0, VLENXS);
    pg1 = svnot_z(pg0, pg2);
  }

  // Nitadori-san's implementation
  inline void set1_at(const int i, svbool_t &pg){
    svbool_t pg0 = svptrue_b32();
    svbool_t pg1 = svwhilelt_b32(0, i);
    svbool_t pg2 = svwhilelt_b32(0, i+1);
    pg = sveor_z(pg0, pg, pg1);
    pg = sveor_z(pg0, pg, pg2);
  }

  inline void rot1_R(uint32_t *u, const int len = VLENXS)
  {
    uint32_t tmp = u[len-1]; // tail
    for(int i = len-1; i >= 1; --i){
      u[i] = u[i-1];
    }
    u[0] = tmp;
  }

  inline void rot1_L(uint32_t *u, const int len = VLENXS)
  {
    uint32_t tmp = u[0]; // head
    for(int i = 0; i < len-1; ++i){
      u[i] = u[i+1];
    }
    u[len-1] = tmp;
  }

  inline void set_idx_predicate_xp_eo(svbool_t  &pg,
                                      svuint32_t &idx,
                                      const int ieo)
  {
    uint32_t u[VLEN];
    for(int i=0; i<VLEN; i++) u[i] = i;

    pg = svpfalse();
    if(0 == ieo){
      // L-shift odd rows
      for(int i = VLENXS; i < VLEN; i += 2*VLENXS){
	set1_at(i, pg); // 3, 11
	rot1_L(u+i, VLENXS);
      }
    }
    if(1 == ieo){
      // L-shift env rows
      for(int i = 0; i < VLEN; i += 2*VLENXS){
	set1_at(i, pg); // 7, 15
	rot1_L(u+i, VLENXS);
      }
    }
    idx = svld1_u32(svptrue_b32(), u);
  }

  inline void set_idx_predicate_xm_eo(svbool_t  &pg,
                                      svuint32_t &idx,
                                      const int ieo)
  {
    uint32_t u[VLEN];
    for(int i = 0; i < VLEN; ++i) u[i] = i;

    pg = svpfalse();
    if(0 == ieo){
      // R-shift evn rows
      for(int i=0; i<VLEN; i+=2*VLENXS){
	set1_at(i + VLENXS-1, pg); // 3, 11
	rot1_R(u+i, VLENXS);
      }
    }
    if(1 == ieo){
      // R-shift odd rows
      for(int i = VLENXS; i < VLEN; i += 2*VLENXS){
	set1_at(i + VLENXS-1, pg); // 7, 15
	rot1_R(u+i, VLENXS);
      }
    }
    idx = svld1_u32(svptrue_b32(), u);
  }

  inline void set_idx_predicate_yp(svbool_t& pg1, svuint32_t &idx)
  {
    pg1 = svwhilelt_b32(0, VLENXS);
    uint32_t u[VLENS];
    for(int i = 0; i < VLENXS*(VLENYS-1); ++i){
      u[i] = i + VLENXS;
    }
    for(int i = 0; i < VLENXS; ++i){
      u[i+VLENXS*(VLENYS-1)] = i;
    }
    idx = svld1_u32(svptrue_b32(), u);
  }

  inline void set_idx_predicate_ym(svbool_t& pg1, svuint32_t &idx)
  {
    svbool_t pg2 = svwhilelt_b32(0, VLENXS*(VLENYS-1));
    pg1 = svnot_z(svptrue_b32(), pg2);
    uint32_t u[VLENS];
    for(int i = 0; i < VLENXS; ++i){
      u[i] = i + VLENXS*(VLENYS-1);
    }
    for(int i = VLENXS; i < VLENS; ++i){
      u[i] = i - VLENXS;
    }
    idx = svld1_u32(svptrue_b32(), u);
  }


  inline void set_vec(svbool_t pg, svfloat32_t& v, float32_t a){
    v  = svdup_f32_m(v, pg, a);
  }

  inline void clear_vec(svbool_t pg, svfloat32_t& v){
    v  = svdup_f32_m(v, pg, 0.0);
  }

  inline void load_vec(svbool_t pg, svfloat32_t& v, const float32_t* vp){
    v = svld1_f32(pg, vp);
  }

  inline void load_add(svbool_t pg, svfloat32_t& v, float32_t* vp){
    svfloat32_t v2;
    v2 = svld1_f32(pg, vp);
    v = svsel_f32(pg, v2, v);
    //    v  = svmul_m(pg, v, 0.0);
    //    v  = svadd_m(pg, v, v2);
  }

  inline void load_svint(svbool_t pg, svint32_t& v, int32_t* vp){
    v = svld1_s32(pg, vp);
  }

  inline void load_svuint(svbool_t pg, svuint32_t& v, uint32_t* vp){
    v = svld1_u32(pg, vp);
  }

  inline void save_vec(svbool_t pg, float* v, svfloat32_t vt){
    svst1(pg, v, vt);
  }

  inline void load_vec_gather(svbool_t pg, svfloat32_t& v, float32_t* vp,
                              svint32_t index){
    v = svld1_gather_s32index_f32(pg, vp, index);
  }

  inline void load_add_gather(svbool_t pg, svfloat32_t& v, float32_t* vp,
                              svint32_t index){
    svfloat32_t v2;
    v2 = svld1_gather_s32index_f32(pg, vp, index);
    v = svsel_f32(pg, v2, v);
    //    v  = svmul_m(pg, v, 0.0);
    //    v  = svadd_m(pg, v, v2);
  }

  inline void flip_sign(svbool_t pg, svfloat32_t& v){
    v = svneg_f32_m(v, pg, v);
  }

  inline void flip_sign(svbool_t pg, svfloat32_t& v1, svfloat32_t &v2){
    v1 = svneg_f32_m(v2, pg, v2);
  }

}

#endif // __ARM_FEATURE_SVE

#endif
