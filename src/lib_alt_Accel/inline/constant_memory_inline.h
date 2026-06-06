#ifndef CONSTANT_MEMORY_INLINE_H
#define CONSTANT_MEMORY_INLINE_H

/*!
      @file    constant_memory_inline.h
      @brief
      @author  Wei-Lun Chen (wlchen)
      @date    $LastChangedDate: 2025-04-08 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "constant_memory.h"

extern __constant__ double const_e_double[NS_MAX];
extern __constant__ double const_f_double[NS_MAX];
extern __constant__ double const_dpinv_double[NS_MAX];
extern __constant__ double const_dm_double[NS_MAX];
extern __constant__ double const_b_double[NS_MAX];
extern __constant__ double const_c_double[NS_MAX];

extern __constant__ float const_e_float[NS_MAX];
extern __constant__ float const_f_float[NS_MAX];
extern __constant__ float const_dpinv_float[NS_MAX];
extern __constant__ float const_dm_float[NS_MAX];
extern __constant__ float const_b_float[NS_MAX];
extern __constant__ float const_c_float[NS_MAX];

// Low words of the float-float coefficient pairs (paired with const_*_float).
extern __constant__ float const_e_ff_lo[NS_MAX];
extern __constant__ float const_f_ff_lo[NS_MAX];
extern __constant__ float const_dpinv_ff_lo[NS_MAX];
extern __constant__ float const_dm_ff_lo[NS_MAX];
extern __constant__ float const_b_ff_lo[NS_MAX];
extern __constant__ float const_c_ff_lo[NS_MAX];

// Triple-word (QTW) coefficient mid/lo arrays. hi reuses const_*_float.
extern __constant__ float const_e_tw_mid[NS_MAX];
extern __constant__ float const_f_tw_mid[NS_MAX];
extern __constant__ float const_dpinv_tw_mid[NS_MAX];
extern __constant__ float const_dm_tw_mid[NS_MAX];
extern __constant__ float const_b_tw_mid[NS_MAX];
extern __constant__ float const_c_tw_mid[NS_MAX];

extern __constant__ float const_e_tw_lo[NS_MAX];
extern __constant__ float const_f_tw_lo[NS_MAX];
extern __constant__ float const_dpinv_tw_lo[NS_MAX];
extern __constant__ float const_dm_tw_lo[NS_MAX];
extern __constant__ float const_b_tw_lo[NS_MAX];
extern __constant__ float const_c_tw_lo[NS_MAX];

// Precomputed coupling a[j] = 0.5*dm[j+1]*dpinv[j] (hi / ff-lo / tw mid+lo).
extern __constant__ double const_a_double[NS_MAX];
extern __constant__ float  const_a_float[NS_MAX];
extern __constant__ float  const_a_ff_lo[NS_MAX];
extern __constant__ float  const_a_tw_mid[NS_MAX];
extern __constant__ float  const_a_tw_lo[NS_MAX];

// Double-build triple-word (double-triple) mid/lo words. hi reuses const_*_double.
extern __constant__ double const_e_td_mid[NS_MAX];
extern __constant__ double const_f_td_mid[NS_MAX];
extern __constant__ double const_dpinv_td_mid[NS_MAX];
extern __constant__ double const_dm_td_mid[NS_MAX];
extern __constant__ double const_b_td_mid[NS_MAX];
extern __constant__ double const_c_td_mid[NS_MAX];
extern __constant__ double const_e_td_lo[NS_MAX];
extern __constant__ double const_f_td_lo[NS_MAX];
extern __constant__ double const_dpinv_td_lo[NS_MAX];
extern __constant__ double const_dm_td_lo[NS_MAX];
extern __constant__ double const_b_td_lo[NS_MAX];
extern __constant__ double const_c_td_lo[NS_MAX];
extern __constant__ double const_a_td_mid[NS_MAX];
extern __constant__ double const_a_td_lo[NS_MAX];

template<typename T> struct ConstantMemoryTraits;

template<> struct ConstantMemoryTraits<double> {
    static __device__ __forceinline__ const double* e() { return const_e_double; }
    static __device__ __forceinline__ const double* f() { return const_f_double; }
    static __device__ __forceinline__ const double* dpinv() { return const_dpinv_double; }
    static __device__ __forceinline__ const double* dm() { return const_dm_double; }
    static __device__ __forceinline__ const double* b() { return const_b_double; }
    static __device__ __forceinline__ const double* c() { return const_c_double; }
    // Double-double (QDW) low words: the DD low word of each coefficient.
    static __device__ __forceinline__ const double* e_lo() { return const_e_td_mid; }
    static __device__ __forceinline__ const double* f_lo() { return const_f_td_mid; }
    static __device__ __forceinline__ const double* dpinv_lo() { return const_dpinv_td_mid; }
    static __device__ __forceinline__ const double* dm_lo() { return const_dm_td_mid; }
    static __device__ __forceinline__ const double* b_lo() { return const_b_td_mid; }
    static __device__ __forceinline__ const double* c_lo() { return const_c_td_mid; }
    // Triple-word (QTW) mid = DD low word; lo = 0 (reserved for a future
    // triple-double host upgrade). Together hi+mid+lo = genuine double-double.
    static __device__ __forceinline__ const double* e_mid() { return const_e_td_mid; }
    static __device__ __forceinline__ const double* f_mid() { return const_f_td_mid; }
    static __device__ __forceinline__ const double* dpinv_mid() { return const_dpinv_td_mid; }
    static __device__ __forceinline__ const double* dm_mid() { return const_dm_td_mid; }
    static __device__ __forceinline__ const double* b_mid() { return const_b_td_mid; }
    static __device__ __forceinline__ const double* c_mid() { return const_c_td_mid; }
    static __device__ __forceinline__ const double* e_lo_tw() { return const_e_td_lo; }
    static __device__ __forceinline__ const double* f_lo_tw() { return const_f_td_lo; }
    static __device__ __forceinline__ const double* dpinv_lo_tw() { return const_dpinv_td_lo; }
    static __device__ __forceinline__ const double* dm_lo_tw() { return const_dm_td_lo; }
    static __device__ __forceinline__ const double* b_lo_tw() { return const_b_td_lo; }
    static __device__ __forceinline__ const double* c_lo_tw() { return const_c_td_lo; }
    // Precomputed coupling a: hi / dd-lo / td mid+lo (now populated for double).
    static __device__ __forceinline__ const double* a() { return const_a_double; }
    static __device__ __forceinline__ const double* a_lo() { return const_a_td_mid; }
    static __device__ __forceinline__ const double* a_mid() { return const_a_td_mid; }
    static __device__ __forceinline__ const double* a_lo_tw() { return const_a_td_lo; }
};

template<> struct ConstantMemoryTraits<float> {
    static __device__ __forceinline__ const float* e() { return const_e_float; }
    static __device__ __forceinline__ const float* f() { return const_f_float; }
    static __device__ __forceinline__ const float* dpinv() { return const_dpinv_float; }
    static __device__ __forceinline__ const float* dm() { return const_dm_float; }
    static __device__ __forceinline__ const float* b() { return const_b_float; }
    static __device__ __forceinline__ const float* c() { return const_c_float; }
    // Low words of the float-float coefficient pairs.
    static __device__ __forceinline__ const float* e_lo() { return const_e_ff_lo; }
    static __device__ __forceinline__ const float* f_lo() { return const_f_ff_lo; }
    static __device__ __forceinline__ const float* dpinv_lo() { return const_dpinv_ff_lo; }
    static __device__ __forceinline__ const float* dm_lo() { return const_dm_ff_lo; }
    static __device__ __forceinline__ const float* b_lo() { return const_b_ff_lo; }
    static __device__ __forceinline__ const float* c_lo() { return const_c_ff_lo; }
    // Triple-word coefficient accessors (h reuses base e()/f()/...).
    static __device__ __forceinline__ const float* e_mid() { return const_e_tw_mid; }
    static __device__ __forceinline__ const float* f_mid() { return const_f_tw_mid; }
    static __device__ __forceinline__ const float* dpinv_mid() { return const_dpinv_tw_mid; }
    static __device__ __forceinline__ const float* dm_mid() { return const_dm_tw_mid; }
    static __device__ __forceinline__ const float* b_mid() { return const_b_tw_mid; }
    static __device__ __forceinline__ const float* c_mid() { return const_c_tw_mid; }
    static __device__ __forceinline__ const float* e_lo_tw() { return const_e_tw_lo; }
    static __device__ __forceinline__ const float* f_lo_tw() { return const_f_tw_lo; }
    static __device__ __forceinline__ const float* dpinv_lo_tw() { return const_dpinv_tw_lo; }
    static __device__ __forceinline__ const float* dm_lo_tw() { return const_dm_tw_lo; }
    static __device__ __forceinline__ const float* b_lo_tw() { return const_b_tw_lo; }
    static __device__ __forceinline__ const float* c_lo_tw() { return const_c_tw_lo; }
    // Precomputed coupling a[j] = 0.5*dm[j+1]*dpinv[j]: hi / ff-lo / tw mid+lo.
    static __device__ __forceinline__ const float* a() { return const_a_float; }
    static __device__ __forceinline__ const float* a_lo() { return const_a_ff_lo; }
    static __device__ __forceinline__ const float* a_mid() { return const_a_tw_mid; }
    static __device__ __forceinline__ const float* a_lo_tw() { return const_a_tw_lo; }
};

#define const_e     (ConstantMemoryTraits<real_t>::e())
#define const_f     (ConstantMemoryTraits<real_t>::f())
#define const_dpinv (ConstantMemoryTraits<real_t>::dpinv())
#define const_dm    (ConstantMemoryTraits<real_t>::dm())
#define const_b     (ConstantMemoryTraits<real_t>::b())
#define const_c     (ConstantMemoryTraits<real_t>::c())

#endif