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

template<typename T> struct ConstantMemoryTraits;

template<> struct ConstantMemoryTraits<double> {
    static __device__ __forceinline__ const double* e() { return const_e_double; }
    static __device__ __forceinline__ const double* f() { return const_f_double; }
    static __device__ __forceinline__ const double* dpinv() { return const_dpinv_double; }
    static __device__ __forceinline__ const double* dm() { return const_dm_double; }
    static __device__ __forceinline__ const double* b() { return const_b_double; }
    static __device__ __forceinline__ const double* c() { return const_c_double; }
};

template<> struct ConstantMemoryTraits<float> {
    static __device__ __forceinline__ const float* e() { return const_e_float; }
    static __device__ __forceinline__ const float* f() { return const_f_float; }
    static __device__ __forceinline__ const float* dpinv() { return const_dpinv_float; }
    static __device__ __forceinline__ const float* dm() { return const_dm_float; }
    static __device__ __forceinline__ const float* b() { return const_b_float; }
    static __device__ __forceinline__ const float* c() { return const_c_float; }
};

#define const_e     (ConstantMemoryTraits<real_t>::e())
#define const_f     (ConstantMemoryTraits<real_t>::f())
#define const_dpinv (ConstantMemoryTraits<real_t>::dpinv())
#define const_dm    (ConstantMemoryTraits<real_t>::dm())
#define const_b     (ConstantMemoryTraits<real_t>::b())
#define const_c     (ConstantMemoryTraits<real_t>::c())

#endif