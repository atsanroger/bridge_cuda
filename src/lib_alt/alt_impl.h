/*!
        @file    alt_impl.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
        @version $LastChangedRevision: 2499 $
*/

#ifndef ALT_IMPL_INCLUDED
#define ALT_IMPL_INCLUDED

enum Impl
{
  CORELIB, SIMD, SIMD2, VECTOR, QXS, ACCEL, OPENACC
};

// alignment for each IMPL
template<Impl IMPL>
constexpr int alignment_size();

#endif
