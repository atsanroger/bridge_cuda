/*!
        @file    complexTraits.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-12-24 19:34:22 #$
        @version $LastChangedRevision: 2562 $
*/

#ifndef COMPLEX_TRAITS_INCLUDED
#define COMPLEX_TRAITS_INCLUDED

#include "bridge_complex.h"

template<typename REALTYPE>
class ComplexTraits;

template<>
class ComplexTraits<double> {
 public:
  typedef dcomplex complex_t;
};

template<>
class ComplexTraits<float> {
 public:
  typedef fcomplex complex_t;
};

#ifdef USE_FP16
template<>
class ComplexTraits<_Float16> {
  //class ComplexTraits<__fp16> {
 public:
  typedef hcomplex complex_t;
};
#endif

#endif
