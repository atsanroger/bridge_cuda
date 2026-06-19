/*!
      @file    aprecond.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#include "lib_alt/Solver/aprecond.h"

#include "lib_alt_SIMD2/Field/afield.h"

// explicit instanciation.
template class APrecond<AField<float,SIMD2> >;
template class APrecond<AField<double,SIMD2> >;

//============================================================END=====
