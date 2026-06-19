/*!
        @file    afprop.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#include "lib_alt/Measurements/Fermion/afprop.h"

#include "lib_alt_SIMD2/Field/afield.h"


// explicit instanciation for AField<double>.

template class AFprop<AField<double,SIMD2> >;

template class AFprop<AField<float,SIMD2> >;

//============================================================END=====
