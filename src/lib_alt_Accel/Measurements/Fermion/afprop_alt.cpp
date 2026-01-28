/*!
        @file    afprop_alt.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#include "lib_alt/Measurements/Fermion/fprop_alt.h"

#include "lib_alt_Accel/Field/afield.h"


// explicit instanciation

template class Fprop_alt<AField<double,ACCEL> >;

template class Fprop_alt<AField<float,ACCEL> >;

//============================================================END=====
