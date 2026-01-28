/*!
        @file    fprop_alt.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
        @version $LastChangedRevision: 2499 $
*/

#include "lib_alt/Measurements/Fermion/fprop_alt.h"

#include "lib_alt_QXS/Field/afield.h"


// explicit instanciation for AField<double>.

template class Fprop_alt<AField<double, QXS> >;

template class Fprop_alt<AField<float, QXS> >;

//============================================================END=====
