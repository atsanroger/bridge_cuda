/*!
        @file    afopr.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
        @version $LastChangedRevision: 2499 $
*/

#include "lib_alt/Fopr/afopr_dd.h"
#include "lib_alt_QXS/Field/afield.h"

// explicit instanciation.
template class AFopr_dd<AField<double, QXS> >;
template class AFopr_dd<AField<float, QXS> >;

//============================================================END=====
