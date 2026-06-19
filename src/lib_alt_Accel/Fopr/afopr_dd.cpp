/*!                                                                            
        @file    afopr.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-04-05 15:35:59 #$
        @version $LastChangedRevision: 2223 $
*/

#include "lib_alt/Fopr/afopr_dd.h"
#include "lib_alt_Accel/Field/afield.h"

// explicit instanciation.
template class AFopr_dd<AField<double,ACCEL> >;
template class AFopr_dd<AField<float,ACCEL> >;

//============================================================END=====
