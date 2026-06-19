/*!
        @file    fopr_Smeared_alt.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-03-07 17:24:38 #$

        @version $LastChangedRevision: 2359 $
*/

#ifndef FOPR_SMEARED_ALT_INCLUDED
#define FOPR_SMEARED_ALT_INCLUDED

#include "lib/Fopr/afopr_Smeared_alt.h"

//! Domain-wall fermion operator.

/*!
    Alternative smeared fermion operatior with AFIELD=Field.
                                    [29 Mar 2023 H.Matsufuru]
 */

typedef AFopr_Smeared_alt<Field> Fopr_Smeared_alt;

#endif
