/*!
        @file    fopr_Overlap_5d.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-03-07 17:24:38 #$

        @version $LastChangedRevision: 2359 $
*/

#ifndef FOPR_OVERLAP_5D_INCLUDED
#define FOPR_OVERLAP_5D_INCLUDED

#include "Fopr/afopr_Overlap_5d.h"
#include "Field/field.h"

//! 5-dimensional overlap fermion operator.

/*!
    This class is the 5-dimensional overlap operator used in
    5D overlap solver.
    After ver.2.0, Fopr_Overlap is an instantiation of
    class template AFopr_Sign<AFIELD> for Field class.
                                      [13 Feb 2022 H.Matsufuru]
 */

typedef AFopr_Overlap_5d<Field> Fopr_Overlap_5d;

#endif
