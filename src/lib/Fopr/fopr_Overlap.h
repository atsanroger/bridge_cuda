/*!
        @file    fopr_Overlap.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/


#ifndef FOPR_OVERLAP_INCLUDED
#define FOPR_OVERLAP_INCLUDED

#include "Fopr/afopr_Overlap.h"
#include "Field/field.h"

//! Overlap fermion operator.

/*!
    This class implements the overlap fermion operator with a given
    kernel fermion operator.
    After ver.2.0, Fopr_Overlap is an instantiation of
    class template AFopr_Sign<AFIELD> for Field class.
                                          [12 Feb 2022 H.Matsufuru]
 */

typedef AFopr_Overlap<Field> Fopr_Overlap;

#endif
