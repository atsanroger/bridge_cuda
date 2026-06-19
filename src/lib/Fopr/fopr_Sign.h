/*!
        @file    fopr_Sign.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-03-07 17:24:38 #$

        @version $LastChangedRevision: 2359 $
*/

#ifndef FOPR_SIGN_INCLUDED
#define FOPR_SIGN_INCLUDED

#include "Fopr/afopr_Sign.h"
#include "Field/field.h"

/*!
    This class implements a sign function of a given fermion
    operator emplying Zolotarev approximation.
    After ver.2.0, Fopr_Sign is an instantiation of
    class template AFopr_Sign<AFIELD> for Field class.
                                    [12 Feb 2022 H.Matsufuru]
 */

typedef AFopr_Sign<Field> Fopr_Sign;

#endif
