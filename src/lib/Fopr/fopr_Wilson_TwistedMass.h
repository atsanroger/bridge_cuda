/*!
        @file    fopr_Wilson_TwistedMass.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-12-17 11:56:41 #$

        @version $LastChangedRevision: 2424 $
*/

#ifndef FOPR_TMWILSON_INCLUDED
#define FOPR_TMWILSON_INCLUDED

#include "Fopr/afopr_Wilson_TwistedMass.h"

//! Twisted-mass Wilson fermion operator.

/*!
    This is a twisted mass fermion operator.
    After ver.2.0, Fopr_Wilson_TwistedMass is an instantiation of
    class template AFopr_Wilson_TwistedMass<AFIELD> for Field class.
                                  [19 Dec 2021 H.Matsufuru]
 */


typedef AFopr_Wilson_TwistedMass<Field> Fopr_Wilson_TwistedMass;


#endif
