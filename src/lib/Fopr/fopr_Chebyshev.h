/*!
        @file    fopr_Chebyshev.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-03-07 17:24:38 #$

        @version $LastChangedRevision: 2359 $
*/


#ifndef FOPR_CHEBYSHEV_INCLUDED
#define FOPR_CHEBYSHEV_INCLUDED

#include "Fopr/afopr_Chebyshev.h"

//! Chebyshev polynomial of fermion operator.

/*!
    This class implements a Chrbyshev polynomial of fermion
    operator assuming a usage in eigenvalue solver.
    After ver.2.0, Fopr_Sign is an instantiation of
    class template AFopr_Chebyshev<AFIELD> for Field class.
                                    [05 Mar 2022 H.Matsufuru]
 */

typedef AFopr_Chebyshev<Field> Fopr_Chebyshev;

#endif
