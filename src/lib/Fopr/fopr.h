/*!
        @file    fopr.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-01-24 20:05:58 #$

        @version $LastChangedRevision: 2452 $
*/


#ifndef FOPR_INCLUDED
#define FOPR_INCLUDED

#include "Fopr/afopr.h"

//! Base class of fermion operator family.

/*!
    This class defines the interface of the fermion operators.
    At present, void functions mult(v,w) and mult_dag(v,w) is
    not purely virtual, because some of subclass have not
    implemented them yet.
                                      [20 Dec 2011 H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks.
                                      [21 Mar 2015 Y.Namekawa]
    parameters factory is introduced. [21 Apr 2015 Y.Namekawa]
    implementation is moved to template class AFopr.
                                      [21 Jan 2023 H.Matsufuru]
*/

typedef AFopr<Field> Fopr;


#endif
