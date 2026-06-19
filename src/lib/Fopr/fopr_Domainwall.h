/*!
        @file    fopr_Domainwall.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-03-07 17:24:38 #$

        @version $LastChangedRevision: 2359 $
*/

#ifndef FOPR_DOMAINWALL_INCLUDED
#define FOPR_DOMAINWALL_INCLUDED

#include "Fopr/afopr_Domainwall.h"

//! Domain-wall fermion operator.

/*!
    This class implements the domain-wall fermion operator.
    The first version was implemented only for the standard
    (Shamir's) form.                [24 Dec 2011 H.Matsufuru]
    Later it was generalized to cover the Mobius form as a
    template class in the alternative branch.
    In ver.2.0, it was merged to trunk by renaming to
    AFopr_Domainwall and Fopr_Domainwall a specialization to
    AFIELD = Field.
                                    [06 Mar 2022 H.Matsufuru]
 */

typedef AFopr_Domainwall<Field> Fopr_Domainwall;

#endif
