/*!
        @file    fopr_Domainwall_eo.h

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-12-17 11:56:41 #$

        @version $LastChangedRevision: 2424 $
*/

#ifndef FOPR_DOMAINWALL_EO_INCLUDED
#define FOPR_DOMAINWALL_EO_INCLUDED

#include "lib/Fopr/afopr_Domainwall_eo.h"

//! Domain-wall (even-odd) fermion operator.

/*!
    This class implements the even-odd domain-wall fermion
    operator.
    It was first implemented in alternative code branch.
    In ver.2.0, it was merged to trunk by renaming to
    AFopr_Domainwall and Fopr_Domainwall a specialization to
    AFIELD = Field.
                                    [06 Mar 2022 H.Matsufuru]
 */

typedef AFopr_Domainwall_eo<Field> Fopr_Domainwall_eo;

#endif
