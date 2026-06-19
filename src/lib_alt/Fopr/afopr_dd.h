/*!
        @file    afopr_dd.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
        @version $LastChangedRevision: 2499 $
*/

#ifndef AFOPR_DD_INCLUDED
#define AFOPR_DD_INCLUDED

#include "lib/Fopr/afopr.h"

//! Base class of fermion operator family.

/*!
    This class defines the interface of the fermion operators
    with domain-decomposition.
                                     [10 Apr 2021 H.Matsufuru]
*/

template<typename AFIELD>
class AFopr_dd : public AFopr<AFIELD>
{
 public:

  // using AFopr<AFIELD>::m_vl;

  virtual ~AFopr_dd() {}

  //! SAP operator
  virtual void mult_sap(AFIELD&, const AFIELD&, const int eo) = 0;

  //! Mult only inside domain
  virtual void mult_dd(AFIELD&, const AFIELD&) = 0;

  //! Upward hopping part of mult
  virtual void mult_dup(AFIELD&, const AFIELD&, const int mu) = 0;

  //! Downward hopping part of mult
  virtual void mult_ddn(AFIELD&, const AFIELD&, const int mu) = 0;

  //! Returns floating operation counts for SAP mult.
  virtual double flop_count_sap() = 0;
};
#endif
