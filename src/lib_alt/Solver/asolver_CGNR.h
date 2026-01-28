/*!
        @file    asolver_CGNR.h
        @brief
        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate: 2015-03-23 23:00:52 #$
        @version $LastChangedRevision: 2499 $
*/

#ifndef ASOLVER_CGNR_INCLUDED
#define ASOLVER_CGNR_INCLUDED

#include "lib_alt/Solver/asolver_CG.h"

//! CGNR solver.

/*!
    solve D x = b
    by
    (Ddag D) x = b'  and  b' = Ddag b.

    delegate to base CG solver.
    NB. iter displayed by vout.detailed does not give #mult of D
        but #mult of DdagD

    Created by copying from core library.
                                       [21 Nov 2018 H.Matsufuru]
 */

template<typename AFIELD>
class ASolver_CGNR : public ASolver_CG<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

  using ASolver_CG<AFIELD>::m_vl;

 private:
  AFopr<AFIELD> *m_fopr;
  AFIELD m_b2;
  std::string m_mode_fopr;
  double m_flop;

 public:

  ASolver_CGNR(AFopr<AFIELD> *fopr)
    : ASolver_CG<AFIELD>(fopr), m_fopr(fopr)
  { init(); }

  ASolver_CGNR(unique_ptr<AFopr<AFIELD> >& fopr)
    : ASolver_CG<AFIELD>(fopr.get()), m_fopr(fopr.get())
  { init(); }

  ~ASolver_CGNR() {}

  void set_parameters(const Parameters& params);

  void solve(AFIELD& solution, const AFIELD& source,
             int& Nconv, real_t& diff);

  AFopr<AFIELD> *get_fopr() { return this->ASolver_CG<AFIELD>::get_fopr(); }

  double flop_count();

 private:
  void init();

#ifdef USE_FACTORY
 private:
  static ASolver<AFIELD> *create_object_with_fopr(AFopr<AFIELD> *fopr)
  {
    return new ASolver_CGNR<AFIELD>(fopr);
  }

 public:
  static bool register_factory()
  {
    return ASolver<AFIELD>::Factory_fopr::Register("CGNR",
                                                   create_object_with_fopr);
  }
#endif
};
#endif
