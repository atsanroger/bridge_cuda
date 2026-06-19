/*!
        @file    solver_CGNR.h

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2015-03-23 23:00:52 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef SOLVER_CGNR_INCLUDED
#define SOLVER_CGNR_INCLUDED

#include "solver_CG.h"

//! CGNR solver.

/*!
    solve D x = b
    by
    (Ddag D) x = b'  and  b' = Ddag b.

    delegate to base CG solver.
    NB. iter displayed by vout.detailed does not give #mult of D
        but #mult of DdagD

    Introduce unique_ptr to avoid memory leaks.
                                   [21 Mar 2015 Y.Namekawa]
    Add flop_count.                [ 8 Aug 2016 Y.Namekawa]
 */

class Solver_CGNR : public Solver_CG
{
 private:
  Field m_b2;

 public:
  static const std::string class_name;

  Solver_CGNR(Fopr *fopr)
    : Solver_CG(fopr) {}

  Solver_CGNR(Fopr *fopr, const Parameters& params)
    : Solver_CG(fopr, params)
  {
    // set_parameters delegate to Solver_CG.
  }

  ~Solver_CGNR() {}

  void set_parameters(const Parameters& params);

  void get_parameters(Parameters& params) const;

  void solve(Field& solution, const Field& source, int& Nconv, double& diff);

  Fopr *get_fopr() { return this->Solver_CG::get_fopr(); }

  double flop_count();

#ifdef USE_FACTORY
 private:
  static Solver *create_object(Fopr *fopr)
  {
    return new Solver_CGNR(fopr);
  }

  static Solver *create_object_with_params(Fopr *fopr, const Parameters& params)
  {
    return new Solver_CGNR(fopr, params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Solver::Factory::Register("CGNR", create_object);
    init &= Solver::Factory_params::Register("CGNR", create_object_with_params);
    return init;
  }
#endif
};
#endif
