/*!
        @file    solver_CGNE.h

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate: 2013-03-21 15:28:34 #$

        @version $LastChangedRevision: 2314 $
*/

#ifndef SOLVER_CGNE_INCLUDED
#define SOLVER_CGNE_INCLUDED

#include "solver_CG.h"

//! CGNE solver.

/*!
    solve D x = b
    by
    (D Ddag) y = b  and  x = Ddag y.

    delegate to base CG solver.
    NB. iter displayed by vout.detailed does not give #mult of D
        but #mult of DdagD

    Introduce unique_ptr to avoid memory leaks.
                               [21 Mar 2015 Y.Namekawa]
    Add flop_count.            [ 8 Aug 2016 Y.Namekawa]
 */

class Solver_CGNE : public Solver_CG
{
 private:
  Field m_y;

 public:
  static const std::string class_name;

  Solver_CGNE(Fopr *fopr)
    : Solver_CG(fopr) {}

  Solver_CGNE(Fopr *fopr, const Parameters& params)
    : Solver_CG(fopr, params)
  {
    // set_parameters delegate to Solver_CG.
  }

  ~Solver_CGNE() {}

  void set_parameters(const Parameters& params);

  void get_parameters(Parameters& params) const;

  void solve(Field& solution, const Field& source, int& Nconv, double& diff);

  Fopr *get_fopr() { return this->Solver_CG::get_fopr(); }

  double flop_count();

#ifdef USE_FACTORY
 private:
  static Solver *create_object(Fopr *fopr)
  {
    return new Solver_CGNE(fopr);
  }

  static Solver *create_object_with_params(Fopr *fopr, const Parameters& params)
  {
    return new Solver_CGNE(fopr, params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Solver::Factory::Register("CGNE", create_object);
    init &= Solver::Factory_params::Register("CGNE", create_object_with_params);
    return init;
  }
#endif
};
#endif
