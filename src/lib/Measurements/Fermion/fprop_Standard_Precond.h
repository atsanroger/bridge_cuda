#ifndef FPROP_STANDARD_PRECOND_INCLUDED
#define FPROP_STANDARD_PRECOND_INCLUDED

#include "Measurements/Fermion/fprop.h"
#include "Field/field.h"
#include "Fopr/fopr.h"
//class Fopr;
class Solver;


//! Solving quark propagator with LU preconditioning.

/*!
    This is temporary implementation.
                                        [5 Mar 2012 H.Matsufuru]
 */
class Fprop_Standard_Precond : public Fprop
{
 private:
  Fopr *m_fopr;
  Solver *m_solver;

 public:
  Fprop_Standard_Precond(Fopr *fopr)
    : m_fopr(fopr)
  {
    init();
  }

  ~Fprop_Standard_Precond()
  {
    tidyup();
  }

  void invert_D(Field&, const Field&, int& Nconv, double& diff);
  void invert_DdagD(Field&, const Field&, int& Nconv, double& diff);

  void set_config(Field *);

  double flop_count();

 private:

  void invert_D_prec(Field&, const Field&, int& Nconv, double& diff);
  void invert_D_plain(Field&, const Field&, int& Nconv, double& diff);

  void invert_DdagD_prec(Field&, const Field&, int& Nconv, double& diff);
  void invert_DdagD_plain(Field&, const Field&, int& Nconv, double& diff);

  void init();
  void tidyup();
};

#endif
