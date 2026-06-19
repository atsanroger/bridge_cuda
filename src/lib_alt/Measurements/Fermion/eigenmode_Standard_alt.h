/*!
        @file    eigenmode_Standard_alt.h
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2022-12-16 15:57:38 #$
        @version $LastChangedRevision: 2422 $
*/

#ifndef EIGENMODE_STANDARD_ALT_INCLUDED
#define EIGENMODE_STANDARD_ALT_INCLUDED

#include "lib/Measurements/Fermion/fprop.h"
#include "lib/Tools/timer.h"
#include "lib/IO/bridgeIO.h"
#include "lib/Eigen/aeigensolver_IRLanczos.h"

#include "lib/Fopr/afopr.h"


//! Detemine the eigenmodes of fermion operator.

/*!
    This is temporary implementation.
                                     [12 Jul 2018 H.Matsufuru]
 */

template<typename AFIELD>
class Eigenmode_Standard_alt
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  AFopr<AFIELD>* m_fopr;
  AEigensolver_IRLanczos<AFIELD, AFopr<AFIELD> >* m_eigensolver;

  int m_Nk, m_Np;

  Timer m_timer;
  double m_flop_count;
  double m_elapsed_time;

  Bridge::VerboseLevel m_vl;

 public:
  Eigenmode_Standard_alt(Parameters& params_fopr,
                         Parameters& params_eigensolver)
   //    : Fprop()
    { init(params_fopr, params_eigensolver); }

  ~Eigenmode_Standard_alt()
    { tidyup(); }

  void set_config(Field *);

  void calc_eigenmode();

  //  void invert_DdagD_prec(Field&, const Field&, int&, double&);

  //  double flop_count();

  //  void reset_performance();

  //  void report_performance();

  //  void mult_performance(const std::string mode, const int Nrepeat);

 private:
  void init(Parameters& params_fopr, Parameters& params_eigensolver);
  void tidyup();

};
#endif
