/*!
      @file    fprop_alt_QWS.cpp
      @brief
      @author  Issaku Kanamori (kanamori)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-03-20 10:52:44 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef FPROP_QWS_INCLUDED
#define FPROP_QWS_INCLUDED

#include "lib_alt/Measurements/Fermion/fprop_alt.h"

#include "lib/Tools/timer.h"

#include "lib/Fopr/afopr.h"
//#include "lib/Smear/director_Smear.h"
//#include "lib_alt/Solver/asolver.h"


//#ifdef USE_QWSLIB
//#include "lib_alt_QXS/extra/qws_lib.h"
//#endif



//====================================================================
template<typename AFIELD>
class Fprop_alt_QWS : public Fprop_alt<AFIELD>
{
 public:

  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;
  using Fprop_alt<AFIELD>::m_vl;
  using Fprop_alt<AFIELD>::m_mode;

 private:
  AFopr<AFIELD> *m_kernel;
  AFopr<AFIELD> *m_fopr;

  Timer m_timer;
  double m_flop_count;
  double m_elapsed_time;
  int m_Niter_d;
  int m_Niter_s;
  int m_nm;
  int m_Nsap;
  int m_Nconv;
  double m_Stop_cond;
  double m_tol_d;
  double m_tol_s;

  Director_Smear *m_dr_smear;

 public:
  Fprop_alt_QWS(const Parameters& params_fopr,
                const Parameters& params_solver)
    : Fprop_alt<AFIELD>()
  { init(params_fopr, params_solver); }

  ~Fprop_alt_QWS()
  { tidyup(); }

  void set_config(Field *);

  void invert(Field&, const Field&, int&, double&);

  void invert_D(Field&, const Field&, int&, double&);

  void invert_DdagD(Field&, const Field&, int&, double&);

  double flop_count();

  void reset_performance();

  void get_performance(double& flop_count, double& elapsed_time);

  void report_performance();

  void mult_performance(const std::string mode, const int Nrepeat);

 private:
  void init(const Parameters& params_fopr,
            const Parameters& params_solver);

  void init(const Parameters& params_fopr,
            const Parameters& params_solver,
            Director_Smear *dr_smear);

  void set_parameters(const Parameters& params_solver);

  void tidyup();
};

#endif


//============================================================END=====
