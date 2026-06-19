/*!
        @file    fprop_Standard_lex.cpp

        @brief

        @author  Satoru Ueda  (ueda)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "fprop_Standard_lex.h"
#include "Tools/timer.h"

const std::string Fprop_Standard_lex::class_name = "Fprop_Standard_lex";

//====================================================================
void Fprop_Standard_lex::set_config(Field *U)
{
  m_solver->get_fopr()->set_config(U);
}


//====================================================================
void Fprop_Standard_lex::invert_D(Field& xq, const Field& b, int& Nconv, double& diff)
{
  m_solver->get_fopr()->set_mode("D");

#pragma omp parallel
  {
    m_solver->solve(xq, b, Nconv, diff);
  }
}


//====================================================================
void Fprop_Standard_lex::invert_DdagD(Field& xq, const Field& b, int& Nconv, double& diff)
{
  m_solver->get_fopr()->set_mode("DdagD");

#pragma omp parallel
  {
    m_solver->solve(xq, b, Nconv, diff);
  }
}


//====================================================================
double Fprop_Standard_lex::flop_count()
{
  return m_solver->flop_count();
}


//====================================================================
void Fprop_Standard_lex::mult_performance(const std::string mode,
                                          const int Nrepeat)
{
  Fopr *fopr = m_solver->get_fopr();

  int nin  = fopr->field_nin();
  int nvol = fopr->field_nvol();
  int nex  = fopr->field_nex();

  Field axq(nin, nvol, nex), abq(nin, nvol, nex);
  abq.set(0.0);
  abq.set(0, 1.0);

  unique_ptr<Timer> timer(new Timer);

  std::string mode_prev = fopr->get_mode();

  fopr->set_mode(mode);

  timer->start();

#pragma omp parallel
  {
    for (int i = 0; i < Nrepeat; ++i) {
      fopr->mult(axq, abq);
      fopr->mult(abq, axq);
    }
  }

  timer->stop();

  double flop_fopr = fopr->flop_count() * 1.e9;
  // note that fopr->flop_count() is in units of GFlop
  double flop_total = flop_fopr * double(2 * Nrepeat);

  double elapsed_time = timer->elapsed_sec();
  double flops        = flop_total / elapsed_time;
  double gflops       = flops * 1.0e-9;

  vout.general(m_vl, "\n");
  vout.general(m_vl, "%s: mult performance:\n", class_name.c_str());
  vout.general(m_vl, "  mult mode = %s\n", mode.c_str());
  vout.general(m_vl, "  Elapsed time = %14.6f sec\n", elapsed_time);
  vout.general(m_vl, "  Flop(Fopr)   = %18.0f\n", flop_fopr);
  vout.general(m_vl, "  Flop(total)  = %18.0f\n", flop_total);
  vout.general(m_vl, "  Performance  = %11.3f GFlops\n", gflops);

  fopr->set_mode(mode_prev);
}


//============================================================END=====
