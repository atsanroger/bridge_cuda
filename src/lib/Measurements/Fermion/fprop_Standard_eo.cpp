/*!
        @file    fprop_Standard_eo.cpp

        @brief

        @author  Satoru Ueda (sueda)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-04-17 13:21:11 #$

        @version $LastChangedRevision: 2514 $
*/

#include "fprop_Standard_eo.h"
#include "Tools/timer.h"

const std::string Fprop_Standard_eo::class_name = "Fprop_Standard_eo";

//====================================================================
void Fprop_Standard_eo::set_config(Field *U)
{
  m_index->convertField(*m_Ueo, *U);

  m_solver->get_fopr()->set_config(U);
}


//====================================================================
void Fprop_Standard_eo::invert_D(Field& xq, const Field& b, int& Nconv, double& diff)
{
  const int Nin  = b.nin();
  const int Nvol = b.nvol();
  const int Nex  = b.nex();

  Field Be(Nin, Nvol / 2, Nex);
  Field bo(Nin, Nvol / 2, Nex);
  Field xe(Nin, Nvol / 2, Nex);

  int Nconv1 = 0;

  Fopr_eo *fopr = dynamic_cast<Fopr_eo*>(m_solver->get_fopr());
  if(fopr==nullptr){
    vout.general(m_vl, "Fprop_Standard_eo: the given Fopr is not Fopr_eo.\n");
    abort();
  }

  //  std::string mode_prev = fopr->get_mode();

  fopr->set_mode("D");
#pragma omp parallel
  {
    fopr->preProp(Be, bo, b);
#pragma omp barrier
    m_solver->solve(xe, Be, Nconv1, diff);
#pragma omp barrier
    fopr->postProp(xq, xe, bo);
  } // omp parallel

  //- NB. #mult is doubled for even-odd
  Nconv = 2 * Nconv1;
}


//====================================================================
void Fprop_Standard_eo::invert_DdagD(Field& xq, const Field& b, int& Nconv, double& diff)
{
  const int Nin  = b.nin();
  const int Nvol = b.nvol();
  const int Nex  = b.nex();

  Field Be(Nin, Nvol / 2, Nex);
  Field bo(Nin, Nvol / 2, Nex);
  Field xe(Nin, Nvol / 2, Nex);

  int    Nconv1 = 0, Nconv2 = 0;
  double diff1 = 1.0, diff2 = 1.0;

  Fopr_eo *fopr = dynamic_cast<Fopr_eo*>(m_solver->get_fopr());
  if(fopr==nullptr){
    vout.crucial(m_vl, "Fprop_Standard_eo: the given Fopr is not Fopr_eo.\n");
    abort();
  }

#pragma omp parallel
  {
    fopr->set_mode("Ddag");
    fopr->preProp(Be, bo, b);
#pragma omp barrier
    m_solver->solve(xe, Be, Nconv1, diff1);
#pragma omp barrier
    fopr->postProp(xq, xe, bo);

    fopr->set_mode("D");
    fopr->preProp(Be, bo, xq);
#pragma omp barrier
    m_solver->solve(xe, Be, Nconv2, diff2); 
#pragma omp barrier
    fopr->postProp(xq, xe, bo);
  }


  //- NB. #mult is doubled for even-odd
  Nconv = 2 * (Nconv1 + Nconv2);

  //- rough estimate of diff
  diff = (diff1 + diff2) / 2.0;
}


//====================================================================
double Fprop_Standard_eo::flop_count()
{
  const int    NPE = CommonParameters::NPE();
  const double eps = CommonParameters::epsilon_criterion();

  //- NB1 Nin = 2 * Nc * Nd, Nex = 1  for field_F
  //- NB2 Nvol/2 for eo
  const int Nin  = 2 * CommonParameters::Nc() * CommonParameters::Nd();
  const int Nvol = CommonParameters::Nvol();
  const int Nex  = 1;

  const double flop_fopr = m_solver->get_fopr()->flop_count();

  if (flop_fopr < eps) {
    vout.crucial(m_vl, "Warning at %s: no fopr->flop_count() is available, setting flop = 0.0.\n", class_name.c_str());
    return 0.0;
  }

  const double flop_axpy = static_cast<double>(Nin * Nex * 2) * (Nvol / 2 * NPE);

  const double flop_preProp  = flop_fopr + flop_axpy;
  const double flop_solver   = m_solver->flop_count();
  const double flop_postProp = flop_fopr + flop_axpy;

  const double flop = flop_preProp + 2 * flop_solver + flop_postProp;

  return flop;
}



//====================================================================
void Fprop_Standard_eo::mult_performance(const std::string mode,
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

  double flop_fopr  = fopr->flop_count() * 1.e9;
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
