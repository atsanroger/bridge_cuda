#include "Measurements/Fermion/fprop_Standard_Precond.h"

#include "Parameters/commonParameters.h"
#include "Communicator/communicator.h"

#include "Fopr/fopr.h"
#include "Solver/solver_CG.h"


//====================================================================
void Fprop_Standard_Precond::set_config(Field *config)
{
  m_fopr->set_config(config);
}


//====================================================================
void Fprop_Standard_Precond::init()
{
  // temporary setup.
  int       Niter     = 500;
  int       Nrestart  = 100;
  double    Stop_cond = 1.0e-28;
  Solver_CG *solver   = new Solver_CG(m_fopr);
  solver->set_parameters(Niter, Nrestart, Stop_cond);

  m_solver = (Solver *)solver;
  //  m_solver->set_parameter_verboselevel(Bridge::DETAILED);
}


//====================================================================
void Fprop_Standard_Precond::tidyup()
{
  delete m_solver;
}


//====================================================================
void Fprop_Standard_Precond::invert_D(Field& xq, const Field& b,
                                      int& Nconv, double& diff)
{
  assert(xq.nin() == b.nin());
  assert(xq.nvol() == b.nvol());
  assert(xq.nex() == b.nex());

  invert_D_prec(xq, b, Nconv, diff);
  //  invert_D_plain(xq,b);
}


//====================================================================
void Fprop_Standard_Precond::invert_DdagD(Field& xq, const Field& b,
                                          int& Nconv, double& diff)
{
  int NinF = b.nin();
  int Nvol = b.nvol();
  int Nex  = b.nex();
  assert(xq.check_size(NinF, Nvol, Nex));

  Field v(NinF, Nvol, Nex);

  invert_DdagD_prec(xq, b, Nconv, diff);
  // invert_DdagD_plain(xq,b);

  m_fopr->set_mode("DdagD");
  m_fopr->mult(v, xq);
  axpy(v, -1.0, b);
  double vv = v.norm2() / b.norm2();
  vout.detailed(m_vl, "  diff2 = %18.8e\n", vv);
}


//====================================================================
void Fprop_Standard_Precond::invert_D_prec(Field& xq, const Field& b,
                                           int& Nconv, double& diff)
{
  int NinF = b.nin();
  int Nvol = b.nvol();
  int Nex  = b.nex();
  assert(xq.check_size(NinF, Nvol, Nex));

  Field v2(NinF, Nvol, Nex);

  m_fopr->set_mode("DdagD_prec");

  //  m_fopr->set_mode("D_prec");
  m_fopr->mult_dag(v2, b, "D_prec");

  m_solver->solve(xq, v2, Nconv, diff);
  vout.detailed(m_vl, "    Solver: Nconv = %6d  diff  = %12.6e\n",
                Nconv, diff);

  //  m_fopr->set_mode("LU_inv");
  copy(v2, xq);
  m_fopr->mult(xq, v2, "Prec");
}


//====================================================================
void Fprop_Standard_Precond::invert_DdagD_prec(Field& xq, const Field& b,
                                               int& Nconv, double& diff)
{
  int NinF = b.nin();
  int Nvol = b.nvol();
  int Nex  = b.nex();
  assert(xq.check_size(NinF, Nvol, Nex));

  Field v2(NinF, Nvol, Nex);

  //  m_fopr->set_mode("LU_inv");
  //  m_fopr->mult_dag(v2, b);
  m_fopr->mult_dag(v2, b, "Prec");

  m_fopr->set_mode("DdagD_prec");
  m_solver->solve(xq, v2, Nconv, diff);
  vout.detailed(m_vl, "    Solver: Nconv = %6d  diff  = %12.6e\n",
                Nconv, diff);

  //  m_fopr->set_mode("LU_inv");
  copy(v2, xq);
  //  m_fopr->mult(xq, v2);
  m_fopr->mult(xq, v2, "Prec");
}


//====================================================================
void Fprop_Standard_Precond::invert_D_plain(Field& xq, const Field& b,
                                            int& Nconv, double& diff)
{
  int NinF = b.nin();
  int Nvol = b.nvol();
  int Nex  = b.nex();
  assert(xq.check_size(NinF, Nvol, Nex));

  Field v2(NinF, Nvol, Nex);

  m_fopr->set_mode("D");
  m_fopr->mult_dag(v2, b);
  m_fopr->set_mode("DdagD");

  m_solver->solve(xq, v2, Nconv, diff);
  vout.detailed(m_vl, "    Solver: Nconv = %6d  diff  = %12.6e\n",
                Nconv, diff);
}


//====================================================================
void Fprop_Standard_Precond::invert_DdagD_plain(Field& xq, const Field& b,
                                                int& Nconv, double& diff)
{
  m_fopr->set_mode("DdagD");

  m_solver->solve(xq, b, Nconv, diff);
  vout.detailed(m_vl, "    Solver: Nconv = %6d  diff  = %12.6e\n",
                Nconv, diff);
}


//====================================================================
double Fprop_Standard_Precond::flop_count()
{
  //  return m_solver->flop_count();

  return 0;  // temporary unsupported.
}


//====================================================================
//============================================================END=====
