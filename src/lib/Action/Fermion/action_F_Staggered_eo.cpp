/*!
        @file    action_F_Staggered_eo.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "action_F_Staggered_eo.h"

#include "Field/field_F_1spinor.h"

#include "Smear/smear_APE.h"
#include "Smear/forceSmear_APE.h"


const std::string Action_F_Staggered_eo::class_name = "Action_F_Staggered_eo";

//====================================================================
void Action_F_Staggered_eo::init()
{
  const std::string str_proj_type = "Stout_SU3";
  m_proj = Projection::New(str_proj_type);

  const std::string str_smear_type = "APE";
  m_smear       = Smear::New(str_smear_type, m_proj);
  m_force_smear = ForceSmear::New(str_smear_type, m_proj);

  const std::string str_solver_type = "BiCGStab_Cmplx";
  m_solver = Solver::New(str_solver_type, m_fopr);
}


//====================================================================
void Action_F_Staggered_eo::tidyup()
{
  delete m_proj;
  delete m_smear;
  delete m_force_smear;
  delete m_solver;
}


//====================================================================
void Action_F_Staggered_eo::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double rho;
  int    Nsmear;

  int err = 0;
  err += params.fetch_double("rho_uniform", rho);
  err += params.fetch_int("number_of_smearing", Nsmear);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(Nsmear, rho);
}


//====================================================================
void Action_F_Staggered_eo::get_parameters(Parameters& params) const
{
  params.set_double("rho_uniform", m_rho);
  params.set_int("number_of_smearing", m_Nsmear);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Action_F_Staggered_eo::set_parameters(const int Nsmear, const double rho)
{
  //- print input parameters
  vout.general(m_vl, "%s:\n", class_name.c_str());
  vout.general(m_vl, "  Nsmear = %4d\n", Nsmear);
  vout.general(m_vl, "  rho    = %8.4f\n", rho);

  //- range check
  // NB. Nsmear,rho == 0 is allowed.

  //- store values
  m_Nsmear = Nsmear;
  m_rho    = rho;

  //- post-process
  const int Nc   = CommonParameters::Nc();
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();
  const int NinG = 2 * Nc * Nc;

#if 1
  Parameters params_smear;

  params_smear.set_double("rho_uniform", rho);
  params_smear.set_string("verbose_level", Bridge::BridgeIO::get_verbose_level(m_vl));

  m_smear->set_parameters(params_smear);
  m_force_smear->set_parameters(params_smear);
#else
  dynamic_cast<Smear_APE *>(m_smear)->set_parameters(rho);
  dynamic_cast<ForceSmear_APE *>(m_force_smear)->set_parameters(rho);
#endif

  m_Usmear.resize(m_Nsmear);
  for (int i_smear = 0; i_smear < m_Nsmear; ++i_smear) {
    m_Usmear[i_smear].reset(NinG, Nvol, Ndim);
  }
}


//====================================================================
double Action_F_Staggered_eo::langevin(RandomNumbers *rand)
{
  const int Nc   = CommonParameters::Nc();
  const int Nd   = CommonParameters::Nd();
  const int Ndim = CommonParameters::Ndim();

  const int Nvol  = CommonParameters::Nvol();
  const int Nvol2 = Nvol / 2;
  const int NPE   = CommonParameters::NPE();

  Field_F_1spinor xi;

  rand->gauss_lex_global(xi);

  m_psf.reset(2 * Nc, Nvol2, 1);

  Field_G Usmear(Nvol, Ndim);
  if (m_Nsmear == 0) {
    Usmear = *m_U;
  }

  for (int i_smear = 0; i_smear < m_Nsmear; ++i_smear) {
    Field_G Uprev(Nvol, Ndim);
    if (i_smear == 0) {
      Uprev = *m_U;
    } else {
      Uprev = (Field_G)m_Usmear[i_smear - 1];
    }
    m_smear->smear(Usmear, Uprev);
    m_Usmear[i_smear] = (Field)Usmear;
  }

  m_fopr->set_config(&Usmear);

  Field_F_1spinor xi_e(Nvol2, 1);
  m_index_eo.convertField(xi_e, xi, 0);

  Field_F_1spinor xi_o(Nvol2, 1);
  m_index_eo.convertField(xi_o, xi, 1);

  m_fopr->mult_dag(m_psf, xi_e, "Dee");
  m_fopr->mult_dag(xi_e, xi_o, "Deo");
  axpy(m_psf, 1.0, xi_e);

  const double H_psf = calcH();

  return H_psf;
}


//====================================================================
double Action_F_Staggered_eo::calcH()
{
  const int Nc   = CommonParameters::Nc();
  const int Ndim = CommonParameters::Ndim();

  const int Nvol = CommonParameters::Nvol();
  const int NPE  = CommonParameters::NPE();

  const int    Niter     = 100;
  const int    Nrestart  = 40;
  const double Stop_cond = 1.0e-24;

  Field_G Usmear(Nvol, Ndim);

  if (m_Nsmear == 0) {
    Usmear = *m_U;
  }

  for (int i_smear = 0; i_smear < m_Nsmear; ++i_smear) {
    Field_G Uprev(Nvol, Ndim);
    if (i_smear == 0) {
      Uprev = *m_U;
    } else {
      Uprev = (Field_G)m_Usmear[i_smear - 1];
    }
    m_smear->smear(Usmear, Uprev);
    m_Usmear[i_smear] = (Field)Usmear;
  }

  m_fopr->set_config(&Usmear);

  m_fopr->set_mode("D");

  Field v1(m_psf);
  assert(v1.nvol() == Nvol / 2);
  int    Nconv;
  double diff;
  m_solver->set_parameters(Niter, Nrestart, Stop_cond);
  m_solver->solve(v1, m_psf, Nconv, diff);

  vout.general(m_vl, "  Nconv = %d\n", Nconv);
  vout.general(m_vl, "  diff  = %.8e\n", diff);

  const double H_psf = dot(v1, m_psf);

  vout.general(m_vl, "H_Fstaggered_eo = %18.8f\n", H_psf);
  vout.general(m_vl, "H_F/dof    = %18.8f\n", H_psf / (Nvol / 2) / NPE / Nc);


  return H_psf;
}


//====================================================================
void Action_F_Staggered_eo::force(Field& force)
{
  const int Nc   = CommonParameters::Nc();
  const int Ndim = CommonParameters::Ndim();

  const int Nin  = m_U->nin();
  const int Nvol = m_U->nvol();
  const int Nex  = m_U->nex();

  assert(force.nin() == Nin);
  assert(force.nvol() == Nvol);
  assert(force.nex() == Nex);

  const int    Niter     = 100;
  const int    Nrestart  = 40;
  const double Stop_cond = 1.0e-24;

  force.set(0.0);

  vout.general(m_vl, "  %s:\n", class_name.c_str());

  Field_G Usmear(Nvol, Ndim);
  if (m_Nsmear == 0) {
    Usmear = *m_U;
  }

  for (int i_smear = 0; i_smear < m_Nsmear; ++i_smear) {
    Field_G Uprev(Nvol, Ndim);
    if (i_smear == 0) {
      Uprev = *m_U;
    } else {
      Uprev = (Field_G)m_Usmear[i_smear - 1];
    }
    m_smear->smear(Usmear, Uprev);
    m_Usmear[i_smear] = (Field)Usmear;
  }

  m_fopr->set_config(&Usmear);
  m_fopr->set_mode("D");

  Field_F_1spinor eta(Nvol / 2, 1);
  int             Nconv;
  double          diff;
  m_solver->set_parameters(Niter, Nrestart, Stop_cond);
  m_solver->solve(eta, m_psf, Nconv, diff);

  vout.general(m_vl, "    Solver: Nconv = %6d  diff  = %12.6e\n", Nconv, diff);


  // force of smeared fermion operator
  m_fopr_force->set_config(&Usmear);

  Field_G force1(Nvol, Nex);
  Field_G force2(Nvol, Nex);

  if (m_Nsmear == 0) {
    m_fopr_force->force_core(force1, eta);
  } else {
    m_fopr_force->force_udiv(force2, eta);

    for (int i_smear = m_Nsmear - 1; i_smear > 0; --i_smear) {
      Usmear = (Field_G)m_Usmear[i_smear - 1];
      m_force_smear->force_udiv(force1, force2, Usmear);
      force1 = force;
    }

    Usmear = *m_U;
    m_force_smear->force_udiv(force1, force2, Usmear);
    force2 = force1;

    for (int mu = 0; mu < Ndim; ++mu) {
      for (int site = 0; site < Nvol; ++site) {
        Mat_SU_N ut(Nc);
        ut = Usmear.mat(site, mu) * force2.mat(site, mu);
        ut.at();
        force1.set_mat(site, mu, ut);
      }
    }
  }

#pragma omp parallel
  {
    axpy(force, -2.0, force1);
#pragma omp barrier
    double Fave, Fmax, Fdev;
    force.stat(Fave, Fmax, Fdev);
    vout.general(m_vl, "  %s:\n", class_name.c_str());
    vout.general(m_vl, "    Fave = %12.6f  Fmax = %12.6f  Fdev = %12.6f\n",
                 Fave, Fmax, Fdev);
  }
}


//====================================================================
//============================================================END=====
