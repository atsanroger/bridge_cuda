/*!
        @file    force_F_Overlap_Nf2.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-04-17 11:32:37 #$

        @version $LastChangedRevision: 2513 $
*/

#include "force_F_Overlap_Nf2.h"

const std::string Force_F_Overlap_Nf2::class_name = "Force_F_Overlap_Nf2";

//====================================================================
void Force_F_Overlap_Nf2::set_parameters(const Parameters& params)
{
  std::string vlevel;
  if (!params.fetch_string("verbose_level", vlevel)) {
    m_vl = vout.set_verbose_level(vlevel);
  }

  //- fetch and check input parameters
  double           mq, M0;
  int              Np;
  double           x_min, x_max;
  int              Niter_ms;
  double           Stop_cond_ms;
  std::vector<int> bc;

  int err = 0;
  err += params.fetch_double("quark_mass", mq);
  err += params.fetch_double("domain_wall_height", M0);
  err += params.fetch_int("number_of_poles", Np);
  err += params.fetch_double("lower_bound", x_min);
  err += params.fetch_double("upper_bound", x_max);
  err += params.fetch_int("maximum_number_of_iteration", Niter_ms);
  err += params.fetch_double("convergence_criterion_squared", Stop_cond_ms);
  err += params.fetch_int_vector("boundary_condition", bc);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }


  set_parameters(mq, M0, Np, x_min, x_max, Niter_ms, Stop_cond_ms, bc);
}


//====================================================================
void Force_F_Overlap_Nf2::get_parameters(Parameters& params) const
{
  params.set_double("quark_mass", m_mq);
  params.set_double("domain_wall_height", m_M0);
  params.set_int("number_of_poles", m_Np);
  params.set_double("lower_bound", m_x_min);
  params.set_double("upper_bound", m_x_max);
  params.set_int("maximum_number_of_iteration", m_Niter_ms);
  params.set_double("convergence_criterion_squared", m_Stop_cond_ms);
  params.set_int_vector("boundary_condition", m_boundary);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));
}


//====================================================================
void Force_F_Overlap_Nf2::set_parameters(const double mq, const double M0,
                                         const int Np, const double x_min, const double x_max,
                                         const int Niter_ms, const double Stop_cond_ms,
                                         const std::vector<int> bc)
{
  const int Ndim = CommonParameters::Ndim();

  //- print input parameters
  vout.general(m_vl, "Paramters of %s:\n", class_name.c_str());
  vout.general(m_vl, "  mq           = %12.8f\n", mq);
  vout.general(m_vl, "  M0           = %12.8f\n", M0);
  vout.general(m_vl, "  Np           = %4d\n", Np);
  vout.general(m_vl, "  x_min        = %12.8f\n", x_min);
  vout.general(m_vl, "  x_max        = %12.8f\n", x_max);
  vout.general(m_vl, "  Niter_ms     = %6d\n", Niter_ms);
  vout.general(m_vl, "  Stop_cond_ms = %8.2e\n", Stop_cond_ms);
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
  }

  //- range check
  int err = 0;
  err += ParameterCheck::non_zero(M0);
  err += ParameterCheck::non_zero(mq);
  err += ParameterCheck::non_zero(Np);
  // NB. x_min,x_max == 0 is allowed.
  err += ParameterCheck::non_zero(Niter_ms);
  err += ParameterCheck::square_non_zero(Stop_cond_ms);

  if (err) {
    vout.crucial(m_vl, "Error at %s: parameter range check failed.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  assert(bc.size() == Ndim);

  //- store values
  m_M0           = M0;
  m_mq           = mq;
  m_Np           = Np;
  m_x_min        = x_min;
  m_x_max        = x_max;
  m_Niter_ms     = Niter_ms;
  m_Stop_cond_ms = Stop_cond_ms;

  m_boundary.resize(Ndim);
  m_boundary = bc;

  //- post-process
  m_kappa = 0.5 / (4.0 - m_M0);

  m_fopr_w = new Fopr_Wilson();
  m_fopr_w->set_parameters(m_kappa, m_boundary);
  m_fopr_w->set_config(m_U);

  m_force_w = new Force_F_Wilson_Nf2;
  m_force_w->set_parameters(m_kappa, m_boundary);
  m_force_w->set_config(m_U);

  //- Zolotarev coefficients
  m_sigma.resize(m_Np);
  m_cl.resize(2 * m_Np);
  m_bl.resize(m_Np);
}


//====================================================================
void Force_F_Overlap_Nf2::force_core(Field& force_, const Field& psi)
{
  //- solving psi = (H^2)^(-1) phi
  Field_F psi5;

  m_fopr_w->mult_gm5(psi5, psi);

  double ff = m_M0 * m_M0 - 0.25 * m_mq * m_mq;
  scal(psi5, ff);

  // force determination
  force_core1(force_, psi, psi5);
  scal(force_, -1.0);
}


//====================================================================
void Force_F_Overlap_Nf2::force_udiv(Field& force, const Field& psi)
{
  vout.crucial(m_vl, "Error at %s: not implemented.\n", __func__);
  exit(EXIT_FAILURE);
}


//====================================================================
void Force_F_Overlap_Nf2::force_core1(Field& force_, const Field& psi_, const Field& psi5_)
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  Field_G force(Nvol, Ndim);
  Field_F psi(psi_);
  Field_F psi5(psi5_);

  force_core1_impl(force, psi, psi5);

  force_ = force;
}


//====================================================================
void Force_F_Overlap_Nf2::force_core1_impl(Field_G& force, const Field_F& psi, const Field_F& psi5)
{
  const int Nc   = CommonParameters::Nc();
  const int Nd   = CommonParameters::Nd();
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();
  const int NinF = 2 * Nc * Nd;

  //- Zolotarev coefficient defined
  double              bmax = m_x_max / m_x_min;
  Math_Sign_Zolotarev sign_func(m_Np, bmax);

  sign_func.get_sign_parameters(m_cl, m_bl);

  for (int i = 0; i < m_Np; ++i) {
    m_sigma[i] = m_cl[2 * i] * m_x_min * m_x_min;
  }

  for (int i = 0; i < m_Np; ++i) {
    vout.general(m_vl, " %3d %12.4e %12.4e %12.4e\n",
                 i, m_cl[i], m_cl[i + m_Np], m_bl[i]);
  }

  //- Shiftsolver
  const int Nshift = m_Np;

  std::vector<Field> psi_l(Nshift);
  for (int i = 0; i < Nshift; ++i) {
    psi_l[i].reset(NinF, Nvol, 1);
  }

  std::vector<Field> psi5_l(Nshift);
  for (int i = 0; i < Nshift; ++i) {
    psi5_l[i].reset(NinF, Nvol, 1);
  }

  vout.general(m_vl, "Shift solver in MD\n");
  vout.general(m_vl, "  Number of shift values = %d\n", m_sigma.size());

  m_fopr_w->set_mode("DdagD");

  {
    Shiftsolver_CG solver(m_fopr_w, m_Niter_ms, m_Stop_cond_ms);
    int            Nconv;
    double         diff;
    solver.solve(psi_l, m_sigma, psi, Nconv, diff);
    solver.solve(psi5_l, m_sigma, psi5, Nconv, diff);
  }

  force.set(0.0);

  Field_F psid;
  psid.set(0.0);

  Field_F psi5d;
  psi5d.set(0.0);

  Field_F w1, w2;
  Field_G force1(Nvol, Ndim);

  for (int ip = 0; ip < m_Np; ++ip) {
    psid.addpart_ex(0, psi_l[ip], 0, m_bl[ip]);
    psi5d.addpart_ex(0, psi5_l[ip], 0, m_bl[ip]);

    m_fopr_w->set_mode("H");
    m_fopr_w->mult(w1, psi_l[ip]);
    m_fopr_w->mult(w2, psi5_l[ip]);

    double coeff_l = (m_cl[2 * ip] - m_cl[2 * m_Np - 1]) * m_bl[ip] * m_x_min;

    //- first term
    Field_F w3;
    m_fopr_w->mult(w3, w2);
    Field_F w4 = (Field_F)psi_l[ip];  //-- convert Field to Field_F

    m_force_w->force_core1(force1, w3, w4);
    axpy(force, coeff_l, force1);

    m_force_w->force_core1(force1, w2, w1);
    axpy(force, coeff_l, force1);

    //- second term
    m_fopr_w->mult(w3, w1);
    w4 = (Field_F)psi5_l[ip];  //-- convert Field to Field_F

    m_force_w->force_core1(force1, w3, w4);
    axpy(force, coeff_l, force1);
    m_force_w->force_core1(force1, w1, w2);
    axpy(force, coeff_l, force1);
  }

  double coeff1 = m_cl[2 * m_Np - 1] * m_x_min * m_x_min;
  double coeff2 = 1.0 / m_x_min;

  //- first term
  m_fopr_w->mult(w1, psid);
  m_fopr_w->mult(w2, w1);
  w2.addpart_ex(0, psid, 0, coeff1);

  m_force_w->force_core1(force1, psi5, w2);
  axpy(force, coeff2, force1);

  //- second term
  m_fopr_w->mult(w1, psi5d);
  m_fopr_w->mult(w2, w1);
  w2.addpart_ex(0, psi5d, 0, coeff1);

  m_force_w->force_core1(force1, psi, w2);
  axpy(force, coeff2, force1);
}


//====================================================================
void Force_F_Overlap_Nf2::force_udiv1(Field& force_, const Field& psi_, const Field& psi5_)
{
  vout.crucial(m_vl, "Error at %s: not implemented.\n", __func__);
  exit(EXIT_FAILURE);
}


//====================================================================
//===========================================================END======
