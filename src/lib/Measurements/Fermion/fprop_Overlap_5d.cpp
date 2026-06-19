/*!
        @file    fprop_Overlap_5d.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-12-16 15:57:38 #$

        @version $LastChangedRevision: 2422 $
*/

#include "fprop_Overlap_5d.h"

const std::string Fprop_Overlap_5d::class_name = "Fprop_Overlap_5d";

//====================================================================
void Fprop_Overlap_5d::set_parameters(const Parameters& params)
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

  std::string repr, kernel_type;

  int err = 0;
  err += params.fetch_double("quark_mass", mq);
  err += params.fetch_double("domain_wall_height", M0);
  err += params.fetch_int("number_of_poles", Np);
  err += params.fetch_double("lower_bound", x_min);
  err += params.fetch_double("upper_bound", x_max);
  err += params.fetch_int("maximum_number_of_iteration", Niter_ms);
  err += params.fetch_double("convergence_criterion_squared", Stop_cond_ms);
  err += params.fetch_int_vector("boundary_condition", bc);

  err += params.fetch_string("kernel_type", kernel_type);
  err += params.fetch_string("gamma_matrix_type", repr);

  if (err) {
    vout.crucial(m_vl, "Error at %s: input parameter not found.\n", class_name.c_str());
    exit(EXIT_FAILURE);
  }

  vout.general(m_vl, "kernel_type input = %s\n", m_kernel_type.c_str());
  m_kernel_type = kernel_type;

  set_parameters(mq, M0, Np, x_min, x_max, Niter_ms, Stop_cond_ms, bc, repr);
}


//====================================================================
void Fprop_Overlap_5d::set_parameters(const Parameters& params_overlap,
                                      const Parameters& params_solver)
{
  const double mq           = params_overlap.get_double("quark_mass");
  const double M0           = params_overlap.get_double("domain_wall_height");
  const int    Np           = params_overlap.get_int("number_of_poles");
  const double x_min        = params_overlap.get_double("lower_bound");
  const double x_max        = params_overlap.get_double("upper_bound");
  const int    Niter_ms     = params_overlap.get_int("maximum_number_of_iteration");
  const double Stop_cond_ms = params_overlap.get_double("convergence_criterion_squared");

  const std::vector<int> bc = params_overlap.get_int_vector("boundary_condition");

  std::string kernel_type = params_overlap.get_string("kernel_type");
  vout.general(m_vl, "kernel_type input = %s\n", m_kernel_type.c_str());
  m_kernel_type = kernel_type;

  const std::string repr = params_overlap.get_string("gamma_matrix_type");

  set_parameters(mq, M0, Np, x_min, x_max, Niter_ms, Stop_cond_ms, bc, repr);
  m_params_solver = params_solver;
}


//====================================================================
void Fprop_Overlap_5d::get_parameters(Parameters& params) const
{
  params.set_double("quark_mass", m_mq);
  params.set_double("domain_wall_height", m_M0);
  params.set_int("number_of_poles", m_Np);
  params.set_double("lower_bound", m_x_min);
  params.set_double("upper_bound", m_x_max);
  params.set_int("maximum_number_of_iteration", m_Niter_ms);
  params.set_double("convergence_criterion_squared", m_Stop_cond_ms);
  params.set_int_vector("boundary_condition", m_boundary);

  params.set_string("gamma_matrix_type", m_repr);

  params.set_string("verbose_level", vout.get_verbose_level(m_vl));

  return;
}


//====================================================================
void Fprop_Overlap_5d::set_parameters(const double mq, const double M0,
                                      const int Np, const double x_min, const double x_max,
                                      const int Niter_ms, const double Stop_cond_ms,
                                      const std::vector<int> bc,
                                      const std::string& repr
                                      )
{
  const int Ndim = CommonParameters::Ndim();

  //- print input parameters
  vout.general(m_vl, "Paramters of %s:\n", class_name.c_str());
  vout.general(m_vl, "  mq           = %10.6f\n", mq);
  vout.general(m_vl, "  M0           = %10.6f\n", M0);
  vout.general(m_vl, "  Np           = %4d\n", Np);
  vout.general(m_vl, "  x_min        = %12.8f\n", x_min);
  vout.general(m_vl, "  x_max        = %12.6f\n", x_max);
  vout.general(m_vl, "  Niter_ms     = %6d\n", Niter_ms);
  vout.general(m_vl, "  Stop_cond_ms = %8.2e\n", Stop_cond_ms);
  for (int mu = 0; mu < Ndim; ++mu) {
    vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
  }
  vout.general(m_vl, "  kernel type  = %s\n", m_kernel_type.c_str());
  vout.general(m_vl, "  gamma matrix = %s\n", repr.c_str());

  //- range check
  int err = 0;
  err += ParameterCheck::non_zero(mq);
  err += ParameterCheck::non_zero(M0);
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
  m_mq           = mq;
  m_M0           = M0;
  m_Np           = Np;
  m_x_min        = x_min;
  m_x_max        = x_max;
  m_Niter_ms     = Niter_ms;
  m_Stop_cond_ms = Stop_cond_ms;

  m_boundary.resize(Ndim);
  m_boundary = bc;

  m_Nsbt = 0;

  m_repr = repr;
}


//====================================================================
void Fprop_Overlap_5d::initialize()
{
  const int Nvol = CommonParameters::Nvol();
  const int Ndim = CommonParameters::Ndim();

  //- setting even-odd envoronement
  m_Ueo = new Field_G(Nvol, Ndim);
  m_index_eo.convertField(*m_Ueo, *m_U);

  //- setting Wilson fermion operator
  const double kappa = 0.5 / (4.0 - m_M0);

  // m_fopr_w_eo = new Fopr_Wilson_eo(m_repr);
  // m_fopr_w_eo->set_parameters(kappa, m_boundary);
  // NB. U is converted to Ueo in fopr_eo->set_config
  //// m_fopr_w_eo->set_config(m_Ueo);
  // m_fopr_w_eo->set_config(m_U);

  //- setting 5D overlap fermion operator
  /*
  m_fopr_ov5d = new Fopr_Overlap_5d(m_fopr_w_eo);
  m_fopr_ov5d->set_parameters(m_mq, m_M0, m_Np, m_x_min, m_x_max,
                              m_Niter_ms, m_Stop_cond_ms, m_boundary);
  m_fopr_ov5d->set_lowmodes(m_Nsbt, m_ev, m_vk);
  */

  Parameters params_ov;
  params_ov.set_string("kernel_type", m_kernel_type);
  params_ov.set_string("gamma_matrix_type", m_repr);
  params_ov.set_double("quark_mass", m_mq);
  params_ov.set_double("domain_wall_height", m_M0);
  params_ov.set_int("number_of_poles", m_Np);
  params_ov.set_double("lower_bound", m_x_min);
  params_ov.set_double("upper_bound", m_x_max);
  params_ov.set_int("maximum_number_of_iteration", m_Niter_ms);
  params_ov.set_double("convergence_criterion_squared", m_Stop_cond_ms);
  params_ov.set_int_vector("boundary_condition", m_boundary);
  params_ov.set_string("verbose_level", vout.get_verbose_level(m_vl));

  //  m_fopr_ov5d = new Fopr_Overlap_5d(m_fopr_w_eo, params_ov);
  m_fopr_ov5d = new Fopr_Overlap_5d(params_ov);
  m_fopr_ov5d->set_config(m_U);
  m_fopr_ov5d->set_lowmodes(m_Nsbt, m_ev, m_vk);


  // ####  object setup  ####
  const string str_solver_type = m_params_solver.get_string("solver_type");
  m_solver = Solver::New(str_solver_type, m_fopr_ov5d);
  m_solver->set_parameters(m_params_solver);

  m_is_initialized = true;
}


//====================================================================
void Fprop_Overlap_5d::finalize()
{
  if (m_is_initialized) {
    //- tidy up
    delete m_solver;
    delete m_fopr_ov5d;
    // delete m_fopr_w_eo;
    delete m_Ueo;
  }
}


//====================================================================
void Fprop_Overlap_5d::invert_D(Field& xq, const Field& b)
{
  if (!m_is_initialized) initialize();

  const double snorm = 1.0 / b.norm2();

  int    Nconv = 0, Nconv1 = 0, Nconv2 = 0;
  double diff = 1.0, diff1 = 1.0, diff2 = 1.0;

  Field v1((Field)b);

  m_fopr_ov->mult_gm5(v1, b);

  solve_overlap_5D_1(xq, v1, Nconv, diff);

  m_fopr_ov->set_mode("D");
  m_fopr_ov->D(v1, xq);
  
  axpy(v1, -1.0, b);
  double vv = dot(v1, v1);
  vv = vv * snorm;

  vout.general(m_vl, "5D overlap solver total:\n");
  vout.general(m_vl, "  Nconv_eo(ov) = %10d\n", Nconv);
  vout.general(m_vl, "  Diff_eo(ov)  = %16.8e\n", diff);
  vout.general(m_vl, "  Diff(ov)     = %16.8e\n", vv);
}


//====================================================================
void Fprop_Overlap_5d::calc_H2inv(Field& xq, const Field& b)
{
  if (!m_is_initialized) initialize();

  const double snorm = 1.0 / b.norm2();

  int    Nconv = 0, Nconv1 = 0, Nconv2 = 0;
  double diff = 1.0, diff1 = 1.0, diff2 = 1.0;

  Field v1(b);

  solve_overlap_5D_1(v1, b, Nconv1, diff1);
  solve_overlap_5D_1(xq, v1, Nconv2, diff2);

  Nconv = Nconv1 + Nconv2;
  diff  = diff1 + diff2;

  m_fopr_ov->set_mode("DdagD");
  m_fopr_ov->mult(v1, xq);

  axpy(v1, -1.0, b);
  double vv = dot(v1, v1);
  vv = vv * snorm;

  vout.general(m_vl, "  5D overlap solver for Q_ov total:\n");
  vout.general(m_vl, "  5D overlap solver total:\n");
  vout.general(m_vl, "  Nconv_eo(ov) = %10d\n", Nconv);
  vout.general(m_vl, "  Diff_eo(ov)  = %16.8e\n", diff);
  vout.general(m_vl, "  Diff(ov)     = %16.8e\n", vv);
}


//====================================================================
void Fprop_Overlap_5d::solve_overlap_5D_1(Field& xq, const Field& b,
                                          int& Nconv, double& diff)
{
  const double snorm = b.norm2();

  // ####  parameter setup  ####
  const int Npl = 2 * m_Np + 1;

  const int Nin   = b.nin();
  const int Nvol  = b.nvol();
  const int Nvol2 = Nvol / 2;

  Field v(Nin, Nvol2, Npl);
  Field s(Nin, Nvol2, Npl);
  Field t(Nin, Nvol2, Npl);
  Field p(Nin, Nvol2, Npl);
  Field x(Nin, Nvol2, Npl);

  Field be(Nin, Nvol2, 1);
  Field bo(Nin, Nvol2, 1);

#pragma omp parallel
 {
  m_index_eo.convertField(be, b, 0);
  m_index_eo.convertField(bo, b, 1);

  v.set(0.0);
  copy(v, 2 * m_Np, be, 0);
  m_fopr_ov5d->LUprecond(t, v, 0);

  v.set(0.0);
  copy(v, 2 * m_Np, bo, 0);
  m_fopr_ov5d->LUprecond(p, v, 1);
  m_fopr_ov5d->Mopr_5d_eo(v, p, 0);
  m_fopr_ov5d->LUprecond(p, v, 0);

  axpy(t, -1.0, p);

  m_fopr_ov5d->DD_5d_eo(s, t, -1);

  m_solver->solve(x, s, Nconv, diff);
  // }

  //- Check
  v.set(0.0);
  copy(v, 2 * m_Np, be, 0);
  m_fopr_ov5d->LUprecond(s, v, 0);

  v.set(0.0);
  copy(v, 2 * m_Np, bo, 0);
  m_fopr_ov5d->LUprecond(p, v, 1);
  m_fopr_ov5d->Mopr_5d_eo(t, p, 0);
  m_fopr_ov5d->LUprecond(p, t, 0);

  axpy(s, -1.0, p);

  m_fopr_ov5d->DD_5d_eo(p, x, 1);
  axpy(p, -1.0, s);

  double diff1 = dot(p, p);

  int ith = ThreadManager::get_thread_id();
  if(ith == 0) diff = diff1 * snorm;
 }

  //- Solution
  Field xqe(Nin, Nvol2, 1);
  Field xqo(Nin, Nvol2, 1);
  Field xt(Nin, Nvol2, 1);

#pragma omp parallel
 {
  copy(xqe, 0, x, 2 * m_Np);

  v.set(0.0);
  copy(v, 2 * m_Np, bo, 0);
  m_fopr_ov5d->LUprecond(p, v, 1);
  copy(xqo, 0, p, 2 * m_Np);

  m_fopr_ov5d->Mopr_5d_eo(p, x, 1);
  m_fopr_ov5d->LUprecond(t, p, 1);
  copy(xt, 0, t, 2 * m_Np);
  axpy(xqo, -1.0, xt);
 }
  m_index_eo.reverseField(xq, xqe, 0);
  m_index_eo.reverseField(xq, xqo, 1);

  vout.general(m_vl, "  solver_eo1(ov):\n");
  vout.general(m_vl, "  Nconv_eo1(ov) = %10d\n", Nconv);
  vout.general(m_vl, "  Diff_eo1(ov)  = %16.8e\n", diff);
}


//====================================================================
double Fprop_Overlap_5d::flop_count()
{
  const int    NPE = CommonParameters::NPE();
  const double eps = CommonParameters::epsilon_criterion();

  //- NB1 Nin = 2 * Nc * Nd, Nex = 1  for field_F
  //- NB2 Nvol/2 for eo
  const int Nin  = 2 * CommonParameters::Nc() * 1;
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
//============================================================END=====
