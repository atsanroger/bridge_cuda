/*!
        @file    fopr_CloverTerm_eo_impl.cpp

        @brief

        @author  UEDA, Satoru (maintained by H.Matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "fopr_CloverTerm_eo_impl.h"

#include "Measurements/Gauge/staple_eo.h"
#include "Field/shiftField_eo.h"
#include "Solver/solver_CG.h"
#include "Tools/gammaMatrixSet.h"

#include "Field/field_thread-inc.h"

namespace Org {
  const std::string Fopr_CloverTerm_eo::class_name
    = "Org::Fopr_CloverTerm_eo";

//====================================================================
  void Fopr_CloverTerm_eo::init(const Parameters& params)
  {
    ThreadManager::assert_single_thread(class_name);

    m_vl = CommonParameters::Vlevel();

    vout.general(m_vl, "%s: construction\n", class_name.c_str());
    vout.increase_indent();

    std::string repr;
    if (!params.fetch_string("gamma_matrix_type", repr)) {
      m_repr = repr;
    } else {
      m_repr = "Dirac"; // default
      vout.general(m_vl, "gamma_matrix_type is not given: defalt = %s\n",
                   m_repr.c_str());
    }
    if ((m_repr != "Dirac") && (m_repr != "Chiral")) {
      vout.crucial("Error in %s: irrelevant mult mode = %s\n",
                   class_name.c_str(), m_repr.c_str());
      exit(EXIT_FAILURE);
    }

    setup();

    set_parameters(params);

    vout.decrease_indent();
    vout.general(m_vl, "%s: construction finished.\n",
                 class_name.c_str());
  }


//====================================================================
  void Fopr_CloverTerm_eo::init(const std::string repr)
  {
    ThreadManager::assert_single_thread(class_name);

    m_vl = CommonParameters::Vlevel();

    vout.general(m_vl, "%s: construction (obsolete)\n",
                 class_name.c_str());
    vout.increase_indent();

    m_repr = repr;
    setup();

    vout.decrease_indent();

    vout.general(m_vl, "%s: construction finished.\n",
                 class_name.c_str());
  }


//====================================================================
  void Fopr_CloverTerm_eo::setup()
  {
    m_Nc    = CommonParameters::Nc();
    m_Nd    = CommonParameters::Nd();
    m_NinF  = 2 * m_Nc * m_Nd;
    m_Nvol  = CommonParameters::Nvol();
    m_Nvol2 = m_Nvol / 2;
    m_Ndim  = CommonParameters::Ndim();

    m_Ndm2 = m_Nd * m_Nd / 2;
    m_T.reset(m_Nvol, m_Ndm2);

    m_boundary.resize(m_Ndim);

    m_Ueo.reset(m_Nvol, m_Ndim);

    m_shift_eo = new ShiftField_eo(2 * m_Nc * m_Nc);
    m_staple   = new Staple_eo;

    m_GM.resize(m_Ndim + 1);
    m_SG.resize(m_Ndim * m_Ndim);

    GammaMatrixSet *gmset = GammaMatrixSet::New(m_repr);

    m_GM[0] = gmset->get_GM(gmset->GAMMA1);
    m_GM[1] = gmset->get_GM(gmset->GAMMA2);
    m_GM[2] = gmset->get_GM(gmset->GAMMA3);
    m_GM[3] = gmset->get_GM(gmset->GAMMA4);
    m_GM[4] = gmset->get_GM(gmset->GAMMA5);

    m_SG[sg_index(0, 1)] = gmset->get_GM(gmset->SIGMA12);
    m_SG[sg_index(1, 2)] = gmset->get_GM(gmset->SIGMA23);
    m_SG[sg_index(2, 0)] = gmset->get_GM(gmset->SIGMA31);
    m_SG[sg_index(3, 0)] = gmset->get_GM(gmset->SIGMA41);
    m_SG[sg_index(3, 1)] = gmset->get_GM(gmset->SIGMA42);
    m_SG[sg_index(3, 2)] = gmset->get_GM(gmset->SIGMA43);

    m_SG[sg_index(1, 0)] = m_SG[sg_index(0, 1)].mult(-1);
    m_SG[sg_index(2, 1)] = m_SG[sg_index(1, 2)].mult(-1);
    m_SG[sg_index(0, 2)] = m_SG[sg_index(2, 0)].mult(-1);
    m_SG[sg_index(0, 3)] = m_SG[sg_index(3, 0)].mult(-1);
    m_SG[sg_index(1, 3)] = m_SG[sg_index(3, 1)].mult(-1);
    m_SG[sg_index(2, 3)] = m_SG[sg_index(3, 2)].mult(-1);

    m_SG[sg_index(0, 0)] = gmset->get_GM(gmset->UNITY);
    m_SG[sg_index(1, 1)] = gmset->get_GM(gmset->UNITY);
    m_SG[sg_index(2, 2)] = gmset->get_GM(gmset->UNITY);
    m_SG[sg_index(3, 3)] = gmset->get_GM(gmset->UNITY);
    // these 4 gamma matrices are actually not used.

    delete gmset;

    m_Fee_inv.reset(m_Nvol2, m_Nc * m_Nd);
    m_Foo_inv.reset(m_Nvol2, m_Nc * m_Nd);

    const int Nfst = 6;
    m_T2.resize(2); // 0: even, 1: odd.
    m_T2[0].reset(m_Nvol2, Nfst);
    m_T2[1].reset(m_Nvol2, Nfst);

    m_v1.reset(m_Nvol2, 1);
    m_v2.reset(m_Nvol2, 1);
    m_v3.reset(m_Nvol2, 1);

    m_w1.reset(m_Nvol2, 1);
    m_w2.reset(m_Nvol2, 1);

    // setup solver
    Parameters params_solver;
    params_solver.set_string("solver_type", "CG");
    params_solver.set_int("maximum_number_of_iteration", 100);
    params_solver.set_int("maximum_number_of_restart", 40);
    params_solver.set_double("convergence_criterion_squared", 1.0e-30);
    params_solver.set_string("use_initial_guess", "false");
    params_solver.set_string("verbose_level", "Crucial");
    // NB. set VerboseLevel to CRUCIAL to suppress frequent messages.

    m_solver = new Solver_CG(this);
    m_solver->set_parameters(params_solver);
  }


//====================================================================
  void Fopr_CloverTerm_eo::tidyup()
  {
    delete m_solver;
    delete m_shift_eo;
    delete m_staple;
  }


//====================================================================
  void Fopr_CloverTerm_eo::set_parameters(const Parameters& params)
  {
    int ith = ThreadManager::get_thread_id();

    std::string vlevel;
    if (!params.fetch_string("verbose_level", vlevel)) {
      if (ith == 0) m_vl = vout.set_verbose_level(vlevel);
    }

    //- fetch and check input parameters
    double           kappa, cSW;
    std::vector<int> bc;

    int err = 0;
    err += params.fetch_double("hopping_parameter", kappa);
    err += params.fetch_double("clover_coefficient", cSW);
    err += params.fetch_int_vector("boundary_condition", bc);

    if (err) {
      vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }

    set_parameters(kappa, cSW, bc);
  }


//====================================================================
  void Fopr_CloverTerm_eo::set_parameters(const double kappa,
                                          const double cSW,
                                          const std::vector<int> bc)
  {
    int ith = ThreadManager::get_thread_id();

    assert(bc.size() == m_Ndim);

    if (ith == 0) {
      m_kappa    = kappa;
      m_cSW      = cSW;
      m_boundary = bc;
    }

    vout.general(m_vl, "%s: input parameters\n", class_name.c_str());
    vout.general(m_vl, "  gamma matrix type = %s\n", m_repr.c_str());
    vout.general(m_vl, "  kappa = %12.8f\n", m_kappa);
    vout.general(m_vl, "  cSW   = %12.8f\n", m_cSW);
    for (int mu = 0; mu < m_Ndim; ++mu) {
      vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
    }
  }


//====================================================================
  void Fopr_CloverTerm_eo::get_parameters(Parameters& params) const
  {
    params.set_double("hopping_parameter", m_kappa);
    params.set_double("clover_coefficient", m_cSW);
    params.set_int_vector("boundary_condition", m_boundary);
    params.set_string("gamma_matrix_type", m_repr);

    params.set_string("verbose_level", vout.get_verbose_level(m_vl));
  }


//====================================================================
  void Fopr_CloverTerm_eo::set_config(Field *U)
  {
    int nth = ThreadManager::get_num_threads();
    vout.detailed(m_vl, "%s: set_config is called: num_threads = %d\n",
                  class_name.c_str(), nth);

    if (nth > 1) {
      set_config_impl(U);
    } else {
      set_config_omp(U);
    }

    vout.detailed(m_vl, "%s: set_config finished\n", class_name.c_str());
  }


//====================================================================
  void Fopr_CloverTerm_eo::set_config_omp(Field *U)
  {
    vout.detailed(m_vl, "  set_config_omp is called.\n");

#pragma omp parallel
    {
      set_config_impl(U);
    }
  }


//====================================================================
  void Fopr_CloverTerm_eo::set_config_impl(Field *U)
  {
    m_index_eo.convertField(m_Ueo, *U);

    set_csw();
    solve_csw_inv();
  }


//====================================================================
  void Fopr_CloverTerm_eo::set_mode(const std::string mode)
  {
#pragma omp barrier

    int ith = ThreadManager::get_thread_id();
    if (ith == 0) m_mode = mode;

#pragma omp barrier
  }


//====================================================================
  void Fopr_CloverTerm_eo::mult(Field& v, const Field& w)
  {
    if (m_mode == "even") {
      D(v, w, 0);
    } else if (m_mode == "odd") {
      D(v, w, 1);
    } else {
      vout.crucial(m_vl, "Error at %s: undefined mode = %s\n",
                   class_name.c_str(), m_mode.c_str());
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_CloverTerm_eo::mult_dag(Field& v, const Field& w)
  {
    if (m_mode == "even") {
      D(v, w, 0);
    } else if (m_mode == "odd") {
      D(v, w, 1);
    } else {
      vout.crucial(m_vl, "Error at %s: undefined mode = %s\n",
                   class_name.c_str(), m_mode.c_str());
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_CloverTerm_eo::D(Field& v, const Field& f, const int ieo)
  {
    // multiplies [ 1 - csw kappa sigma_{mu nu} F_{mu nu} ]

#pragma omp barrier

    copy(m_v1, f); //  m_v1 = (Field_F)f;
    copy(m_v3, f);

    const double coeff = -m_kappa * m_cSW;

    // i s_23 F_23
    mult_iGM(m_v2, m_SG[sg_index(1, 2)], m_v1);
    multadd_Field_Gn(m_v3, 0, m_T2[ieo], 0, m_v2, 0, coeff);

    // i s_31 F_31
    mult_iGM(m_v2, m_SG[sg_index(2, 0)], m_v1);
    multadd_Field_Gn(m_v3, 0, m_T2[ieo], 1, m_v2, 0, coeff);

    // i s_12 F_12
    mult_iGM(m_v2, m_SG[sg_index(0, 1)], m_v1);
    multadd_Field_Gn(m_v3, 0, m_T2[ieo], 2, m_v2, 0, coeff);

    // i s_41 F_41
    mult_iGM(m_v2, m_SG[sg_index(3, 0)], m_v1);
    multadd_Field_Gn(m_v3, 0, m_T2[ieo], 3, m_v2, 0, coeff);

    // i s_42 F_42
    mult_iGM(m_v2, m_SG[sg_index(3, 1)], m_v1);
    multadd_Field_Gn(m_v3, 0, m_T2[ieo], 4, m_v2, 0, coeff);

    // i s_43 F_43
    mult_iGM(m_v2, m_SG[sg_index(3, 2)], m_v1);
    multadd_Field_Gn(m_v3, 0, m_T2[ieo], 5, m_v2, 0, coeff);

    copy(v, m_v3);

#pragma omp barrier
  }


//====================================================================
  void Fopr_CloverTerm_eo::mult_csw_inv(Field& v,
                                        const Field& f, const int ieo)
  {
    // multiplies [ 1 - csw kappa sigma_{mu nu} F_{mu nu} ]^{-1}

#pragma omp barrier

    copy(m_v1, f);
    m_v2.set(0.0);

    Field_F *csw_inv;
    if (ieo == 0) {
      csw_inv = &m_Fee_inv;
    } else if (ieo == 1) {
      csw_inv = &m_Foo_inv;
    } else {
      vout.crucial(m_vl, "Error at %s: wrong parameter, ieo=%d.\n",
                   class_name.c_str(), ieo);
      exit(EXIT_FAILURE);
    }

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol2);

    for (int isite = is; isite < ns; ++isite) {
      for (int id = 0; id < m_Nd; ++id) {
        for (int ic = 0; ic < m_Nc; ++ic) {
          double re = 0.0;
          double im = 0.0;
          for (int jd = 0; jd < m_Nd; ++jd) {
            for (int jc = 0; jc < m_Nc; ++jc) {
              int icd = jc + m_Nc * jd;
              // Hermiticity of clover term is used here.
              re += csw_inv->cmp_r(ic, id, isite, icd)
                    * m_v1.cmp_r(jc, jd, isite, 0);
              re += csw_inv->cmp_i(ic, id, isite, icd)
                    * m_v1.cmp_i(jc, jd, isite, 0);

              im += csw_inv->cmp_r(ic, id, isite, icd)
                    * m_v1.cmp_i(jc, jd, isite, 0);
              im -= csw_inv->cmp_i(ic, id, isite, icd)
                    * m_v1.cmp_r(jc, jd, isite, 0);
            }
          }

          m_v2.set_ri(ic, id, isite, 0, re, im);
        }
      }
    }

    copy(v, m_v2);

#pragma omp barrier
  }


//====================================================================
  void Fopr_CloverTerm_eo::mult_isigma(Field_F& v, const Field_F& w,
                                       const int mu, const int nu)
  {
#pragma omp barrier

    assert(mu != nu);
    mult_iGM(v, m_SG[sg_index(mu, nu)], w);

#pragma omp barrier
  }


//====================================================================
  void Fopr_CloverTerm_eo::set_csw()
  {
#pragma omp barrier

    m_T.set(0.0);

    set_fieldstrength(m_Ft, 1, 2); // F_23
    copy(m_T, 0, m_Ft, 0);

    set_fieldstrength(m_Ft, 2, 0); // F_31
    copy(m_T, 1, m_Ft, 0);

    set_fieldstrength(m_Ft, 0, 1); // F_12
    copy(m_T, 2, m_Ft, 0);

    set_fieldstrength(m_Ft, 3, 0); // F_41
    copy(m_T, 3, m_Ft, 0);

    set_fieldstrength(m_Ft, 3, 1); // F_42
    copy(m_T, 4, m_Ft, 0);

    set_fieldstrength(m_Ft, 3, 2); // F_43
    copy(m_T, 5, m_Ft, 0);

#pragma omp barrier

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol2);

    const int Nfst = 6;
    const int NinG = m_T2[0].nin();

    for (int ieo = 0; ieo < 2; ++ieo) {
      for (int ex = 0; ex < Nfst; ++ex) {
        for (int isite = is; isite < ns; ++isite) {
          for (int in = 0; in < NinG; ++in) {
            double cmp = m_T.cmp(in, m_index_eo.site(isite, ieo), ex);
            m_T2[ieo].set(in, isite, ex, cmp);
          }
        }
      }
    }

#pragma omp barrier
  }


//====================================================================
  void Fopr_CloverTerm_eo::solve_csw_inv()
  {
    const double eps2 = CommonParameters::epsilon_criterion2();

#pragma omp barrier

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol2);

    for (int id = 0; id < m_Nd; ++id) {
      for (int ic = 0; ic < m_Nc; ++ic) {
        int icd = ic + m_Nc * id;

        m_w1.set(0.0);
        for (int isite = is; isite < ns; ++isite) {
          m_w1.set_ri(ic, id, isite, 0, 1, 0);
        }
#pragma omp barrier

        if (m_cSW * m_cSW < eps2) {
          copy(m_Fee_inv, icd, m_w1, 0);
          copy(m_Foo_inv, icd, m_w1, 0);
        } else {
          int    Nconv;
          double diff;
          set_mode("even");
          m_solver->solve(m_w2, m_w1, Nconv, diff);
          copy(m_Fee_inv, icd, m_w2, 0);

          vout.detailed(m_vl, "  Nconv,diff = %d %12.4e\n", Nconv, diff);

          set_mode("odd");
          m_solver->solve(m_w2, m_w1, Nconv, diff);
          copy(m_Foo_inv, icd, m_w2, 0);

          vout.detailed(m_vl, "  Nconv,diff = %d %12.4e\n", Nconv, diff);
        }
      }
    }
#pragma omp barrier

    // redefine the inverse matrix with its dagger.
    for (int ics = 0; ics < m_Nc * m_Nd; ++ics) {
      for (int site = is; site < ns; ++site) {
        for (int id = 0; id < m_Nd; ++id) {
          for (int ic = 0; ic < m_Nc; ++ic) {
            double re = m_Foo_inv.cmp_r(ic, id, site, ics);
            double im = m_Foo_inv.cmp_i(ic, id, site, ics);

            m_Foo_inv.set_ri(ic, id, site, ics, re, -im);

            re = m_Fee_inv.cmp_r(ic, id, site, ics);
            im = m_Fee_inv.cmp_i(ic, id, site, ics);

            m_Fee_inv.set_ri(ic, id, site, ics, re, -im);
          }
        }
      }
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_CloverTerm_eo::set_fieldstrength(Field_G& Fst,
                                             const int mu, const int nu)
  {
    copy(m_Umu, 0, m_Ueo, mu);
    m_staple->upper(m_Cup, m_Ueo, mu, nu);
    m_staple->lower(m_Cdn, m_Ueo, mu, nu);

    mult_Field_Gnd(Fst, 0, m_Umu, 0, m_Cup, 0);
    multadd_Field_Gnd(Fst, 0, m_Umu, 0, m_Cdn, 0, -1.0);

    mult_Field_Gdn(m_u1, 0, m_Cup, 0, m_Umu, 0);
    multadd_Field_Gdn(m_u1, 0, m_Cdn, 0, m_Umu, 0, -1.0);

    m_shift_eo->forward(m_u2, m_u1, mu);
    axpy(Fst, 1.0, m_u2);

    ah_Field_G(Fst, 0);
    scal(Fst, 0.25);
  }


//====================================================================
  double Fopr_CloverTerm_eo::flop_count()
  {
    // Counting of floating point operations in giga unit.
    // The following counting explicitly depends on the implementation
    // and to be recalculated when the code is modified.
    // Present counting is based on rev.1107. [24 Aug 2014 H.Matsufuru]

    const int Nvol = CommonParameters::Nvol();
    const int NPE  = CommonParameters::NPE();

    const int flop_site = 8 * m_Nc * m_Nc * m_Nd * m_Nd;

    const double gflop = flop_site * ((Nvol / 2) * (NPE / 1.0e+9));

    return gflop;
  }


//====================================================================
}
//============================================================END=====
