/*!
        @file    fopr_CloverTerm_eo_impl.cpp

        @brief

        @author  UEDA, Satoru (maintained by H.Matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "fopr_CloverTerm_eo_impl.h"

#include "Field/shiftField_eo.h"
#include "Measurements/Gauge/staple_eo.h"
#include "Solver/solver_CG.h"
#include "Tools/gammaMatrixSet.h"

#include "Field/field_thread-inc.h"
// because of no communication in mult, threading for Field is applied.

#if defined USE_GROUP_SU3
#include "fopr_Wilson_impl_SU3-inc.h"
#elif defined USE_GROUP_SU2
#include "fopr_Wilson_impl_SU2-inc.h"
#elif defined USE_GROUP_SU_N
#include "fopr_Wilson_impl_SU_N-inc.h"
#endif

namespace Imp {
  const std::string Fopr_CloverTerm_eo::class_name
    = "Imp::Fopr_CloverTerm_eo";

//====================================================================
  void Fopr_CloverTerm_eo::init(const Parameters& params)
  {
    ThreadManager::assert_single_thread(class_name);

    m_vl = CommonParameters::Vlevel();

    vout.general(m_vl, "%s: construction\n", class_name.c_str());
    vout.increase_indent();

    setup();

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

    setup_gamma_matrices();

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

    setup();

    m_repr = repr;

    setup_gamma_matrices();

    vout.decrease_indent();

    vout.general(m_vl, "%s: construction finished.\n",
                 class_name.c_str());
  }


//====================================================================
  void Fopr_CloverTerm_eo::setup()
  {
    m_Nc    = CommonParameters::Nc();
    m_Nd    = CommonParameters::Nd();
    m_Ndim  = CommonParameters::Ndim();
    m_NinF  = 2 * m_Nc * m_Nd;
    m_Nvol  = CommonParameters::Nvol();
    m_Nvol2 = m_Nvol / 2;

    m_boundary.resize(m_Ndim);

    m_Ueo.reset(m_Nvol, m_Ndim);

    m_Ndm2 = m_Nd * m_Nd / 2;
    m_T.reset(m_Nvol, m_Ndm2);

    m_GM.resize(m_Ndim + 1);
    m_SG.resize(m_Ndim * m_Ndim);

    m_Fee_inv.reset(m_Nvol2, m_Nc * m_Nd);
    m_Foo_inv.reset(m_Nvol2, m_Nc * m_Nd);

    int Ndf = 2 * m_Nc * m_Nc;
    m_shift_eo = new ShiftField_eo(Ndf);

    m_staple = new Staple_eo;

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
  void Fopr_CloverTerm_eo::setup_gamma_matrices()
  {
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
    m_idx.convertField(m_Ueo, *U);

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
    // multiplies 1-csw kappa sigma_{mu nu} F_{mu nu}
    if (m_mode == "even") {
      D(v, w, 0);
    } else if (m_mode == "odd") {
      D(v, w, 1);
    } else {
      vout.crucial("Error at %s: undefined mode = %s\n",
                   class_name.c_str(), m_mode.c_str());
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_CloverTerm_eo::mult_dag(Field& v, const Field& w)
  {
    // multiplies 1-csw kappa sigma_{mu nu} F_{mu nu}
    mult(v, w);
  }


//====================================================================
  void Fopr_CloverTerm_eo::solve_csw_inv()
  {
    const double eps2 = CommonParameters::epsilon_criterion2();

    for (int ispin = 0; ispin < m_Nd; ++ispin) {
      for (int icolor = 0; icolor < m_Nc; ++icolor) {
        int spin_color = icolor + m_Nc * ispin;

        m_w1.set(0.0);
#pragma omp barrier

        for (int isite = 0; isite < m_Nvol2; ++isite) {
          m_w1.set_ri(icolor, ispin, isite, 0, 1, 0);
        }
#pragma omp barrier

        if (m_cSW * m_cSW < eps2) {
          copy(m_Fee_inv, spin_color, m_w1, 0);
          copy(m_Foo_inv, spin_color, m_w1, 0);
        } else {
          int    Nconv;
          double diff;
          set_mode("even");
#pragma omp barrier
          m_solver->solve(m_w2, m_w1, Nconv, diff);
          copy(m_Fee_inv, spin_color, m_w2, 0);
#pragma omp barrier
          vout.detailed(m_vl, "  Nconv,diff = %d %12.4e\n", Nconv, diff);

          set_mode("odd");
#pragma omp barrier
          m_solver->solve(m_w2, m_w1, Nconv, diff);
          copy(m_Foo_inv, spin_color, m_w2, 0);
#pragma omp barrier
          vout.detailed(m_vl, "  Nconv,diff = %d %12.4e\n", Nconv, diff);
        }
      }
    }

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol2);

    // redefine the inverse matrix with its dagger.
    for (int ics = 0; ics < m_Nc * m_Nd; ++ics) {
      for (int site = is; site < ns; ++site) {
        for (int id = 0; id < m_Nd; ++id) {
          for (int ic = 0; ic < m_Nc; ++ic) {
            double re, im;
            re = m_Foo_inv.cmp_r(ic, id, site, ics);
            im = m_Foo_inv.cmp_i(ic, id, site, ics);
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
  void Fopr_CloverTerm_eo::mult_csw_inv(Field& v,
                                        const Field& w, const int ieo)
  {
    // multiplies [ 1 - csw kappa sigma_{mu nu} F_{mu nu} ]^{-1}
    if (m_repr == "Dirac") {
      mult_csw_inv_dirac(v, w, ieo);
    } else if (m_repr == "Chiral") {
      mult_csw_inv_chiral(v, w, ieo);
    }
  }


//====================================================================
  void Fopr_CloverTerm_eo::mult_csw_inv_dirac(Field& v, const Field& w,
                                              const int ieo)
  {
    // multiplies [ 1 - csw kappa sigma_{mu nu} F_{mu nu} ]^{-1}
    const int Nvc = 2 * m_Nc;

    const double *v1 = w.ptr(0);
    double       *v2 = v.ptr(0);

    double *csw_inv = 0;
    if (ieo == 0) {
      csw_inv = m_Fee_inv.ptr(0);
    } else if (ieo == 1) {
      csw_inv = m_Foo_inv.ptr(0);
    }

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol2);

#pragma omp barrier

    const int Nd2 = m_Nd / 2;
    for (int site = is; site < ns; ++site) {
      for (int icd = 0; icd < m_Nc * Nd2; ++icd) {
        int iv2 = 2 * icd + m_NinF * site;
        v2[iv2]     = 0.0;
        v2[iv2 + 1] = 0.0;
        for (int jd = 0; jd < m_Nd; ++jd) {
          int jcd = Nvc * jd;
          int iv  = jcd + m_NinF * site;

          int ig = jcd + m_NinF * (site + m_Nvol2 * icd);
          v2[iv2]     += mult_uv_r(&csw_inv[ig], &v1[iv], m_Nc);
          v2[iv2 + 1] += mult_uv_i(&csw_inv[ig], &v1[iv], m_Nc);
        }
      }

      for (int icd = 0; icd < m_Nc * Nd2; ++icd) {
        int iv2 = 2 * (icd + m_Nc * Nd2) + m_NinF * site;
        v2[iv2]     = 0.0;
        v2[iv2 + 1] = 0.0;
        for (int jd = 0; jd < m_Nd; ++jd) {
          int jd2 = (jd + Nd2) % m_Nd;
          int iv  = Nvc * jd + m_NinF * site;
          int ig  = Nvc * jd2 + m_NinF * (site + m_Nvol2 * icd);

          v2[iv2]     += mult_uv_r(&csw_inv[ig], &v1[iv], m_Nc);
          v2[iv2 + 1] += mult_uv_i(&csw_inv[ig], &v1[iv], m_Nc);
        }
      }
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_CloverTerm_eo::mult_csw_inv_chiral(Field& v, const Field& w,
                                               const int ieo)
  {
    // multiplies [ 1 - csw kappa sigma_{mu nu} F_{mu nu} ]^{-1}
    const int Nvc = 2 * m_Nc;

    const double *v1 = w.ptr(0);
    double       *v2 = v.ptr(0);

    double *csw_inv = 0;
    if (ieo == 0) {
      csw_inv = m_Fee_inv.ptr(0);
    } else if (ieo == 1) {
      csw_inv = m_Foo_inv.ptr(0);
    }

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol2);

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      for (int icd = 0; icd < m_Nc * m_Nd / 2; ++icd) {
        int iv2 = 2 * icd + m_NinF * site;
        v2[iv2]     = 0.0;
        v2[iv2 + 1] = 0.0;

        for (int jd = 0; jd < m_Nd / 2; ++jd) {
          int jcd = Nvc * jd;
          int iv  = jcd + m_NinF * site;
          int ig  = jcd + m_NinF * (site + m_Nvol2 * icd);

          v2[iv2]     += mult_uv_r(&csw_inv[ig], &v1[iv], m_Nc);
          v2[iv2 + 1] += mult_uv_i(&csw_inv[ig], &v1[iv], m_Nc);
        }
      }

      for (int icd = m_Nc * m_Nd / 2; icd < m_Nc * m_Nd; ++icd) {
        int iv2 = 2 * icd + m_NinF * site;
        v2[iv2]     = 0.0;
        v2[iv2 + 1] = 0.0;

        for (int jd = m_Nd / 2; jd < m_Nd; ++jd) {
          int jcd = Nvc * jd;
          int iv  = jcd + m_NinF * site;
          int ig  = jcd + m_NinF * (site + m_Nvol2 * icd);

          v2[iv2]     += mult_uv_r(&csw_inv[ig], &v1[iv], m_Nc);
          v2[iv2 + 1] += mult_uv_i(&csw_inv[ig], &v1[iv], m_Nc);
        }
      }
    }

#pragma omp barrier
  }


//====================================================================
  void Fopr_CloverTerm_eo::D(Field& v, const Field& w, const int ieo)
  {
    // multiplies [ 1 - csw kappa sigma_{mu nu} F_{mu nu} ]

    assert(w.check_size(m_NinF, m_Nvol2, 1));
    assert(v.check_size(m_NinF, m_Nvol2, 1));

    if (m_repr == "Dirac") {
      D_dirac(v, w, ieo);
    } else if (m_repr == "Chiral") {
      D_chiral(v, w, ieo);
    }
  }


//====================================================================
  void Fopr_CloverTerm_eo::D_dirac(Field& v, const Field& w, const int ieo)
  {
    // multiplies [ 1 - csw kappa sigma_{mu nu} F_{mu nu} ]

    const double *wp = w.ptr(0);
    double       *vp = v.ptr(0);
    double       *tp = m_T.ptr(0, m_Nvol2 * ieo, 0);

    const int Nvc  = 2 * m_Nc;
    const int Nd2  = m_Nd / 2;
    const int NinF = 2 * m_Nc * m_Nd;
    const int NinG = 2 * m_Nc * m_Nc;

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol2);

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      for (int id = 0; id < Nd2; ++id) {
        for (int ic = 0; ic < m_Nc; ++ic) {
          int icd = ic + m_Nc * id;

          int iv2 = 2 * icd + NinF * site;
          vp[iv2]     = 0.0;
          vp[iv2 + 1] = 0.0;
          for (int jd = 0; jd < m_Nd; ++jd) {
            int iv = Nvc * jd + NinF * site;
            int ig = Nvc * ic + NinG * (site + m_Nvol * (id * m_Nd + jd));

            vp[iv2]     += mult_uv_r(&tp[ig], &wp[iv], m_Nc);
            vp[iv2 + 1] += mult_uv_i(&tp[ig], &wp[iv], m_Nc);
          }

          iv2        += Nvc * Nd2;
          vp[iv2]     = 0.0;
          vp[iv2 + 1] = 0.0;
          for (int jd = 0; jd < m_Nd; ++jd) {
            int jd2 = (2 + jd) % m_Nd;
            int iv  = Nvc * jd + NinF * site;
            int ig  = Nvc * ic + NinG * (site + m_Nvol * (id * m_Nd + jd2));

            vp[iv2]     += mult_uv_r(&tp[ig], &wp[iv], m_Nc);
            vp[iv2 + 1] += mult_uv_i(&tp[ig], &wp[iv], m_Nc);
          }
        }
      }
    }

#pragma omp barrier
  }


//====================================================================
  void Fopr_CloverTerm_eo::D_chiral(Field& v,
                                    const Field& w, const int ieo)
  {
    // multiplies [ 1 - csw kappa sigma_{mu nu} F_{mu nu} ]

    const double *wp = w.ptr(0);
    double       *vp = v.ptr(0);
    double       *tp = m_T.ptr(0, m_Nvol2 * ieo, 0);

    const int Nvc  = 2 * m_Nc;
    const int Nd2  = m_Nd / 2;
    const int NinF = 2 * m_Nc * m_Nd;
    const int NinG = 2 * m_Nc * m_Nc;

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol2);

#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      for (int id = 0; id < Nd2; ++id) {
        for (int ic = 0; ic < m_Nc; ++ic) {
          int icd = ic + m_Nc * id;

          double *vp2 = &vp[NinF * site];
          vp2[2 * icd]     = 0.0;
          vp2[2 * icd + 1] = 0.0;
          for (int jd = 0; jd < Nd2; ++jd) {
            int iv = Nvc * jd + NinF * site;
            int ig = Nvc * ic + NinG * (site + m_Nvol * (id * Nd2 + jd));
            vp2[2 * icd]     += mult_uv_r(&tp[ig], &wp[iv], m_Nc);
            vp2[2 * icd + 1] += mult_uv_i(&tp[ig], &wp[iv], m_Nc);
          }

          vp2              = &vp[Nvc * Nd2 + NinF * site];
          vp2[2 * icd]     = 0.0;
          vp2[2 * icd + 1] = 0.0;
          for (int jd = 0; jd < Nd2; ++jd) {
            int iv = Nvc * (Nd2 + jd) + NinF * site;
            int ig = Nvc * ic + NinG * (site + m_Nvol * (m_Nd + id * Nd2 + jd));
            vp2[2 * icd]     += mult_uv_r(&tp[ig], &wp[iv], m_Nc);
            vp2[2 * icd + 1] += mult_uv_i(&tp[ig], &wp[iv], m_Nc);
          }
        }
      }
    }


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
    if (m_repr == "Dirac") {
      set_csw_dirac();
    } else if (m_repr == "Chiral") {
      set_csw_chiral();
    } else {
      vout.crucial(m_vl, "Error at %s: irrelevant gamma_matrix_type: %s\n",
                   class_name.c_str(), m_repr.c_str());
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_CloverTerm_eo::set_csw_dirac()
  {
    // The clover term in the Dirac representation is as spin-space
    // matrix
    //   [ P Q ]
    //   [ Q P ],
    // where P and Q are 2x2 block matrices as
    //   P =  [          iF(1,2)   F(3,1) + iF(2,3) ]
    //        [-F(3,1) + iF(2,3)          - iF(1,2) ]
    // and
    //   Q =  [        - iF(4,3)  -F(4,2) - iF(4,1) ]
    //        [ F(4,2) - iF(4,1)            iF(4,3) ]
    // up to the coefficient.
    // in the following what defined is
    // [ P Q ] = [ T(0) T(1)  T(2) T(3) ]
    //           [ T(4) T(5)  T(6) T(7) ].

#pragma omp barrier

    m_T.set(0.0);
#pragma omp barrier

    //- sigma23
    //  Field_G F;
    set_fieldstrength(m_Ft, 1, 2);
    m_Ft.xI();
#pragma omp barrier
    axpy(m_T, 1, 1.0, m_Ft, 0);
    axpy(m_T, 4, 1.0, m_Ft, 0);

    //- sigma31
    set_fieldstrength(m_Ft, 2, 0);
    axpy(m_T, 1, 1.0, m_Ft, 0);
    axpy(m_T, 4, -1.0, m_Ft, 0);

    //- sigma12
    set_fieldstrength(m_Ft, 0, 1);
    m_Ft.xI();
#pragma omp barrier
    axpy(m_T, 0, 1.0, m_Ft, 0);
    axpy(m_T, 5, -1.0, m_Ft, 0);

    //- sigma41
    set_fieldstrength(m_Ft, 3, 0);
    m_Ft.xI();
#pragma omp barrier
    axpy(m_T, 3, -1.0, m_Ft, 0);
    axpy(m_T, 6, -1.0, m_Ft, 0);

    //- sigma42
    set_fieldstrength(m_Ft, 3, 1);
    axpy(m_T, 3, -1.0, m_Ft, 0);
    axpy(m_T, 6, 1.0, m_Ft, 0);

    //- sigma43
    set_fieldstrength(m_Ft, 3, 2);
    m_Ft.xI();
#pragma omp barrier
    axpy(m_T, 2, -1.0, m_Ft, 0);
    axpy(m_T, 7, 1.0, m_Ft, 0);

    scal(m_T, -m_kappa * m_cSW);

    //    Field_G Unity;
    m_Ft.set_unit();
#pragma omp barrier
    axpy(m_T, 0, 1.0, m_Ft, 0);
    axpy(m_T, 5, 1.0, m_Ft, 0);

#pragma omp barrier
  }


//====================================================================
  void Fopr_CloverTerm_eo::set_csw_chiral()
  {
    // The clover term in the Dirac representation is
    // as spin-space matrix
    //  [ P+Q  0  ]
    //  [ 0   P-Q ],
    // where P and Q are 2x2 block matrices as
    //        [          iF(1,2) |  F(3,1) + iF(2,3) ]
    //   P =  [ -----------------+------------------ ]
    //        [-F(3,1) + iF(2,3) |         - iF(1,2) ]
    // and
    //        [        - iF(4,3) | -F(4,2) - iF(4,1) ]
    //   Q =  [ -----------------+------------------ ]
    //        [ F(4,2) - iF(4,1) |           iF(4,3) ]
    // up to the coefficient.
    // in the following what defined is
    //        [ T(0) | T(1) ]          [ T(4) | T(5) ]
    //  P+Q = [ -----+----- ]  P - Q = [ -----+----- ]
    //        [ T(2) | T(3) ]          [ T(6) | T(7) ]

#pragma omp barrier

    m_T.set(0.0);

#pragma omp barrier

    //- sigma23
    set_fieldstrength(m_Ft, 1, 2);
    m_Ft.xI();
#pragma omp barrier
    axpy(m_T, 1, 1.0, m_Ft, 0);
    axpy(m_T, 2, 1.0, m_Ft, 0);
    axpy(m_T, 5, 1.0, m_Ft, 0);
    axpy(m_T, 6, 1.0, m_Ft, 0);

    //- sigma31
    set_fieldstrength(m_Ft, 2, 0);
    axpy(m_T, 1, 1.0, m_Ft, 0);
    axpy(m_T, 2, -1.0, m_Ft, 0);
    axpy(m_T, 5, 1.0, m_Ft, 0);
    axpy(m_T, 6, -1.0, m_Ft, 0);

    //- sigma12
    set_fieldstrength(m_Ft, 0, 1);
    m_Ft.xI();
#pragma omp barrier
    axpy(m_T, 0, 1.0, m_Ft, 0);
    axpy(m_T, 3, -1.0, m_Ft, 0);
    axpy(m_T, 4, 1.0, m_Ft, 0);
    axpy(m_T, 7, -1.0, m_Ft, 0);

    //- sigma41
    set_fieldstrength(m_Ft, 3, 0);
    m_Ft.xI();
#pragma omp barrier
    axpy(m_T, 1, -1.0, m_Ft, 0);
    axpy(m_T, 2, -1.0, m_Ft, 0);
    axpy(m_T, 5, 1.0, m_Ft, 0);
    axpy(m_T, 6, 1.0, m_Ft, 0);

    //- sigma42
    set_fieldstrength(m_Ft, 3, 1);
    axpy(m_T, 1, -1.0, m_Ft, 0);
    axpy(m_T, 2, 1.0, m_Ft, 0);
    axpy(m_T, 5, 1.0, m_Ft, 0);
    axpy(m_T, 6, -1.0, m_Ft, 0);

    //- sigma43
    set_fieldstrength(m_Ft, 3, 2);
    m_Ft.xI();
#pragma omp barrier
    axpy(m_T, 0, -1.0, m_Ft, 0);
    axpy(m_T, 3, 1.0, m_Ft, 0);
    axpy(m_T, 4, 1.0, m_Ft, 0);
    axpy(m_T, 7, -1.0, m_Ft, 0);
#pragma omp barrier

    scal(m_T, -m_kappa * m_cSW);

#pragma omp barrier

    m_Ft.set_unit();
#pragma omp barrier
    axpy(m_T, 0, 1.0, m_Ft, 0);
    axpy(m_T, 3, 1.0, m_Ft, 0);
    axpy(m_T, 4, 1.0, m_Ft, 0);
    axpy(m_T, 7, 1.0, m_Ft, 0);

#pragma omp barrier
  }


//====================================================================
  void Fopr_CloverTerm_eo::set_fieldstrength(Field_G& Fst,
                                             const int mu, const int nu)
  {
    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol);

    m_Cup.set(0.0);
#pragma omp barrier

    m_staple->upper(m_Cup, m_Ueo, mu, nu);

    m_Cdn.set(0.0);
    m_staple->lower(m_Cdn, m_Ueo, mu, nu);

    copy(m_Umu, 0, m_Ueo, mu);
#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      m_ut3.set_mat(site, 0, m_Umu.mat(site) * m_Cup.mat_dag(site));
    }

    for (int site = is; site < ns; ++site) {
      m_ut2.set_mat(site, 0, m_Umu.mat(site) * m_Cdn.mat_dag(site));
    }
#pragma omp barrier

    axpy(m_ut3, -1.0, m_ut2);
#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      m_ut1.set_mat(site, 0, m_Cup.mat_dag(site) * m_Umu.mat(site));
    }

    for (int site = is; site < ns; ++site) {
      m_ut2.set_mat(site, 0, m_Cdn.mat_dag(site) * m_Umu.mat(site));
    }
#pragma omp barrier

    axpy(m_ut1, -1.0, m_ut2);
#pragma omp barrier

    m_shift_eo->forward(m_ut2, m_ut1, mu);

    axpy(m_ut3, 1.0, m_ut2);
#pragma omp barrier

    for (int site = is; site < ns; ++site) {
      Fst.set_mat(site, 0, m_ut3.mat(site).ah());
    }
#pragma omp barrier

    scal(Fst, 0.25);

#pragma omp barrier
  }


//====================================================================
  double Fopr_CloverTerm_eo::flop_count()
  {
    // Counting of floating point operations in giga unit.
    // The following counting explicitly depends on the implementation
    // and to be recalculated when the code is modified.
    // Present counting is based on rev.1107. [24 Aug 2014 H.Matsufuru]

    const int NPE = CommonParameters::NPE();

    int flop_site = 0; // superficial initialization

    if (m_repr == "Dirac") {
      flop_site = 8 * m_Nc * m_Nc * m_Nd * m_Nd;
    } else if (m_repr == "Chiral") {
      flop_site = 4 * m_Nc * m_Nc * m_Nd * m_Nd;
    }

    const double gflop = flop_site * ((m_Nvol / 2) * (NPE / 1.0e+9));

    return gflop;
  }


//====================================================================
}
//============================================================END=====
