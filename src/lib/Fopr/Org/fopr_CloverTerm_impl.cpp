/*!
        @file    fopr_CloverTerm_impl.cpp

        @brief

        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "fopr_CloverTerm_impl.h"

namespace Org {
  const std::string Fopr_CloverTerm::class_name = "Org::Fopr_CloverTerm";

//====================================================================
  void Fopr_CloverTerm::init(const Parameters& params)
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
  void Fopr_CloverTerm::init(const std::string repr)
  {
    ThreadManager::assert_single_thread(class_name);

    m_vl = CommonParameters::Vlevel();

    vout.general(m_vl, "%s: construction (obsolete)\n",
                 class_name.c_str());

    m_repr = repr;

    setup();

    vout.general(m_vl, "%s: construction finished.\n",
                 class_name.c_str());
  }


//====================================================================
  void Fopr_CloverTerm::setup()
  {
    m_Nvol = CommonParameters::Nvol();
    m_Ndim = CommonParameters::Ndim();
    m_Nc   = CommonParameters::Nc();
    m_Nd   = CommonParameters::Nd();
    m_NinF = 2 * m_Nc * m_Nd;

    m_U = 0;

    m_boundary.resize(m_Ndim);

    m_shift = new ShiftField_lex(2 * m_Nc * m_Nc);

    m_SG.resize(m_Ndim * m_Ndim);

    GammaMatrixSet *gmset = GammaMatrixSet::New(m_repr);

    m_GM5 = gmset->get_GM(gmset->GAMMA5);

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
  void Fopr_CloverTerm::tidyup()
  {
    delete m_shift;
  }


//====================================================================
  void Fopr_CloverTerm::set_parameters(const Parameters& params)
  {
#pragma omp barrier

    int ith = ThreadManager::get_thread_id();

    std::string vlevel;
    if (!params.fetch_string("verbose_level", vlevel)) {
      if (ith == 0) m_vl = vout.set_verbose_level(vlevel);
    }
#pragma omp barrier

    double           kappa, cSW;
    std::vector<int> bc;

    int err = 0;
    err += params.fetch_double("hopping_parameter", kappa);
    err += params.fetch_double("clover_coefficient", cSW);
    err += params.fetch_int_vector("boundary_condition", bc);

    if (err) {
      vout.crucial(m_vl, "Error at %s: input parameter not found\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }

    set_parameters(kappa, cSW, bc);
  }


//====================================================================
  void Fopr_CloverTerm::set_parameters(const double kappa,
                                       const double cSW,
                                       const std::vector<int> bc)
  {
    assert(bc.size() == m_Ndim);

#pragma omp barrier

    int ith = ThreadManager::get_thread_id();

    if (ith == 0) {
      m_kappa    = kappa;
      m_cSW      = cSW;
      m_boundary = bc;
    }

#pragma omp barrier

    vout.general(m_vl, "%s: input parameters\n", class_name.c_str());
    vout.general(m_vl, "  gamma matrix type = %s\n", m_repr.c_str());
    vout.general(m_vl, "  kappa = %12.8f\n", kappa);
    vout.general(m_vl, "  cSW   = %12.8f\n", cSW);
    for (int mu = 0; mu < m_Ndim; ++mu) {
      vout.general(m_vl, "  boundary[%d] = %2d\n", mu, bc[mu]);
    }
  }


//====================================================================
  void Fopr_CloverTerm::get_parameters(Parameters& params) const
  {
    params.set_double("hopping_parameter", m_kappa);
    params.set_double("clover_coefficient", m_cSW);
    params.set_int_vector("boundary_condition", m_boundary);
    params.set_string("gamma_matrix_type", m_repr);

    params.set_string("verbose_level", vout.get_verbose_level(m_vl));
  }


//====================================================================
  void Fopr_CloverTerm::set_config(Field *U)
  {
    int nth = ThreadManager::get_num_threads();
    vout.detailed(m_vl, "%s: set_config is called: num_threads = %d\n",
                  class_name.c_str(), nth);

    if (nth > 1) {
      set_config_impl(U);
    } else {
      set_config_omp(U);
    }

    vout.detailed(m_vl, "%s: set_config finished.\n", class_name.c_str());
  }


//====================================================================
  void Fopr_CloverTerm::set_config_omp(Field *U)
  {
#pragma omp parallel
    {
      set_config_impl(U);
    }
  }


//====================================================================
  void Fopr_CloverTerm::set_config_impl(Field *U)
  {
#pragma omp barrier

    int ith = ThreadManager::get_thread_id();
    if (ith == 0) m_U = (Field_G *)U;

#pragma omp barrier

    set_csw();

#pragma omp barrier
  }


//====================================================================
  void Fopr_CloverTerm::set_mode(const std::string mode)
  {
#pragma omp barrier
    int ith = ThreadManager::get_thread_id();
    if (ith == 0) m_mode = mode;
#pragma omp barrier
  }


//====================================================================
  void Fopr_CloverTerm::mult(Field& v, const Field& w)
  {
    // csw kappa sigma_{mu nu} F_{mu nu}

    if ((m_mode == "D") || (m_mode == "F")) {
      mult_sigmaF(v, w);
    } else {
      vout.crucial(m_vl, "Error at %s: undefined mode = %s\n",
                   class_name.c_str(), m_mode.c_str());
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_CloverTerm::mult_dag(Field& v, const Field& w)
  {
    // csw kappa sigma_{mu nu} F_{mu nu}

    if ((m_mode == "D") || (m_mode == "F")) {
      mult_sigmaF(v, w);
    } else {
      vout.crucial(m_vl, "Error at %s: undefined mode = %s\n",
                   class_name.c_str(), m_mode.c_str());
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_CloverTerm::mult_gm5(Field& v, const Field& w)
  {
    assert(w.check_size(m_NinF, m_Nvol, 1));
    assert(v.check_size(m_NinF, m_Nvol, 1));

#pragma omp barrier

    copy(m_v1, w);
    mult_GM(m_v2, m_GM5, m_v1);
    copy(v, m_v2);

#pragma omp barrier
  }


//====================================================================
  void Fopr_CloverTerm::mult_isigma(Field_F& v, const Field_F& w,
                                    const int mu, const int nu)
  {
    assert(mu != nu);

#pragma omp barrier

    mult_iGM(v, m_SG[sg_index(mu, nu)], w);

#pragma omp barrier
  }


//====================================================================
  void Fopr_CloverTerm::mult_sigmaF(Field& v, const Field& f)
  {
    mult_csw(v, f);
  }


//====================================================================
  void Fopr_CloverTerm::mult_csw(Field& v, const Field& w)
  {
    // multiplies csw kappa sigma_{mu nu} F_{mu nu}
    // NOTE: this is NOT  1 - csw kappa sigma_{mu nu} F_{mu nu}

#pragma omp barrier

    copy(m_v1, w);

    m_v3.set(0.0);

    mult_iGM(m_v2, m_SG[sg_index(1, 2)], m_v1);
    multadd_Field_Gn(m_v3, 0, m_Bx, 0, m_v2, 0, 1.0);

    mult_iGM(m_v2, m_SG[sg_index(2, 0)], m_v1);
    multadd_Field_Gn(m_v3, 0, m_By, 0, m_v2, 0, 1.0);

    mult_iGM(m_v2, m_SG[sg_index(0, 1)], m_v1);
    multadd_Field_Gn(m_v3, 0, m_Bz, 0, m_v2, 0, 1.0);

    mult_iGM(m_v2, m_SG[sg_index(3, 0)], m_v1);
    multadd_Field_Gn(m_v3, 0, m_Ex, 0, m_v2, 0, 1.0);

    mult_iGM(m_v2, m_SG[sg_index(3, 1)], m_v1);
    multadd_Field_Gn(m_v3, 0, m_Ey, 0, m_v2, 0, 1.0);

    mult_iGM(m_v2, m_SG[sg_index(3, 2)], m_v1);
    multadd_Field_Gn(m_v3, 0, m_Ez, 0, m_v2, 0, 1.0);

    scal(m_v3, m_kappa * m_cSW);

    copy(v, m_v3);

#pragma omp barrier
  }


//====================================================================
  void Fopr_CloverTerm::set_csw()
  {
    set_fieldstrength(m_Bx, 1, 2);
    set_fieldstrength(m_By, 2, 0);
    set_fieldstrength(m_Bz, 0, 1);
    set_fieldstrength(m_Ex, 3, 0);
    set_fieldstrength(m_Ey, 3, 1);
    set_fieldstrength(m_Ez, 3, 2);
  }


//====================================================================
  void Fopr_CloverTerm::set_fieldstrength(Field_G& Fst,
                                          const int mu, const int nu)
  {
    copy(m_Umu, 0, *m_U, mu);
    m_staple.upper(m_Cup, *m_U, mu, nu);
    m_staple.lower(m_Cdn, *m_U, mu, nu);

    mult_Field_Gnd(Fst, 0, m_Umu, 0, m_Cup, 0);
    multadd_Field_Gnd(Fst, 0, m_Umu, 0, m_Cdn, 0, -1.0);

    mult_Field_Gdn(m_u1, 0, m_Cup, 0, m_Umu, 0);
    multadd_Field_Gdn(m_u1, 0, m_Cdn, 0, m_Umu, 0, -1.0);

    m_shift->forward(m_u2, m_u1, mu);
    axpy(Fst, 1.0, m_u2);

    ah_Field_G(Fst, 0);
    scal(Fst, 0.25);
  }


//====================================================================
  double Fopr_CloverTerm::flop_count()
  {
    // Counting of floating point operations in giga unit.
    // The following counting explicitly depends on the implementation
    // and to be recalculated when the code is modified.
    // Present counting is based on rev.1131. [10 Sep 2014 H.Matsufuru]

    const int Nvol = CommonParameters::Nvol();
    const int NPE  = CommonParameters::NPE();

    const int flop_site = m_Nc * m_Nd * (2 + 12 + 48 * m_Nc);

    const double gflop = flop_site * (Nvol * (NPE / 1.0e+9));

    return gflop;
  }


//====================================================================
}
//============================================================END=====
