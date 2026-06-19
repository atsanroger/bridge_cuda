/*!
        @file    fopr_Wilson_impl.cpp

        @brief

        @author  Hideo Matsufuru  (matufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "Fopr/Org/fopr_Wilson_impl.h"

namespace Org {
#ifdef USE_FACTORY_AUTOREGISTER
  namespace {
    bool init = Fopr_Wilson::register_factory();
  }
#endif

  const std::string Fopr_Wilson::class_name = "Org::Fopr_Wilson";

//====================================================================
  void Fopr_Wilson::init(const Parameters& params)
  {
    ThreadManager::assert_single_thread(class_name);

    m_vl = CommonParameters::Vlevel();

    vout.general(m_vl, "%s: construction\n", class_name.c_str());
    vout.increase_indent();

    std::string repr;
    if (!params.fetch_string("gamma_matrix_type", repr)) {
      m_repr = repr;
    } else {
      m_repr = "Dirac"; // default gamma-matrix type
      vout.general(m_vl, "gamma_matrix_type is not given: defalt = %s\n",
                   m_repr.c_str());
    }

    if ((m_repr != "Dirac") && (m_repr != "Chiral")) {
      vout.crucial("Error at %s: unsupported gamma-matrix type: %s\n",
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
  void Fopr_Wilson::init(const std::string repr)
  {
    ThreadManager::assert_single_thread(class_name);

    m_vl = CommonParameters::Vlevel();

    vout.general(m_vl, "%s: construction (obsolete)\n",
                 class_name.c_str());
    vout.increase_indent();

    m_repr = repr;
    vout.general(m_vl, "gamma-matrix type is set to %s\n",
                 m_repr.c_str());

    setup();

    vout.decrease_indent();
    vout.general(m_vl, "%s: construction finished.\n",
                 class_name.c_str());
  }


//====================================================================
  void Fopr_Wilson::setup()
  {
    m_Nc = CommonParameters::Nc();
    m_Nd = CommonParameters::Nd();
    int Nin = 2 * m_Nc * m_Nd;

    m_Nvol = CommonParameters::Nvol();
    m_Ndim = CommonParameters::Ndim();
    m_boundary.resize(m_Ndim);

    m_U = 0;

    m_shift = new ShiftField_lex(Nin);

    m_GM.resize(m_Ndim + 1);

    GammaMatrixSet *gmset = GammaMatrixSet::New(m_repr);

    m_GM[0] = gmset->get_GM(GammaMatrixSet::GAMMA1);
    m_GM[1] = gmset->get_GM(GammaMatrixSet::GAMMA2);
    m_GM[2] = gmset->get_GM(GammaMatrixSet::GAMMA3);
    m_GM[3] = gmset->get_GM(GammaMatrixSet::GAMMA4);
    m_GM[4] = gmset->get_GM(GammaMatrixSet::GAMMA5);

    delete gmset;

    m_w1.reset(Nin, m_Nvol, 1);
    m_w2.reset(Nin, m_Nvol, 1);
    m_w3.reset(Nin, m_Nvol, 1);

    m_v1.reset(m_Nvol, 1);  // Field_F
    m_v2.reset(m_Nvol, 1);  // Field_F
  }


//====================================================================
  void Fopr_Wilson::tidyup()
  {
    delete m_shift;
  }


//====================================================================
  void Fopr_Wilson::set_parameters(const Parameters& params)
  {
    std::string vlevel;
    if (!params.fetch_string("verbose_level", vlevel)) {
      m_vl = vout.set_verbose_level(vlevel);
    }

    //- fetch and check input parameters
    double           kappa;
    std::vector<int> bc;

    int err = 0;
    err += params.fetch_double("hopping_parameter", kappa);
    err += params.fetch_int_vector("boundary_condition", bc);

    if (err) {
      vout.crucial(m_vl, "Error at %s: input parameter not found.\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }

    set_parameters(kappa, bc);
  }


//====================================================================
  void Fopr_Wilson::set_parameters(const double kappa,
                                   const std::vector<int> bc)
  {
    //- range check
    // NB. kappa = 0 is allowed.
    assert(bc.size() == m_Ndim);

    int ith = ThreadManager::get_thread_id();
    if (ith == 0) {
      m_kappa    = kappa;
      m_boundary = bc;
    }

    //- print input parameters
    vout.general(m_vl, "%s: input parameters\n", class_name.c_str());
    vout.general(m_vl, "  gamma matrix type = %s\n", m_repr.c_str());
    vout.general(m_vl, "  kappa = %12.8f\n", m_kappa);
    for (int mu = 0; mu < m_Ndim; ++mu) {
      vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
    }
  }


//====================================================================
  void Fopr_Wilson::get_parameters(Parameters& params) const
  {
    params.set_double("hopping_parameter", m_kappa);
    params.set_int_vector("boundary_condition", m_boundary);
    params.set_string("gamma_matrix_type", m_repr);

    params.set_string("verbose_level", vout.get_verbose_level(m_vl));
  }


//====================================================================
  void Fopr_Wilson::set_config(Field *U)
  {
#pragma omp barrier

    int ith = ThreadManager::get_thread_id();
    if (ith == 0) m_U = (Field_G *)U;

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson::set_mode(const std::string mode)
  {
#pragma omp barrier
    int ith = ThreadManager::get_thread_id();
    if (ith == 0) m_mode = mode;
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson::mult(Field& v, const Field& w)
  {
    if (m_mode == "D") {
      D(v, w);
    } else if (m_mode == "Ddag") {
      Ddag(v, w);
    } else if (m_mode == "DdagD") {
      DdagD(v, w);
    } else if (m_mode == "DDdag") {
      DDdag(v, w);
    } else if (m_mode == "H") {
      H(v, w);
    } else {
      vout.crucial(m_vl, "Error at %s: input mode is undefined.\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_Wilson::mult_dag(Field& v, const Field& w)
  {
    if (m_mode == "D") {
      Ddag(v, w);
    } else if (m_mode == "Ddag") {
      D(v, w);
    } else if (m_mode == "DdagD") {
      DdagD(v, w);
    } else if (m_mode == "DDdag") {
      DDdag(v, w);
    } else if (m_mode == "H") {
      H(v, w);
    } else {
      vout.crucial(m_vl, "Error at %s: input mode is undefined.\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_Wilson::mult(Field& v, const Field& w,
                         const std::string mode)
  {
    if (mode == "D") {
      D(v, w);
    } else if (mode == "Ddag") {
      Ddag(v, w);
    } else if (mode == "DdagD") {
      DdagD(v, w);
    } else if (mode == "DDdag") {
      DDdag(v, w);
    } else if (mode == "H") {
      H(v, w);
    } else {
      vout.crucial(m_vl, "Error at %s: undefined mode: %s\n",
                   class_name.c_str(), mode.c_str());
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_Wilson::mult_dag(Field& v, const Field& w,
                             const std::string mode)
  {
    if (mode == "D") {
      Ddag(v, w);
    } else if (mode == "Ddag") {
      D(v, w);
    } else if (mode == "DdagD") {
      DdagD(v, w);
    } else if (mode == "DDdag") {
      DDdag(v, w);
    } else if (mode == "H") {
      H(v, w);
    } else {
      vout.crucial(m_vl, "Error at %s: undefined mode: %s\n",
                   class_name.c_str(), mode.c_str());
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_Wilson::DdagD(Field& v, const Field& w)
  {
    D(m_w1, w);
    mult_gm5(v, m_w1);
    D(m_w1, v);
    mult_gm5(v, m_w1);
  }


//====================================================================
  void Fopr_Wilson::Ddag(Field& v, const Field& w)
  {
    mult_gm5(v, w);
    D(m_w1, v);
    mult_gm5(v, m_w1);
  }


//====================================================================
  void Fopr_Wilson::DDdag(Field& v, const Field& w)
  {
    mult_gm5(m_w1, w);
    D(v, m_w1);
    mult_gm5(m_w1, v);
    D(v, m_w1);
  }


//====================================================================
  void Fopr_Wilson::H(Field& v, const Field& w)
  {
    D(m_w1, w);
    mult_gm5(v, m_w1);
  }


//====================================================================
  void Fopr_Wilson::D(Field& v, const Field& w)
  {
#pragma omp barrier

    v.set(0.0);
    for (int mu = 0; mu < m_Ndim; ++mu) {
      mult_up(mu, v, w);
      mult_dn(mu, v, w);
    }
    aypx(-m_kappa, v, w);

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson::D_ex(Field& v, const int ex1,
                         const Field& w, const int ex2)
  {
#pragma omp barrier

    copy(m_w2, 0, w, ex2);
    m_w3.set(0.0);

    for (int mu = 0; mu < m_Ndim; ++mu) {
      mult_up(mu, m_w3, m_w2);
      mult_dn(mu, m_w3, m_w2);
    }

    aypx(-m_kappa, m_w3, m_w2);
    copy(v, ex1, m_w3, 0);

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson::mult_gm5(Field& v, const Field& w)
  {
    int Nin = 2 * m_Nc * m_Nd;
    assert(w.check_size(Nin, m_Nvol, 1));
    assert(v.check_size(Nin, m_Nvol, 1));

#pragma omp barrier

    copy(m_v1, w);
    mult_GM(m_v2, m_GM[4], m_v1);
    copy(v, m_v2);

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson::proj_chiral(Field& v, const int ex1,
                                const Field& w, const int ex2,
                                const int ipm)
  {
    assert(ipm == 1 || ipm == -1);

    double fpm = 1.0;
    if (ipm == -1) fpm = -1.0;

    copy(m_w1, 0, w, ex2);
    mult_gm5(m_w2, m_w1);
    axpy(m_w1, ex1, fpm, m_w2, 0);
    copy(v, ex1, m_w1, 0);

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson::mult_gm5p(const int mu, Field& v, const Field& w)
  {
    assert(mu >= 0 && mu < m_Ndim);

    m_w1.set(0.0);
    mult_up(mu, m_w1, w);
    mult_gm5(v, m_w1);
  }


//====================================================================
  void Fopr_Wilson::mult_up(const int mu, Field& v, const Field& w)
  {
#pragma omp barrier

    for (int ex = 0; ex < w.nex(); ++ex) {
      copy(m_v2, 0, w, ex);
      m_shift->backward(m_v1, m_v2, m_boundary[mu], mu);
      mult_Field_Gn(m_v2, 0, *m_U, mu, m_v1, 0);
      mult_GMproj2(m_v1, -1, m_GM[mu], m_v2);
      axpy(v, ex, 1.0, m_v1, 0);
    }
#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson::mult_dn(const int mu, Field& v, const Field& w)
  {
#pragma omp barrier

    for (int ex = 0; ex < w.nex(); ++ex) {
      copy(m_v1, 0, w, ex);
      mult_Field_Gd(m_v2, 0, *m_U, mu, m_v1, ex);
      m_shift->forward(m_v1, m_v2, m_boundary[mu], mu);
      mult_GMproj2(m_v2, 1, m_GM[mu], m_v1);
      axpy(v, ex, 1.0, m_v2, 0);
    }
#pragma omp barrier
  }


//====================================================================
  double Fopr_Wilson::flop_count()
  {
    // This counting is based on the Org-implementation of ver.1.2.0.
    // Flop count of mult_GMproj2 is different for upward and downward
    // directions due to the implemetation in Field_F.cpp.
    // The present counting is based on rev.1130. [10 Sep 2014 H.Matsufuru]

    const int Nvol = CommonParameters::Nvol();
    const int NPE  = CommonParameters::NPE();

    const int Nc = m_Nc;
    const int Nd = m_Nd;

    int flop_per_site = Nc * Nd * 2 * 8 * (4 * Nc - 1); // #(mult_Field_Gn/d)

    flop_per_site += Nc * Nd * 2 * (4 * 3 + 4 * 2);     // #(mult_GMproj2)
    flop_per_site += Nc * Nd * 2 * 8;                   // #(addpart_ex)
    flop_per_site += Nc * Nd * 2 * 2;                   // #(aypx(kappa))

    double gflop = flop_per_site * (Nvol * (NPE / 1.0e+9));

    if ((m_mode == "DdagD") || (m_mode == "DDdag")) {
      gflop *= 2;
    }

    return gflop;
  }


//====================================================================
}
//============================================================END=====
