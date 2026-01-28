/*!
        @file    fopr_Wilson_eo_impl.cpp

        @brief

        @author  Satoru UEDA (maintailed by H.Matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "Fopr/Org/fopr_Wilson_eo_impl.h"

#include "Field/field_thread-inc.h"

namespace Org {
#ifdef USE_FACTORY_AUTOREGISTER
  namespace {
    bool init = Fopr_Wilson_eo::register_factory();
  }
#endif

  const std::string Fopr_Wilson_eo::class_name = "Org::Fopr_Wilson_eo";

//====================================================================
  void Fopr_Wilson_eo::init(const Parameters& params)
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
  void Fopr_Wilson_eo::init(const std::string repr)
  {
    ThreadManager::assert_single_thread(class_name);

    m_vl = CommonParameters::Vlevel();

    vout.general(m_vl, "%s: construction (obsolete)\n",
                 class_name.c_str());
    vout.increase_indent();

    m_repr = repr;
    if ((m_repr != "Dirac") && (m_repr != "Chiral")) {
      vout.crucial(m_vl, "Error at %s: input repr is undefined.\n",
                   class_name.c_str());
      exit(EXIT_FAILURE);
    }

    vout.general(m_vl, "gamma-matrix type is set to %s\n",
                 m_repr.c_str());

    vout.decrease_indent();
    vout.general(m_vl, "%s: construction finished.\n",
                 class_name.c_str());

    setup();
  }


//====================================================================
  void Fopr_Wilson_eo::setup()
  {
    m_Nc    = CommonParameters::Nc();
    m_Nd    = CommonParameters::Nd();
    m_Nvol  = CommonParameters::Nvol();
    m_Nvol2 = m_Nvol / 2;
    m_Ndim  = CommonParameters::Ndim();

    m_boundary.resize(m_Ndim);
    m_GM.resize(m_Ndim + 1);

    GammaMatrixSet *gmset = GammaMatrixSet::New(m_repr);

    m_GM[0] = gmset->get_GM(gmset->GAMMA1);
    m_GM[1] = gmset->get_GM(gmset->GAMMA2);
    m_GM[2] = gmset->get_GM(gmset->GAMMA3);
    m_GM[3] = gmset->get_GM(gmset->GAMMA4);
    m_GM[4] = gmset->get_GM(gmset->GAMMA5);

    delete gmset;

    int NinF = 2 * m_Nc * m_Nd;
    m_shift = new ShiftField_eo(NinF);

    m_Ueo.reset(m_Nvol, m_Ndim);

    m_t1.reset(m_Nvol2, 1);
    m_t2.reset(m_Nvol2, 1);
    m_w1.reset(m_Nvol2, 1);
    m_w2.reset(m_Nvol2, 1);
    m_v1.reset(m_Nvol2, 1);
    m_v2.reset(m_Nvol2, 1);
    m_v3.reset(m_Nvol2, 1);
  }


//====================================================================
  void Fopr_Wilson_eo::tidyup()
  {
    delete m_shift;
  }


//====================================================================
  void Fopr_Wilson_eo::set_parameters(const Parameters& params)
  {
#pragma omp barrier
    int         ith = ThreadManager::get_thread_id();
    std::string vlevel;
    if (!params.fetch_string("verbose_level", vlevel)) {
      if (ith == 0) m_vl = vout.set_verbose_level(vlevel);
    }
#pragma omp barrier

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
  void Fopr_Wilson_eo::set_parameters(const double kappa,
                                      const std::vector<int> bc)
  {
    assert(bc.size() == m_Ndim);

#pragma omp barrier

    int ith = ThreadManager::get_thread_id();
    if (ith == 0) {
      m_kappa    = kappa;
      m_boundary = bc;
    }

#pragma omp barrier

    //- print input parameters
    vout.general(m_vl, "%s: input parameters\n", class_name.c_str());
    vout.general(m_vl, "  gamma matrix type = %s\n", m_repr.c_str());
    vout.general(m_vl, "  kappa = %12.8f\n", m_kappa);
    for (int mu = 0; mu < m_Ndim; ++mu) {
      vout.general(m_vl, "  boundary[%d] = %2d\n", mu, m_boundary[mu]);
    }
  }


//====================================================================
  void Fopr_Wilson_eo::get_parameters(Parameters& params) const
  {
    params.set_double("hopping_parameter", m_kappa);
    params.set_int_vector("boundary_condition", m_boundary);
    params.set_string("gamma_matrix_type", m_repr);

    params.set_string("verbose_level", vout.get_verbose_level(m_vl));
  }


//====================================================================
  void Fopr_Wilson_eo::set_config(Field *U)
  {
#pragma omp barrier

    m_index.convertField(m_Ueo, *U);

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::set_mode(const std::string mode)
  {
#pragma omp barrier

    int ith = ThreadManager::get_thread_id();
    if (ith == 0) m_mode = mode;

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::mult(Field& v, const Field& w)
  {
    if (m_mode == "D") {
      D(v, w);
    } else if (m_mode == "Ddag") {
      Ddag(v, w);
    } else if (m_mode == "DdagD") {
      DdagD(v, w);
    } else if (m_mode == "DDdag") {
      DDdag(v, w);
    } else {
      vout.crucial("Error at %s: irrelevant mult mode = %s.\n",
                   class_name.c_str(), m_mode.c_str());
      exit(EXIT_FAILURE);
    }
  }


//====================================================================
  void Fopr_Wilson_eo::mult_dag(Field& v, const Field& w)
  {
    if (m_mode == "D") {
      Ddag(v, w);
    } else if (m_mode == "Ddag") {
      D(v, w);
    } else if (m_mode == "DdagD") {
      DdagD(v, w);
    } else if (m_mode == "DDdag") {
      DDdag(v, w);
    } else {
      vout.crucial("Error at %s: irrelevant mult mode = %s.\n",
                   class_name.c_str(), m_mode.c_str());
      exit(EXIT_FAILURE);
    }
  }


//=====================================================================
  void Fopr_Wilson_eo::mult(Field& v, const Field& w,
                            const std::string mode)
  {
    if (mode == "Deo") {
      Meo(v, w, 0);
    } else if (mode == "Doe") {
      Meo(v, w, 1);
    } else {
      vout.crucial("Error at %s: irrelevant mult mode = %s.\n",
                   class_name.c_str(), mode.c_str());
      exit(EXIT_FAILURE);
    }
  }


//=====================================================================
  void Fopr_Wilson_eo::mult_dag(Field& v, const Field& w,
                                const std::string mode)
  {
    if (mode == "Deo") {
      Meo(v, w, 1);               // "Doe" for mult
    } else if (mode == "Doe") {
      Meo(v, w, 0);               // "Deo" for mult
    } else {
      vout.crucial("Error at %s: irrelevant mult mode = %s.\n",
                   class_name.c_str(), mode.c_str());
      exit(EXIT_FAILURE);
    }
  }


//=====================================================================
  void Fopr_Wilson_eo::mult_gm5(Field& v, const Field& w)
  {
#pragma omp barrier

    copy(m_w1, w);
    mult_GM(m_w2, m_GM[4], m_w1);
    copy(v, m_w2);

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::preProp(Field& Be, Field& bo, const Field& b)
  {
    int Nin = 2 * m_Nc * m_Nd;
    assert(Be.check_size(Nin, m_Nvol2, 1));
    assert(bo.check_size(Nin, m_Nvol2, 1));
    assert(b.check_size(Nin, m_Nvol, 1));

#pragma omp barrier

    m_index.convertField(m_v1, b, 0);
    m_index.convertField(bo, b, 1);

    copy(Be, m_v1);

    if (m_mode == "D") {
      Meo(m_v1, bo, 0);
    } else if (m_mode == "Ddag") {
      Mdageo(m_v1, bo, 0);
    } else {
      vout.crucial("Error at %s: irrelevant preProp mode = %s.\n",
                   class_name.c_str(), m_mode.c_str());
      exit(EXIT_FAILURE);
    }

    axpy(Be, -1.0, m_v1);

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::postProp(Field& x, const Field& xe, const Field& bo)
  {
    int Nin = 2 * m_Nc * m_Nd;
    assert(bo.check_size(Nin, m_Nvol2, 1));
    assert(xe.check_size(Nin, m_Nvol2, 1));
    assert(x.check_size(Nin, m_Nvol, 1));

#pragma omp barrier

    copy(m_v2, bo);

    if (m_mode == "D") {
      Meo(m_v1, xe, 1);
    } else if (m_mode == "Ddag") {
      Mdageo(m_v1, xe, 1);
    } else {
      vout.crucial("Error at %s: irrelevant postProp mode = %s.\n",
                   class_name.c_str(), m_mode.c_str());
      exit(EXIT_FAILURE);
    }

    axpy(m_v2, -1.0, m_v1);

    m_index.reverseField(x, xe, 0);
    m_index.reverseField(x, m_v2, 1);

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::DdagD(Field& v, const Field& w)
  {
    D(m_v3, w);
    Ddag(v, m_v3);
  }


//====================================================================
  void Fopr_Wilson_eo::DDdag(Field& v, const Field& w)
  {
    Ddag(m_v3, w);
    D(v, m_v3);
  }


//====================================================================
  void Fopr_Wilson_eo::H(Field& v, const Field& w)
  {
    D(m_v3, w);
    mult_gm5(v, m_v3);
  }


//====================================================================
  void Fopr_Wilson_eo::Ddag(Field& v, const Field& w)
  {
    assert(w.check_size(2 * m_Nc * m_Nd, m_Nvol2, 1));
    assert(v.check_size(2 * m_Nc * m_Nd, m_Nvol2, 1));

#pragma omp barrier

    copy(v, w);

    Mdageo(m_v1, w, 1);

    Mdageo(m_v2, m_v1, 0);

    axpy(v, -1.0, m_v2);

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::D(Field& v, const Field& w)
  {
    assert(w.check_size(2 * m_Nc * m_Nd, m_Nvol2, 1));
    assert(v.check_size(2 * m_Nc * m_Nd, m_Nvol2, 1));

#pragma omp barrier

    copy(v, w);

    Meo(m_v1, w, 1);

    Meo(m_v2, m_v1, 0);

    axpy(v, -1.0, m_v2);

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::Meo(Field& v, const Field& w, const int ieo)
  {
#pragma omp barrier

    copy(m_w1, w);
    m_w2.set(0.0);

    for (int mu = 0; mu < m_Ndim; ++mu) {
      mult_up(mu, m_w2, m_w1, ieo);
      mult_dn(mu, m_w2, m_w1, ieo);
    }

    scal(m_w2, -m_kappa);

    copy(v, m_w2);

#pragma omp barrier
  }


//====================================================================
  void Fopr_Wilson_eo::Mdageo(Field& v, const Field& w, const int ieo)
  {
#pragma omp barrier

    copy(m_w2, w);
    mult_GM(m_w1, m_GM[4], m_w2);

    m_w2.set(0.0);

    for (int mu = 0; mu < m_Ndim; ++mu) {
      mult_up(mu, m_w2, m_w1, ieo);
      mult_dn(mu, m_w2, m_w1, ieo);
    }

    scal(m_w2, -m_kappa);

    mult_GM(m_w1, m_GM[4], m_w2);

    copy(v, m_w1);

#pragma omp barrier
  }


//=====================================================================
  void Fopr_Wilson_eo::gm5p(const int mu, Field& v, const Field& w)
  {
    vout.crucial(m_vl, "Error at %s: gm5p is undefined.\n",
                 class_name.c_str());
    exit(EXIT_FAILURE);
  }


//=====================================================================
  void Fopr_Wilson_eo::mult_up(const int mu, Field_F& v, const Field_F& w,
                               const int ieo)
  {
    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol2);

    int Nex = w.nex();
    for (int ex = 0; ex < Nex; ++ex) {
      copy(m_t2, 0, w, ex);

      m_shift->backward_h(m_t1, m_t2, m_boundary[mu], mu, ieo);

      for (int site = is; site < ns; ++site) {
        for (int id = 0; id < m_Nd; ++id) {
          Vec_SU_N vec(m_Nc);
          vec = m_Ueo.mat(site + ieo * m_Nvol2, mu) * m_t1.vec(id, site, 0);
          m_t2.set_vec(id, site, 0, vec);
        }
      }

      mult_GMproj2(m_t1, -1, m_GM[mu], m_t2);

      axpy(v, ex, 1.0, m_t1, 0);
    }
  }


//=====================================================================
  void Fopr_Wilson_eo::mult_dn(const int mu, Field_F& v, const Field_F& w,
                               const int ieo)
  {
    int jeo = 1 - ieo;

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, m_Nvol2);

    int Nex = w.nex();
    for (int ex = 0; ex < Nex; ++ex) {
      for (int site = is; site < ns; ++site) {
        for (int id = 0; id < m_Nd; ++id) {
          Vec_SU_N vec(m_Nc);
          vec = m_Ueo.mat_dag(site + jeo * m_Nvol2, mu)
                * w.vec(id, site, ex);
          m_t1.set_vec(id, site, 0, vec);
        }
      }
      m_shift->forward_h(m_t2, m_t1, m_boundary[mu], mu, ieo);
      mult_GMproj2(m_t1, 1, m_GM[mu], m_t2);
      axpy(v, ex, 1.0, m_t1, 0);
    }
  }


//====================================================================
  double Fopr_Wilson_eo::flop_count()
  {
    // Counting of floating point operations in giga unit.
    // this counting is based on the Org-implementation of ver.1.2.0.
    // Flop count of mult_GMproj2 is different for upward and downward
    // directions due to the implemetation in Field_F.cpp.
    // Present counting is based on rev.1130. [10 Sep 2014 H.Matsufuru]

    const int Nvol = CommonParameters::Nvol();
    const int NPE  = CommonParameters::NPE();

    const int Nc = m_Nc;
    const int Nd = m_Nd;

    int flop_Meo = Nc * Nd * 2 * 8 * (4 * Nc - 1); // #(mult_Field_Gn/d)

    flop_Meo += Nc * Nd * 2 * (4 * 3 + 4 * 2);     // #(mult_GMproj2)
    flop_Meo += Nc * Nd * 2 * 8;                   // #(addpart_ex)
    flop_Meo += Nc * Nd * 2;                       // #(scal(kappa))

    const int flop_per_site = 2 * flop_Meo + Nc * Nd * 2 * 2;
    // #(2*Meo + axpy)

    double gflop = flop_per_site * ((Nvol / 2) * (NPE / 1.0e+9));

    if ((m_mode == "DdagD") || (m_mode == "DDdag")) {
      gflop *= 2;
    }

    return gflop;
  }


//====================================================================
} // namespace Org
//============================================================END=====
