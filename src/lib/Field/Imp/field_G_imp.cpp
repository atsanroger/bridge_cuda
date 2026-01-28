/*!
        @file    field_G_imp.cpp

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "Field/field_G.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

#include "Field/field_thread-inc.h"

#if defined USE_GROUP_SU3
#include "field_G_imp_SU3-inc.h"
#elif defined USE_GROUP_SU2
#include "field_G_imp_SU2-inc.h"
#elif defined USE_GROUP_SU_N
#include "field_G_imp_SU_N-inc.h"
#endif

//====================================================================
void Field_G::check()
{
  assert(NC == CommonParameters::Nc());

  check_Nc();
}


//====================================================================
void Field_G::set_unit()
{
  Mat_SU_N ut(nc());

  ut.unit();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

  for (int mu = 0; mu < m_Nex; ++mu) {
    for (int site = is; site < ns; ++site) {
      set_mat(site, mu, ut);
    }
  }
}


//====================================================================
void Field_G::set_random(RandomNumbers *rand)
{
  rand->gauss_lex_global(*this);

  Mat_SU_N ut(m_Nc);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

  for (int mu = 0; mu < m_Nex; ++mu) {
    for (int site = is; site < m_Nvol; ++site) {
      this->mat(ut, site, mu);
      ut.reunit();
      this->set_mat(site, mu, ut);
    }
  }
}


//====================================================================
void Field_G::reunit()
{
  Mat_SU_N ut(nc());

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, m_Nvol);

  for (int mu = 0; mu < m_Nex; ++mu) {
    for (int site = is; site < ns; ++site) {
      mat(ut, site, mu);
      ut.reunit();
      set_mat(site, mu, ut);
    }
  }
}


//====================================================================
void mult_Field_Gnn(Field_G& W, const int ex,
                    const Field_G& U1, const int ex1,
                    const Field_G& U2, const int ex2)
{
  int Nvol = W.nvol();

  assert(ex < W.nex());
  assert(ex1 < U1.nex());
  assert(ex2 < U2.nex());
  assert(U1.nvol() == W.nvol());
  assert(U2.nvol() == W.nvol());

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  double       *g  = W.ptr(0, 0, ex);
  const double *g1 = U1.ptr(0, 0, ex1);
  const double *g2 = U2.ptr(0, 0, ex2);

  const int Nc  = W.nc();
  const int Nc2 = 2 * Nc;
  const int Ndf = Nc2 * Nc;

  for (int site = is; site < ns; ++site) {
    int ig = Ndf * site;

    for (int ic2 = 0; ic2 < Nc; ++ic2) {
      for (int ic1 = 0; ic1 < Nc; ++ic1) {
        int jg2 = ic2 * 2 + ig;
        int jg1 = ic1 * Nc2 + ig;
        g[2 * ic2 + Nc2 * ic1 + ig]     = mult_Gnn_r(&g1[jg1], &g2[jg2], Nc);
        g[2 * ic2 + 1 + Nc2 * ic1 + ig] = mult_Gnn_i(&g1[jg1], &g2[jg2], Nc);
      }
    }
  }
}


//====================================================================
void mult_Field_Gdn(Field_G& W, const int ex,
                    const Field_G& U1, const int ex1,
                    const Field_G& U2, const int ex2)
{
  int Nvol = W.nvol();

  assert(ex < W.nex());
  assert(ex1 < U1.nex());
  assert(ex2 < U2.nex());
  assert(U1.nvol() == W.nvol());
  assert(U2.nvol() == W.nvol());

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  double       *g  = W.ptr(0, 0, ex);
  const double *g1 = U1.ptr(0, 0, ex1);
  const double *g2 = U2.ptr(0, 0, ex2);

  const int Nc  = W.nc();
  const int Nc2 = 2 * Nc;
  const int Ndf = Nc2 * Nc;

  for (int site = is; site < ns; ++site) {
    int ig = Ndf * site;

    for (int ic2 = 0; ic2 < Nc; ++ic2) {
      for (int ic1 = 0; ic1 < Nc; ++ic1) {
        int jg2 = ic2 * 2 + ig;
        int jg1 = ic1 * 2 + ig;
        g[2 * ic2 + Nc2 * ic1 + ig]     = mult_Gdn_r(&g1[jg1], &g2[jg2], Nc);
        g[2 * ic2 + 1 + Nc2 * ic1 + ig] = mult_Gdn_i(&g1[jg1], &g2[jg2], Nc);
      }
    }
  }
}


//====================================================================
void mult_Field_Gnd(Field_G& W, const int ex,
                    const Field_G& U1, const int ex1,
                    const Field_G& U2, const int ex2)
{
  int Nvol = W.nvol();

  assert(ex < W.nex());
  assert(ex1 < U1.nex());
  assert(ex2 < U2.nex());
  assert(U1.nvol() == W.nvol());
  assert(U2.nvol() == W.nvol());

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  double       *g  = W.ptr(0, 0, ex);
  const double *g1 = U1.ptr(0, 0, ex1);
  const double *g2 = U2.ptr(0, 0, ex2);

  const int Nc  = W.nc();
  const int Nc2 = 2 * Nc;
  const int Ndf = Nc2 * Nc;

  for (int site = is; site < ns; ++site) {
    int ig = Ndf * site;

    for (int ic2 = 0; ic2 < Nc; ++ic2) {
      for (int ic1 = 0; ic1 < Nc; ++ic1) {
        int jg2 = ic2 * Nc2 + ig;
        int jg1 = ic1 * Nc2 + ig;
        g[2 * ic2 + Nc2 * ic1 + ig]     = mult_Gnd_r(&g1[jg1], &g2[jg2], Nc);
        g[2 * ic2 + 1 + Nc2 * ic1 + ig] = mult_Gnd_i(&g1[jg1], &g2[jg2], Nc);
      }
    }
  }
}


//====================================================================
void mult_Field_Gdd(Field_G& W, const int ex,
                    const Field_G& U1, const int ex1,
                    const Field_G& U2, const int ex2)
{
  int Nvol = W.nvol();

  assert(ex < W.nex());
  assert(ex1 < U1.nex());
  assert(ex2 < U2.nex());
  assert(U1.nvol() == W.nvol());
  assert(U2.nvol() == W.nvol());

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  double       *g  = W.ptr(0, 0, ex);
  const double *g1 = U1.ptr(0, 0, ex1);
  const double *g2 = U2.ptr(0, 0, ex2);

  const int Nc  = W.nc();
  const int Nc2 = 2 * Nc;
  const int Ndf = Nc2 * Nc;

  for (int site = is; site < ns; ++site) {
    int ig = Ndf * site;

    for (int ic2 = 0; ic2 < Nc; ++ic2) {
      for (int ic1 = 0; ic1 < Nc; ++ic1) {
        int jg2 = ic2 * Nc2 + ig;
        int jg1 = ic1 * 2 + ig;
        g[2 * ic2 + Nc2 * ic1 + ig]     = mult_Gdd_r(&g1[jg1], &g2[jg2], Nc);
        g[2 * ic2 + 1 + Nc2 * ic1 + ig] = mult_Gdd_i(&g1[jg1], &g2[jg2], Nc);
      }
    }
  }
}


//====================================================================
void multadd_Field_Gnn(Field_G& W, const int ex,
                       const Field_G& U1, const int ex1,
                       const Field_G& U2, const int ex2,
                       const double ff)
{
  int Nvol = W.nvol();

  assert(ex < W.nex());
  assert(ex1 < U1.nex());
  assert(ex2 < U2.nex());
  assert(U1.nvol() == W.nvol());
  assert(U2.nvol() == W.nvol());

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  double       *g  = W.ptr(0, 0, ex);
  const double *g1 = U1.ptr(0, 0, ex1);
  const double *g2 = U2.ptr(0, 0, ex2);

  const int Nc  = W.nc();
  const int Nc2 = 2 * Nc;
  const int Ndf = Nc2 * Nc;

  for (int site = is; site < ns; ++site) {
    int ig = Ndf * site;

    for (int ic2 = 0; ic2 < Nc; ++ic2) {
      for (int ic1 = 0; ic1 < Nc; ++ic1) {
        int jg2 = ic2 * 2 + ig;
        int jg1 = ic1 * Nc2 + ig;
        g[2 * ic2 + Nc2 * ic1 + ig] +=
          ff * mult_Gnn_r(&g1[jg1], &g2[jg2], Nc);
        g[2 * ic2 + 1 + Nc2 * ic1 + ig] +=
          ff * mult_Gnn_i(&g1[jg1], &g2[jg2], Nc);
      }
    }
  }
}


//====================================================================
void multadd_Field_Gdn(Field_G& W, const int ex,
                       const Field_G& U1, const int ex1,
                       const Field_G& U2, const int ex2,
                       const double ff)
{
  int Nvol = W.nvol();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  assert(ex < W.nex());
  assert(ex1 < U1.nex());
  assert(ex2 < U2.nex());
  assert(U1.nvol() == W.nvol());
  assert(U2.nvol() == W.nvol());

  double       *g  = W.ptr(0, 0, ex);
  const double *g1 = U1.ptr(0, 0, ex1);
  const double *g2 = U2.ptr(0, 0, ex2);

  const int Nc  = W.nc();
  const int Nc2 = 2 * Nc;
  const int Ndf = Nc2 * Nc;

  for (int site = is; site < ns; ++site) {
    int ig = Ndf * site;

    for (int ic2 = 0; ic2 < Nc; ++ic2) {
      for (int ic1 = 0; ic1 < Nc; ++ic1) {
        int jg2 = ic2 * 2 + ig;
        int jg1 = ic1 * 2 + ig;
        g[2 * ic2 + Nc2 * ic1 + ig] +=
          ff * mult_Gdn_r(&g1[jg1], &g2[jg2], Nc);
        g[2 * ic2 + 1 + Nc2 * ic1 + ig] +=
          ff * mult_Gdn_i(&g1[jg1], &g2[jg2], Nc);
      }
    }
  }
}


//====================================================================
void multadd_Field_Gnd(Field_G& W, const int ex,
                       const Field_G& U1, const int ex1,
                       const Field_G& U2, const int ex2,
                       const double ff)
{
  int Nvol = W.nvol();

  assert(ex < W.nex());
  assert(ex1 < U1.nex());
  assert(ex2 < U2.nex());
  assert(U1.nvol() == W.nvol());
  assert(U2.nvol() == W.nvol());

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  double       *g  = W.ptr(0, 0, ex);
  const double *g1 = U1.ptr(0, 0, ex1);
  const double *g2 = U2.ptr(0, 0, ex2);

  const int Nc  = W.nc();
  const int Nc2 = 2 * Nc;
  const int Ndf = Nc2 * Nc;

  for (int site = is; site < ns; ++site) {
    int ig = Ndf * site;

    for (int ic2 = 0; ic2 < Nc; ++ic2) {
      for (int ic1 = 0; ic1 < Nc; ++ic1) {
        int jg2 = ic2 * Nc2 + ig;
        int jg1 = ic1 * Nc2 + ig;
        g[2 * ic2 + Nc2 * ic1 + ig] +=
          ff * mult_Gnd_r(&g1[jg1], &g2[jg2], Nc);
        g[2 * ic2 + 1 + Nc2 * ic1 + ig]
          += ff * mult_Gnd_i(&g1[jg1], &g2[jg2], Nc);
      }
    }
  }
}


//====================================================================
void multadd_Field_Gdd(Field_G& W, const int ex,
                       const Field_G& U1, const int ex1,
                       const Field_G& U2, const int ex2,
                       const double ff)
{
  int Nvol = W.nvol();

  assert(ex < W.nex());
  assert(ex1 < U1.nex());
  assert(ex2 < U2.nex());
  assert(U1.nvol() == W.nvol());
  assert(U2.nvol() == W.nvol());

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  double       *g  = W.ptr(0, 0, ex);
  const double *g1 = U1.ptr(0, 0, ex1);
  const double *g2 = U2.ptr(0, 0, ex2);

  const int Nc  = W.nc();
  const int Nc2 = 2 * Nc;
  const int Ndf = Nc2 * Nc;

  for (int site = is; site < ns; ++site) {
    int ig = Ndf * site;

    for (int ic2 = 0; ic2 < Nc; ++ic2) {
      for (int ic1 = 0; ic1 < Nc; ++ic1) {
        int jg2 = ic2 * Nc2 + ig;
        int jg1 = ic1 * 2 + ig;
        g[2 * ic2 + Nc2 * ic1 + ig] +=
          ff * mult_Gdd_r(&g1[jg1], &g2[jg2], Nc);
        g[2 * ic2 + 1 + Nc2 * ic1 + ig] +=
          ff * mult_Gdd_i(&g1[jg1], &g2[jg2], Nc);
      }
    }
  }
}


//====================================================================
void at_Field_G(Field_G& W, const int ex)
{
  int Nvol = W.nvol();

  assert(ex < W.nex());

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  double *g = W.ptr(0, 0, ex);

  const int Nc  = W.nc();
  const int Ndf = 2 * Nc * Nc;

  for (int site = is; site < ns; ++site) {
    int ig = Ndf * site;

    for (int a = 0; a < Nc; ++a) {
      for (int b = a + 1; b < Nc; ++b) {
        double re = g[2 * (Nc * a + b) + ig] - g[2 * (Nc * b + a) + ig];
        double im = g[2 * (Nc * a + b) + 1 + ig] + g[2 * (Nc * b + a) + 1 + ig];

        g[2 * (Nc * a + b) + ig]     = 0.5 * re;
        g[2 * (Nc * a + b) + 1 + ig] = 0.5 * im;

        g[2 * (Nc * b + a) + ig]     = -0.5 * re;
        g[2 * (Nc * b + a) + 1 + ig] = 0.5 * im;
      }
    }
    double tr = 0.0;
    for (int cc = 0; cc < Nc; ++cc) {
      tr += g[2 * (Nc * cc + cc) + 1 + ig];
    }
    tr = tr / Nc;
    for (int cc = 0; cc < Nc; ++cc) {
      g[2 * (Nc * cc + cc) + ig]      = 0.0;
      g[2 * (Nc * cc + cc) + 1 + ig] -= tr;
    }
  }
}


//====================================================================
void ah_Field_G(Field_G& W, const int ex)
{
  int Nvol = W.nvol();

  assert(ex < W.nex());

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  double *g = W.ptr(0, 0, ex);

  const int Nc  = W.nc();
  const int Ndf = 2 * Nc * Nc;

  for (int site = is; site < ns; ++site) {
    int ig = Ndf * site;

    for (int a = 0; a < Nc; ++a) {
      for (int b = a; b < Nc; ++b) {
        double re = g[2 * (Nc * a + b) + ig] - g[2 * (Nc * b + a) + ig];
        double im = g[2 * (Nc * a + b) + 1 + ig] + g[2 * (Nc * b + a) + 1 + ig];

        g[2 * (Nc * a + b) + ig]     = 0.5 * re;
        g[2 * (Nc * a + b) + 1 + ig] = 0.5 * im;

        g[2 * (Nc * b + a) + ig]     = -0.5 * re;
        g[2 * (Nc * b + a) + 1 + ig] = 0.5 * im;
      }
    }
  }
}


//====================================================================
void mult_exp_Field_G(Field_G& W,
                      const double alpha, const Field_G& iP,
                      const Field_G& U, const int Nprec)
{
  // W = exp(alpha * iP) * U
  //   = (U + alpha * iP * (U + alpha/2 * iP * (...(U + alpha/n * iP * U)...)

  const int Nvol = U.nvol();
  const int Nex  = U.nex();
  const int Nc   = CommonParameters::Nc();

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  //- not imp version, yet.

  Mat_SU_N u0(Nc), v0(Nc), w0(Nc);

  for (int ex = 0; ex < Nex; ++ex) {
    for (int site = is; site < ns; ++site) {
      u0 = U.mat(site, ex);
      v0 = iP.mat(site, ex);
      w0 = SU_N::mat_exp(alpha, v0, u0, Nprec);
      W.set_mat(site, ex, w0);
    }
  }
}


//============================================================END=====
