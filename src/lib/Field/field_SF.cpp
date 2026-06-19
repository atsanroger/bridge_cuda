/*!
        @file    field_SF.cpp

        @brief

        @author  "Taniguchi Yusuke"  (tanigchi)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#include "field_SF.h"
#include "lib/Field/field_thread-inc.h"


//free function version

namespace Field_SF {
//====================================================================

/*!
  Set the boundary spatial link \f$U_k(t=0)\f$ to its proper Dirichlet bounday \f$W_k\f$.
  <ul>
  <li>Supposed to be used for u_mu given by u_mu.setpart_ex(0, *g,mu) when mu is a spatial.
  <li>Lorentz index is set to mn=0.
  <li>Field is given as an argument.
  <li>Introduced by Yusuke Taniguchi for SF bc.
  </ul>
*/
  void set_boundary_wk(Field_G& u, const Mat_SU_N& wk)
  {
    assert(u.nex() == 1);

    const int mn   = 0;
    int       Nt   = CommonParameters::Nt();
    int       Nvol = u.nvol();
    int       Svol = Nvol / Nt;

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Svol);

    if (Communicator::ipe(3) == 0) {
      for (int site = is; site < ns; ++site) {
        u.set_mat(site, mn, wk);
      }
    }
  }


//====================================================================

/*!
  Set the boundary spatial link \f$U_k(t=N_t-1)\f$ to its proper Dirichlet bounday \f$W_k'\f$.
  <ul>
  <li>Supposed to be used for v given by shift.backward(v,u_nu,mu) when mu is a temporal.
  <li>Lorentz index is set to mn=0.
  <li>Field is given as an argument.
  <li>Introduced by Yusuke Taniguchi for SF bc.
  </ul>
*/
  void set_boundary_wkpr(Field_G& u, const Mat_SU_N& wkpr)
  {
    assert(u.nex() == 1);

    const int mn   = 0;
    int       Nt   = CommonParameters::Nt();
    int       NPEt = CommonParameters::NPEt();
    int       Nvol = u.nvol();
    int       Svol = Nvol / Nt;

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Svol);

    if (Communicator::ipe(3) == NPEt - 1) {
      int site0 = Nvol - Svol;
      for (int site2 = is; site2 < ns; ++site2) {
        int site = site2 + site0;
        u.set_mat(site, mn, wkpr);
      }
    }
  }


//====================================================================

/*!
  Set the boundary matrix to zero.
  <ul>
  <li>Supposed to be used for a staple for the boundary spatial plaquette in order to keep the corresponding force to zero.
  <li>Lorentz index is set to mn=0.
  <li>Introduced by Yusuke Taniguchi for SF bc.
  </ul>
*/
  void set_boundary_zero(Field_G& u)
  {
    const int mn   = 0;
    int       Nc   = CommonParameters::Nc();
    int       Nt   = CommonParameters::Nt();
    int       Nvol = u.nvol();
    int       Svol = Nvol / Nt;

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Svol);

    if (Communicator::ipe(3) == 0) {
      for (int site = is; site < ns; ++site) {
        for (int cc = 0; cc < 2 * Nc * Nc; ++cc) {
          u.set(cc, site, mn, 0.0);
        }
      }
    }
  }


//====================================================================
  void set_boundary_zero(Field& f)
  {
#pragma omp barrier

    int Nt   = CommonParameters::Nt();
    int Nvol = f.nvol();
    int Svol = Nvol / Nt;
    int Nin  = f.nin();

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Svol);

    if (Communicator::ipe(3) == 0) {
      for (int site = is; site < ns; ++site) {
        for (int in = 0; in < Nin; ++in) {
          f.set(in, site, 0, 0.0);
        }
      }
    }
#pragma omp barrier
  }


//====================================================================

/*!
  Set the boundary spatial link to zero.
  <ul>
  <li>Supposed to be used for a force for the boundary spatial link in order to keep it zero.
  <li>Lorentz index mn runs for spatial value.
  <li>Introduced by Yusuke Taniguchi for SF bc.
  </ul>
*/
  void set_boundary_spatial_link_zero(Field_G& u)
  {
    int Nc   = CommonParameters::Nc();
    int Nt   = CommonParameters::Nt();
    int Nvol = u.nvol();
    int Svol = Nvol / Nt;

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Svol);

    if (Communicator::ipe(3) == 0) {
      for (int mn = 0; mn < 3; ++mn) {
        for (int site = is; site < ns; ++site) {
          for (int cc = 0; cc < 2 * Nc * Nc; ++cc) {
            u.set(cc, site, mn, 0.0);
          }
        }
      }
    }
  }


//====================================================================

/*!
  Multiply the boundary improvement factor ct or ctr to an SU(N) matrix object attached to the boundary t=0 or t=Nt.
  <ul>
  <li>Supposed to be used for a matrix given by upper(g,mu,nu) or lower(g,mu,nu) when it is attached to the t=0 boundary.
  <li>The SU(N) matrix object may belong to a site at t=0 or t=1 or t=Nt.
  <li>This function performs no check if the object is a proper one at the boundary.
  <li>Lorentz index is set to mn=0.
  <li>Introduced by Yusuke Taniguchi for SF bc.
  </ul>
*/
  void mult_ct_boundary(Field_G& u, const int t, const double ct)
  {
    int Nc   = CommonParameters::Nc();
    int Nx   = CommonParameters::Nx();
    int Ny   = CommonParameters::Ny();
    int Nz   = CommonParameters::Nz();
    int Nt   = CommonParameters::Nt();
    int NPEt = CommonParameters::NPEt();
    int Svol = Nx * Ny * Nz;

    const int mn  = 0;
    const int ini = Nx * Ny * Nz * t;
    // const int fin = Nx * Ny * Nz * (t + 1);

    int ith, nth, is, ns;
    set_threadtask(ith, nth, is, ns, Svol);

    if (((Communicator::ipe(3) == 0) && ((t == 0) || (t == 1))) ||
        ((Communicator::ipe(3) == NPEt - 1) && (t == Nt - 1))) {
      for (int site2 = is; site2 < ns; ++site2) {
        int site = site2 + ini;
        for (int cc = 0; cc < 2 * Nc * Nc; ++cc) {
          double ut = u.cmp(cc, site, mn);
          u.set(cc, site, mn, ut * ct);
        }
      }
    }
  }


//====================================================================
  void set_boundary_matrix(Mat_SU_N& wk, const std::vector<double>& phi)
  {
    const int    Lx     = CommonParameters::Lx();
    const double Lx_inv = 1.0 / double(Lx);

    double c0r = cos(phi[0] * Lx_inv);
    double c0i = sin(phi[0] * Lx_inv);
    double c1r = cos(phi[1] * Lx_inv);
    double c1i = sin(phi[1] * Lx_inv);
    double c2r = cos(phi[2] * Lx_inv);
    double c2i = sin(phi[2] * Lx_inv);

    wk.zero();
    wk.set(0, 0, c0r, c0i);
    wk.set(1, 1, c1r, c1i);
    wk.set(2, 2, c2r, c2i);
  }
} // namespace
//============================================================END=====
