/*!
        @file    aforce_F_Smeared.cpp
        @brief
        @author  Hideo Matsufuru (matufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
        @version $LastChangedRevision: 2160 $
*/

#include "lib/Force/Fermion/aforce_F_Rational.h"

#include<cassert>

// include files in core library
#include "lib/ResourceManager/threadManager_OpenMP.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

// include files in alt-code dorectories
#include "lib_alt_SIMD2/Field/afield.h"
#include "lib_alt_SIMD2/Field/index_lex_alt.h"
#include "lib_alt_SIMD2/Field/afield-inc.h"

#include "lib/Force/Fermion/aforce_F_Rational-tmpl.h"

//====================================================================
// explicit instanciation for AField<double,SIMD2>.

template<>
const std::string AForce_F_Rational<AField<double,SIMD2> >::class_name
                               = "AForce_F_Rational<Afield<double,SIMD2> >";

//====================================================================
template<>
void AForce_F_Rational<AField<double,SIMD2> >::force_udiv_impl(
                  AField<double,SIMD2>& force, const AField<double,SIMD2>& eta)
{
  int Nvol = CommonParameters::Nvol();
  int Ndim = CommonParameters::Ndim();

  int NinF  = eta.nin();
  int NvolF = eta.nvol();
  int NexF  = eta.nex();

  //- Shiftsolver
  const int Nshift = m_Np;

  std::vector<AField<double,SIMD2> > psi(Nshift);
  for (int i = 0; i < Nshift; ++i) {
    psi[i].reset(NinF, NvolF, NexF);
  }

  vout.general(m_vl, "    Shift solver in force calculation\n");
  vout.general(m_vl, "      Number of shift values = %d\n", m_cl.size());
  m_fopr->set_mode("DdagD");

  AShiftsolver_CG<AField<double,SIMD2>, AFopr<AField<double,SIMD2> > >
     solver(m_fopr, m_Niter, m_Stop_cond);

  int    Nconv;
  double diff;
  solver.solve(psi, m_cl, eta, Nconv, diff);
  vout.general(m_vl, "      diff(max) = %22.15e  \n", diff);

  int NinG = force.nin();
  AField<double,SIMD2> force1(NinG, Nvol, Ndim);

  force.set(0.0);
  for (int i = 0; i < Nshift; ++i) {
    m_force->force_udiv(force1, psi[i]);
    scal(force1, m_bl[i]);
    axpy(force, 1.0, force1);
  }

}

//====================================================================
template<>
void AForce_F_Rational<AField<double,SIMD2> >::force_udiv(
                  AField<double,SIMD2>& force, const AField<double,SIMD2>& eta)
{
  force_udiv_impl(force, eta);

}

//====================================================================

template class AForce_F_Rational<AField<double,SIMD2> >;

//============================================================END=====
