/*!
        @file    aindex_eo-inc.h
        @brief
        @author  <Hideo Matsufuru> hideo.matsufuru@kek.jp(matsufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2024-01-10 20:11:05 #$
        @version $LastChangedRevision: 2567 $
*/

#include "lib_alt_Accel/Field/aindex_eo.h"

#include "lib_alt_Accel/inline/define_index.h"
#include "lib_alt_Accel/BridgeACC/bridgeACC_Index_eo_alt.h"

#include <assert.h>

//====================================================================
template<typename REALTYPE>
template <typename AFIELD>
void AIndex_eo<REALTYPE,ACCEL>::split(AFIELD& field_e,
                                      AFIELD& field_o,
                                      const AFIELD& field_lex)
{
  typedef REALTYPE real_t;

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    int Nin = field_lex.nin();
    int Nex = field_lex.nex();
    int Nvol  = field_lex.nvol();
    int Nvol2 = Nvol/2;

    int Nvol_pad  = field_lex.nvol_pad();
    int Nvol2_pad = field_e.nvol_pad();
    assert(field_o.nvol_pad() == Nvol2_pad);

    assert(field_e.check_size(Nin, Nvol2, Nex));
    assert(field_o.check_size(Nin, Nvol2, Nex));

    int ieo_org = ieo_origin();

    for(int ex = 0; ex < Nex; ++ex){
      int idx_lex = IDX2(Nin, 0, Nvol_pad  * ex);
      int idx_eo  = IDX2(Nin, 0, Nvol2_pad * ex);

      real_t *w  = const_cast<AFIELD*>(&field_lex)->ptr(idx_lex);
      real_t *ve = field_e.ptr(idx_eo);
      real_t *vo = field_o.ptr(idx_eo);

      BridgeACC::split(ve, vo, w,
                       ieo_org, Nin, m_Nsize, Nvol2_pad, Nvol_pad);
    }

  }
#pragma omp barrier
}

//====================================================================
template<typename REALTYPE>
template <typename AFIELD>
void AIndex_eo<REALTYPE,ACCEL>::split(
                               AFIELD& field_e, const int ex_e,
                               AFIELD& field_o, const int ex_o,
                               const AFIELD& field_lex, const int ex)
{
  typedef REALTYPE real_t;

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    int Nin   = field_lex.nin();
    int Nvol  = field_lex.nvol();
    int Nvol2 = Nvol/2;

    int ieo_org   = ieo_origin();
    int Nvol_pad  = field_lex.nvol_pad();
    int Nvol2_pad = field_e.nvol_pad();
    assert(field_o.nvol_pad() == Nvol2_pad);

    AIndex_lex<REALTYPE,ACCEL> index_lex;
    int idx_lex = index_lex.idx(0, Nin, 0, ex);

    int idx_e = idxh(0, Nin, 0, ex_e);
    int idx_o = idxh(0, Nin, 0, ex_o);

    real_t *w  = const_cast<AFIELD*>(&field_lex)->ptr(idx_lex);
    real_t *ve = field_e.ptr(idx_e);
    real_t *vo = field_o.ptr(idx_o);

    BridgeACC::split(ve, vo, w,
                     ieo_org, Nin, m_Nsize, Nvol2_pad, Nvol_pad);
  }

#pragma omp barrier
}

//====================================================================
template<typename REALTYPE>
template <typename AFIELD>
void AIndex_eo<REALTYPE,ACCEL>::split_gauge(
                                    AFIELD& Ueo, const AFIELD& Ulex)
{
  typedef REALTYPE real_t;

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    int Ndf = NDF;
    int Ndim = NDIM;

    AIndex_lex<real_t,AFIELD::IMPL> index_lex;

    int nvol_pad  = Ulex.nvol_pad();
    int nvol2_pad = Ueo.nvol_pad();

    for(int mu = 0; mu < Ndim; ++mu){
      real_t *ue   = Ueo.ptr(idxh(0, Ndf, 0, 2*mu));
      real_t *uo   = Ueo.ptr(idxh(0, Ndf, 0, 2*mu+1));
      real_t *ulex = const_cast<AFIELD*>(&Ulex)->ptr(
                                       index_lex.idx(0, Ndf, 0, mu));

      BridgeACC::split(ue, uo, ulex, ieo_origin(),
                       Ndf, m_Nsize, nvol2_pad, nvol_pad);
    }
  }

#pragma omp barrier
}

//====================================================================
template<typename REALTYPE>
template <typename AFIELD>
void AIndex_eo<REALTYPE,ACCEL>::merge(AFIELD& field_lex,
                         const AFIELD& field_e, const AFIELD& field_o)
{
  typedef REALTYPE real_t;

#pragma omp barrier

  int ith = ThreadManager::get_thread_id();
  if (ith == 0) {
    int Nin = field_lex.nin();
    int Nex = field_lex.nex();
    int Nvol = field_lex.nvol();
    int Nvol2 = Nvol/2;

    int Nvol2_pad = field_e.nvol_pad();
    int Nvol_pad  = field_lex.nvol_pad();
    assert(field_o.nvol_pad() == Nvol2_pad);

    assert(field_e.check_size(Nin, Nvol2, Nex));
    assert(field_o.check_size(Nin, Nvol2, Nex));

    int ieo_org   = ieo_origin();

    for(int ex = 0; ex < Nex; ++ex){
      int idx_lex = IDX2(Nin, 0, Nvol_pad  * ex);
      int idx_eo  = IDX2(Nin, 0, Nvol2_pad * ex);

      real_t *v  = field_lex.ptr(idx_lex);
      real_t *we = const_cast<AFIELD*>(&field_e)->ptr(idx_eo);
      real_t *wo = const_cast<AFIELD*>(&field_o)->ptr(idx_eo);

      BridgeACC::merge(v, we, wo,
                       ieo_org, Nin, m_Nsize, Nvol2_pad, Nvol_pad);
    }
  }

#pragma omp barrier
}

//============================================================END=====
