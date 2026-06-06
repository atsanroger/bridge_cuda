/*!
        @file    aindex_eo-inc.h
        @brief
        @author  <Hideo Matsufuru> hideo.matsufuru@kek.jp(matsufuru)
                 $LastChangedBy: matufuru $
        @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
        @version $LastChangedRevision: 2543 $
*/

#include "lib_alt_Vector/Field/aindex_eo.h"
#include "lib_alt_Vector/BridgeVec/bridgeVec_aindex_eo.h"

#include <assert.h>

//====================================================================
template<typename REALTYPE>
template <typename AFIELD>
void AIndex_eo<REALTYPE,VECTOR>::split(AFIELD& field_e,
                                           AFIELD& field_o,
                                           const AFIELD& field_lex,
                                           bool /*qdw*/)
{
  typedef REALTYPE real_t;

  int Nin = field_lex.nin();
  int Nex = field_lex.nex();
  int Nvol = field_lex.nvol();
  int Nvol2 = Nvol/2;

  assert(field_e.check_size(Nin, Nvol2, Nex));
  assert(field_o.check_size(Nin, Nvol2, Nex));

  real_t *w  = const_cast<AFIELD*>(&field_lex)->ptr(0);
  real_t *ve = field_e.ptr(0);
  real_t *vo = field_o.ptr(0);

  int ieo_org   = ieo_origin();
  int nvol2_pad = field_e.nvol_pad();
  assert(field_o.nvol_pad() == nvol2_pad);
  int nvol_pad  = field_lex.nvol_pad();

  BridgeVec::split_dev(ve, vo, w,
                       ieo_org, Nin, m_Nsize, nvol2_pad, nvol_pad);

  /*
  AIndex_lex<REALTYPE,VECTOR> index_lex;

  real_t *fp = const_cast<AFIELD*>(&field_lex)->ptr(0);
  int size = Nin * Nvol * Nex;
  copy_from_device(fp, size);

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for(int ex = 0; ex < Nex; ++ex){
    for(int ist = is; ist < ns; ++ist){
      int ix = ist % Nx;
      int iyzt = ist/Nx;
      int ist2 = ist/2;
      int ieo = (ix + Leo[iyzt]) % 2;
      if(ieo == 0){
        for(int in = 0; in < Nin; ++in){
          int index1 = index_lex.idx(in, Nin, ist, ex);
          int index2 = idxh(in, Nin, ist2, ex);
          field_e.set_host(index2, field_lex.cmp_host(index1));
        }
      }else{
        for(int in = 0; in < Nin; ++in){
          int index1 = index_lex.idx(in, Nin, ist, ex);
          int index2 = idxh(in, Nin, ist2, ex);
          field_o.set_host(index2, field_lex.cmp_host(index1));
        }
      }
    }
  }

  real_t *fpe = field_e.ptr(0);
  real_t *fpo = field_o.ptr(0);
  int size_h = Nin * Nvol2 * Nex;
  copy_to_device(fpe, size_h);
  copy_to_device(fpo, size_h);
  */

}

//====================================================================
template<typename REALTYPE>
template <typename AFIELD>
void AIndex_eo<REALTYPE,VECTOR>::split_gauge(
                                    AFIELD& Ueo, const AFIELD& Ulex)
{
  typedef REALTYPE real_t;

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
    BridgeVec::split_dev(ue, uo, ulex, ieo_origin(),
                         Ndf, m_Nsize, nvol2_pad, nvol_pad);
  }

}

//====================================================================
template<typename REALTYPE>
template <typename AFIELD>
void AIndex_eo<REALTYPE,VECTOR>::merge(AFIELD& field_lex,
                         const AFIELD& field_e, const AFIELD& field_o,
                         bool /*qdw*/)
{
  typedef REALTYPE real_t;

  int Nin = field_lex.nin();
  int Nex = field_lex.nex();
  int Nvol = field_lex.nvol();
  int Nvol2 = Nvol/2;

  assert(field_e.check_size(Nin, Nvol2, Nex));
  assert(field_o.check_size(Nin, Nvol2, Nex));

  real_t *v  = field_lex.ptr(0);
  real_t *we = const_cast<AFIELD*>(&field_e)->ptr(0);
  real_t *wo = const_cast<AFIELD*>(&field_o)->ptr(0);

  int ieo_org   = ieo_origin();
  int nvol2_pad = field_e.nvol_pad();
  assert(field_o.nvol_pad() == nvol2_pad);
  int nvol_pad  = field_lex.nvol_pad();

  BridgeVec::merge_dev(v, we, wo, ieo_org, Nin, m_Nsize, nvol2_pad, nvol_pad);

  /*
  real_t *fpe = const_cast<AFIELD*>(&field_e)->ptr(0);
  real_t *fpo = const_cast<AFIELD*>(&field_o)->ptr(0);
  int size_h = Nin * Nvol2 * Nex;
  copy_from_device(fpe, size_h);
  copy_from_device(fpo, size_h);

  AIndex_lex<REALTYPE,VECTOR> index_lex;

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, Nvol);

  for(int ex = 0; ex < Nex; ++ex){
    for(int ist = is; ist < ns; ++ist){
      int ix = ist % Nx;
      int iyzt = ist/Nx;
      int ist2 = ist/2;
      int ieo = (ix + Leo[iyzt]) % 2;
      if(ieo == 0){
        for(int in = 0; in < Nin; ++in){
          int index1 = index_lex.idx(in, Nin, ist, ex);
          int index2 = idxh(in, Nin, ist2, ex);
          field_lex.set_host(index1, field_e.cmp_host(index2));
        }
      }else{
        for(int in = 0; in < Nin; ++in){
          int index1 = index_lex.idx(in, Nin, ist, ex);
          int index2 = idxh(in, Nin, ist2, ex);
          field_lex.set_host(index1, field_o.cmp_host(index2));
        }
      }
    }
  }

  real_t *fp = field_lex.ptr(0);
  int size = Nin * Nvol * Nex;
  copy_to_device(fp, size);
  */

}

//============================================================END=====
