/*!
        @file    $Id:: index_eo_alt-inc.#$

        @brief

        @author  <Hideo Matsufuru> hideo.matsufuru@kek.jp(matsufuru)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2021-02-14 15:05:06 #$

        @version $LastChangedRevision: 2160 $
*/


#include "index_eo_alt.h"

#include <assert.h>

#include "afield.h"

//====================================================================
template<typename REALTYPE>
template <typename AFIELD>
void Index_eo_alt<REALTYPE,SIMD2>::split(AFIELD& field_e, AFIELD& field_o,
                         const AFIELD& field_lex)
{
  int Nin = field_lex.nin();
  int Nex = field_lex.nex();
  int Nvol = field_lex.nvol();
  int Nvol2 = Nvol/2;

  assert(field_e.check_size(Nin, Nvol2, Nex));
  assert(field_o.check_size(Nin, Nvol2, Nex));

  Index_lex_alt<REALTYPE,AFIELD::IMPL> index_lex;

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
          field_e.set(index2, field_lex.cmp(index1));
        }
      }else{
        for(int in = 0; in < Nin; ++in){
          int index1 = index_lex.idx(in, Nin, ist, ex);
          int index2 = idxh(in, Nin, ist2, ex);
          field_o.set(index2, field_lex.cmp(index1));
        }
      }
    }
  }

}

//====================================================================
template<typename REALTYPE>
template <typename AFIELD>
void Index_eo_alt<REALTYPE,SIMD2>::merge(AFIELD& field_lex,
                         const AFIELD& field_e, const AFIELD& field_o)
{
  int Nin = field_lex.nin();
  int Nex = field_lex.nex();
  int Nvol = field_lex.nvol();
  int Nvol2 = Nvol/2;

  assert(field_e.check_size(Nin, Nvol2, Nex));
  assert(field_o.check_size(Nin, Nvol2, Nex));

  Index_lex_alt<REALTYPE,AFIELD::IMPL> index_lex;

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
          field_lex.set(index1, field_e.cmp(index2));
        }
      }else{
        for(int in = 0; in < Nin; ++in){
          int index1 = index_lex.idx(in, Nin, ist, ex);
          int index2 = idxh(in, Nin, ist2, ex);
          field_lex.set(index1, field_o.cmp(index2));
        }
      }
    }
  }

}

//====================================================================
//============================================================END=====
