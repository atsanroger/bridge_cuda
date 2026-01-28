/*!
      @file    define_index.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2543 $
*/

#ifndef DEFINE_INDEX_INCLUDED
#define DEFINE_INDEX_INCLUDED

// Index functions for lib_alt_Vector implementation.
// Note that this file must be included after define_params_SU3.h.
//                                         [21 Jul 2019 H.Matsufuru]

// the following definitions explicitly assume the SU(3) gauge group.

#define  ID1     0
#define  ID2     6
#define  ID3    12
#define  ID4    18

// index of corelib
#define IDX_CORE(nin, in, ist)  (in + nin * ist)

// general index
#define IDXV(nin, nvol, in, ist, ex)  (ist + (nvol+VPAD)*(in + nin * ex))

// general index
#define IDXV_BARE(nin, nvol, in, ist, ex)  (ist + nvol*(in + nin * ex))

// for 4 component spinor (Wilson)
#define IDXV_SP_R(nvol, ic, id, ist)   (ist + (nvol+VPAD)*(  2*ic + NVC * id))
#define IDXV_SP_I(nvol, ic, id, ist)   (ist + (nvol+VPAD)*(1+2*ic + NVC * id))
#define IDXV_SP( nvol, ivc, id, ist)   (ist + (nvol+VPAD)*(   ivc + NVC * id))

// for 2 component spinor (communication in Wilson)
#define IDXV_2SP_R(nvol, ic, id, ist)  (ist + nvol*(  2*ic + NVC * id))
#define IDXV_2SP_I(nvol, ic, id, ist)  (ist + nvol*(1+2*ic + NVC * id))
#define IDXV_2SP( nvol, ivc, id, ist)  (ist + nvol*(   ivc + NVC * id))

// for 1 component spinor (staggered)
#define IDXV_1SP_R(nvol, ic, ist)   (ist + (nvol+VPAD)*(  2*ic))
#define IDXV_1SP_I(nvol, ic, ist)   (ist + (nvol+VPAD)*(1+2*ic))
#define IDXV_1SP( nvol, ivc, ist)   (ist + (nvol+VPAD)*ivc)

// for gauge field
#define IDXV_G_R(nvol, ic1, ic2, ist, ex)  (ist + (nvol+VPAD)*(  2*ic1 + NVC*(ic2 + NC * ex)))
#define IDXV_G_I(nvol, ic1, ic2, ist, ex)  (ist + (nvol+VPAD)*(1+2*ic1 + NVC*(ic2 + NC * ex)))


#endif
