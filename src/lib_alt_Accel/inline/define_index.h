/*!
      @file    define_index.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2613 $
*/

#ifndef ACCEL_DEFINE_INDEX_INCLUDED
#define ACCEL_DEFINE_INDEX_INCLUDED

// Index functions for lib_alt_Accel implementation.
// Note that this file must be included after define_params_SU3.h.
//                                         [21 Jul 2019 H.Matsufuru]

// the following definitions explicitly assume the SU(3) gauge group.
#define  ID1     0
#define  ID2     6
#define  ID3    12
#define  ID4    18

#define  NDD     8   // assuming ND = 4

// index of corelib
#define IDX_CORE(nin, in, ist)  ((in) + (nin) * (ist))

// general index
#define IDX2(nin, in, ist)      (((ist)%NWP) + NWP*((in) + (nin)*((ist)/NWP)))

// for Wilson-type fermion
#define IDX2_SP_R(ic, id, ist)  (((ist)%NWP) + NWP*(2*(ic)   + NVC*((id) + ND*((ist)/NWP))))
#define IDX2_SP_I(ic, id, ist)  (((ist)%NWP) + NWP*(2*(ic)+1 + NVC*((id) + ND*((ist)/NWP))))
#define IDX2_SP(ivc, id, ist)   (((ist)%NWP) + NWP*((ivc)    + NVC*((id) + ND*((ist)/NWP))))

// for Domainwall-type fermion
#define IDX2_SP_5D_R(ic, id, is, Ns, ist)  (((ist)%NWP) + NWP*(2*(ic)   + NVC*((id) + ND*( (is) + (Ns) * ((ist)/NWP)))))
#define IDX2_SP_5D_I(ic, id, is, Ns, ist)  (((ist)%NWP) + NWP*(2*(ic)+1 + NVC*((id) + ND*( (is) + (Ns) * ((ist)/NWP)))))


// for communication buffer of Wilson-type fermion
#define IDXBF_R(ic, id, ist)   (((ist)%NWP) + NWP*(2*(ic)   + NVC*((id) + ND2*((ist)/NWP))))
#define IDXBF_I(ic, id, ist)   (((ist)%NWP) + NWP*(2*(ic)+1 + NVC*((id) + ND2*((ist)/NWP))))
#define IDXBF(ivc, id, ist)    (((ist)%NWP) + NWP*((ivc)    + NVC*((id) + ND2*((ist)/NWP))))

// for 1 component spinor (staggered)
#define IDX2_1SP_R(ic, ist)   (((ist)%NWP) + NWP*(  2*((ic) + NC*((ist)/NWP))))
#define IDX2_1SP_I(ic, ist)   (((ist)%NWP) + NWP*(1+2*((ic) + NC*((ist)/NWP))))
#define IDX2_1SP(ivc, ist)    (((ist)%NWP) + NWP*((ivc) + NVC*((ist)/NWP)))

// for gauge field
#define IDX2_G_R(ic1, ic2, ist)  (((ist)%NWP) + NWP*(  2*((ic1) + NC*((ic2) + NC*((ist)/NWP)))))
#define IDX2_G_I(ic1, ic2, ist)  (((ist)%NWP) + NWP*(1+2*((ic1) + NC*((ic2) + NC*((ist)/NWP)))))

// for clover term
#define IDX_CT_R(ic1, ic2, idd, ist)  (((ist)%NWP) + NWP*(  2*((ic1) + NC*((ic2) + NC*((idd) + NDD*((ist)/NWP))))))
#define IDX_CT_I(ic1, ic2, idd, ist)  (((ist)%NWP) + NWP*(1+2*((ic1) + NC*((ic2) + NC*((idd) + NDD*((ist)/NWP))))))

// for Domainwall-type fermion
#define IDX2_5D(nin, in, is, Ns, ist)  (((ist) % NWP) + NWP * ((in) + (nin) * (is + (Ns) * ((ist) / NWP))))

#define IDX2_SP_5D_R(ic, id, is, Ns, ist)  (((ist)%NWP) + NWP*(2*(ic)   + NVC*((id) + ND*( (is) + (Ns) * ((ist)/NWP)))))
#define IDX2_SP_5D_I(ic, id, is, Ns, ist)  (((ist)%NWP) + NWP*(2*(ic)+1 + NVC*((id) + ND*( (is) + (Ns) * ((ist)/NWP)))))

#endif
