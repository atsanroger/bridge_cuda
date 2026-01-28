/*!
      @file    afield_Gauge-inc.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
      @version $LastChangedRevision: 2543 $
*/

#ifndef AFIELD_GAUGE_INC_INCLUDED
#define AFIELD_GAUGE_INC_INCLUDED

//#include <cstdlib> 

#include "afield_Gauge_vector.h"

//====================================================================
template <typename REALTYPE>
void multadd_Gnn(AField<REALTYPE,VECTOR>& u, const int exu,
                 const AField<REALTYPE,VECTOR>& v, const int exv,
                 const AField<REALTYPE,VECTOR>& w, const int exw,
                 const REALTYPE a)
{
  typedef AField<REALTYPE,VECTOR> AFIELD;
  typedef REALTYPE real_t;

  int Nst = v.nvol();

  real_t* up = u.ptr(0);
  real_t* vp = const_cast<AFIELD*>(&v)->ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  multadd_Gnn(up, exu, vp, exv, wp, exw, a, Nst);

}

//====================================================================
template <typename REALTYPE>
void mult_Gnn(AField<REALTYPE,VECTOR>& u, const int exu,
              const AField<REALTYPE,VECTOR>& v, const int exv,
              const AField<REALTYPE,VECTOR>& w, const int exw)
{
  typedef AField<REALTYPE,VECTOR> AFIELD;
  typedef REALTYPE real_t;

  int Nst = v.nvol();

  real_t* up = u.ptr(0);
  real_t* vp = const_cast<AFIELD*>(&v)->ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  mult_Gnn(up, exu, vp, exv, wp, exw, Nst);

}

//====================================================================
template <typename REALTYPE>
void multadd_Gnd(AField<REALTYPE,VECTOR>& u, const int exu,
                 const AField<REALTYPE,VECTOR>& v, const int exv,
                 const AField<REALTYPE,VECTOR>& w, const int exw,
                 const REALTYPE a)
{
  typedef AField<REALTYPE,VECTOR> AFIELD;
  typedef REALTYPE real_t;

  int Nst = v.nvol();

  real_t* up = u.ptr(0);
  real_t* vp = const_cast<AFIELD*>(&v)->ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  multadd_Gnd(up, exu, vp, exv, wp, exw, a, Nst);

}

//====================================================================
template <typename REALTYPE>
void mult_Gnd(AField<REALTYPE,VECTOR>& u, const int exu,
              const AField<REALTYPE,VECTOR>& v, const int exv,
              const AField<REALTYPE,VECTOR>& w, const int exw)
{
  typedef AField<REALTYPE,VECTOR> AFIELD;
  typedef REALTYPE real_t;

  int Nst = v.nvol();

  real_t* up = u.ptr(0);
  real_t* vp = const_cast<AFIELD*>(&v)->ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  mult_Gnd(up, exu, vp, exv, wp, exw, Nst);

}

//====================================================================
template <typename REALTYPE>
void multadd_Gdn(AField<REALTYPE,VECTOR>& u, const int exu,
                 const AField<REALTYPE,VECTOR>& v, const int exv,
                 const AField<REALTYPE,VECTOR>& w, const int exw,
                 const REALTYPE a)
{
  typedef AField<REALTYPE,VECTOR> AFIELD;
  typedef REALTYPE real_t;

  int Nst = v.nvol();

  real_t* up = u.ptr(0);
  real_t* vp = const_cast<AFIELD*>(&v)->ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  multadd_Gdn(up, exu, vp, exv, wp, exw, a, Nst);

}

//====================================================================
template <typename REALTYPE>
void mult_Gdn(AField<REALTYPE,VECTOR>& u, const int exu,
              const AField<REALTYPE,VECTOR>& v, const int exv,
              const AField<REALTYPE,VECTOR>& w, const int exw)
{
  typedef AField<REALTYPE,VECTOR> AFIELD;
  typedef REALTYPE real_t;

  int Nst = v.nvol();

  real_t* up = u.ptr(0);
  real_t* vp = const_cast<AFIELD*>(&v)->ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  mult_Gdn(up, exu, vp, exv, wp, exw, Nst);

}

//====================================================================
template <typename REALTYPE>
void mult_Gdd(AField<REALTYPE,VECTOR>& u, const int exu,
              const AField<REALTYPE,VECTOR>& v, const int exv,
              const AField<REALTYPE,VECTOR>& w, const int exw)
{
  typedef AField<REALTYPE,VECTOR> AFIELD;
  typedef REALTYPE real_t;

  int Nst = v.nvol();

  real_t* up = u.ptr(0);
  real_t* vp = const_cast<AFIELD*>(&v)->ptr(0);
  real_t* wp = const_cast<AFIELD*>(&w)->ptr(0);

  mult_Gdd(up, exu, vp, exv, wp, exw, Nst);

}

//====================================================================
// anti-hermitian
template <typename REALTYPE>
void ah_G(AField<REALTYPE,VECTOR>& u, const int ex)
{
  typedef REALTYPE real_t;

  int Nst = u.nvol();

  real_t* up = u.ptr(0);

  ah_G(up, ex, Nst);

}

//====================================================================
// anti-hermitian traceless
template <typename REALTYPE>
void at_G(AField<REALTYPE,VECTOR>& u, const int ex)
{
  typedef REALTYPE real_t;

  int Nst = u.nvol();

  real_t* up = u.ptr(0);

  at_G(up, ex, Nst);

}

//====================================================================
template <typename REALTYPE>
void add_unit(AField<REALTYPE,VECTOR>& u, const int ex, REALTYPE a)
{         // u = u + a * I (I: unit matrix)

  typedef REALTYPE real_t;

  int Nst = u.nvol();

  real_t* up = u.ptr(0);

  add_unit(up, ex, a, Nst);

}

//============================================================END=====
#endif
