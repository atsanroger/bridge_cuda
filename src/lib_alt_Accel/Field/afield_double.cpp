/*!
      @file    afield_double.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-08-20 14:25:12 #$
      @version $LastChangedRevision: 2535 $
*/

#include "lib_alt_Accel/Field/afield.h"

#include <cassert>

#include "lib/ResourceManager/threadManager.h"

#include "lib_alt_Accel/Field/aindex_lex.h"

typedef double real_t;

#include "lib_alt_Accel/inline/afield_th-inc.h"

// library
#include "lib_alt_Accel/BridgeACC/bridgeACC.h"
#include "lib_alt_Accel/BridgeACC/bridgeACC_AField.h"
#include "lib_alt_Accel/BridgeACC/bridgeACC_AField_Gauge.h"

// class template file
#include "lib_alt_Accel/Field/afield-tmpl.h"


template<>
const std::string AField<double, ACCEL>::class_name = "AField<double, ACCEL>";


// explicit instanciation.
template class AField<double, ACCEL>;

// copy between fields with different precisions.
//====================================================================
void copy(AField<double, ACCEL>& v, const AField<float, ACCEL> &w)
{
  assert(v.check_size(w.nin(), w.nvol(), w.nex()));

  float *wp = const_cast<AField<float, ACCEL>*>(&w)->ptr(0);
  BridgeACC::copy_from_device(wp, w.size_pad());

  int ith, nth, is, ns;
  set_threadtask(ith, nth, is, ns, v.size_pad());

  for (int i = is; i < ns; ++i) {
    v.set_host(i, double(w.cmp(i)) );
  }

  BridgeACC::copy_to_device(v.ptr(0), v.size_pad());

}

//====================================================================
/*
void axpy(AField<double, ACCEL>& v, const dcomplex& a, const AField<double, ACCEL>& w)
{
  double ar = real(a);
  double ai = imag(a);
  v.axpy(ar, ai, w);
}

//====================================================================
void aypx(const dcomplex& a, AField<double, ACCEL>& v, const AField<double, ACCEL>& w)
{
  double ar = real(a);
  double ai = imag(a);
  v.aypx(ar, ai, w);
}

//====================================================================
dcomplex dotc(const AField<double, ACCEL>& v, const AField<double, ACCEL>& w)
{
  double ar, ai;
  v.dotc(ar, ai, w);
  dcomplex a = cmplx(ar, ai);
  return a;
}
*/
//============================================================END=====
