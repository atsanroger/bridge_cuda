/*!
      @file    afield.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2021-02-14 15:05:06 #$
      @version $LastChangedRevision: 2160 $
*/

#ifndef AFIELD_SIMD2_INCLUDED
#define AFIELD_SIMD2_INCLUDED

#include <cstddef>
#include <cstdio>
#include <string>

#include  "lib_alt/Field/afield_base.h"  // primary template

#include  "bridge_defs.h"
#include  "bridge_complex.h"
#include  "complexTraits.h"
#include  "lib/Parameters/commonParameters.h"
#include  "lib/Communicator/communicator.h"
#include  "lib/Field/field.h"
#include  "lib/IO/bridgeIO.h"
using Bridge::vout;

template<typename REALTYPE>
class AField<REALTYPE, SIMD2>
{

public:
    typedef REALTYPE real_t;
    typedef typename ComplexTraits<REALTYPE>::complex_t complex_t;
    typedef Element_type::type element_type;  // see bridge_defs.h

    static const Impl IMPL = SIMD2;
    static const std::string class_name;

protected:
    int m_nin;
    int m_nvol;
    int m_nex;
    element_type m_element_type;
    Bridge::VerboseLevel m_vl;

    std::size_t m_nsize;
    real_t* m_field;

public:
    //! constructor without argument
    AField()
      : m_element_type(Element_type::COMPLEX)
    { init(0, 0, 0); }

    //! standard constructor
    AField(const int nin, const int nvol, const int nex,
           const element_type cmpl = Element_type::COMPLEX)
      : m_element_type(cmpl)
      { init(nin, nvol, nex); }

    //! copy constructor
    AField(const AField<real_t,SIMD2>& w)
      : m_element_type(w.field_element_type())
    {
      init(w.nin(), w.nvol(), w.nex());
      if(m_nsize > 0) copy(w);
    }

    //! destructor
    ~AField(){ tidyup(); }

 private:
    //! initializer
    void init(const int nin, const int nvol, const int nex);

    //! final clean-up
    void tidyup();

 public:

    // resetting object.
    void reset(const int nin, const int nvol, const int nex,
               const element_type cmpl = Element_type::COMPLEX);

    //! alias of reset (without element_type): use reset.
    void resize(const int nin, const int nvol, const int nex)
    { reset(nin, nvol, nex); }

    //! returning size of inner (on site) d.o.f.
    int nin() const { return m_nin; }

    //! returning size of site d.o.f.
    int nvol() const { return m_nvol; }

    //! returning size of extra d.o.f.
    int nex() const { return m_nex; }

    //! returning element_type (real or complex).
    element_type field_element_type() const { return m_element_type; }

    //! checking size parameters.
    bool check_size(const int nin, const int nvol, const int nex) const
    {
      bool chk = true;
      if(m_nin != nin || m_nvol != nvol || m_nex != nex) chk = false;
      return chk;
    }

    bool check_size(const AField<REALTYPE,SIMD2>& w) const
    {
      bool chk = true;
      if(m_nin != w.nin() || m_nvol != w.nvol() || m_nex != w.nex())
        chk = false;
      return chk;
    }

    //! reference of data element.  to be discarded.
    inline real_t& e(const int index)
    { return m_field[index]; }

    //! return the array size
    inline int size(void) const { return m_nsize; }

    //! return the array size
    inline int ntot(void) const { return m_nsize; }

    real_t* ptr(int i){ return &m_field[i]; }

    //! reference of data element
    real_t cmp(const int index) const
    { return m_field[index]; }

    void set(const int index, const real_t a)
    { m_field[index] = a; }

    void set(const real_t a);

    void copy(const AField<real_t,SIMD2>& w);

    void copy(const int ex,
              const AField<real_t,SIMD2>& w, const int ex_w);

    // real axpy
    void axpy(const real_t, const AField<real_t,SIMD2>&);

    // real axpy with ex indices
    void axpy(const int ex, const real_t a, 
              const AField<real_t,SIMD2>& w, const int ex_w);

    //! complex version of axpy: real and imaginary parts in order.
    void axpy(const real_t a_r, const real_t a_i,
              const AField<real_t,SIMD2>&);

    //! complex version of axpy: real and imaginary parts in order.
    void axpy(const int ex, const real_t a_r, const real_t a_i, 
              const AField<real_t,SIMD2>& w, const int ex_w);

    //! complex version of axpy: complex_t is set by ComplexTraits.
    void axpy(const complex_t a,
              const AField<real_t, SIMD2>& w);

    //! complex version of axpy: complex_t is set by ComplexTraits.
    void axpy(const int ex, const complex_t a,
              const AField<real_t, SIMD2>& w, const int ex_w);

    // real aypx (y = ay + x)
    void aypx(const real_t, const AField<real_t,SIMD2>&);

    // real aypx (y = ay + x) with ex indices
    void aypx(const int ex, const real_t,
              const AField<real_t,SIMD2>& w, const int ex_w);

    //! complex version of aypx: real and imaginary parts in order.
    void aypx(const real_t a_r, const real_t a_i,
              const AField<real_t,SIMD2>&);

    //! complex version of aypx: real and imaginary parts in order.
    void aypx(const int ex, const real_t a_r, const real_t a_i,
              const AField<real_t, SIMD2>& w, const int ex_w);

    //! complex version of aypx: complex_t is set by ComplexTraits.
    void aypx(const complex_t a,
              const AField<real_t, SIMD2>& w);

    //! complex version of aypx: complex_t is set by ComplexTraits.
    void aypx(const int ex, const complex_t a,
              const AField<real_t, SIMD2>& w, const int ex_w);


    void scal(const real_t);

    void scal(const real_t, const int ex);

    void scal(const complex_t);

    void scal(const complex_t, const int ex);

    real_t dot(const AField<real_t,SIMD2>&);

    //! complex inner-product: real and imaginary parts in order.
    void dotc(real_t&, real_t&, const AField<real_t,SIMD2>&) const;


    real_t norm2(void) const;

    void xI();

    void conjg();

};

#endif
