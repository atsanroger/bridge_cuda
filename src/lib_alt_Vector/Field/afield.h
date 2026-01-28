/*!
      @file    afield.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-10-02 13:50:31 #$
      @version $LastChangedRevision: 2543 $
*/

#ifndef AFIELD_VECTOR_INCLUDED
#define AFIELD_VECTOR_INCLUDED

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

#include  "lib_alt_Vector/inline/define_params.h"

template<typename REALTYPE>
class AField<REALTYPE, VECTOR>
{

public:
    typedef REALTYPE real_t;
    typedef typename ComplexTraits<REALTYPE>::complex_t complex_t;
    typedef Element_type::type element_type;  // see bridge_defs.h

    static const Impl IMPL = VECTOR;
    static const std::string class_name;

protected:
    int m_nin;
    int m_nvol;
    int m_nex;
    element_type m_element_type;
    Bridge::VerboseLevel m_vl;

    real_t* m_field_base;
    real_t* m_field;

    int m_nvol_pad;           //!< after padding.
    std::size_t m_nsize;
    std::size_t m_nsize_pad;  //!< after padding.

public:
    //! constructor without argument
    AField()
    { init(0, 0, 0, Element_type::COMPLEX); }

    //! constructor
    AField(const int nin, const int nvol, const int nex,
           const element_type cmpl = Element_type::COMPLEX){
      init(nin, nvol, nex, cmpl);
    }

    //! copy constructor
    AField(const AField<real_t, VECTOR>& w){
      init(w.nin(), w.nvol(), w.nex(), w.field_element_type());
      copy(w);
    }

    //! copy constructor
    AField(const Field& w){
      init(w.nin(), w.nvol(), w.nex(), w.field_element_type());
      copy(w);
    }

    //! destructor
    ~AField(){ tidyup(); };

 private:

    void init(const int nin, const int nvol, const int nex,
              const element_type cmpl);

    void tidyup();

 public:

    // resetting object.
    void reset(const int nin, const int nvol, const int nex,
               const element_type cmpl = Element_type::COMPLEX);

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

    bool check_size(const AField<REALTYPE, VECTOR>& w) const
    {
      bool chk = true;
      if(m_nin != w.nin() || m_nvol != w.nvol() || m_nex != w.nex())
        chk = false;
      return chk;
    }

    //! return the array size
    inline int size(void) const { return m_nsize; }

    //! return the address of data at given index.
    real_t* ptr(int i){ return &m_field[i]; }

    //! reference of data element
    real_t cmp(const int index) const;

    void set(const int index, const real_t a);

    void set(const real_t a);

    //###  BLAS-like methods  ### 

    void clear();

    void copy(const Field &w);

    void copy(const AField<real_t, VECTOR>& w);

    void copy(const int ex,
              const AField<real_t, VECTOR>& w, const int ex_w);

    void axpy(const real_t, const AField<real_t, VECTOR>&);

    void axpy(const int ex, const real_t a, 
              const AField<real_t, VECTOR>& w, const int ex_w);

    //! complex version of axpy: real and imaginary parts in order.
    void axpy(const real_t a_r, const real_t a_i,
              const AField<real_t, VECTOR>&);

    //! complex version of axpy: real and imaginary parts in order.
    void axpy(const int ex, const real_t a_r, const real_t a_i, 
              const AField<real_t, VECTOR>& w, const int ex_w);

    //! complex version of axpy: complex_t is set by ComplexTraits.
    void axpy(const int ex, const complex_t a,
              const AField<real_t, VECTOR>& w, const int ex_w);

    //! complex version of axpy: complex_t is set by ComplexTraits.
    void axpy(const complex_t a,
              const AField<real_t, VECTOR>& w);

    void aypx(const real_t, const AField<real_t, VECTOR>&);

    void aypx(const int ex, const real_t,
              const AField<real_t, VECTOR>& w, const int ex_w);

    //! complex version of aypx: real and imaginary parts in order.
    void aypx(const real_t a_r, const real_t a_i,
              const AField<real_t, VECTOR>&);

    //! complex version of aypx: real and imaginary parts in order.
    void aypx(const int ex, const real_t a_r, const real_t a_i,
              const AField<real_t, VECTOR>& w, const int ex_w);

    //! complex version of aypx: complex_t is set by ComplexTraits.
    void aypx(const complex_t a,
              const AField<real_t, VECTOR>& w);

    //! complex version of aypx: complex_t is set by ComplexTraits.
    void aypx(const int ex, const complex_t a,
              const AField<real_t, VECTOR>& w, const int ex_w);

    void scal(const real_t);

    void scal(const real_t, const int ex);

    real_t dot(const AField<real_t, VECTOR>&);

    //! complex inner-product: real and imaginary parts in order.
    void dotc(real_t&, real_t&, const AField<real_t, VECTOR>&) const;

    real_t norm2(void) const;

    void xI();

    void conjg();

    //###  Methods specific to Vector implementation  ###

    void set_host(const real_t a);

    real_t cmp_host(const int index) const
    { return m_field[index]; }

    void set_host(const int index, const real_t a)
    { m_field[index] = a; }

    void update_host();

    void update_device();

    //! volume size ceiled with the multiple of NWP.
    //    int ceil_nwp(const int nst);

    //! returning volume size after padding.
    int nvol_pad() const { return m_nvol_pad; }

    //! returning volume size after padding.
    int size_pad() const { return m_nsize_pad; }

};

// copy between fields with different precisions.

void copy(AField<double, VECTOR>&, const AField<float, VECTOR>&);
void copy(AField<float, VECTOR>&, const AField<double, VECTOR>&);


#endif // AFIELD_INCLUDED
