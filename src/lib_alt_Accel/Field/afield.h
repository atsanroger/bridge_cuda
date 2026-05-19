/*!
      @file    afield.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2023-08-20 14:25:12 #$
      @version $LastChangedRevision: 2535 $
*/

#ifndef ACCEL_AFIELD_INCLUDED
#define ACCEL_AFIELD_INCLUDED

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

#include  "lib_alt_Accel/inline/define_params.h"

template<typename REALTYPE>
class AField<REALTYPE, ACCEL>
{

public:
    typedef REALTYPE real_t;
    typedef typename ComplexTraits<REALTYPE>::complex_t complex_t;
    typedef Element_type::type element_type;  // see bridge_defs.h

    static const Impl IMPL = ACCEL;
    static const std::string class_name;

    static real_t* m_red1;
    static real_t* m_red2;
    static int m_num_instance;

protected:
    int m_nin;
    int m_nvol;
    int m_nex;
    element_type m_element_type;
    Bridge::VerboseLevel m_vl;

    int m_nvol_pad;           //!< after padding.
    std::size_t m_nsize;
    std::size_t m_nsize_pad;  //!< after padding.

    real_t* m_field;

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
    AField(const AField<real_t, ACCEL>& w){
      init(w.nin(), w.nvol(), w.nex(), w.field_element_type());
      if(m_nsize > 0) copy(w);
    }

    //! copy constructor
    AField(const Field& w){
      init(w.nin(), w.nvol(), w.nex(), w.field_element_type());
      if(m_nsize > 0) copy(w);
    }

    //! operator overload
    AField& operator=(const Field& w) {
      assert(check_size(w.nin(), w.nvol(), w.nex()));
      if(m_nsize > 0) copy(w);
      return *this;
    }

    //! operator overload
    AField& operator=(const AField& w) {
      assert(check_size(w));
      if(m_nsize > 0) copy(w);
      return *this;
    }

    //! destructor
    ~AField(){ tidyup(); };

 private:

    //! initial setup
    void init(const int nin, const int nvol, const int nex,
              const element_type cmpl);

    //! final cleanup
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

    //! checking size parameters
    bool check_size(const int nin, const int nvol, const int nex) const
    {
      bool chk = true;
      if(m_nin != nin || m_nvol != nvol || m_nex != nex) chk = false;
      return chk;
    }

    //! checking size parameters
    bool check_size(const AField<REALTYPE, ACCEL>& w) const
    {
      bool chk = true;
      if(m_nvol != w.nvol() || m_nex != w.nex()) chk = false;
      if(m_nin != w.nin() && m_nin != 2 * w.nin() && m_nin * 2 != w.nin()) chk = false;
      return chk;
    }

    //! return the array size
    inline int size(void) const { return m_nsize; }

    //! return the address of data at given index.
    real_t* ptr(int i){ return &m_field[i]; }

    //complex_t* ptr_complex(int i){ return (complex_t*)&m_field[2*i]; }
    complex_t* ptr_complex(int i){ return (complex_t*)&m_field[i]; }

    //! reference of data element
    real_t cmp(const int index) const;

    void set(const int index, const real_t a);

    void set(const real_t a);

    // ###  BLAS-like methods  ### 

    void copy(const Field &w);

    void copy(const AField<real_t, ACCEL>& w);

    void copy(const int ex,
              const AField<real_t, ACCEL>& w, const int ex_w);

    void axpy(const real_t, const AField<real_t, ACCEL>&, const int use_qdw = 0);

    void axpy(const int ex, const real_t a, 
              const AField<real_t, ACCEL>& w, const int ex_w, const int use_qdw = 0);

    //! complex version of axpy: real and imaginary parts in order.
    void axpy(const real_t a_r, const real_t a_i,
              const AField<real_t, ACCEL>&);

    //! complex version of axpy: real and imaginary parts in order.
    void axpy(const int ex, const real_t a_r, const real_t a_i, 
              const AField<real_t, ACCEL>& w, const int ex_w);

    //! complex version of axpy: complex_t is set by ComplexTraits.
    void axpy(const int ex, const complex_t a,
              const AField<real_t, ACCEL>& w, const int ex_w);

    //! complex version of axpy: complex_t is set by ComplexTraits.
    void axpy(const complex_t a,
              const AField<real_t, ACCEL>& w);

    void aypx(const real_t, const AField<real_t, ACCEL>&, const int use_qdw = 0);

    void aypx(const int ex, const real_t,
              const AField<real_t, ACCEL>& w, const int ex_w, const int use_qdw = 0);

    //! complex version of aypx: real and imaginary parts in order.
    void aypx(const real_t a_r, const real_t a_i,
              const AField<real_t, ACCEL>&, const int use_qdw = 0);

    //! complex version of aypx: real and imaginary parts in order.
    void aypx(const int ex, const real_t a_r, const real_t a_i,
              const AField<real_t, ACCEL>& w, const int ex_w, const int use_qdw = 0);

    //! complex version of aypx: complex_t is set by ComplexTraits.
    void aypx(const complex_t a,
              const AField<real_t, ACCEL>& w, const int use_qdw = 0);

    //! complex version of aypx: complex_t is set by ComplexTraits.
    void aypx(const int ex, const complex_t a,
              const AField<real_t, ACCEL>& w, const int ex_w, const int use_qdw = 0);

    void scal(const real_t);

    void scal(const real_t, const int ex);

    void scal(const complex_t);

    void scal(const complex_t, const int ex);

    real_t dot(const AField<real_t, ACCEL>&, const int use_qdw = 0);

    //! complex inner-product: real and imaginary parts in order.
    void dotc(real_t&, real_t&, const AField<real_t, ACCEL>&, const int use_qdw = 0) const;

    //! square norm squared (|v|^2).
    real_t norm2(const int use_qdw = 0) const;

    //! square norm squared (|v|^2).
    real_t norm2_host(void) const;

    //! normalize QDW/QTW drift
    void normalize(const int use_qdw = 0);

    //! multiplying imaginary unit.
    void xI();

    //! taking complex conjugate.
    void conjg();

    //###  Methods specific to Accel implementation  ###

    void set_host(const real_t a);

    real_t cmp_host(const int index) const
    { return m_field[index]; }

    void set_host(const int index, const real_t a)
    { m_field[index] = a; }

    void update_host() const;

    void update_device();

    //! returning volume size after padding.
    int nvol_pad() const { return m_nvol_pad; }

    //! returning volume size after padding.
    int size_pad() const { return m_nsize_pad; }

};

// copy between fields with different precisions.

void copy(AField<double, ACCEL>&, const AField<float, ACCEL>&);
void copy(AField<float, ACCEL>&, const AField<double, ACCEL>&);


#endif // ACCEL_AFIELD_INCLUDED
