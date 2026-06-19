/*!
      @file    afopr_Staggered_eo.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef QXS_AFOPR_STAGGERED_EO_INCLUDED
#define QXS_AFOPR_STAGGERED_EO_INCLUDED

#include <cstdio>
#include <cstdlib>

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "lib/Fopr/afopr.h"
#include "lib/Communicator/communicator.h"
#include "lib/Communicator/communicator_impl.h"
#include "lib/Parameters/commonParameters.h"
#include "lib/Parameters/parameters.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt_QXS/Field/aindex_eo.h"
#include "lib_alt_QXS/Field/shiftAField_eo.h"

//! Staggered_eo fermion operator.

/*!
  Implementation of the standard staggered fermion operator
  (even-odd site index) for the A64FX architecture.
                                        [26 Sep 2021 H.Matsufuru]
 */

template<typename AFIELD>
class AFopr_Staggered_eo : public AFopr<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 private:
  int m_Nc, m_Nvc, m_Ndf, m_Ndim;
  int m_Nx, m_Ny, m_Nz, m_Nt, m_Nst;
  int m_Nx2, m_Nst2;
  int m_Nx2v, m_Nyv, m_Nst2v;  //!< for SIMD arrays

  real_t m_mq;                 //!< quark mass.
  std::vector<int> m_boundary; //!< boundary conditions.
  Bridge::VerboseLevel m_vl;   //!< verbose level

  std::vector<int> m_Leo;      //!< Leo = 0 (even site) or 1 (odd site).

  AFIELD m_stg_phase;          //!< staggered phase
  AFIELD m_parity;             //!< site parity for multiplying gamma_5

  AFIELD m_Ueo;                //!< gauge field multiplied by staggered phase
  AFIELD m_Ulex;               //!< converted lexical gauge field

  ShiftAField_eo<AFIELD> *m_shift;

  std::string m_mode;

  AFIELD m_w1, m_w2; //!< working vectors
  AFIELD m_v1, m_v2; //!< working vectors

  int do_comm[4];    //!< switchs of communication (4=Ndim): (0: n, 1: y).
  int do_comm_any;   //!< switchs of communication (if any): (0: n, 1: y).

  std::vector<int> m_Nbdsize;
  using allocator_t = typename AFIELD::template aligned_allocator<char>;
  using Channel     = Channel_impl<allocator_t>;
  std::vector<Channel> chsend_up, chrecv_up, chsend_dn, chrecv_dn;
  ChannelSet chset_send, chset_recv;

  int m_Nsize[4]; //!< lattice size in units of SIMD vector

 public:
  //! constructor.
  AFopr_Staggered_eo(const Parameters& params) : AFopr<AFIELD>()
  { init(params); }

  //! destructor.
  ~AFopr_Staggered_eo() { tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters by values.
  void set_parameters(const real_t mq, const std::vector<int> bc);

  //! getting parameters via a Parameters object.
  void get_parameters(Parameters& params) const;

  //! setting gauge configuration (common interface).
  void set_config(Field *U);

  bool needs_convert() { return false; }

  void set_mode(std::string mode);

  std::string get_mode() const { return m_mode; }

  void mult(AFIELD&, const AFIELD&);
  void mult_dag(AFIELD&, const AFIELD&);
  void mult_gm5(AFIELD&, const AFIELD&);

  void mult(AFIELD&, const AFIELD&, const std::string mode);
  void mult_dag(AFIELD&, const AFIELD&, const std::string mode);

  void mult_gm5(AFIELD&);

  void normalize_fprop(AFIELD& v);

  void normalize_fopr(AFIELD& v);

  int field_nvol() { return m_Nst2; }
  int field_nin()  { return m_Nvc; }
  int field_nex()  { return 1; }

  //! returns floating operation counts.
  double flop_count();

  //! returns floating operation counts.
  double flop_count(const std::string mode);

 private:
  void init(const Parameters& params);
  void tidyup();

  void set_staggered_phase();

  void setup_channels();

  //! setting gauge configuration (setting omp parallel).
  void set_config_omp(Field *u);

  //! setting gauge configuration (implementation).
  void set_config_impl(Field *u);

  void D(AFIELD&, const AFIELD&);
  void Ddag(AFIELD&, const AFIELD&);
  void DdagD(AFIELD&, const AFIELD&);
  void H(AFIELD&, const AFIELD&);

  void Meo(AFIELD&, const AFIELD&, int ieo);
  void Meo(AFIELD&, const AFIELD&, const AFIELD&,
           int ieo, int iflag);

  void mult_Meo_qxs(AFIELD&, const AFIELD&, const AFIELD&,
                    int ieo, int iflag);
  void mult_Meo_alt(AFIELD&, const AFIELD&, const AFIELD&,
                    int ieo, int iflag);

  void clear(real_t *);

  void axpby(real_t, real_t *, real_t, real_t *);

  void mult_up(int mu, AFIELD&, const AFIELD&);
  void mult_dn(int mu, AFIELD&, const AFIELD&);

  void mult_up(int mu, AFIELD&, const AFIELD&, int ieo);
  void mult_dn(int mu, AFIELD&, const AFIELD&, int ieo);


#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object_with_params(const Parameters& params)
  { return new AFopr_Staggered_eo(params); }

 public:
  static bool register_factory()
  {
    bool init1 = AFopr<AFIELD>::Factory_params::Register("Staggered_eo",
                                                         create_object_with_params);
    return init1;
  }
#endif
};
#endif
