/*!
      @file    afopr_Wilson.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef QXS_AFOPR_WILSON_INCLUDED
#define QXS_AFOPR_WILSON_INCLUDED

#include <cstdio>
#include <cstdlib>

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "lib/Fopr/afopr.h"
#include "lib/Parameters/commonParameters.h"
#include "lib/Communicator/communicator.h"
#include "lib/Communicator/communicator_impl.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/aindex_lex.h"

class Field;

/*!
  Implementation of the Wilson fermion operator (lexical index)
  in the QXS branch.
                                      [24 Dec 2022 H.Matsufuru]
*/
template<typename AFIELD>
class AFopr_Wilson : public AFopr<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 protected:
  int m_Nc, m_Nd, m_Nvc, m_Ndf, m_Ndim;
  int m_Nx, m_Ny, m_Nz, m_Nt, m_Nst;
  int m_Nxv, m_Nyv, m_Nstv;

  real_t m_CKs;                //!< hopping parameter.
  std::vector<int> m_boundary; //!< pointer to boundary condition
  std::string m_repr;          //!< gamma matrix representation

  Bridge::VerboseLevel m_vl;   //!< verbose level

  Field *m_conf;               //!< original gauge config.
  AFIELD m_U;                  //!< copied gauge config. with boundary conditions.

  std::string m_mode;          //!< mult mode

  AFIELD m_v2;

  int m_Nsize[4];  // lattice sizes (Nxv in x-direction)

  int do_comm[4];  // switchs of communication (4=Ndim): (0: n, 1: y).
  int do_comm_any; // switchs of communication (if any): (0: n, 1: y).

  std::vector<int> m_bdsize;
  using allocator_t = typename AFIELD::template aligned_allocator<char>;
  using Channel     = Channel_impl<allocator_t>;
  std::vector<Channel> chsend_up, chrecv_up, chsend_dn, chrecv_dn;
  ChannelSet chset_send, chset_recv;

 public:
  //! constructor.
  AFopr_Wilson(const Parameters& params) : AFopr<AFIELD>()
  { init(params); }

  //! destructor.
  ~AFopr_Wilson() { tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters by values.
  void set_parameters(real_t CKs, std::vector<int> bc);

  //! get parameters via a Parameter object
  void get_parameters(Parameters& params) const;

  //! setting gauge configuration (common interface).
  void set_config(Field *u);

  //! QXS version requires convert of spinor field.
  bool needs_convert() { return true; }

  //! convert of spinor field.
  void convert(AFIELD& v, const Field& w);

  //! reverse of spinor field.
  void reverse(Field& v, const AFIELD& w);

  //! setting mult mode.
  void set_mode(std::string mode);

  //! returns mult mode.
  std::string get_mode() const;

  void mult(AFIELD&, const AFIELD&);
  void mult_dag(AFIELD&, const AFIELD&);
  void mult_gm5(AFIELD&, const AFIELD&);

  void mult_up(int mu, AFIELD&, const AFIELD&);
  void mult_dn(int mu, AFIELD&, const AFIELD&);

  //! returns inner size parameter.
  int field_nin() { return 2 * m_Nc * m_Nd; }

  //! returns local volume size parameter.
  int field_nvol() { return m_Nst; }

  //! returns external size parameter.
  int field_nex() { return 1; }

  //! returns floating operation counts.
  double flop_count() { return flop_count(m_mode); }

  //! returns floating operation counts for given mode.
  double flop_count(const std::string mode);

 private:
  //! initial setup.
  void init(const Parameters& params);

  //! final tidy-up.
  void tidyup();

  //! setup channels for communication.
  void setup_channels();

  //! setting gauge configuration (setting omp parallel).
  void set_config_omp(Field *u);

  //! setting gauge configuration (implementation).
  void set_config_impl(Field *u);

  void DdagD(AFIELD&, const AFIELD&);
  void Ddag(AFIELD&, const AFIELD&);
  void H(AFIELD&, const AFIELD&);
  void D(AFIELD&, const AFIELD&);

  void mult_gm4(AFIELD&, const AFIELD&);

  //! standard D mult.
  void mult_D(AFIELD&, const AFIELD&);

  //! D mult using mult_xp, etc.
  void mult_D_alt(AFIELD&, const AFIELD&);

  void mult_xp(real_t *, real_t *);
  void mult_xm(real_t *, real_t *);
  void mult_yp(real_t *, real_t *);
  void mult_ym(real_t *, real_t *);
  void mult_zp(real_t *, real_t *);
  void mult_zm(real_t *, real_t *);
  void mult_tp(real_t *, real_t *);
  void mult_tm(real_t *, real_t *);

  void clear(real_t *);

  void aypx(real_t, real_t *, real_t *);
  void gm5_aypx(real_t, real_t *, real_t *);


#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object_with_params(const Parameters& params)
  { return new AFopr_Wilson(params); }

 public:
  static bool register_factory()
  {
    bool init1 = AFopr<AFIELD>::Factory_params::Register("Wilson",
                                                         create_object_with_params);
    return init1;
  }
#endif
};

#endif
