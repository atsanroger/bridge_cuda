/*!
      @file    afopr_Clover_QWS_dd.h
      @brief
      @author  Issaku Kanamori (kanamori)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2499 $
*/

#ifndef QXS_AFOPR_CLOVER_QWS_DD_INCLUDED
#define QXS_AFOPR_CLOVER_QWS_DD_INCLUDED

#include <cstdio>
#include <cstdlib>

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "lib_alt/Fopr/afopr_dd.h"

#include "lib/Parameters/commonParameters.h"
#include "lib/Parameters/parameters.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib/Communicator/communicator.h"
#include "lib/Communicator/communicator_impl.h"

#include "lib/Fopr/fopr_CloverTerm.h"

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/aindex_lex.h"

class Field;

template<typename AFIELD>
class AFopr_Clover_QWS_dd : public AFopr_dd<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 protected:
  int m_Nc, m_Nd, m_Nvc, m_Ndf, m_Ndim;
  int m_Nx, m_Ny, m_Nz, m_Nt, m_Nst;
  int m_Nxv, m_Nyv, m_Nstv;

  real_t m_CKs;                       //!< hopping parameter.
  std::vector<int> m_boundary;        //!< pointer to boundary condition
  std::string m_repr;                 //!< gamma matrix representation

  Bridge::VerboseLevel m_vl;          //!< verbose level

  std::vector<int> m_block_size;      //!< block size

  int m_Ieo;                          //!< even-odd label of origin in units of block

  Field *m_conf;                      //!< original gauge config.
  AFIELD m_U;                         //!< copied gauge config. with boundary conditions.
  AFIELD m_Ublock;                    //!< copied gauge config. with block condition.

  std::string m_mode;                 //!< mult mode

  Fopr_CloverTerm *m_fopr_csw;        //!< clover term (corelib)
  Fopr_CloverTerm *m_fopr_csw_chiral; //!< clover term (corelib)

  AFIELD m_v1, m_v2;

  AFIELD m_T;           //!< clover term
  AFIELD m_T_qws;       //!< clover term

  int m_Nsize[4];       //!< lattice sizes (Nxv in x-direction)
  int m_block_sizev[4]; //!< block size in units of SIMD vector

  int do_comm[4];       // switchs of communication (4=Ndim): (0: n, 1: y).
  int do_comm_any;      // switchs of communication (if any): (0: n, 1: y).

  std::vector<int> m_bdsize;
  using allocator_t = typename AFIELD::template aligned_allocator<char>;
  using Channel     = Channel_impl<allocator_t>;
  std::vector<Channel> chsend_up, chrecv_up, chsend_dn, chrecv_dn;
  ChannelSet chset_send, chset_recv;

 public:
  //! constructor.
  AFopr_Clover_QWS_dd(const Parameters& params) : AFopr_dd<AFIELD>()
  {
    init();
    set_parameters(params);
  }

  //! destructor.
  ~AFopr_Clover_QWS_dd() { tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters by values.
  void set_parameters(real_t CKs, std::vector<int> bc,
                      std::vector<int> block_size);

  //! setting gauge configuration.
  void set_config(Field *u);

  //! QXS version requires convert of spinor field.
  bool needs_convert() { return true; }

  //! convert of spinor field (to QWS specific).
  void convert(AFIELD& v, const Field& w);

  //! convert of spinor field (from QXS to QWS specific).
  void convert(AFIELD& v, const AFIELD& w);

  //! reverse of spinor field (from QWS specific).
  void reverse(Field& v, const AFIELD& w);

  //! reverse of spinor field (from QWS specfic to QXS).
  void reverse(AFIELD& v, const AFIELD& w);

  //! returns the pointer to gauge configuration.
  inline Field *get_conf(void) { return m_conf; }

  //! setting mult mode.
  //void set_mode(std::string mode){ m_mode = mode; }
  void set_mode(std::string mode);

  //! returns mult mode.
  std::string get_mode() const;

  void mult(AFIELD&, const AFIELD&);
  void mult_dag(AFIELD&, const AFIELD&);
  void mult_gm5(AFIELD&, const AFIELD&);

  void project_chiral(AFIELD&, const AFIELD&, int ch);

  void mult_up(int mu, AFIELD&, const AFIELD&);
  void mult_dn(int mu, AFIELD&, const AFIELD&);

  void mult_sap(AFIELD&, const AFIELD&, const int ieo);
  void mult_dd(AFIELD&, const AFIELD&);
  void mult_dup(AFIELD&, const AFIELD&, const int mu);
  void mult_ddn(AFIELD&, const AFIELD&, const int mu);

  void mult_clv_inv(AFIELD& v);
  void mult_clv(AFIELD& v);

  // to be discarded
  void mult_block_hop(AFIELD& v, const AFIELD& w, const int mu)
  { mult_dup(v, w, mu); }

  //! returns inner size parameter.
  int field_nin() { return 2 * m_Nc * m_Nd; }

  //! returns local volume size parameter.
  int field_nvol() { return m_Nst; }

  //! returns external size parameter.
  int field_nex() { return 1; }

  //! returns floating operation counts.
  double flop_count() { return flop_count(m_mode); }

  //! returns floating operation counts of mult_sap.
  double flop_count_sap();

  //! returns floating operation counts for given mode.
  double flop_count(const std::string mode);

 private:
  //! initial setup.
  void init();

  //! final tidy-up.
  void tidyup();

  //! setup channels for communication.
  void setup_channels();

  //! inpose the boundary condition to link variable.
  void set_boundary();

  //! partially inpose the boundary condition to link variable.
  void set_boundary_config(AFIELD&, const int);

  //! inpose the block condition to link variable.
  void set_block_config(AFIELD&);

  //! set_csw now assumes Dirac repr.
  void set_csw();

  //! set_csw now assumes Dirac repr. with internally using Chiral repr.
  void set_csw_chrot();

  //! set_csw to use qws implemenation
  void set_csw_qws();

  void solve_csw_qws_inv(Field&, const Field&);

  //! set_csw now assumes Dirac repr.
  void mult_csw(real_t *, real_t *);

  void DdagD(AFIELD&, const AFIELD&);
  void Ddag(AFIELD&, const AFIELD&);
  void H(AFIELD&, const AFIELD&);
  void D(AFIELD&, const AFIELD&);

  void mult_gm4(AFIELD&, const AFIELD&);

  //! standard D mult.
  void mult_D(AFIELD&, const AFIELD&);

  //! D mult using QWS library.
  void mult_D_qws(AFIELD&, const AFIELD&);

  //! D mult using mult_xp, etc.
  void mult_D_alt(AFIELD&, const AFIELD&);

  void mult_xp(real_t *, real_t *, int);
  void mult_xm(real_t *, real_t *, int);
  void mult_yp(real_t *, real_t *, int);
  void mult_ym(real_t *, real_t *, int);
  void mult_zp(real_t *, real_t *, int);
  void mult_zm(real_t *, real_t *, int);
  void mult_tp(real_t *, real_t *, int);
  void mult_tm(real_t *, real_t *, int);

  void mult_gm5(real_t *, real_t *);

  void clear(real_t *);

  void scal_local(real_t *, real_t);
  void aypx(real_t, real_t *, real_t *);
  void gm5_aypx(real_t, real_t *, real_t *);


#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object_with_params(const Parameters& params)
  { return new AFopr_Clover_QWS_dd(params); }

 public:
  static bool register_factory()
  {
    bool init1 = AFopr<AFIELD>::Factory_params::Register("Clover_QWS_dd",
                                                         create_object_with_params);
    return init1;
  }
#endif
};

#endif  // AFOPR_CLOVER_DD_INCLUDED
