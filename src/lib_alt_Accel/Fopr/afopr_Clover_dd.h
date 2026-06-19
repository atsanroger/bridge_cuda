/*!
      @file    afopr_Clover_dd.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2537 $
*/

#ifndef ACCEL_AFOPR_CLOVER_DD_INCLUDED
#define ACCEL_AFOPR_CLOVER_DD_INCLUDED

#include <cstdio>
#include <cstdlib>

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "lib_alt/Fopr/afopr_dd.h"

#include "lib/Parameters/commonParameters.h"
#include "lib/Parameters/parameters.h"
#include "lib/Communicator/communicator.h"
#include "lib/Communicator/communicator_impl.h"
#include "lib/Tools/timer.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt_Accel/Fopr/afopr_CloverTerm.h"

class Field;

template<typename AFIELD>
class AFopr_Clover_dd : public AFopr_dd<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 protected:
  enum Gamma_repr {DIRAC, CHIRAL} m_repr;  //!< gamma matrix representation

  int m_Nc, m_Nd, m_Ndim;
  int m_Nvc, m_Ndf, m_Ndm2;
  int m_Nx, m_Ny, m_Nz, m_Nt, m_Nst;
  //  int m_Nxv, m_Nstv;

  real_t  m_CKs;               //!< hopping parameter.
  real_t  m_cSW;               //!< clover coefficient.
  std::vector<int> m_boundary; //!< pointer to boundary condition
  //  std::string  m_repr;         //!< gamma matrix representation
  Bridge::VerboseLevel m_vl;   //!< verbose level

  Field *m_conf;        //!< original gauge config.
  AFIELD m_U;   //!< copied gauge config. with boundary conditions.

  std::string  m_mode;  //!< mult mode

  AFopr_CloverTerm<AFIELD>  *m_fopr_ct;

  AFIELD m_T;  //!< clover term

  bool m_memory_saved;

  AFIELD m_v1, m_v2;    //!< working spinor field.

  int do_comm[4];  // switchs of communication (4=Ndim): (0: n, 1: y).
  int do_comm_any; // switchs of communication (if any): (0: n, 1: y).

  int m_Nsize[4];
  int m_bc[4];
  int m_bc2[4];

  std::vector<int> m_Nbdsize;
  using Channel = Channel_impl<std::allocator<real_t> >;
  std::vector<Channel> chsend_up, chrecv_up, chsend_dn, chrecv_dn;
  ChannelSet chset_send, chset_recv;

  Timer m_timer;

  std::vector<int> m_block_size;  //!< block size

  int m_Ieo;   //!< even-odd label of origin in units of block

  AFIELD m_Ublock;  //!< copied gauge config. with block condition.
  AFIELD m_Uhop;    //!< = m_U - m_Ublock

public:
  //! constructor (to be discarded).
  //AFopr_Clover_dd(){ init(); }

  //! constructor.
  AFopr_Clover_dd(const Parameters& params){
    init(params);
    //set_parameters(params);
  }

  //! destructor.
  ~AFopr_Clover_dd(){ tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters by values.
  void set_parameters(real_t CKs, real_t cSW, std::vector<int> bc,
                      std::vector<int> block_size);

  void get_parameters(Parameters& params) const;

  //! setting gauge configuration.
  void set_config(Field* u);

  //  void set_config(unique_ptr<Field_G>& u)
  //  { set_config(u.get()); }

  //! returns the pointer to gauge configuration.
  inline Field* get_conf(void){return m_conf;};

  //! setting mult mode.
  void set_mode(std::string mode);

  //! returns mult mode.
  std::string get_mode() const;

  void mult(AFIELD&, const AFIELD&);
  void mult_dag(AFIELD&, const AFIELD&);
  void mult_gm5(AFIELD&, const AFIELD&);
  void mult_gm4(AFIELD&, const AFIELD&);

  void mult(AFIELD&, const AFIELD&,
            const std::string mode);

  void mult_dag(AFIELD&, const AFIELD&,
                const std::string mode);

  //! upward nearest neighbor hopping term.
  virtual void mult_up(int mu, AFIELD&, const AFIELD&);

  //! downward nearest neighbor hopping term.
  virtual void mult_dn(int mu, AFIELD&, const AFIELD&);

  void DdagD(AFIELD&, const AFIELD&);
  void Ddag(AFIELD&, const AFIELD&);
  void H(AFIELD&, const AFIELD&);
  void D(AFIELD&, const AFIELD&);
  void gm5(AFIELD&, const AFIELD&);

  void H_alt(AFIELD&, const AFIELD&);
  void D_alt(AFIELD&, const AFIELD&);

  void mult_sap(AFIELD&, const AFIELD&, const int ieo);

  void mult_dd(AFIELD&, const AFIELD&);

  void mult_dup(AFIELD&, const AFIELD&, const int mu);

  void mult_ddn(AFIELD&, const AFIELD&, const int mu);

  // to be discarded
  void mult_block_hop(AFIELD &v, const AFIELD &w, const int mu)
  { mult_dup(v, w, mu); }

  //! returns inner size parameter.
  int field_nin() { return 2 * m_Nc * m_Nd; }

  //! returns local volume size parameter.
  int field_nvol() { return m_Nst; }

  //! returns external size parameter.
  int field_nex() { return 1; }

  //! returns floating operation counts.
  double flop_count(){ return flop_count(m_mode); }

  //! returns floating operation counts of mult_sap.
  double flop_count_sap();

  //! returns floating operation counts for given mode.
  double flop_count(const std::string mode);

 private:
  //! initial setup.
  //  void init(const Parameters& params);
  void init(const Parameters& params);

  //! final tidy-up.
  void tidyup();

  //! setup channels for communication.
  void setup_channels();

  void set_config_omp(Field* u);
  void set_config_impl(Field* u);
  void set_config_saved(Field* u);

  //! inpose the block condition to link variable.
  void set_block_config(AFIELD&);

  void mult_xp(real_t*, real_t*, int);
  void mult_xm(real_t*, real_t*, int);
  void mult_yp(real_t*, real_t*, int);
  void mult_ym(real_t*, real_t*, int);
  void mult_zp(real_t*, real_t*, int);
  void mult_zm(real_t*, real_t*, int);
  void mult_tp(real_t*, real_t*, int);
  void mult_tm(real_t*, real_t*, int);

  void scal_local(real_t *v, real_t a);

  //! multiply diagonal term.
  void multadd_csw(real_t* v2, real_t* v1);

  //! multiply diagonal term at each site (Dirac repr).
  void multadd_csw_dirac(real_t* v2, real_t* v1, int site);

  //! multiply diagonal term at each site (Chiral repr).
  void multadd_csw_chiral(real_t* v2, real_t* v1, int site);

  void mult_gm5(real_t*, real_t*);
  void mult_gm4(real_t*, real_t*);

  // void clear(real_t*);
  //  void aypx(real_t, real_t*, real_t*);
  //  void axpy(real_t*, real_t, real_t*);
  //  void gm5_aypx(real_t, real_t*, real_t*);

  //  void clear(real_t*, const int, const int);
  //  void aypx(real_t, real_t*, real_t*, const int, const int);
  //  void gm5_aypx(real_t, real_t*, real_t*, const int, const int);


#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object_with_params(const Parameters& params)
  { return new AFopr_Clover_dd(params); }

 public:
  static bool register_factory()
  {
    bool init1 = AFopr<AFIELD>::Factory_params::Register("Clover_dd",
                                              create_object_with_params);
    return init1;
  }
#endif

};

#endif  // AFOPR_CLOVER_DD_H
