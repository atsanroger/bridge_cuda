/*!
      @file    afopr_Clover.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2160 $
*/

#ifndef AFOPR_CLOVER_H
#define AFOPR_CLOVER_H

#include <cstdio>
#include <cstdlib>

#include <string>
using std::string;

#include <vector>
using std::vector;

#include "lib/Fopr/afopr.h"

#include "lib/Parameters/commonParameters.h"
#include "lib/Parameters/parameters.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib/Communicator/communicator.h"
#include "lib/Tools/timer.h"

#if defined USE_MPI
#include "lib/Communicator/MPI/communicator_mpi.h"
#include "lib/Communicator/MPI/channel.h"
#else
#include "lib/Communicator/Single/channel.h"
#endif

#include "lib_alt_SIMD2/Fopr/afopr_CloverTerm.h"

//template<typename AFIELD> class AStaple_lex;
class Field;

template<typename AFIELD>
class AFopr_Clover : public AFopr<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 protected:
  int m_Nc, m_Nd, m_Ndim;
  int m_Nvc, m_Ndf, m_Ndm2;
  int m_Nx, m_Ny, m_Nz, m_Nt, m_Nst;
  int m_Nxv, m_Nstv;

  real_t  m_CKs;               //!< hopping parameter.
  real_t  m_cSW;               //!< clover coefficient.
  std::vector<int> m_boundary; //!< pointer to boundary condition
  std::string  m_repr;         //!< gamma matrix representation
  Bridge::VerboseLevel m_vl;   //!< verbose level

  Field *m_conf;        //!< original gauge config.
  AFIELD m_U;   //!< copied gauge config. with boundary conditions.

  std::string  m_mode;  //!< mult mode

  AFopr_CloverTerm<AFIELD>  *m_fopr_ct;

  AFIELD m_T;  //!< clover term

  AFIELD m_v1, m_v2;    //!< working spinor field.

  int do_comm[4];  // switchs of communication (4=Ndim): (0: n, 1: y).
  int do_comm_any; // switchs of communication (if any): (0: n, 1: y).

  std::vector<int> m_Nbdsize;
  std::vector<Channel*> chsend_up, chrecv_up, chsend_dn, chrecv_dn;
  ChannelSet chset_send, chset_recv;

  Timer m_timer;

public:
  //! constructor (to be discarded).
  AFopr_Clover(){ init(); }

  //! constructor.
  AFopr_Clover(const Parameters& params){
    init();
    set_parameters(params);
  }

  //! destructor.
  ~AFopr_Clover(){ tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters by values.
  void set_parameters(real_t CKs, real_t csW, std::vector<int> bc);

  //! setting gauge configuration.
  void set_config(Field* u);

  void set_config(unique_ptr<Field_G>& u)
  { set_config(u.get()); }

  //! setting fermion boundary condition to gauge field.
  void set_boundary_config(AFIELD& U, const int mu);
    
  //! returns the pointer to gauge configuration.
  inline Field* get_conf(void){return m_conf;};

  //! setting mult mode.
  void set_mode(std::string mode){ m_mode = mode; }

  //! returns mult mode.
  std::string get_mode() const;

  void mult(AFIELD&, const AFIELD&);
  void mult_dag(AFIELD&, const AFIELD&);
  void mult_gm5(AFIELD&, const AFIELD&);
  void mult_gm4(AFIELD&, const AFIELD&);

  void mult(AFIELD&, const AFIELD&,
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

  //! returns inner size parameter.
  int field_nin() { return 2 * m_Nc * m_Nd; }

  //! returns local volume size parameter.
  int field_nvol() { return m_Nst; }

  //! returns external size parameter.
  int field_nex() { return 1; }

  //! returns floating operation counts.
  double flop_count();

 private:
  //! initial setup.
  //  void init(const Parameters& params);
  void init();

  //! final tidy-up.
  void tidyup();

  //! setup channels for communication.
  void setup_channels();

  void mult_xp(real_t*, real_t*);
  void mult_xm(real_t*, real_t*);
  void mult_yp(real_t*, real_t*);
  void mult_ym(real_t*, real_t*);
  void mult_zp(real_t*, real_t*);
  void mult_zm(real_t*, real_t*);
  void mult_tp(real_t*, real_t*);
  void mult_tm(real_t*, real_t*);

  //! multiply diagonal term.
  void multadd_csw(real_t* v2, real_t* v1);

  //! multiply diagonal term at each site (Dirac repr).
  void multadd_csw_dirac(real_t* v2, real_t* v1, int site);

  //! multiply diagonal term at each site (Chiral repr).
  void multadd_csw_chiral(real_t* v2, real_t* v1, int site);

  void mult_gm5(real_t*, real_t*);
  void mult_gm4(real_t*, real_t*);

  void clear(real_t*);
  void aypx(real_t, real_t*, real_t*);
  void gm5_aypx(real_t, real_t*, real_t*);

  void clear(real_t*, const int, const int);
  void aypx(real_t, real_t*, real_t*, const int, const int);
  void gm5_aypx(real_t, real_t*, real_t*, const int, const int);


#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object_with_params(const Parameters& params)
  { return new AFopr_Clover(params); }

 public:
  static bool register_factory()
  {
    bool init1 = AFopr<AFIELD>::Factory_params::Register("Clover",
                                              create_object_with_params);
    return init1;
  }
#endif

};

#endif  // AFOPR_CLOVER_H
