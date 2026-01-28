/*!
      @file    afopr_Wilson_eo.h
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2517 $
*/

#ifndef QXS_AFOPR_WILSON_EO_INCLUDED
#define QXS_AFOPR_WILSON_EO_INCLUDED

#include <cstdio>
#include <cstdlib>

#include <string>
using std::string;
#include <vector>
using std::vector;

#include "lib/Fopr/afopr.h"
#include "lib/Parameters/commonParameters.h"
#include "lib/Communicator/communicator_impl.h"
#include "lib/IO/bridgeIO.h"
using Bridge::vout;

#include "lib_alt_QXS/Field/afield.h"
#include "lib_alt_QXS/Field/aindex_eo.h"

class Field;

/*!
  Implementation of the Wilson fermion operator (even-odd index)
  in the QXS branch.
                                      [24 Dec 2022 H.Matsufuru]
 */
template<typename AFIELD>
class AFopr_Wilson_eo : public AFopr<AFIELD>
{
 public:
  typedef typename AFIELD::real_t real_t;
  static const std::string class_name;

 protected:
  int m_Nc, m_Nd, m_Ndf, m_Nvc;
  int m_Nst, m_Nx, m_Ny, m_Nz, m_Nt;
  int m_Ndim;
  int m_Nx2, m_Nst2;
  int m_Nx2v, m_Nst2v, m_Nyv;

  real_t m_CKs;                //!< hopping parameter.
  std::vector<int> m_boundary; //!< pointer to boundary condition
  std::string m_repr;          //!< gamma matrix representation
  Bridge::VerboseLevel m_vl;   //!< verbose level

  std::vector<int> m_Leo;      //!< Leo = 0 (even site) or 1 (odd site).

  AIndex_eo<real_t, AFIELD::IMPL> m_index;
  Field *m_conf;      //!< original gauge configuration
  AFIELD m_Ulex;      //!< converted gauge config with boundary conditions
  AFIELD m_Ueo;       //!< gauge config in even-odd index

  std::string m_mode; //!< mult mode

  AFIELD m_v1, m_v2;  //!< working field.
  AFIELD m_z1;        //!< lexical field: used in convert/reverse.

  int m_Nsize[4];     // lattice sizes (Nx2v in x-direction)

  int do_comm[4];     // switchs of communication (4=Ndim): (0: n, 1: y).
  int do_comm_any;    // switchs of communication (if any): (0: n, 1: y).

  std::vector<int> m_bdsize;
  using allocator_t = typename AFIELD::template aligned_allocator<char>;
  using Channel     = Channel_impl<allocator_t>;
  std::vector<Channel> chsend_up, chrecv_up, chsend_dn, chrecv_dn;
  ChannelSet chset_send, chset_recv;

 public:
  //! constructor.
  AFopr_Wilson_eo(const Parameters& params)
  {
    init(params);
  }

  //! destructor.
  ~AFopr_Wilson_eo() { tidyup(); }

  //! setting parameters by a Parameter object.
  void set_parameters(const Parameters& params);

  //! setting parameters by values.
  void set_parameters(real_t CKs, std::vector<int> bc);

  //! get parameters via a Parameter object
  void get_parameters(Parameters& params) const;

  //! setting gauge configuration (common interface).
  void set_config(Field *u);

  //! setting gauge configuration (setting omp parallel).
  void set_config_omp(Field *u);

  //! setting gauge configuration (implementation).
  void set_config_impl(Field *u);

  //! setting mult mode.
  void set_mode(std::string mode);

  //! returns mult mode.
  std::string get_mode() const { return m_mode; }

  void mult(AFIELD&, const AFIELD&);
  void mult_dag(AFIELD&, const AFIELD&);
  void mult_gm5(AFIELD&, const AFIELD&);

  void mult_gm4(AFIELD&, const AFIELD&);

  void mult(AFIELD&, const AFIELD&,
            const std::string mode);

  void DdagD(AFIELD&, const AFIELD&);
  void Ddag(AFIELD&, const AFIELD&);
  void D(AFIELD&, const AFIELD&);
  void H(AFIELD&, const AFIELD&);

  void aypx(real_t, AFIELD&, const AFIELD&);

  //! Fermion matrix with ieo = 0: even <-- odd, 1: odd <-- even.
  void Meo(AFIELD&, const AFIELD&, const int ieo);

  //! Fermion matrix with ieo = 0: even <-- odd, 1: odd <-- even.
  void Meo(AFIELD&, const AFIELD&, const AFIELD&,
           const int ieo, const int iflag);

  //! Meo implementation: standard
  void mult_Meo(AFIELD&, const AFIELD&, const AFIELD&,
                const int ieo, const int iflag);

  //! Meo implementation: using mult_xp etc.
  void mult_Meo_alt(AFIELD&, const AFIELD&, const AFIELD&,
                    const int ieo, const int iflag);

  //! Meo implementation: using qxs library
  void mult_Meo_qxs(AFIELD&, const AFIELD&, const AFIELD&,
                    const int ieo, const int iflag);

  //! returns inner size parameter.
  int field_nin() { return 2 * m_Nc * m_Nd; }

  //! returns local volume size parameter.
  int field_nvol() { return m_Nst2; }

  //! returns external size parameter.
  int field_nex() { return 1; }

  //! returns floating operation counts.
  double flop_count() { return flop_count(m_mode); }

  //! returns floating operation counts.
  double flop_count(const std::string mode);

  //! field convert is necessary in this implementation.
  bool needs_convert() { return true; }

  //! convert Field to AField for this class.
  void convert(AFIELD&, const Field&);

  //! reverse AField to Field.
  void reverse(Field&, const AFIELD&);

 private:
  //! initial setup.
  void init(const Parameters& params);

  //! final tidy-up.
  void tidyup();

  //! setup channels for communication.
  void setup_channels();

  void mult_gm4(real_t *, real_t *);

  void clear(real_t *);

  void aypx(real_t, real_t *, real_t *);
  void scal(real_t *, const real_t);

  void clear(real_t *, const int, const int);

  void aypx(real_t, real_t *, real_t *, const int, const int);
  void scal(real_t *, const real_t, const int, const int);


#ifdef USE_FACTORY
 private:
  static AFopr<AFIELD> *create_object_with_params(const Parameters& params)
  { return new AFopr_Wilson_eo(params); }

 public:
  static bool register_factory()
  {
    bool init1 = AFopr<AFIELD>::Factory_params::Register("Wilson_eo",
                                                         create_object_with_params);
    return init1;
  }
#endif
};

#endif
