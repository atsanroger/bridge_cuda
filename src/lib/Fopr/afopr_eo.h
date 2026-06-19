/*!
        @file    afopr_eo.h

        @brief

        @author  Satoru Ueda  (sueda)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-03-07 17:24:38 #$

        @version $LastChangedRevision: 2359 $
*/

#ifndef AFOPR_EO_INCLUDED
#define AFOPR_EO_INCLUDED

#include "Fopr/afopr.h"

//! Base class of fermion operator family.

/*!
    In Bridge-OF, the fermion operator implies an operator which
    transforms a field to other field, irrespective of physical
    formulation.
    This class defines the interface of the fermion operators.
    At present, void functions mult(v,w) and mult_dag(v,w) is
    not purely virtual, because some of subclass have not
    implemented them yet.
                          [20 Dec 2011  H.Matsufuru]
    unique_ptr is introduced to avoid memory leaks
                          [21 Mar 2015 Y.Namekawa]
*/

template<typename AFIELD>
class AFopr_eo : public AFopr<AFIELD>
{
 protected:
  static const std::string class_name;

 public:
  virtual ~AFopr_eo() {}

  // methods for even odd fermion operator
  virtual void preProp(AFIELD&, AFIELD&, const AFIELD&)        = 0;
  virtual void postProp(AFIELD&, const AFIELD&, const AFIELD&) = 0;

  //! setting pointer to the gauge configuration.
  virtual void set_config(Field *) = 0;

  //! \brief multiplies fermion operator to a given field (2nd argument)
  //   and set the resultant field to the 1st argument.
  virtual void mult(AFIELD&, const AFIELD&) {}

  //! hermitian conjugate of mult(Field&, const Field&).
  virtual void mult_dag(AFIELD&, const AFIELD&) {}

  virtual void mult(AFIELD&, const AFIELD&, const std::string) {}

  virtual void mult_dag(AFIELD&, const AFIELD&, const std::string) {}

  //! \brief setting the mode of multiplication if necessary.
  //!  Default implementation here is just to avoid irrelevant call.
  virtual void set_mode(const std::string mode)
  {
    vout.general("AFopr_eo: set_mode not implemented.\n");
  }

  std::string get_mode() const
  {
    vout.general("AFopr_eo: get_mode not implemented.\n");
    return std::string();
  }

  //! returns the volume for which the fermion operator is defined.
  virtual int field_nvol() = 0;

  //! returns the on-site d.o.f. for which the fermion operator is defined.
  virtual int field_nin() = 0;

  //! returns the external d.o.f. for which the fermion operator is defined.
  virtual int field_nex() = 0;
};
#endif
