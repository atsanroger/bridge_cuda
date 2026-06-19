/*!
        @file    director.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/


#ifndef DIRECTOR_INCLUDED
#define DIRECTOR_INCLUDED

#include "bridge_defs.h"
#include "Parameters/commonParameters.h"
#include "Parameters/parameters.h"

#include "Field/field_G.h"

#include "IO/bridgeIO.h"

//! Manager of commonly used data object in HMC.

/*!
    Director-type class manages data which requires memory cost
    and/or computational cost to be held in multiple objects.
    Examples are smeared configurations and eigenvectors.
    This mechanism is useful mainly in HMC, while also in
    defining such as smeared fermion operators.
                                    [28 Dec 2011 H.Matsufuru]
 */

class Director
{
 public:
  Director() {}
  virtual ~Director() {}

 private:
  // non-copyable
  Director(const Director&);
  Director& operator=(const Director&);

 public:
  // To be called when link variable is updated.
  virtual void notify_linkv() = 0;

  virtual void set_parameters(const Parameters& params) = 0;

  virtual void get_parameters(Parameters& params) const = 0;

  virtual void set_config(Field *U) = 0;
};
#endif
