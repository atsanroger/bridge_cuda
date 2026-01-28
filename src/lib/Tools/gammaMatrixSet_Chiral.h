/*!
        @file    gammaMatrixSet_Chiral.h

        @brief

        @author  Hideo Matsufuru  (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2019-01-21 17:01:33 #$

        @version $LastChangedRevision: 1928 $
*/

#ifndef GAMMAMATRIXSET_CHIRAL_INCLUDED
#define GAMMAMATRIXSET_CHIRAL_INCLUDED

#include "gammaMatrixSet.h"


//! Set of Gamma Matrix: chiral representation.

/*!
                                        [4 Feb 2012 H.Matsufuru]
 */

class GammaMatrixSet_Chiral : public GammaMatrixSet
{
 public:
  static const std::string class_name;

 public:
  GammaMatrixSet_Chiral()
  {
    init_GM();
  }

  void print();

 private:
  void init_GM();

#ifdef USE_FACTORY
 private:
  static GammaMatrixSet *create_object()
  {
    return new GammaMatrixSet_Chiral();
  }

 public:
  static bool register_factory()
  {
    return GammaMatrixSet::Factory::Register("Chiral", create_object);
  }
#endif
};
#endif
