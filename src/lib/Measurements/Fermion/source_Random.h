/*!
        @file    source_Random.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef SOURCE_RANDOM_INCLUDED
#define SOURCE_RANDOM_INCLUDED

#include "source.h"

#include "Field/index_lex.h"
#include "Tools/randomNumberManager.h"

//! Random noise source in the space at a given timeslice

/*!
    Random noise source in the space at a given timeslice
      src = \eta_{i_color,i_spin}^{i_noise}(x),
    where \eta is a random number field satisfying
      \sum_{i_noise} \eta^{i_noise} = 0
      \sum \eta_i^{i_noise,\dagger} \eta_j^{i_noise} = \delta_{ij}
                                     [10 Jan 2017 Y.Namekawa]
  Add two functions set_all_space_time(ic), set_all_space_time(ic,is) intended to be used in VEV calculation with fermion_flow.
  [21 December 2017 Y.Taniguchi]
 */

class Source_Random : public Source
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  RandomNumbers *m_rand;
  Index_lex m_index;
  std::vector<int> m_source_position;
  std::vector<int> m_source_momentum;
  std::string m_str_noise_type;
  bool m_in_node;

 public:
  Source_Random()
    : m_vl(CommonParameters::Vlevel()),
    m_rand(RandomNumberManager::getInstance())
  {
  }

  Source_Random(const Parameters& params)
    : m_vl(CommonParameters::Vlevel()),
    m_rand(RandomNumberManager::getInstance())
  {
    set_parameters(params);
  }

  /* Source_Random(RandomNumbers *rand) */
  /*   : Source(), m_rand(rand) */
  /*   {} */

  void set_parameters(const Parameters& params);
  void set_parameters(const std::vector<int>& source_position,
                      const std::vector<int>& source_momentum,
                      const std::string noise_type);

  void get_parameters(Parameters& params) const;

  void set(Field& src, const int idx);
  void set(Field& src, const int i_color, const int i_spin);
  void set_all_color(Field& src, const int i_spin);
  void set_all_color_spin(Field& src);

  //! Setting a noise vector. Filling all the sites and spin indices for color index "ic". The same random number is set for all the spin indices at one site.
  void set_all_space_time(Field& src, const int ic);

  //! Setting a noise vector. Filling all the sites for spin-color index "is" and "ic".
  void set_all_space_time(Field& src, const int ic, const int is);

#ifdef USE_FACTORY
 private:
  static Source *create_object()
  {
    return new Source_Random();
  }

  // static Source *create_object_with_arg(RandomNumbers *rand)
  // {
  //   return new Source_Random(rand);
  // }

  static Source *create_object_with_params(const Parameters& params)
  {
    return new Source_Random(params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Source::Factory::Register("Random", create_object);
    // init &= Source::Factory_rand::Register("Random", create_object_with_arg);
    init &= Source::Factory_params::Register("Random", create_object_with_params);
    return init;
  }
#endif
};
#endif
