/*!
      @file    aindex_block_lex.cpp
      @brief
      @author  Hideo Matsufuru (matufuru)
               $LastChangedBy: matufuru $
      @date    $LastChangedDate:: 2021-04-05 15:35:59 #$
      @version $LastChangedRevision: 2223 $
*/

#include "lib_alt_Accel/Field/aindex_block_lex.h"

template<typename REALTYPE>
const std::string AIndex_block_lex<REALTYPE,ACCEL>::class_name
                           = "AIndex_block_lex<REALTYPE,ACCEL>";

//====================================================================
template<typename REALTYPE>
void AIndex_block_lex<REALTYPE,ACCEL>::setup_buffer()
{
  // buffer field assumes complex variables on sites.
  buf_red.reset(2, m_fine_nvol, 1);
  buf_out.reset(2, m_coarse_nvol, 1);
}

//====================================================================
// explicit instanciation.

template<>
const std::string AIndex_block_lex<double, ACCEL>::class_name
                             = "AIndex_block_lex<double,ACCEL>";

template class AIndex_block_lex<double, ACCEL>;

template<>
const std::string AIndex_block_lex<float, ACCEL>::class_name
                             = "AIndex_block_lex<float,ACCEL>";

template class AIndex_block_lex<float, ACCEL>;

//============================================================END=====
