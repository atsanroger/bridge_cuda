/*!
        @file    gaugeConfig_SF.h

        @brief

        @author  Yusuke Taniguchi
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2022-10-21 12:25:27 #$

        @version $LastChangedRevision: 2408 $
*/


#ifndef GAUGECONFIG_SF_INCLUDED
#define GAUGECONFIG_SF_INCLUDED

#include <string>
using std::string;

#include "gaugeConfig.h"
#include "Field/index_lex.h"

#include "bridgeIO.h"
using Bridge::vout;

/*!
  Setup the tree level gauge configuration for the SF boundary condition, which minimize the classical action.
  \f[
  U_3=1
  \f]
  \f[
  U_k(x)=\exp\left(aB_k(x)\right),\quad
  B_k(x)=\frac{1}{T}\left(x_0C_k'+(T-x_0)C_k\right)
  \f]
  \f[
  C_k=\frac{i}{L}\pmatrix{\phi_1\cr&\phi_2\cr&&\phi_3\cr},\quad
  C_k'=\frac{i}{L}\pmatrix{\phi'_1\cr&\phi'_2\cr&&\phi'_3\cr}
  \f]
 */
class GaugeConfig_SF : public GaugeConfig
{
 public:
  static const std::string class_name;

 public:
  GaugeConfig_SF(const string& type) : GaugeConfig(type) {}

  void set_cold_SF(Field_G& U, double phi[3], double phipr[3]);
};
#endif
