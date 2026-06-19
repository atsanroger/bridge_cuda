/*!
      @file    bridgeACC_ShiftAField_lex.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2539 $
*/

#ifndef BRIDGEACC_SHIFTAFIELD_LEX_INCLUDED
#define BRIDGEACC_SHIFTAFIELD_LEX_INCLUDED

namespace BridgeACC {

  // real_t = double

  void shift_lex_xp1(real_t* buf, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_xp2(real_t* v2, real_t* buf, 
                     int nin, int *Nsize, int *bc);

  void shift_lex_xpb(real_t* v2, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_xm1(real_t* buf, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_xm2(real_t* v2, real_t* buf, 
                     int nin, int *Nsize, int *bc);

  void shift_lex_xmb(real_t* v2, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_yp1(real_t* buf, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_yp2(real_t* v2, real_t* buf, 
                     int nin, int *Nsize, int *bc);

  void shift_lex_ypb(real_t* v2, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_ym1(real_t* buf, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_ym2(real_t* v2, real_t* buf, 
                     int nin, int *Nsize, int *bc);

  void shift_lex_ymb(real_t* v2, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_zp1(real_t* buf, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_zp2(real_t* v2, real_t* buf, 
                     int nin, int *Nsize, int *bc);

  void shift_lex_zpb(real_t* v2, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_zm1(real_t* buf, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_zm2(real_t* v2, real_t* buf, 
                     int nin, int *Nsize, int *bc);

  void shift_lex_zmb(real_t* v2, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_tp1(real_t* buf, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_tp2(real_t* v2, real_t* buf, 
                     int nin, int *Nsize, int *bc);

  void shift_lex_tpb(real_t* v2, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_tm1(real_t* buf, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_tm2(real_t* v2, real_t* buf, 
                     int nin, int *Nsize, int *bc);

  void shift_lex_tmb(real_t* v2, real_t* v1,
                     int nin, int *Nsize, int *bc);


  // real_t = double

  void shift_lex_xp1(real_t* buf, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_xp2(real_t* v2, real_t* buf, 
                     int nin, int *Nsize, int *bc);

  void shift_lex_xpb(real_t* v2, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_xm1(real_t* buf, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_xm2(real_t* v2, real_t* buf, 
                     int nin, int *Nsize, int *bc);

  void shift_lex_xmb(real_t* v2, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_yp1(real_t* buf, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_yp2(real_t* v2, real_t* buf, 
                     int nin, int *Nsize, int *bc);

  void shift_lex_ypb(real_t* v2, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_ym1(real_t* buf, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_ym2(real_t* v2, real_t* buf, 
                     int nin, int *Nsize, int *bc);

  void shift_lex_ymb(real_t* v2, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_zp1(real_t* buf, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_zp2(real_t* v2, real_t* buf, 
                     int nin, int *Nsize, int *bc);

  void shift_lex_zpb(real_t* v2, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_zm1(real_t* buf, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_zm2(real_t* v2, real_t* buf, 
                     int nin, int *Nsize, int *bc);

  void shift_lex_zmb(real_t* v2, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_tp1(real_t* buf, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_tp2(real_t* v2, real_t* buf, 
                     int nin, int *Nsize, int *bc);

  void shift_lex_tpb(real_t* v2, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_tm1(real_t* buf, real_t* v1,
                     int nin, int *Nsize, int *bc);

  void shift_lex_tm2(real_t* v2, real_t* buf, 
                     int nin, int *Nsize, int *bc);

  void shift_lex_tmb(real_t* v2, real_t* v1,
                     int nin, int *Nsize, int *bc);

}

#endif
