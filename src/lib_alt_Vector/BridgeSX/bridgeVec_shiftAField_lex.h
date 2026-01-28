/*!
      @file    bridgeVec_shiftAField_lex.h
      @brief
      @author  Hideo Matsufuru (matufuru)
      @date    $LastChangedDate: 2013-01-22 13:51:53 #$
      @version $LastChangedRevision: 2543 $
*/

#ifndef BRIDGEVEC_SHIFTFIELD_LEX_INCLUDED
#define BRIDGEVEC_SHIFTFIELD_LEX_INCLUDED

namespace BridgeVec {

//====================================================================

// real_t = double

void shift_lex_xp1(double *RESTRICT buf, double *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_xp2(double *RESTRICT v2, double *RESTRICT buf, 
                   int nin, int *Nsize, int *bc);

void shift_lex_xpb(double *RESTRICT v2, double *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_xm1(double *RESTRICT buf, double *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_xm2(double *RESTRICT v2, double *RESTRICT buf, 
                   int nin, int *Nsize, int *bc);

void shift_lex_xmb(double *RESTRICT v2, double *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_yp1(double *RESTRICT buf, double *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_yp2(double *RESTRICT v2, double *RESTRICT buf, 
                   int nin, int *Nsize, int *bc);

void shift_lex_ypb(double *RESTRICT v2, double *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_ym1(double *RESTRICT buf, double *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_ym2(double *RESTRICT v2, double *RESTRICT buf, 
                   int nin, int *Nsize, int *bc);

void shift_lex_ymb(double *RESTRICT v2, double *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_zp1(double *RESTRICT buf, double *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_zp2(double *RESTRICT v2, double *RESTRICT buf, 
                   int nin, int *Nsize, int *bc);

void shift_lex_zpb(double *RESTRICT v2, double *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_zm1(double *RESTRICT buf, double *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_zm2(double *RESTRICT v2, double *RESTRICT buf, 
                   int nin, int *Nsize, int *bc);

void shift_lex_zmb(double *RESTRICT v2, double *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_tp1(double *RESTRICT buf, double *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_tp2(double *RESTRICT v2, double *RESTRICT buf, 
                   int nin, int *Nsize, int *bc);

void shift_lex_tpb(double *RESTRICT v2, double *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_tm1(double *RESTRICT buf, double *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_tm2(double *RESTRICT v2, double *RESTRICT buf, 
                   int nin, int *Nsize, int *bc);

void shift_lex_tmb(double *RESTRICT v2, double *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

//====================================================================

// real_t = float

void shift_lex_xp1(float *RESTRICT buf, float *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_xp2(float *RESTRICT v2, float *RESTRICT buf, 
                   int nin, int *Nsize, int *bc);

void shift_lex_xpb(float *RESTRICT v2, float *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_xm1(float *RESTRICT buf, float *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_xm2(float *RESTRICT v2, float *RESTRICT buf, 
                   int nin, int *Nsize, int *bc);

void shift_lex_xmb(float *RESTRICT v2, float *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_yp1(float *RESTRICT buf, float *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_yp2(float *RESTRICT v2, float *RESTRICT buf, 
                   int nin, int *Nsize, int *bc);

void shift_lex_ypb(float *RESTRICT v2, float *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_ym1(float *RESTRICT buf, float *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_ym2(float *RESTRICT v2, float *RESTRICT buf, 
                   int nin, int *Nsize, int *bc);

void shift_lex_ymb(float *RESTRICT v2, float *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_zp1(float *RESTRICT buf, float *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_zp2(float *RESTRICT v2, float *RESTRICT buf, 
                   int nin, int *Nsize, int *bc);

void shift_lex_zpb(float *RESTRICT v2, float *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_zm1(float *RESTRICT buf, float *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_zm2(float *RESTRICT v2, float *RESTRICT buf, 
                   int nin, int *Nsize, int *bc);

void shift_lex_zmb(float *RESTRICT v2, float *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_tp1(float *RESTRICT buf, float *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_tp2(float *RESTRICT v2, float *RESTRICT buf, 
                   int nin, int *Nsize, int *bc);

void shift_lex_tpb(float *RESTRICT v2, float *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_tm1(float *RESTRICT buf, float *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

void shift_lex_tm2(float *RESTRICT v2, float *RESTRICT buf, 
                   int nin, int *Nsize, int *bc);

void shift_lex_tmb(float *RESTRICT v2, float *RESTRICT v1,
                   int nin, int *Nsize, int *bc);

}

#endif

//====================================================================

