/*
        @file    $Id: run_test_afopr_coarse.cpp #$

        @brief

        @author  Tatsumi Aoyama <aoym@post.kek.jp> (aoym)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate: 2013-04-08 18:00:27 #$

        @version $LastChangedRevision: 2533 $
*/

#include <stdlib.h>

#ifdef USE_ALT_QXS
#include "lib_alt_QXS/bridge_alt_qxs.h"
#endif

#ifdef USE_ALT_ACCEL
#include "lib_alt_Accel/bridge_alt_accel.h"
#endif

#include "testlist_MG.h"

//====================================================================
int run_test_afopr_coarse()
{
  //  test_afopr_coarse::test();
  //test_afopr_coarse::test_parallel();
  vout.general("hoge!  run_test_afopr_coarse finised.\n");

  return EXIT_SUCCESS;
}
