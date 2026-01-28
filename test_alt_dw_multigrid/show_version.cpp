/*
        @file    $Id: show_version.cpp #$

        @brief   utility to print the version

        @author  Issaku Kanamori
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:  #$

        @version $LastChangedRevision: 2578 $
*/
#include "IO/bridgeIO.h"

namespace Bridge {
  void show_version()
  {
#define GET_CHAR_FROM_MACRO_(str)    # str
#define GET_CHAR_FROM_MACRO(str)     GET_CHAR_FROM_MACRO_(str)

    vout.general("version info:\n");
#ifdef _BRIDGE_VERSION
    vout.general("  Bridge: %s\n", GET_CHAR_FROM_MACRO(_BRIDGE_VERSION));
#endif
#ifdef _COMMIT_ID
    vout.general("  %s\n", GET_CHAR_FROM_MACRO(_COMMIT_ID));
#endif
#ifdef _COMMIT_AUTHOR
    vout.general("  %s\n", GET_CHAR_FROM_MACRO(_COMMIT_AUTHOR));
#endif
#ifdef _COMMIT_DATE
    vout.general("  %s\n", GET_CHAR_FROM_MACRO(_COMMIT_DATE));
#endif
    vout.general("\n");

#undef GET_CHAR_FROM_MACRO
  }
}
