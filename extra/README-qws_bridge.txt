Note on the extra library qws
------------------------------------

The directory qws_20210322_0f95_bridge contains a modified version of
qws library, which is developed by
  Yoshifumi Nakamura, Yuta Mukai, Ken-Ichi Ishikawa and Issaku Kanamori,
  https://github.com/RIKEN-LQCD/qws/
The version included in the Bridge++ is based on
  commit 0f958b3ae063252199c30b32234f4cfbd2797248  (Mon Mar 22 21:55:47 2021 +0900)
and the following files are modified:
    modified:   Makefile
    modified:   main.cc
    modified:   qws.cc
    modified:   rankmap_lib_utofu.c
    new file:   version_git

The modifications are to suppress generating redundant log files and
set the default compile option suitable to be link to Bridge++.
