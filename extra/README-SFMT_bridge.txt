Note on the extra library SFMT-1.4.1_bridge
------------------------------------

The directory SFMT-1.4.1_bridge contains the original source code of 
SFMT (SIMD-oriented Fast Mersenne Twister) random number generator 
developed by Mutsuo Saito and Makoto Matsumoto, 
    http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/
and the jump-ahead feature presented in 
    http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/SFMT/JUMP/index.html
under the subdirectory jump/ .

In addition, the modified version of the jump-ahead routines are 
included in the library package for the use with the Bridge++ lattice 
QCD code set, which are:
  SFMT-jump-alt.{h, cpp}, SFMT-jump-params.h, and characteristic.[n].h 
where [n] are the Mersenne exponents. Makefile and SFMT.h are also 
modified. 

To use the jump-ahead feature, NTL (A Library for doing Number Theory) 
is also required.
http://shoup.net/ntl/

