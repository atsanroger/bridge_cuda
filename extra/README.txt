a) LIME
The latest version of LIME is available at http://usqcd-software.github.io/c-lime/
lime-1.3.2 is verified for Bridge++

b) SFMT (Mersenne Twister)
See README-SFMT_bridge.txt

c) NTL (used with SFMT)
The latest version of NTL is available at http://www.shoup.net/ntl/
ntl-11.5.1 is verified for Bridge++

d) FFTW
The latest version of FFTW is available at http://www.fftw.org/
fftw-3.3.10 is verified for Bridge++

e) yaml-cpp
The latest version of yaml-cpp is available at https://github.com/jbeder/yaml-cpp/
yaml-cpp-release-0.7.0 is verified for Bridge++

f) TinyXML-2
The latest version of TinyXML-2 is available at https://github.com/leethomason/tinyxml2/
tinyxml2-9.0.0 is verified for Bridge++.
unpack the release package under extra/tinyxml2, and build using cmake as follows:
  mkdir build && cd build
  cmake .. -DCMAKE_INSTALL_PREFIX=.. -DBUILD_SHARED_LIBS=OFF && make && make install

