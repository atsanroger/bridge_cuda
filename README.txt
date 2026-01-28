Release Note for Bridge++ ver.2.0.0                 01 Mar 2023
====================================

   *************************************************************
     This package is free software available under
     "the GNU General Public License" published by the FSF.
   *************************************************************

Overview:
---------

Bridge++ is a code set for lattice gauge theory simulations
on Linux/Unix workstations and supercomputers (Fugaku,
FX1000:Fujitsu, SX-Aurora TSUBASA:NEC, etc) using "C/C++"
standard language with MPI for communication and OpenMP for
multi-threading.


Software Compatibility:
-----------------------

This code set follows the standard specification of C++.
The code has been tested with GNU C/C++ 9.4, Clang C++ 14.0,
Intel C++ 2020.4 and 2023.0, NVIDIA HPC SDK C++ 22.7,
Fujitsu C++ (Clang mode) 1.2.36, NEC SX C++ 3.5.1.


Installation:
-------------

Edit Makefile appropriately to your environment, and execute
"make" (or in some cases, "gmake").

See also the file "INSTALL.txt".


Using the code set:
--------------------

After building the executable by make, executing "cd build/tests/"
and "./bridge.elf" in interactive environment will run the program
using test files. "./bridge.elf -h" shows the help messages.

  usage: ./bridge.elf
    (without args) -- run tests interactively.
    test_names     -- run specified tests.
    -a | --all     -- run all tests.
    -l | --list    -- list registered tests.
    -h | --help    -- show this message.

Sample directories are also prepared, such as sample_spectrum/.
The details are explained in the "first-step-guide", which is
available on our wiki website.


Getting Help:
-------------

Visit our website, http://bridge.kek.jp/Lattice-code/ for a
contact information. See also our first step guide and the
doxygen document.

Bug reports are welcome. Please send them to the contact
address on the above website.


Acknowledging Bridge++:
-----------------------

If you use this software for your research, please cite
J.Phys.Conf.Ser. 523 (2014) 012046
J.Phys.Conf.Ser. 2207 (2022) 012053
as well as our website.

For example, "This work is in part based on Bridge++ code
(http://bridge.kek.jp/Lattice-code/)~\cite{Ueda:2014rya,Akahoshi:2021gvk}".


Project Member:
---------------
Tatsumi   Aoyama    (Univ. Tokyo)
Issaku    Kanamori  (RIKEN R-CCS)
Kazuyuki  Kanaya    (Univ. of Tsukuba)
Hideo     Matsufuru (KEK)
Yusuke    Namekawa  (Hiroshima Univ.)
Hidekatsu Nemura    (Kyoto Univ.)
Yusuke    Taniguchi (Univ. of Tsukuba)

Contributed by:
---------------
Sinya     Aoki      (Kyoto Univ.)
Takumi    Doi       (RIKEN)
Shoji     Hashimoto (KEK)
Noriyoshi Ishii     (Osaka Univ.)
Ken-ichi  Ishikawa  (Hiroshima Univ.)
Takashi   Kaneko    (KEK)
Yoshinobu Kuramashi (Univ. of Tsukuba)
Keigo     Nitadori  (RIKEN R-CCS)
Kenji     Sasaki    (Osaka Univ.)
Naoya     Ukita     (Univ. of Tsukuba)
Tomoteru  Yoshie    (Univ. of Tsukuba)

Former Member:
---------------
Yutaro    Akahoshi
Guido     Cossu
Takaya    Miyamoto
Shinji    Motoki
Jun-Ichi  Noaki
Kenji     Ogawa
Hana      Saito
Satoru    Ueda


This work is supported by the following grants, research programs,
and organizations:
- Grant-in-Aid for Scientific Research on Innovative Areas (of
  the Japanese Ministry of Education, Culture, Sports, Science
  and Technology) (Nos.20105001, 20105005)
- JSPS KAKENHI Grant Numbers 25400284, 15K05068, 16K05340,
  16H03988, 20K03961.
- HPCI Strategic Program Field 5 ("The origin of matter and
  the universe"),
- Priority Issue on Post-K computer (Elucidation of the
  Fundamental Laws and Evolution of the Universe)
- Joint Institute for Computational Fundamental Science (JICFuS)
- Program for Promoting Researches on the Supercomputer Fugaku
  "Simulation for basic science: from fundamental laws of
   particles to creation of nuclei"
- Multidisciplinary Cooperative Research Program, Center for
  Computational Science, University of Tsukuba (15a33, 16a42,
  17a41, 18-14, 19-43, 20-79, 21-42, 22-60)
- KEK Large-scale Simulation Program (Ohgata 09-22, 09/10-23,
  10-18, (T)11-09, 12-18, 12/13-15)
- Particle, Nuclear, and Astro Physics Simulation Program,
  Institute of Particle and Nuclear Studies, KEK (2019-T003,
  2019-004, 2020-002, 2021-006, 2022-002)
- KEK Computing Research Center
- Research Center for Nuclear Physics, Osaka University
- Yukawa Institute for Theoretical Physics, Kyoto University
