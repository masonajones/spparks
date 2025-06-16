/* ----------------------------------------------------------------------
   SPPARKS - Stochastic Parallel PARticle Kinetic Simulator
   http://www.cs.sandia.gov/~sjplimp/spparks.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2008) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level SPPARKS directory.
------------------------------------------------------------------------- */
#ifndef SPK_RANDOM_FAST_H
#define SPK_RANDOM_FAST_H

#define AVX2

#include "spktype.h"
#include <vector>

#ifdef AVX2
//#ifdef __AVX2__
#include "xoshiro256_avx.h"

namespace SPPARKS_NS {

class RandomFast {
 public:
  uint64_t seed;

  RandomFast(int);
  RandomFast(double);
  ~RandomFast() {}
  
  void reset(double, int, int);
  void tagreset(double, tagint, int);
  int irandom(int);
  void init_bulkRand();
  double uniform();
  tagint tagrandom(tagint);
  bigint bigrandom(bigint);

 private:
   static constexpr int numRand = 8000; 
   std::vector<double> bulkRand;    
   int iter;
   Xoshiro256AVX2 prng;
};

}

#else

namespace SPPARKS_NS {

class RandomFast {
 public:
  uint64_t seed;

  RandomFast(int);
  RandomFast(double);
  ~RandomFast() {}
  
  void reset(double, int, int);
  void tagreset(double, tagint, int);
  int irandom(int);
  void init_bulkRand();
  double uniform();
  tagint tagrandom(tagint);
  bigint bigrandom(bigint);
              
 private:
   static constexpr int numRand = 4000; 
   std::vector<double> bulkRand; 
   int iter;
   double uniform_slow();

};

}

#endif


#endif
