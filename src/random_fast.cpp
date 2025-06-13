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
#define AVX2
#include "spktype.h"
#include "math.h"
#include "random_fast.h"

using namespace SPPARKS_NS;


#define IM 2147483647

#include <vector>


#ifdef AVX2
#include "xoshiro256_avx.h"

#else
#define IA 16807
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define A 48271

#endif

/* ---------------------------------------------------------------------- 
   Buffered Xoshiro256 RNG
   Falls back to Park if AVX2 not enabled, but maintains buffer
   assume iseed is a positive int
------------------------------------------------------------------------ */

RandomFast::RandomFast(int iseed) : seed(static_cast<uint64_t>(iseed)), bulkRand(numRand), iter(numRand+1), prng(static_cast<uint64_t>(iseed))
{
  static_assert(numRand%4 == 0, "RandomFast: Buffer size (numRand) must be multiple of 4 when AVX2 is enabled");
  //seed = iseed;
}

/* ---------------------------------------------------------------------- 
   set seed to positive int
   assume 0.0 <= rseed < 1.0
------------------------------------------------------------------------ */

RandomFast::RandomFast(double rseed) : seed(static_cast<uint64_t> (rseed*IM)), bulkRand(numRand), iter(numRand+1), prng(static_cast<uint64_t>(rseed * IM))
{
  //seed = static_cast<int> (rseed*IM);
  if (seed == 0) seed = 1;
}

// RandomFast::~RandomFast() 
// {

// }


/* ---------------------------------------------------------------------- 
   reset seed to a positive int based on rseed and offset
   assume 0.0 <= rseed < 1.0 and offset is an int >= 0
   fmod() insures no overflow when static cast to int
   warmup the new RNG if requested
   typically used to setup one RN generator per proc or site or particle
------------------------------------------------------------------------ */
#ifdef AVX2
// Xoshiro RNG uses internal state variables, need to set them differently
void RandomFast::reset(double rseed, int offset, int warmup)
{
  seed = static_cast<int> (fmod(rseed*IM+offset,IM));
  if (seed < 0) seed = -seed;
  if (seed == 0) seed = 1;
  //prng = Xoshiro256AVX2(static_cast<uint64_t>(seed));
  prng.reset(static_cast<uint64_t>(seed));
  for (int i = 0; i < warmup; i++) uniform();
}

void RandomFast::tagreset(double rseed, tagint offset, int warmup)
{
  seed = static_cast<int> (fmod(rseed*IM+offset,IM));
  if (seed < 0) seed = -seed;
  if (seed == 0) seed = 1;
  //prng = Xoshiro256AVX2(static_cast<uint64_t>(seed));
  prng.reset(static_cast<uint64_t>(seed));
  for (int i = 0; i < warmup; i++) uniform();
}

#else
void RandomFast::reset(double rseed, int offset, int warmup)
{
  seed = static_cast<int> (fmod(rseed*IM+offset,IM));
  if (seed < 0) seed = -seed;
  if (seed == 0) seed = 1;
  for (int i = 0; i < warmup; i++) uniform();
}

void RandomFast::tagreset(double rseed, tagint offset, int warmup)
{
  seed = static_cast<int> (fmod(rseed*IM+offset,IM));
  if (seed < 0) seed = -seed;
  if (seed == 0) seed = 1;
  for (int i = 0; i < warmup; i++) uniform();
}

#endif


/* ----------------------------------------------------------------------
   integer RN between 1 and N inclusive
------------------------------------------------------------------------- */

int RandomFast::irandom(int n)
{
  int i = (int) (uniform()*(n-1)) + 1;
  //if (i > n) i = n;
  return i;
}


/* ----------------------------------------------------------------------
   uniform RN 
------------------------------------------------------------------------- */

double RandomFast::uniform()
{
  if(iter<numRand){
    iter++;
    return bulkRand[iter-1];
  }
  else
  {
    init_bulkRand();
    iter = 1;
    return bulkRand[0];
  }
  
}

/* ----------------------------------------------------------------------
   M values RN between 1 and N inclusive
------------------------------------------------------------------------- */
#ifdef AVX2
void RandomFast::init_bulkRand()
{
  int i = 0;
  while (i + 4 <= numRand) {
    __m256d rand_vals = prng.next_batch();
    _mm256_storeu_pd(&bulkRand[i], rand_vals);
    i += 4;
  }
}

#else
/* ----------------------------------------------------------------------
  Use Park RNG if AVX not enabled 
  uniform RN 
------------------------------------------------------------------------- */

double RandomFast::uniform_slow()
{
  int k = seed/IQ;
  seed = IA*(seed-k*IQ) - IR*k;
  if (seed < 0) seed += IM;
  double ans = AM*seed;
  return ans;
}

/* ----------------------------------------------------------------------
   M values RN between 1 and N inclusive
------------------------------------------------------------------------- */

void RandomFast::init_bulkRand()
{
  int i = 0;
  while (i<= numRand) {
    bulkRand[i] = uniform_slow();
    i ++;
  }
}
#endif


/* ----------------------------------------------------------------------
   tagint RN between 1 and N inclusive
------------------------------------------------------------------------- */

tagint RandomFast::tagrandom(tagint n)
{
  tagint i = (tagint) (uniform()*(n-1)) + 1;
  //if (i > n) i = n;
  return i;
}

/* ----------------------------------------------------------------------
   bigint RN between 1 and N inclusive
------------------------------------------------------------------------- */

bigint RandomFast::bigrandom(bigint n)
{
  bigint i = (bigint) (uniform()*(n-1)) + 1;
  //if (i > n) i = n;
  return i;
}
