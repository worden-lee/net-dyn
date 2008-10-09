#ifndef POINT_H
#define POINT_H

#include "LParameters.h"
#include <limits.h>
#include <iostream>
#include <sstream>

#include <boost/random/linear_congruential.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_int.hpp>

/*
 * Point class embodies a list of equal-size blocks of bits
 *
 * The BlockFitnessLandscape class is well prepared to evaluate fitness
 *  of these objects
 *
 * WARNING:  Make sure there aren't any stray bits in the padding
 *  section of the block structs, as they may get into the fitness
 *  function, which uses DES encryption with a granularity of (I think)
 *  bytes
 */

// -----------------------------------------------------------
// NBLOCKS and BLOCKSIZE should be the only things to change,
// to alter the landscape structure.
//  multiple blocks == correlation.
//  BLOCKSIZE is number of bits per block (not bytes). Can be large.

// -----------------------------------------------------------

#if UINT_MAX == 4294967295U
#  define BITS_PER_WORD 32
#  define BLOCK_WORD_TYPE unsigned int
#else
#  error Point.h has to be rewritten for different size int!
#endif

#define WORDS_PER_BLOCK (BLOCKSIZE / BITS_PER_WORD)
#define N_EXTRA_BITS (BLOCKSIZE - WORDS_PER_BLOCK * BITS_PER_WORD)

#if N_EXTRA_BITS > 0
#  define USE_EXTRA_BITS 1
#else
#  define USE_EXTRA_BITS 0
#endif

class Point
{
public:
  // class stuff
  const static int nBlocks = NBLOCKS, blockSize = BLOCKSIZE;
  const static int wordsPerBlock = WORDS_PER_BLOCK;
  const static int bitsPerWord = BITS_PER_WORD;
  const static int nExtraBits = N_EXTRA_BITS;
  typedef struct block_str
  { BLOCK_WORD_TYPE words[WORDS_PER_BLOCK]; 
#if USE_EXTRA_BITS
    BLOCK_WORD_TYPE extraBits:N_EXTRA_BITS;
#endif
  } block;

  // instance stuff
  block genome[NBLOCKS];

  // functions
  Point(void);
    // arguments must be ((wordsPerBlock + (useExtraBits?1:0)) * nBlocks)
    //  number of BLOCK_WORD_TYPEs
  Point(BLOCK_WORD_TYPE b0, ...);

  unsigned int hammingDistance(const Point&) const;

  // mutate *this at random
  template<typename rng_t>
  void mutate(rng_t&rng);
  // mutate at bit 'i'
  void mutate(int i);
  // mutate at block 'bl', bit 'bi'
  void mutate(int bl, int bi);
  // store mutation of *this in *destination
  template<typename rng_t>
  void mutate(Point *destination, rng_t&rng) const;
  void mutate(Point *destination, int bl, int bi) const;

  template<typename rng_t>
  void randomize(rng_t&rng);

  static Point *wildType(void);

  Point &operator++(void);

  Point &operator=(const Point&other);

  bool operator == (const Point&other) const;
  bool operator != (const Point&other) const { return !(*this==other); }
  bool operator < (const Point&other) const; // for sorting

  const char *hexString(void) const;
  friend ostream& operator<< (ostream &o, const Point &comm);
  friend void print_object(const Point&p, ostream&o)
  { o << p; }
 protected:
};

template<typename rng_t>
void Point::mutate(rng_t&rng)
{
  mutate(*this,rng);
}

template<typename rng_t>
void Point::mutate(Point *mutant, rng_t&rng) const
{
  typedef boost::variate_generator<rng_t&, boost::uniform_int<BLOCK_WORD_TYPE> > ui_t;
  ui_t ui_bl(rng,boost::uniform_int<BLOCK_WORD_TYPE>(0,Point::nBlocks));
  ui_t ui_bi(rng,boost::uniform_int<BLOCK_WORD_TYPE>(0,Point::blockSize));
  int mbl = ui_bl();
  int mbi = ui_bi();
  mutate(mutant, mbl, mbi);
}

template<typename rng_t>
void Point::randomize(rng_t&rng)
{
  typedef boost::variate_generator<rng_t&, boost::uniform_int<BLOCK_WORD_TYPE> > ui_t;
  ui_t ui_word(rng,boost::uniform_int<BLOCK_WORD_TYPE>(0,(1UL<<Point::bitsPerWord-1)));
#if USE_EXTRA_BITS
  // careful with this 1U if sizeof extraBits gets big
  ui_t ui_extra(rng,boost::uniform_int<BLOCK_WORD_TYPE>(0,(1UL<<Point::nExtraBits-1)));
#endif
    // initialize random genotype
  for (int j = 0; j < Point::nBlocks; j++)
  {
    for ( int k = 0; k < Point::wordsPerBlock; k++ )
      genome[j].words[k] = ui_word();
#if USE_EXTRA_BITS
    genome[j].extraBits = ui_extra();
#endif
  }
}

// === mutations for points ===

class PointMutationGenerator
{ 
public:
  class mutation
  { 
    int i;
  public:
    mutation(int _i=0) : i(_i) {}
    void operator()(Point&p) const
    { p.mutate(i); }
    string asString() const
    { ostringstream os;
      os << i;
      return os.str();
    }
    static mutation &undo(mutation&_m)
    { return _m; }
  };
  
  class mutation_iterator
  { int i;
  public:
    mutation_iterator(int _i) : i(_i) {}
    mutation operator*() const
    { return mutation(i); }
    mutation_iterator&operator++()
    { ++i; return *this; }
    bool operator!=(const mutation_iterator&other) const
    { return i != other.i; }
  };
  PointMutationGenerator() {}
  mutation_iterator begin(Point&p)
  { return mutation_iterator(0); }
  mutation_iterator end()
  { return mutation_iterator(NBLOCKS*BLOCKSIZE); }
};

#endif//POINT_H
