/* -*- c++ -*-
 */
#include "BooleanLandscape.h"

template<unsigned D, unsigned B>
class BlockFitness
{
protected:

public:
  BlockFitness() {}
  typedef BooleanLandscape<D>::Point point_t;
  double operator()(const point_t &x)
  {}
  class Block
  {
  protected:
    unsigned bit0, bitn; // will work with bit0 .. (bitn-1)
    // cryptographic key/salt/what have you to go here
  public:
    // Block(init_data)
    double operator()(const BlockFitness<D,B>::point_t)
    { // do it
    }
  };

protected:
  // this is local storage, not the list of what blocks we use
  // some may be stored here, some may be stored by other objects
  vector<Block> block_storage;
  vector<Block*> blocks;
};
