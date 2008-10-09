/* -*- c++ -*-
 */

#include "FitnessLandscape.h"
#include "Point.h"
#include <openssl/des.h>

class BlockFitnessLandscape : public FitnessLandscape<Point>
{
public:
  BlockFitnessLandscape(string deskey);
  double blockFitness(int *blockp, int blockno);
  double fitness(const Point&x);
private:
  DES_key_schedule keyschedule; // internal form of encryption key
};
