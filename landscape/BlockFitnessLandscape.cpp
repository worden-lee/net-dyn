#include "BlockFitnessLandscape.h"

/* BlockFitnessLandscape produces uncorrelated fitness function for each
   block by using DES encryption
*/
BlockFitnessLandscape::BlockFitnessLandscape(string deskey)
{
  des_cblock keyblock;
  des_string_to_key(deskey.c_str(),&keyblock);
  des_check_key = true;
  if ( des_set_key(&keyblock,keyschedule) )
    cerr << "DES Error setting key!!\n";
}

/* returns fitness from 0 to maxFecundity for each block
   thus these values should be averaged, not summed
*/
double BlockFitnessLandscape::blockFitness(int *block, int blockno) 
{
  des_cblock outblock;
  // ivec is a kind of seed for the encryption
  //  include blockno so each block's landscape is different
  des_cblock ivec = {0xfe,0xdc,0xba,0x98,0x76,0x54,0x32,0x10+blockno};
  // compute des checksum from the block value
  des_cbc_cksum((unsigned char*)block, &outblock, (Point::blockSize+7)/8,
		keyschedule, &ivec);
  // reduce checksum into a single int using xor
  // this code will be in trouble if cblock stops being the size of 2 ints
  unsigned int *cksump = (unsigned int*)&outblock;
  unsigned int cksum = cksump[0] ^ cksump[1];

  static double factor = 1 / (1 + (double)UINT_MAX);
  return cksum * factor;
}

double BlockFitnessLandscape::fitness(const Point&x)
{
  double fitness = 0;
  for ( int i = 0; i < Point::nBlocks; i++ )
    fitness += blockFitness((int*)&x.genome[i], i);
  fitness /= Point::nBlocks;
  return fitness;
}
