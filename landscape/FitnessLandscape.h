/* -*- c++ -*-
 */
#include <vector>
#include <sys/types.h>

template<class _point_t>
class FitnessLandscape
{
public:
  typedef _point_t point_t;
  FitnessLandscape() {}
  virtual double fitness(const point_t&x) = 0;
};
