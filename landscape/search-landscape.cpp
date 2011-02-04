/* -*- c++ -*- 
 *
 */
#include "BlockFitnessLandscape.h"
#include "LParameters.h"
#define SKIP_SIM_ANNEAL
#include "../network/search.h"
#include <boost/random/linear_congruential.hpp>
#include <boost/random/variate_generator.hpp>

template<class fitness_landscape_t>
class fitnessindicator
{ fitness_landscape_t&ls;
public:
  fitnessindicator(fitness_landscape_t&_ls) : ls(_ls) {}
  double operator()(typename fitness_landscape_t::point_t&p)
  { return ls.fitness(p); }
  friend bool indicator_is_stochastic(fitnessindicator&i) { return false; }
};

template<typename params_t>
void write_snapshot(Point&p, double t, params_t&params)
{}

int main(int argc, char**argv)
{
  // ===== initialize =====

  LParameters parameters;
  parameters.handleArgs(argc,argv,"settings/default.settings");
  parameters.afterSetting();

  boost::minstd_rand rng;
  rng.seed(parameters.randSeed());
  // rand() is used in random_shuffle()
  srand(parameters.randSeed());

  // ===== create landscape, do statistics =====

  BlockFitnessLandscape landscape(parameters.desKey());
  fitnessindicator<BlockFitnessLandscape> indicator(landscape);

  Point p;
  p.randomize(rng);
  PointMutationGenerator mutagen;
  hill_climb(p, indicator, mutagen, &parameters, &parameters);
}
