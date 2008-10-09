//  -*- C++ -*-
//=======================================================================
// measure the egt sensitivity (influence) statistics for a given network
//=======================================================================
#include <boost/config.hpp>
#include <vector>
#include <string>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/wavefront.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/lambda/lambda.hpp>
#include <math.h>

#include "NetworkExperimentParameters.h"
typedef NetworkExperimentParameters ParamsClass;

#include "indicators-network.h"
#include "networks.h"
#include "search.h"
#include "sensitivity-simulation.h"
// undefine DISPLAY to make an executable without any of the display
// classes compiled in
#ifdef DISPLAY
#include "BoostDotDisplay.h"
#include "IndicatorsDisplay.h"
#else
template<typename arg_t>
class DisplayFarm{};
#endif//DISPLAY

using namespace boost;
using namespace boost::lambda;
using namespace std;

rng_t rng;

// ------------------------------------------------------------------------
//  main
// ------------------------------------------------------------------------

ParamsClass parameters;

typedef adjacency_list<setS,vecS,bidirectionalS> network_t;

int
main(int argc, char **argv)
{  
  // ===== initialize =====

  parameters.loadDefaults("settings/defaults-network.settings");
  parameters.handleArgs(argc,argv);
  parameters.afterSetting();

  rng.seed(parameters.randSeed());
  // rand() is used in random_shuffle()
  srand(parameters.randSeed());
  
  // ===== create a random or custom graph =====

  network_t n(parameters.n_vertices());
  // network created empty, now initialize it  
  //n.inheritParametersFrom(parameters);
  construct_network(n,parameters,rng);

  // ===== do the work =====

  sensitivity_simulation<network_t,rng_t,ParamsClass> ss(rng,&parameters);
  ss.do_simulation(n);
  return 0;
}
