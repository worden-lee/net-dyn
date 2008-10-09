//  -*- C++ -*-
//=======================================================================
// evaluate fixation probability for a given network
//=======================================================================
#include <boost/config.hpp>
#include <boost/tuple/tuple.hpp>
//#include <boost/graph/wavefront.hpp>
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
//#include "search.h"
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

  network_selection_indicator<network_t,rng_t,ParamsClass> sind(rng);
  sind.inheritParametersFrom(parameters);
  mean_fixation_time_indicator<typeof sind> ft_ind(sind);
  mean_extinction_time_indicator<typeof sind> et_ind(sind);
  mean_absorption_time_indicator<typeof sind> at_ind(sind);
  mean_peak_before_extinction_indicator<typeof sind> mp_ind(sind);

  if (parameters.print_stuff())
  { cout << "network:\n";
    print_object(n,cout);
    cout << "density: " << density(n) << endl;
    cout << "moran probability: "
         << moran_probability(parameters.mutantFitness()
                              /parameters.residentFitness(),
                              parameters.n_vertices()) << '\n';
  }
  double fixp = sind(n);
  network_selection_indicator<network_t,rng_t,ParamsClass>::fixation_stats
    &cache = sind.provide_stats(n);
  if (parameters.print_stuff())
  { cout << "fixation probability: " << fixp << '\n';
    cout << "n of trials: " << cache.n_trials << '\n';
    if (fixp > 0)
    {
      cout << "mean fixation time: "   << ft_ind(n) << '\n';
      cout << "mean extinction time: " << et_ind(n) << '\n';
      cout << "mean absorption time: " << at_ind(n) << '\n';
      cout << "mean peak before extinction: " << mp_ind(n) << '\n';
      cout << "biggest peak before extinction: "
           << cache.max_peakbeforeextinction << '\n';
    }
  }

  // record the results
  CSVDisplay network_csv(parameters.outputDirectory()+"/network.csv");
  // record other measures like mu_{-1} ?
  network_csv << "nv" << "in_exp" << "out_exp" << "density" << "mutuality"
              << "p" << "accuracy" << "moran_probability";
  network_csv.newRow();
  network_csv << parameters.n_vertices() << parameters.pl_in_exp()
              << parameters.pl_out_exp() << density(n) << mutuality_ind(n)
              << fixp << cache.fixation_accuracy()
              << moran_probability(parameters.mutantFitness()
                                   /parameters.residentFitness(),
                                   parameters.n_vertices());
  network_csv.newRow();

  return 0;
}
