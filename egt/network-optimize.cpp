//  -*- C++ -*-
//=======================================================================
// how to optimize or improve a given network property
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

  parameters.handleArgs(argc,argv,"settings/defaults-network.settings");
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

  if (parameters.print_stuff())
  { cout << "moran probability: "
         << moran_probability(parameters.mutantFitness()
                              /parameters.residentFitness(),
                              parameters.n_vertices()) << '\n';
    cout << "initial network\n";
    print_object(n,cout);
    cout << "initial density: " << density(n) << endl;
  }
  
  // adjust it to (near) specified density
  //improveNetwork(n,density_indicator(0.3342847893));

  //   analytic_selection_indicator_fixed_origin<network_t> aind(1);
  //   selection_indicator_fixed_origin<network_t,rng_t> sind(1,rng);

  //   deterministic_indicator_caching<network_t,double,
  //bad_at_selection_but_strongly_connected_indicator<network_t> ind;
  //    deterministic_indicator_caching<network_t,double,
  //      probability_to_amplification< analytic_selection_indicator<network_t> > >
  //        pind;
  typedef network_selection_indicator<network_t,rng_t,ParamsClass> sind_t;
  sind_t sind(rng,&parameters);

  typedef analytic_selection_indicator<network_t,ParamsClass> aind_t;
  aind_t aind(&parameters);

  typedef deterministic_indicator_caching<network_t,double,aind_t> pind_t;
  pind_t pind(aind);

  typedef moran_indicator mind_t;
  mind_t mind(parameters.mutantFitness()/parameters.residentFitness());

  //   relative_indicator_t<pind_t,constant_indicator> rind(pind,mind);
  //   relative_indicator_t< network_selection_indicator<network_t,rng_t>,
  //     constant_indicator > srind(sind,mind);

  typedef comparing_indicator< pind_t, sind_t > cind_t;
  cind_t cind(pind,sind);

  //     typedef relative_indicator_t<pind_t, mind_t> rpind_t;
  //     rpind_t rpind(pind,mind);
  // temp: no caching, it interferes with influence
  typedef relative_indicator_t<aind_t, mind_t> rpind_t;
  rpind_t rpind(aind,mind);

  typedef relative_indicator_t<cind_t, mind_t> rcind_t;
  rcind_t rcind(cind,mind);

  typedef relative_indicator_t<sind_t, mind_t> rsind_t;
  rsind_t rsind(sind,mind);

  num_vertices_indicator<network_t> nvind;
  switching_indicator_t<rpind_t,rsind_t,num_vertices_indicator<network_t> >
    swind(rpind, nvind, (double)parameters.max_n_vertices_for_matrix(),
          rsind);

  //double(*ind)(const network_t&) = zero_indicator;

  // here, for convenience, which indicator we are to optimize
  generic_indicator<network_t> sg(rsind);
  generic_indicator<network_t> pg(rpind);
  bool use_montecarlo = true;
  //  = parameters.n_vertices() > parameters.max_n_vertices_for_matrix();
  generic_indicator<network_t> swg(swind);

  generic_indicator<network_t> &indicator
    = swg;
  //= (use_montecarlo ? sg : pg);
  //parameters.setindicator_is_stochastic(use_montecarlo);
  //     if (parameters.trialsForLocalOptimum() == 0)
  //       parameters.settrialsForLocalOptimum(use_montecarlo ? 200 : 1);

#ifdef DISPLAY
  DisplayFarm<network_t,ParamsClass> displayfarm;
  displayfarm.inheritParametersFrom(parameters);

  BoostDotGraphController<network_t,ParamsClass> opt_animation;
  if (parameters.n_vertices() < 40)
  { displayfarm.installController(opt_animation);
    //       opt_animation.add_vertex_property("color",
    //         make_color_temperature(n));
    //       opt_animation.add_vertex_property("fontcolor",
    //         make_color_influence(n,indicator,parameters));
    opt_animation.add_vertex_property("color",
      make_color_influence(n,indicator,parameters));
    opt_animation.add_edge_property("style",
      make_bool_string_property("bold","",belongs_to_a_2_loop<network_t>,n));
    opt_animation.params.setdisplayToScreen(
      displayfarm.params.animate_optimization());
  }

  IndicatorsDisplayController<network_t,ParamsClass> ind_display;
  displayfarm.installController(ind_display);
  ind_display.params.setdisplayToScreen(
    displayfarm.params.animate_optimization());
  ind_display.installIndicator(indicator,"indicator");
  ind_display.installIndicator(&density<network_t>,"density");
  ind_display.installIndicator(&mutuality_ind<network_t>,"mutuality");
  ind_display.installIndicator(&transitivity_ind<network_t>,"transitivity");
  ind_display.installIndicator(&median_temperature<network_t>,
                               "median temperature");
  ind_display.installIndicator(&median_out_density<network_t>,
                               "median out density");
  ind_display.installIndicator(&median_in_density<network_t>,
                               "median in density");
  ind_display.installIndicator(&relative_diameter_ind<network_t>,
                               "relative diameter");
  ind_display.installIndicator(&relative_mean_path_length<network_t>,
                               "relative mean path length");
  //   ind_display.installIndicator(&relative_mean_wavefront_ind<network_t>,
  // 			       "relative mean wavefront");
  ind_display.installIndicator(&central_point_dominance_ind<network_t>,
                               "central point dominance");

  IndicatorsDisplayController<network_t,ParamsClass> mft_display;
  displayfarm.installController(mft_display);
  mft_display.params.setdisplayToScreen(
    displayfarm.params.animate_optimization());
  mean_fixation_time_indicator<typeof sind> ft_sind(sind);
  mft_display.installIndicator(ft_sind,"mean fixation time");
  mean_extinction_time_indicator<typeof sind> et_sind(sind);
  mft_display.installIndicator(et_sind,"mean extinction time");
  mean_absorption_time_indicator<typeof sind> at_sind(sind);
  mft_display.installIndicator(at_sind,"mean absorption time");
  expected_absorption_time_indicator<typeof aind> at_aind(aind);
  //     mft_display.installIndicator(at_aind,"expected absorption time");
  
  //improveNetwork(n,indicator,&parameters,&displayfarm);
  network_mutation_generator<network_t,rng_t,ParamsClass> mutator(rng);
  mutator.inheritParametersFrom(parameters);
  //    hill_climb(n,indicator,mutator,&parameters,&displayfarm);
  parameters.setannealing_rate_factor(num_vertices(n)*(num_vertices(n)-1));
  simulated_annealing(n,indicator,mutator,&parameters,rng,&displayfarm);
#else
  network_mutation_generator<network_t,rng_t> mutator(rng);
  //    hill_climb(n,indicator,mutator,&parameters,&displayfarm);
  simulated_annealing(n,indicator,mutator,&parameters,rng);
#endif//DISPLAY

  if (parameters.print_stuff())
  { cout << "final network\n";
    print_object(n,cout);
    //  cout << num_edges(n) << " edges" << endl;
    cout << "density " << density(n) << endl;
  }

  return 0;
}
