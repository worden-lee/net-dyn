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
// typedef adjacency_list<setS,vecS,bidirectionalS> base_network_t;
// class network_with_params :
//   public base_network_t, public ObjectWithParameters<ParamsClass>
// {
// public:
//   network_with_params(base_network_t::vertices_size_type nv)
//     : base_network_t(nv) {}
//   network_with_params() {}
//   using base_network_t::clear;
// };

// string canonical(const network_with_params&n)
// { string can = canonical((base_network_t&)n);
//   const string *pvs = n.get("perturb_r_at_vertex");
//   if (pvs)
//     can += " " + *pvs;
//   return can;
// }

// typedef network_with_params network_t;

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

//   double moranprob =
//     (1 - 1/parameters.mutantFitness()) /
//     (1 - 1/pow(parameters.mutantFitness(),
//                (double)parameters.n_vertices()));

  // ===== do the work =====

  if (parameters.experiment() == "OPTIMIZE")
  { if (parameters.print_stuff())
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
  }
  else if (parameters.experiment() == "SAMPLE")
  {
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
      double fixp = sind(n);
      cout << "fixation probability: " << fixp << '\n';
      network_selection_indicator<network_t,rng_t,ParamsClass>::fixation_stats
        &cache = sind.provide_stats(n);
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
  }
  else if (parameters.experiment() == "SENSITIVITY")
    // sensitivity of fixation probability to fitness at each vertex
  { if (parameters.n_vertices() < 0) // skip it for now
    { typedef analytic_selection_indicator<network_t,ParamsClass> aind_t;
      aind_t aind(&parameters);
      typedef network_selection_indicator<network_t,rng_t&,ParamsClass> sind_t;
      sind_t sind(rng,&parameters);
      num_vertices_indicator<network_t> nvind;
      switching_indicator_t<aind_t,sind_t,num_vertices_indicator<network_t> >
        indicator(aind,nvind,parameters.max_n_vertices_for_matrix(),sind);
      // if this is requested, perturbed_mutant_fitness should be set
      double dr
        = parameters.perturbed_mutant_fitness() - parameters.mutantFitness();
      sind.setmonte_carlo_selection_accuracy(dr); // a good estimate
      double p = indicator(n);
      cout << "base probability of fixation: " << p << '\n';
      vector<double> sensitivities(num_vertices(n));
      cout << "sensitivities:\n";
      graph_traits<network_t>::vertex_iterator vi, vend;
      for (tie(vi,vend) = vertices(n); vi != vend; ++vi)
      { parameters.setperturb_r_at_vertex(*vi);
        sensitivities[*vi] = (indicator(n) - p) / dr;
        cout << sensitivities[*vi] << endl;
      }
    }
    else // more sophisticated monte carlo version
    { sensitivity_simulation<network_t,rng_t,ParamsClass> ss(rng,&parameters);
      ss.do_simulation(n);
    }
  }
  return 0;
}
