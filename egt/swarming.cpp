// try the swarm-fixation code
#include "SwarmingExperimentParameters.h"
typedef SwarmingExperimentParameters ParamsClass;

#include "swarming.h"
#include "indicators-swarm.h"
#include "IndicatorsDisplay.h"
#include "BoostDotDisplay.h"
#include "search.h"

ParamsClass parameters;

template<typename swp_t, typename ind_t, typename params_t, typename RNG_t>
void do_optimize(swp_t&swp,ind_t&ind,params_t&params,RNG_t&rng)
{ // can't do, except for lattice type
}

// it gets passed the swp, but it actually displays the graph that
// object contains
template<typename params_t>
class LatticeDisplayFarm
  : public DisplayFarm<static_lattice_swarming_pattern::graph_t, params_t>
{public:
  void update(double t, static_lattice_swarming_pattern&swp)
  { DisplayFarm<static_lattice_swarming_pattern::graph_t,params_t>
      ::update(t, swp.lattice.graph);
  }
};

// confusing: swp is a params object, when to use it vs. params?
template<typename ind_t, typename params_t, typename RNG_t>
void do_optimize(static_lattice_swarming_pattern&swp, ind_t&ind,
                 params_t&params, RNG_t&rng)
{ lattice_moves_generator<static_lattice_swarming_pattern::Dimension,RNG_t>
    mutator(rng);

  LatticeDisplayFarm<params_t> displayfarm;
  displayfarm.inheritParametersFrom(params);

  BoostDotGraphController<static_lattice_swarming_pattern::graph_t,
                          params_t,true>
    opt_animation;
  opt_animation.add_vertex_property("color",
    make_color_temperature(swp.lattice.graph));
  opt_animation.add_vertex_property("pos",
                                    //make_vertex_id(swp.lattice.graph,80));
//     vertex_location_map<static_lattice_swarming_pattern::Dimension,
//                         static_lattice_swarming_pattern::graph_t>
//                                     (swp.lattice.graph,80,2));
       place_vertex<static_lattice_swarming_pattern::Dimension,
                    static_lattice_swarming_pattern::graph_t>
                                    (swp.lattice.graph,80));

// this is silly for an undirected graph
//   opt_animation.add_edge_property("style",
//     make_bool_string_property("bold","",
//       belongs_to_a_2_loop<static_lattice_swarming_pattern::graph_t>,
//       swp.lattice.graph));

  if (swp.n_individuals() < 40)
  { displayfarm.installController(opt_animation);
    opt_animation.params.setdisplayToScreen(
      displayfarm.params.animate_optimization());
    //opt_animation.params.setdisplayEvery(50);
  }

  //hill_climb(swp,ind,mutator,&params,&displayfarm);
  swp.setannealing_rate_factor(swp.n_individuals());
  simulated_annealing(swp,ind,mutator,&params,rng,&displayfarm);
}

template<typename swp_t, typename sind_t, typename rng_t>
void do_experiment(swp_t&swp, sind_t&sind, rng_t&rng)
{
  constant_indicator mind(
    moran_probability(swp.mutantFitness()/swp.residentFitness(),
                      swp.n_individuals()));
  relative_indicator_t<sind_t, constant_indicator> rind(sind,mind);

  if (parameters.experiment() == "SWEEP")
  { IndicatorsDisplayController<swp_t,ParamsClass>
      ind_display;
    ind_display.inheritParametersFrom(swp);
    ind_display.setrecordEvery(0);
    ind_display.setdisplayEvery(0);
    ind_display.installIndicator(rind,"fixation/Moran prob");

//     //   cout << "fixation / Moran probability is " << rind(swp) << endl;
//     for (double rsr = 0.001; rsr <= 1000; rsr *= 1.1/*1.01*/)
//       //   for (swp.room_switch_rate = 1000; swp.room_switch_rate >= 0.001;
//       //        swp.room_switch_rate /= 1.01)
//     { double p_switch = rsr / (1 + rsr);
//       swp.setroom_switch_rate(rsr);
//       ind_display.update(p_switch,swp);
//     }

    for (double mp = 0.01; mp < 1; mp += 0.01)
    { swp.setlattice_move_probability(mp);
      ind_display.update(mp,swp);
    }
  }
  else if (parameters.experiment() == "SAMPLE")
  { cout << "swarming pattern:\n" << canonical(swp);
    double pfix = sind(swp);
    cout << "fixation probability: " << pfix << '\n';
    cout << "Moran probability: " << mind(swp) << '\n';
  }
  else if (parameters.experiment() == "OPTIMIZE")
  { cout << "Moran probability: " << mind(swp) << '\n';
    do_optimize(swp,sind,swp,rng);
  }
}

int main(int argc, char **argv)
{
  parameters.handleArgs(argc,argv,"settings/defaults-swarming.settings");
  parameters.afterSetting();

  rng_t rng;
  rng.seed(parameters.randSeed());
  // rand() used in random_shuffle()
  srand(parameters.randSeed());

  // set up the (initial) swarming pattern and do the experiment
  if (parameters.swarming_pattern() == "NETWORK")
  { graph_swarming_pattern swp;
    swp.inheritParametersFrom(parameters);
    const unsigned nr = swp.n_rooms();
    swp.fixed.set_size(nr);
    swp.n = graph_swarming_pattern::network_t(nr);
    construct_network(swp.n,swp,rng);
    int fixed = (swp.n_individuals() - swp.n_mobile());
    int fpr = fixed / nr;
    for (int i = 0; i < nr; ++i)
    { swp.fixed[i] = fpr;
      fixed -= fpr;
    }
    if (fixed > 0)
      swp.fixed[0] += fixed;
    swarm_selection_indicator<graph_swarming_pattern,rng_t,ParamsClass>
      sind(rng);
    do_experiment(swp,sind,rng);
  }
  else if (parameters.swarming_pattern() == "STATIC_LATTICE")
  { static_lattice_swarming_pattern swp;
    swp.inheritParametersFrom(parameters);
    swp.initialize(rng);
    static_lattice_swarm_selection_indicator<rng_t,ParamsClass>
      sind(rng);
    do_experiment(swp,sind,rng);
  }
  else if (parameters.swarming_pattern() == "DYNAMIC_LATTICE")
  { dynamic_lattice_swarming_pattern dsp;
    dsp.inheritParametersFrom(parameters);
    dynamic_lattice_swarm_selection_indicator<rng_t,ParamsClass>
      sind(rng);
    do_experiment(dsp,sind,rng);
  }
}
