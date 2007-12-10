// -*- C++ -*-
#ifndef __network_optimize_h__
#define __network_optimize_h__

#include <boost/random.hpp>

// ------ parameters for experiment -------

//  --- global stuff

// whether to print information to the terminal
bool print_stuff = true;

//  --- for the main function

// start with a random graph with this many vertices and edges
int init_n_vertices = 13;//5;//12;//15
double init_density = 0;

//  --- for the optimizer

// to control whether it displays the graphs while optimizing
bool animate_optimization = true;//false;
#define DISPLAY 1

// whether the 'fitness' needs to be evaluated multiple times
// simpler if not
bool indicator_is_stochastic = false;//true

// if no changes improve the network, but some are neutral, take this
// many neutral steps before quitting
int maxNeutralSteps = 10000;

// when we think we have an optimum, we try this many times before
//  we're sure.  >1 is only useful when the indicator is stochastic.
int trialsForLocalOptimum = indicator_is_stochastic ? 200 : 1;//10;

//  --- for the fixation process

// if true, the contact process is:
//   choose a parent weighted by fitness
//   choose a child among the parent's out-neighbors, weighed by edge weight
//   (this is the evolutionary graph theory process)
// if false, instead it's:
//   choose a child at random
//   choose a parent among the child's in-neighbors, weighted by
//     fitness and edge weight
// (though so far, edge weights are assumed to be 1)
bool choose_parent_first = true;//false;

// fitness at the 'normal' vertices
double residentFitness = 1;

// fitness of the mutant type, which might go to fixation
double mutantFitness = 1.2;//1.1;
//double mutantFitness = 0.9; // less than resident to select for drift

// this animation isn't implemented
bool animate_fixation = false;

// simulate introduction and possible fixation this many times at a
//  go, to estimate probability of fixation
// but note, these results are stored and more trials are added in at
//  later calls, so this is only the minimum sample size.
int timesToTestFixation = 300;

// if no fixation after this many reproduction events, quit trying
int maxStepsToFixation = 1000;

//  --- for the swarm process

bool mobilityIncreasesWithFitness = true;

// -----------------------------------------

#if 0
// this is like the plain network object, but it's a regular network
// plus a diff
template<typename network_t>
class variation
{
public:
  typedef typename graph_traits<network_t>::vertex_descriptor vertex;
  typedef pair<vertex,vertex> possible_edge_t;
  typedef enum{ subtract, add } edge_operation_t;
  typedef pair<edge_operation_t, possible_edge_t> diff_t;
  network_t &n;
  diff_t _diff;
  variation(network_t&_n,diff_t _d) : n(_n), _diff(_d)
  {}
  network_t &network()
  { return n;
  }
  diff_t &diff()
  { return _diff;
  }  
};
#endif

// -----------------------------------------

// global random number generator, used throughout
typedef boost::minstd_rand rng_t;
extern rng_t rng;

// use these types to actually get random numbers
typedef boost::variate_generator<rng_t&, boost::uniform_int<> >  ui_t;
typedef boost::variate_generator<rng_t&, boost::uniform_real<> > ur_t;

#endif // __network_optimize_h__
