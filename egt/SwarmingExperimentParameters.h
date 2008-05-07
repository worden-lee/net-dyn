// -*- C++ -*-
#ifndef __swarming_experiment_parameters_h__
#define __swarming_experiment_parameters_h__
#include "Parameters.h"
#include <boost/random.hpp>

class SwarmingExperimentParameters : public Parameters
{
public:
  // ===== declarations of parameters =====

  //  --- global stuff

  // what kind of swarming pattern
  //  NETWORK (the rooms model), STATIC_LATTICE, DYNAMIC_LATTICE
  DECLARE_PARAM(string,swarming_pattern)

  //  --- for the general swarming process

  // how many people
  DECLARE_PARAM(unsigned, n_individuals)

  //  --- for the rooms swarming process

  // use a graph of rooms with this many vertices
  DECLARE_PARAM(unsigned, n_rooms)

  // how many of them are mobile
  DECLARE_PARAM(unsigned, n_mobile)

  // this controls the rate of room switching relative to
  // the rate of infection.
  DECLARE_PARAM(double, room_switch_rate)

  // if true, like it says, infected mobiles move more often
  //  than uninfected ones
  DECLARE_PARAM(bool, mobilityIncreasesWithFitness)

  //  --- for lattice swarming

  // # of dimensions of the lattice
  // there's a hard-coded maximum, which is 10 as of this writing
  DECLARE_PARAM(unsigned, lattice_dimensions)

  // size on each dimension of the lattice
  // lattice_dim0, lattice_dim1, etc can be used, if not found this
  // one is used
  DECLARE_PARAM(unsigned, lattice_dim_generic)

  // how close is considered a neighbor on the lattice
  DECLARE_PARAM(unsigned, lattice_neighbor_distance)

  // whether coordinates wrap around edges
  DECLARE_PARAM(bool, lattice_is_torus)

  //  --- for dynamic lattice swarming

  // prob. a given event will be a move rather than an infection
  DECLARE_PARAM(double, lattice_move_probability)

  //  --- for the optimizer

  // to control whether it displays moving graphics of the optimization
  // process
  DECLARE_PARAM(bool,animate_optimization)

  // if this isn't defined, the Gnuplot and Graphviz interfaces won't
  // even be compiled into the program.  useful for running batches on
  // clusters.
#define DISPLAY 1

};

// -----------------------------------------

// random number generator, used throughout
typedef boost::minstd_rand rng_t;
//extern rng_t rng;

// use these types to actually get random numbers
typedef boost::variate_generator<rng_t&, boost::uniform_int<> >  ui_t;
typedef boost::variate_generator<rng_t&, boost::uniform_real<> > ur_t;

#endif // __network_experiment_parameters_h__
