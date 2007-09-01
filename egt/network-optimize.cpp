//=======================================================================
// how to optimize or improve a given network property
//=======================================================================
#include <boost/config.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <math.h>
#include "BoostDotDisplay.h"

using namespace boost;
using namespace std;

// ------ parameters for experiment -------

// start with a random graph with this many vertices and edges
const int init_n_vertices = 5;
const double init_density = 0.1;

// set these to control whether it displays graphs while computing
const bool animate_optimization = true;
const bool animate_fixation = false;
// if no changes improve the network, but some are neutral, take this
// many neutral steps before quitting
const int maxNeutralSteps = 100;

typedef adjacency_list<> network_t;
typedef graph_traits<network_t>::vertex_iterator vertex_iterator;
typedef graph_traits<network_t>::vertex_descriptor vertex_descriptor;
typedef graph_traits<network_t>::edge_descriptor edge_descriptor;
typedef graph_traits<network_t>::adjacency_iterator adjacency_iterator_t;

// global random number generator, used throughout
typedef boost::minstd_rand rng_t;
rng_t rng;

// use these types to actually get random numbers
typedef variate_generator<rng_t&, uniform_int<> >  ui_t;
typedef variate_generator<rng_t&, uniform_real<> > ur_t;

// you have to provide one of these
class indicator_t
{ public: virtual double operator()(const network_t&) const = 0;
};

// ------------------------------------------------------------------------
//  network search algorithms
// ------------------------------------------------------------------------

void print_network(network_t &n)
{
  vertex_iterator vi,vend;
  for (tie(vi,vend) = vertices(n); vi != vend; ++vi)
  { adjacency_iterator_t ai,aend;
    cout << *vi << ":";
    for (tie(ai,aend) = adjacent_vertices(*vi,n); ai != aend; ++ai)
      cout << ' ' << *ai;
    cout << endl;
  }
}

// little fn to turn an edge on/off
void toggle_edge(vertex_descriptor f, vertex_descriptor t,
		 network_t &n)
{ edge_descriptor ed;
  bool ed_found, zb;
  tie(ed,ed_found) = edge(f,t,n);
  if (ed_found)
    remove_edge(ed,n);
  else
    add_edge(f,t,n);
}

// improve the network until the indicator can't be improved any more.
// args:
//   n is the network
//   indicator is a function object that assigns a value to the network.
// returns: void, but modifies n.
void improveNetwork(network_t &n, const indicator_t &indicator)
{ static BoostDotGraphController<network_t> opt_animation;
  double currentIndicator = indicator(n);
  cout << "indicator: " << currentIndicator << '\n';
  
  bool foundImprovement = true;
  bool foundNeutralNeighbor = false;
  pair<vertex_descriptor,vertex_descriptor> neutral_neighbor;
  int neutralSteps = 0;
  while(foundImprovement || neutralSteps < maxNeutralSteps)
  {
//     cout << "current network\n";
//     print_network(n);
    opt_animation.update(n,"network");
    
    // create a randomly shuffled list of 'from' vertices
    vector<vertex_descriptor> from_vertices;
    vertex_iterator vb,ve;
    tie(vb,ve) = vertices(n);
    copy(vb,ve,back_inserter(from_vertices));
    random_shuffle(from_vertices.begin(),from_vertices.end());

    foundImprovement = false;
    // for each 'from' vertex
    for (vector<vertex_descriptor>::iterator fvi = from_vertices.begin();
	 !foundImprovement && fvi != from_vertices.end(); ++fvi)
    { // create a randomly shuffled list of 'to' vertices
      vector<vertex_descriptor> to_vertices;
      tie(vb,ve) = vertices(n);
      copy(vb,ve,back_inserter(to_vertices));
      random_shuffle(to_vertices.begin(),to_vertices.end());      
      // for each 'to' vertex
      for (vector<vertex_descriptor>::iterator tvi = to_vertices.begin();
	   tvi != to_vertices.end(); ++tvi)
      { toggle_edge(*fvi,*tvi,n);
	// see if the indicator is better
	double newIndicator = indicator(n);
	cout << "toggle (" << *fvi << ',' << *tvi << "): "
	     << newIndicator << endl;
	if (newIndicator > currentIndicator)
	{ currentIndicator = newIndicator;
	  foundImprovement = true;
	  cout << "found improvement\n";
	  break;
	}
	else if (newIndicator == currentIndicator && !foundNeutralNeighbor)
	{ neutral_neighbor = make_pair(*fvi,*tvi);
	  foundNeutralNeighbor = true;
	}
	// if we arrive here, we didn't find an improvement, put it back
	toggle_edge(*fvi,*tvi,n);
      }
    }
    if (foundImprovement)
      neutralSteps = 0;
    else if (foundNeutralNeighbor && neutralSteps++ < maxNeutralSteps)
      toggle_edge(neutral_neighbor.first,neutral_neighbor.second,n);
  }
  if (neutralSteps >= maxNeutralSteps)
    cout << "reached max neutral steps\n";
  else
    cout << "at local optimum\n";
  opt_animation.update(n,"network");
}

// ------------------------------------------------------------------------
//  indicators
// ------------------------------------------------------------------------

// density_indicator(d) is an indicator that evaluates how close your
//  network's density is to d.
// i.e.
// { density_indicator di(0.3);
//   double ind = di(n);
// }
// or
//  double ind = density_indicator(0.3)(n);
class density_indicator : public indicator_t
{
public:
  density_indicator(double _target) : target(_target) {}
  double operator()(const network_t&n) const
  { network_t::vertices_size_type nv = num_vertices(n);
    double dens = double(num_edges(n)) / (nv*nv);
    return 1 - fabs(dens-target);
  }
protected:
  double target;
};

class selection_indicator : public indicator_t
{
protected:
  const static int n_trials;
  const static double resident_fitness;
  const static double mutant_fitness;
public:
  double operator()(const network_t &n) const
  {
    // reusable state vector for fixation simulation
    vector<double>state(num_vertices(n));
    int n_fixations = 0;
    cout << "testing for fixation: ";
    for (int i = 0; i < n_trials; ++i)
    { bool fixed = fixation_once(state,n);
      cout << (fixed ? 'y':'n') << flush;
      if (fixed)
	++n_fixations;
    }
    cout << '\n';
    return n_fixations / double(n_trials);
  }
protected:
  // assign one random mutant in a resident population, simulate
  //  whether it fixates (within a given threshold of patience).
  bool fixation_once(vector<double>&state, const network_t &n) const
  { const int max_steps = 10000;
    network_t::vertices_size_type nv = num_vertices(n);
    if (nv <= 0)
      return false;
    // initial state is all residents but one
    state.resize(nv);
    state.assign(nv,resident_fitness);
    // generates random ints from [0,nv)
    ui_t choose_a_node(rng,uniform_int<>(0,nv-1));
    // generates random doubles from [0,1)
    ur_t choose_01(rng,uniform_real<>());
    state[choose_a_node()] = mutant_fitness;

    network_t::vertices_size_type n_infected = 1;
    int n_steps = 0;
    // do infection process until fixation reached
    while (0 < n_infected && n_infected < nv && n_steps++ < max_steps)
    { //cout << n_infected << ' ';
      vertex_descriptor parent, child;
      // trying to select a parent and child
      // select a parent vertex
      double accum =
	choose_01() * (nv*resident_fitness +
		n_infected*(mutant_fitness-resident_fitness));
      vertex_iterator vi,vend;
      for (tie(vi,vend) = vertices(n); vi != vend; ++vi)
	if ((accum -= state[*vi]) <= 0)
	  break;
      if (vi == vend)
      { cerr << "Error finding a parent!" << endl; continue; }
      parent = *vi;
      // now choose a child vertex
      network_t::degree_size_type n_children = out_degree(parent,n);
      adjacency_iterator_t ci,cend;
      if (n_children <= 0)
	continue; // no out edges, start over
      else if (n_children == 1)
      { tie(ci,cend) = adjacent_vertices(parent,n);
        child = *ci;
      }
      else
      { ui_t choose_a_child(rng,uniform_int<>(0,n_children-1));
	int child_selection = choose_a_child();
	for(tie(ci,cend) = adjacent_vertices(parent,n); ci != cend; ++ci)
	  if (child_selection-- == 0)
	    break;
	if (child_selection >= 0)
	{ cerr << "Error finding a child!" << endl; continue; } 
	child = *ci;
      }
      // if we get to here, we found a parent and child.
      //cout << parent << ',' << child << ' ';
      // now copy the parent's state to the child.
      if (state[child] == mutant_fitness)
      { if (state[parent] == resident_fitness)
	  --n_infected;
        state[child] = state[parent];
      }
      else
      { if (state[parent] == mutant_fitness)
	  ++n_infected;
        state[child] = state[parent];
      }
    }// while n_infected
    //cout << endl;
    return (n_infected == nv);
  }// end of function
};
const int selection_indicator::n_trials = 300;
const double selection_indicator::resident_fitness = 1;
const double selection_indicator::mutant_fitness = 1.1;

// ------------------------------------------------------------------------
//  main
// ------------------------------------------------------------------------

// Erdos-Renyi graph generator
typedef boost::erdos_renyi_iterator<boost::minstd_rand, network_t> ERGen;

int
main()
{  
  rng.seed((unsigned)time(0));
  
  // create network object with given size and density
  network_t n(ERGen(rng, init_n_vertices, init_density),
	      ERGen(), init_n_vertices);

  cout << "initial random network\n";
  print_network(n);

  // see what we can do to improve the indicator
  //improveNetwork(n,density_indicator(0.3342847893));
  improveNetwork(n,selection_indicator());

  cout << "final random network\n";
  print_network(n);
  cout << num_edges(n) << " edges" << endl;

  return 0;
}
