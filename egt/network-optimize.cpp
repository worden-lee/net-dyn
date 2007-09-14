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
#include "indicators.h"

using namespace boost;
using namespace std;

// setS is to disallow parallel edges
typedef adjacency_list<setS,vecS,bidirectionalS> network_t;
typedef graph_traits<network_t>::vertex_iterator    vertex_iterator;
typedef graph_traits<network_t>::vertex_descriptor  vertex_descriptor;
typedef graph_traits<network_t>::edge_iterator      edge_iterator;
typedef graph_traits<network_t>::edge_descriptor    edge_descriptor;
typedef graph_traits<network_t>::adjacency_iterator adjacency_iterator_t;

rng_t rng;

// ------------------------------------------------------------------------
//  network search algorithms
// ------------------------------------------------------------------------

template<typename network_t>
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

typedef pair<vertex_descriptor,vertex_descriptor> possible_edge;

// little fn to turn an edge on/off
void toggle_edge(possible_edge ft, network_t &n)
{ edge_descriptor ed;
  bool ed_found, ed_added;
  tie(ed,ed_found) = edge(ft.first,ft.second,n);
  if (ed_found)
    remove_edge(ed,n);
  else
    tie(ed,ed_added) = add_edge(ft.first,ft.second,n);
}

// improve the network until the indicator can't be improved any more.
// args:
//   n is the network
//   indicator is a function object that assigns a value to the network.
// returns: void, but modifies n.
template<typename indicator_t>
void improveNetwork(network_t &n, indicator_t &indicator)
{ static BoostDotGraphController<network_t>
    opt_animation(animate_optimization?1:0, 1);
  int nSteps = 0; 
  int consecutiveTimesFoundLocalOptimum = 0;
  double currentIndicator;

  // create a list of all edges to toggle on/off
  typename network_t::vertices_size_type nv = num_vertices(n);
  vector<possible_edge> edges_to_try;
  edges_to_try.reserve(nv*nv);
  vertex_iterator fb,fe,tb,te;
  for (tie(fb,fe) = vertices(n); fb!=fe; ++fb)
    for (tie(tb,te) = vertices(n); tb!=te; ++tb)
      if (*tb!=*fb)
	edges_to_try.push_back(make_pair(*fb,*tb));
  
  while (consecutiveTimesFoundLocalOptimum < trialsForLocalOptimum)
  { currentIndicator = indicator(n);
    if (print_stuff)
      cout << "indicator: " << currentIndicator << '\n';
  
    bool foundImprovement = true;
    bool foundNeutralNeighbor = false;
    possible_edge neutral_neighbor;
    int neutralSteps = 0;
    while(foundImprovement ||
	  (foundNeutralNeighbor && neutralSteps < maxNeutralSteps))
    {
//     cout << "current network\n";
//     print_network(n);
      //if (animate_optimization)
	opt_animation.update(nSteps,n,"network");
      
      // shuffle the candidate edges.  this matters because the first
      //  one found that's an improvement will be used, or if none,
      //  the first one found that's neutral.
      random_shuffle(edges_to_try.begin(),edges_to_try.end());
    
      foundImprovement = foundNeutralNeighbor = false;
      // for each possible edge
      for (vector<possible_edge>::iterator ei = edges_to_try.begin();
	   !foundImprovement && ei != edges_to_try.end(); ++ei)
      { edge_descriptor ed;
        bool ed_found, ed_added;
	tie(ed,ed_found) = edge(ei->first,ei->second,n);
	if (ed_found)
	  remove_edge(ed,n);
	else
	  tie(ed,ed_added) = add_edge(ei->first,ei->second,n);
	//toggle_edge(*ei,n);
        // see if the indicator is better
        double newIndicator = indicator(n);
	// for stochastic indicators: if it's big, call it again
	//  and make sure it's really big
	// (with the caching I've got in there, the more you call it
	//  the more accurate it is)
	if (indicator_is_stochastic &&
	    newIndicator >= currentIndicator && newIndicator > 0)
	{ newIndicator = indicator(n);
	} 
// 	cout << "toggle (" << ei->first << ',' << ei->second << "): "
// 	     << newIndicator << endl;
	if (newIndicator >= currentIndicator &&
	    newIndicator < currentIndicator + 1e-7)
	{ if (!foundNeutralNeighbor)
  	  { neutral_neighbor = *ei;
 	    foundNeutralNeighbor = true;
	  }
	}
	else if (newIndicator > currentIndicator)
	{ currentIndicator = newIndicator;
          foundImprovement = true;
	  if (print_stuff)
	  { cout << "found improvement " << (ed_found?'-':'+') << '('
		 << ei->first << ',' << ei->second << "): "
		 << currentIndicator << '\n';
  	  //print_network(n);
	  }
	  break;
	}

	// if we arrive here, we didn't find an improvement, put it back
	if (ed_found)
	  add_edge(ei->first,ei->second,n);
	else if (ed_added)
	  remove_edge(ed,n);
	//toggle_edge(*ei,n);
      }
      if (foundImprovement)
      { ++nSteps;
	neutralSteps = 0;
	consecutiveTimesFoundLocalOptimum = 0;
      }
      else if (foundNeutralNeighbor && neutralSteps++ < maxNeutralSteps)
      { if (print_stuff)
  	  cout << "taking neutral step ("
	       << neutral_neighbor.first << ','
	       << neutral_neighbor.second << ")\n";      
        toggle_edge(neutral_neighbor,n);
	if (print_stuff)
	{ cout << "indicator: " << currentIndicator << '\n';
	  //print_network(n);
	}
        ++nSteps;
	consecutiveTimesFoundLocalOptimum = 0;
      }
    }
    if (neutralSteps >= maxNeutralSteps)
    { if (print_stuff)
	cout << "reached max neutral steps\n";
      break;
    }
    else
    { if (print_stuff)
        cout << "at local optimum\n";
      ++consecutiveTimesFoundLocalOptimum;
    }
    //if (animate_optimization)
      opt_animation.update(nSteps,n,"network");
  }
  if (print_stuff)
    cout << "indicator: " << currentIndicator << '\n';
}

// ------------------------------------------------------------------------
//  main
// ------------------------------------------------------------------------

// Erdos-Renyi graph generator
typedef boost::erdos_renyi_iterator<boost::minstd_rand, network_t> ERGen;

int
main()
{  
  rng.seed((unsigned)time(0));
  // rand() used in random_shuffle()
  srand((unsigned)time(0));
  
  ur_t choose_01(rng,uniform_real<>());
  double random_density = choose_01(); // rather than fixed init. dens.
  
  // create network object with given size and density
  network_t n(ERGen(rng, init_n_vertices, random_density),
	      ERGen(), init_n_vertices);

  if (print_stuff)
  { cout << "initial random network\n";
    print_network(n);
    cout << "initial density: " << random_density << endl;
  }
  
  // see what we can do to improve the indicator
  //improveNetwork(n,density_indicator(0.3342847893));
  //selection_indicator<network_t,rng_t> ind(rng);
  //selection_indicator_fixed_origin sfi(1);
  //comparing_indicator<network_t,rng_t> ind(rng);
  analytic_selection_indicator<network_t> ind;
  //bad_at_selection_but_strongly_connected_indicator<network_t> ind;
  improveNetwork(n,ind);

  if (print_stuff)
  { cout << "final network\n";
    print_network(n);
    //  cout << num_edges(n) << " edges" << endl;
    cout << "density " << density(n) << endl;
  }
  return 0;
}
