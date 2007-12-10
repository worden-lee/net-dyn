//=======================================================================
// how to optimize or improve a given network property
//=======================================================================
#include <boost/config.hpp>
#include <iostream>
#include <vector>
#include <deque>
#include <string>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/wavefront.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/graph/reverse_graph.hpp>
#include <boost/lambda/lambda.hpp>
#include <math.h>
#include "BoostDotDisplay.h"
#include "IndicatorsDisplay.h"
#include "indicators-network.h"

using namespace boost;
using namespace boost::lambda;
using namespace std;

rng_t rng;

template<typename network_t>
struct custom_types
{
  typedef pair< typename graph_traits<network_t>::vertex_descriptor,
		typename graph_traits<network_t>::vertex_descriptor >
    possible_edge;

  // for creating random graphs
  typedef boost::erdos_renyi_iterator<boost::minstd_rand, network_t>
    ER_graph_generator;
};

// ------------------------------------------------------------------------
//  network search algorithms
// ------------------------------------------------------------------------

template<typename network_t>
void print_network(network_t &n)
{
  typename graph_traits<network_t>::vertex_iterator ii,iend,ji,jend;
//   for (tie(vi,vend) = vertices(n); vi != vend; ++vi)
//   { typename graph_traits<network_t>::adjacency_iterator ai,aend;
//     cout << *vi << ":";
//     for (tie(ai,aend) = adjacent_vertices(*vi,n); ai != aend; ++ai)
//       cout << ' ' << *ai;
//     cout << endl;
//   }
  for (tie(ii,iend) = vertices(n); ii != iend; ++ii)
  { for (tie(ji,jend) = vertices(n); ji != jend; ++ji)
    { typename graph_traits<network_t>::edge_descriptor ed;
      bool ed_found;
      tie(ed,ed_found) = edge(*ii,*ji,n);
      cout << (ed_found ? '1':'0') << ' ';
    }
    cout << endl;
  }
}

// little fn to turn an edge on/off
template<typename network_t>
void toggle_edge(typename custom_types<network_t>::possible_edge ft,
		 network_t &n)
{ typename graph_traits<network_t>::edge_descriptor ed;
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
template<typename network_t, typename indicator_t>
void improveNetwork(network_t &n, indicator_t &indicator,
                    DisplayFarm<network_t>*df=0)
{
  int nSteps = 0; 
  int consecutiveTimesFoundLocalOptimum = 0;
  double currentIndicator;

  // create a list of all edges to toggle on/off
  typename network_t::vertices_size_type nv = num_vertices(n);
  typedef vector< typename custom_types<network_t>::possible_edge >
    edges_to_try_t;
  edges_to_try_t edges_to_try;
  edges_to_try.reserve(nv*nv);
  typename graph_traits<network_t>::vertex_iterator fb,fe,tb,te;
  for (tie(fb,fe) = vertices(n); fb!=fe; ++fb)
    for (tie(tb,te) = vertices(n); tb!=te; ++tb)
      if (*tb!=*fb)
	edges_to_try.push_back(make_pair(*fb,*tb));
  
  currentIndicator = indicator(n);
  while (consecutiveTimesFoundLocalOptimum < trialsForLocalOptimum)
  { if (print_stuff)
      cout << "indicator: " << currentIndicator << '\n';
  
    bool foundImprovement = true;
    bool foundNeutralNeighbor = false;
    bool nadd;
    typename custom_types<network_t>::possible_edge neutral_neighbor;
    int neutralSteps = 0;
    while(foundImprovement ||
	  (foundNeutralNeighbor && neutralSteps < maxNeutralSteps))
    {
//     cout << "current network\n";
//     print_network(n);
      //if (animate_optimization)
#ifdef DISPLAY
      if(df)
        df->update(nSteps,n);
#endif//DISPLAY
      
      // shuffle the candidate edges.  this matters because the first
      //  one found that's an improvement will be used, or if none,
      //  the first one found that's neutral.
      random_shuffle(edges_to_try.begin(),edges_to_try.end());
    
      foundImprovement = foundNeutralNeighbor = false;
      // for each possible edge
      for (typename edges_to_try_t::iterator ei = edges_to_try.begin();
	   !foundImprovement && ei != edges_to_try.end(); ++ei)
      { typename graph_traits<network_t>::edge_descriptor ed;
        bool ed_found, ed_added;
	tie(ed,ed_found) = edge(ei->first,ei->second,n);
	if (ed_found)
	  remove_edge(ed,n);
	else
	  tie(ed,ed_added) = add_edge(ei->first,ei->second,n);
	//toggle_edge(*ei,n);
        // see if the indicator is better
        double newIndicator = indicator(n);
// 	if (print_stuff)
// 	  cout << newIndicator << endl;

	// for stochastic indicators: if it's big, call it again
	//  and make sure it's really big
	// (with the caching I've got in there, the more you call it
	//  the more accurate it is)
	if (indicator_is_stochastic &&
	    newIndicator >= currentIndicator)
	  // && newIndicator > 0) // what's this for?
	{ newIndicator = indicator(n);
//   	  if (print_stuff)
// 	    cout << newIndicator << endl;
	} 
// 	cout << "toggle (" << ei->first << ',' << ei->second << "): "
// 	     << newIndicator << endl;
	if (newIndicator >= currentIndicator &&
	    // need <= in case of infinity
	    newIndicator <= currentIndicator + 1e-7)
	{ if (!foundNeutralNeighbor)
  	  { neutral_neighbor = *ei;
            nadd = !ed_found;
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
            print_network(n);
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
      { toggle_edge(neutral_neighbor,n);
        if (print_stuff)
  	{ cout << "taking neutral step " << (nadd?'+':'-') << '('
	       << neutral_neighbor.first << ','
	       << neutral_neighbor.second << "): "
               << currentIndicator << '\n';
	  print_network(n);
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
#ifdef DISPLAY
    if (df)
      df->update(nSteps,n);
#endif//DISPLAY
  }
  if (print_stuff)
    cout << "indicator: " << currentIndicator << '\n';
}

// ------------------------------------------------------------------------
//  main
// ------------------------------------------------------------------------

// this edge_inserter is modeled on stl's front_inserter, etc.
template<typename graph_t>
class edge_insert_iterator
{
  graph_t&graph;
public:
  typedef typename graph_traits<graph_t>::vertices_size_type
    vertices_size_type;

  typedef typename graph_traits<graph_t>::edge_descriptor edge_type;
  typedef pair<vertices_size_type,vertices_size_type>     other_edge_type;

  typedef std::output_iterator_tag                    iterator_category;
  typedef edge_type                                   value_type;
  typedef void                                        difference_type;
  typedef void                                        pointer;
  typedef void                                        reference;

  edge_insert_iterator(graph_t&_g) : graph(_g) {}
  typedef edge_insert_iterator self_t;
  self_t &operator*() 
  { return *this;
  }
  self_t &operator=(const edge_type &e)
  { add_edge(source(e,graph),target(e,graph),graph);
    return *this;
  }
  self_t &operator=(const other_edge_type &e)
  { add_edge(source(e,graph),target(e,graph),graph);
    return *this;
  }
  self_t &operator++()
  { return *this;
  }
};

// factory function for convenience
template<typename graph_t>
edge_insert_iterator<graph_t> edge_inserter(graph_t&graph)
{ return edge_insert_iterator<graph_t>(graph);
}
  
int
main()
{  
  rng.seed((unsigned)time(0));
  // rand() used in random_shuffle()
  srand((unsigned)time(0));
  
  ur_t choose_01(rng,uniform_real<>());
  //double init_density = choose_01(); // rather than fixed init. dens.
  
  typedef adjacency_list<setS,vecS,bidirectionalS> network_t;
  network_t n(init_n_vertices);
  // network created empty, now initialize it  
  switch(0)
  {
  case 0: // Erdos-Renyi random graph with given size and density
    typedef custom_types<network_t>::ER_graph_generator ERgen;
    copy(ERgen(rng, init_n_vertices, init_density),
	 ERgen(), edge_inserter(n));
    break;
  case 1: // complete graph
    { graph_traits<network_t>::vertex_iterator si,send,ti,tend;
      for (tie(si,send) = vertices(n); si != send; ++si)
	for (tie(ti,tend) = vertices(n); ti != tend; ++ti)
	  if (*si != *ti)
	    add_edge(*si,*ti,n);
      break;
    }
  case 2: // star
    { graph_traits<network_t>::vertex_iterator ti,tend;
      for (tie(ti,tend) = vertices(n); ti != tend; ++ti)
	if (*ti > 0)
	{ add_edge(0,*ti,n);
  	  add_edge(*ti,0,n);
	}
      break;
    }
  case 3: // superstar
    { graph_traits<network_t>::vertex_iterator vi,vend;
      typedef graph_traits<network_t>::vertex_descriptor vertex;
      deque<vertex> q;
      tie(vi,vend) = vertices(n);
      copy(vi,vend,back_inserter(q));
      vertex center = q.front(); q.pop_front();
      // m is # of sub-centers, and # of feeders per sub-center
      // (though there may be a few missing at the end)
      unsigned m = static_cast<unsigned>(sqrt(num_vertices(n)));
      while(!q.empty())
      { vertex sub = q.front(); q.pop_front();
        add_edge(sub,center,n);
	int i;
	for (i = 0; i < m && !q.empty(); ++i)
	{ vertex feeder = q.front(); q.pop_front();
	  add_edge(feeder,sub,n);
	  add_edge(center,feeder,n);
	}
	if (i == 0)
	  add_edge(center,sub,n);
      }
      break;
    }
  default: // empty graph, I guess    
    break;
  }
  if (false) // reverse all the edges
  { graph_traits<network_t>::edge_iterator ei, eend;
    tie(ei,eend) = edges(n);
    typedef deque< graph_traits<network_t>::edge_descriptor > eq_t;
    eq_t eq;
    copy(ei,eend,back_inserter(eq));
    n.clear();
    // see <boost/lambda.h>
//     typedef graph_traits<network_t>::edge_descriptor ed;
//     transform(eq.begin(),eq.end(),edge_inserter(n),
//               edge(target(ret<const ed&>(*_1),n),source(*_1,n),n));
    while(!eq.empty())
    { graph_traits<network_t>::edge_descriptor ed = eq.front(); eq.pop_front();
      add_edge(target(ed,n),source(ed,n),n);
    }
    // does this work?
    cout << "tried reversing edges:\n";
    print_network(n);
  }
  if (false) // perturb one edge
  { int v1 = random_vertex(n,rng), v2;
    do v2 = random_vertex(n,rng); while (v2 == v1);
    toggle_edge(make_pair(v1,v2),n);
  }
  
  double moranprob =
      (1 - 1/mutantFitness) /
      (1 - 1/pow(mutantFitness,(double)init_n_vertices));
   
  if (print_stuff)
  { cout << "moran probability = " << moranprob << '\n';
    cout << "initial random network\n";
    print_network(n);
    cout << "initial density: " << density(n) << endl;
  }
  
  // adjust it to (near) specified density
  //improveNetwork(n,density_indicator(0.3342847893));

//   analytic_selection_indicator<network_t> aind;
  network_selection_indicator<network_t,rng_t> sind(rng);
//   analytic_selection_indicator_fixed_origin<network_t> aind(1);
//   selection_indicator_fixed_origin<network_t,rng_t> sind(1,rng);

//   comparing_indicator< analytic_selection_indicator<network_t>,
//     selection_indicator<network_t,rng_t> > ind(aind,sind);

//   deterministic_indicator_caching<network_t,double,
  //bad_at_selection_but_strongly_connected_indicator<network_t> ind;
//    deterministic_indicator_caching<network_t,double,
//      probability_to_amplification< analytic_selection_indicator<network_t> > >
//        pind;
  typedef deterministic_indicator_caching<network_t,double,
    analytic_selection_indicator<network_t> > pind_t;
  pind_t pind;
   
  constant_indicator mind(moranprob);

//   relative_indicator_t<pind_t,constant_indicator> rind(pind,mind);
//   relative_indicator_t< network_selection_indicator<network_t,rng_t>,
//     constant_indicator > srind(sind,mind);

  typedef comparing_indicator< pind_t,
    network_selection_indicator<network_t,rng_t> >
    cind_t;
  cind_t cind(pind,sind);

  relative_indicator_t<cind_t, constant_indicator> rind(cind,mind);

  //double(*ind)(const network_t&) = zero_indicator;

  // here, declare which indicator we are to optimize
  typeof(rind) &indicator(rind);

#ifdef DISPLAY
  DisplayFarm<network_t> displayfarm;
  BoostDotGraphController<network_t>
    opt_animation(animate_optimization?1:0, 1);
  displayfarm.installController(opt_animation);

  IndicatorsDisplayController<network_t>
    ind_display(animate_optimization?1:0);
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
  ind_display.installIndicator(
    &central_point_dominance_ind<network_t>,"central point dominance");
  displayfarm.installController(ind_display);

  IndicatorsDisplayController<network_t>
    mft_display(animate_optimization?1:0);
  mean_fixation_time_indicator<typeof sind> mft_sind(sind);
  mft_display.installIndicator(mft_sind,"mean fixation time");
  displayfarm.installController(mft_display);
  
  improveNetwork(n,indicator,&displayfarm);
#else
  improveNetwork(n,indicator,0);
#endif//DISPLAY

  if (print_stuff)
  { cout << "final network\n";
    print_network(n);
    //  cout << num_edges(n) << " edges" << endl;
    cout << "density " << density(n) << endl;
  }
  return 0;
}
