/* -*- C++ -*-
 * 
 * indicator classes for network optimization
 *
 * an indicator is anything that you can do this to:
 *   double value = my_indicator(my_network)
 * in other words it needs to have an operator() that takes a network
 *  and returns a double.
 *
 * you give an indicator to hill_climb() or simulated_annealing() or
 *  what have you, and it looks for a network that maximizes the
 *  indicator.  probability of fixation, for instance.
 *
 * indicators are also used to gather statistics on the networks,
 *  since each statistical measurement is an indicator: density, mean
 *  path length, etc.
*/
#ifndef __indicators_network_h__
#define __indicators_network_h__
#include "indicators-general.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/random.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/vector_property_map.hpp>
#include <boost/lambda/lambda.hpp>
#include <numeric>
#include <sstream>
#include <algorithm>
#ifdef DISPLAY
#include "IndicatorsDisplay.h"
#include "HistogramDisplay.h"
#endif
using namespace boost;
using namespace std;
  
// density of a graph is not found in the boost headers
template<typename network_t>
double density(const network_t &n)
{ typename network_t::vertices_size_type nv = num_vertices(n);
  return double(num_edges(n)) / (nv*(nv-1));
}

// having some trouble with template function types, decided
// to just make a class to give number of vertices of a graph
template<typename network_t>
class num_vertices_indicator
{public:
  double operator()(const network_t &n)
  { return num_vertices(n); 
  }
};

// the optimizing algorithms (search.h) call
// indicator_is_stochastic(indicator) to judge whether to call it more
// than once, to improve the estimate.  though, since adding the
// indicator_value type, this could possibly be done away with.
// meanwhile, it needs to be provided for every indicator class.
// indicators that are functions are all taken to be deterministic,
// for better or worse - there's a template function to that effect in
// indicators-general.h 
template<typename network_t>
bool indicator_is_stochastic(num_vertices_indicator<network_t>&i)
{ return false; }
  
// density_indicator(d) is an indicator that evaluates how close your
//  network's density is to d.
// i.e.
// { density_indicator di(0.3);
//   double ind = di(n);
// }
// or
//  double ind = density_indicator(0.3)(n);
class density_indicator
{
protected:
  double target;
public:
  density_indicator(double _target) : target(_target) {}
  template<typename network_t>
  double operator()(const network_t&n) const
  { return 1 - fabs(density(n)-target);
  }
};

// Given an unweighted directed graph, for the fixation process,
// strictly each edge has weight such that all edges leaving a vertex
// are equal and sum to 1.  You need to use these weights to compute
// temperature of a vertex (sum of its in edge weights).  This thing
// provides a 'weight property map' where one is needed, giving those
// weights.
template<typename network_t>
class stochastic_edge_weight_map
  : public put_get_helper<double,
			  stochastic_edge_weight_map<network_t> >
{
public:
  typedef typename graph_traits<network_t>::edge_descriptor key_type;
  typedef double value_type;
  typedef value_type reference;
  typedef boost::readable_property_map_tag category;

  const network_t&n;
  stochastic_edge_weight_map(const network_t &_n) : n(_n) {}
  inline value_type operator[](const key_type& e) const
  { typename graph_traits<network_t>::vertex_descriptor
      v = source(e,n);
    typename graph_traits<network_t>::degree_size_type
      outdeg = out_degree(v,n);
    return 1.0/outdeg;
  }
};

#if 0
//#include <boost/numeric/ublas/matrix.hpp>
#include <vnl/vnl_matrix.h>
#include <boost/graph/floyd_warshall_shortest.hpp>
template<typename network_t>
vnl_matrix<int> shortest_paths_matrix_raw(const network_t&n)
{ typename graph_traits<network_t>::vertices_size_type nv = num_vertices(n);
  typedef vnl_matrix<int> matrix_t;
  matrix_t distance_matrix(nv,nv);
  typedef constant_property_map<
    typename graph_traits<network_t>::edge_descriptor, int > wm_t;
  wm_t weight_map(1);
  floyd_warshall_all_pairs_shortest_paths(n,distance_matrix,
					  boost::weight_map(weight_map));
  return distance_matrix;
}

template<typename network_t>
vnl_matrix<int> shortest_paths_matrix(const network_t&n)
{ vnl_matrix<int> (*spmr)(const network_t&)
    = &shortest_paths_matrix_raw<network_t>;
  static deterministic_indicator_caching<network_t, vnl_matrix<int> >
    shortest_paths_caching(spmr);
  return shortest_paths_caching(n);
}

template<typename network_t>
double diameter_ind(const network_t&n)
{ vnl_matrix<int> distance_matrix = shortest_paths_matrix(n);
  return distance_matrix.array_inf_norm();
}

template<typename network_t>
double relative_diameter_ind(const network_t&n)
{ double rd = diameter_ind(n) / (num_vertices(n) - 1);
  return rd < 1 ? rd : 1;
}

template<typename network_t>
double mean_path_length(const network_t&n)
{ vnl_matrix<int> distance_matrix = shortest_paths_matrix(n);
  double sum = 0;
  for (int i = 0; i < distance_matrix.rows(); ++i)
    for (int j = 0; j < distance_matrix.cols(); ++j)
      if (i != j)
	sum += distance_matrix(i,j);
  typename graph_traits<network_t>::vertices_size_type nv = num_vertices(n);
  double mpl = sum / (double)(nv*(nv-1));
  return mpl;
}

template<typename network_t>
double relative_mean_path_length(const network_t&n)
{ double mpl = mean_path_length(n) / (num_vertices(n) - 1);
  return mpl < 1 ? mpl : 1;
}

template<typename network_t>
double mutuality_ind(const network_t&n)
{ typename graph_traits<network_t>::edge_iterator ei,eend;
  unsigned mut_count = 0; 
  for(tie(ei,eend) = edges(n); ei != eend; ++ei)
  { typename graph_traits<network_t>::edge_descriptor reverse;
    bool found;
    tie(reverse,found) = edge(target(*ei,n),source(*ei,n),n);
    if (found)
      ++mut_count;
  }
  return mut_count/(double)num_edges(n);
}

template<typename network_t>
double transitivity_ind(const network_t&n)
{ typename graph_traits<network_t>::edge_iterator e1i,e1end;
  typename graph_traits<network_t>::out_edge_iterator e2i,e2end;
  unsigned pair_count = 0, hypot_count = 0;
  for (tie(e1i,e1end) = edges(n); e1i != e1end; ++e1i)
    for (tie(e2i,e2end) = out_edges(target(*e1i,n),n); e2i != e2end; ++e2i)
      if (target(*e2i,n) != source(*e1i,n))
      { ++pair_count;
        typename graph_traits<network_t>::edge_descriptor hypot;
	bool found;
	tie(hypot,found) = edge(target(*e2i,n),source(*e1i,n),n);
	if (found)
	  ++hypot_count;
      }
  return hypot_count / (double)pair_count;
}

template<typename network_t>
double relative_mean_wavefront_ind(const network_t&n)
{ return aver_wavefront(n) / (double)num_vertices(n);
}
#endif//0

// this gives the median of a range of sorted values.
// I think it makes unstated assumptions about the iterator type.
template<typename iterator_t>
typename iterator_t::value_type median(iterator_t begin, unsigned n)
{ for (unsigned i = 0; 2*i < n; ++i)
    ++begin;
  typename iterator_t::value_type med = *begin;
  if (n % 2 == 1)
    med = (med + *++begin)/2;
  return med;
}

#if 0
template<typename network_t>
double median_out_density(const network_t&n)
{ typedef multiset<typename network_t::degree_size_type> distrib_t;
  distrib_t distrib; 
  typename network_t::vertices_size_type nv = num_vertices(n);
  typename graph_traits<network_t>::vertex_iterator vi,vend;
  for(tie(vi,vend) = vertices(n); vi != vend; ++vi)
    distrib.insert(out_degree(*vi,n));
  unsigned sz = distrib.size();
  typename distrib_t::iterator di = distrib.begin();
  for (unsigned i = 0; 2*i < sz; ++i)
    ++di;
  double med = *di;
  if (sz % 2 == 1)
    med = (med + *++di)/2.0;
  return med / nv;
}

template<typename network_t>
double median_in_density(const network_t&n)
{ typedef multiset<typename network_t::degree_size_type> distrib_t;
  distrib_t distrib; 
  typename network_t::vertices_size_type nv = num_vertices(n);
  typename graph_traits<network_t>::vertex_iterator vi,vend;
  for(tie(vi,vend) = vertices(n); vi != vend; ++vi)
    distrib.insert(in_degree(*vi,n));
  return median(distrib.begin(),distrib.size()) / nv;
}

// temperature of a vertex is sum_j w_ji
template<typename network_t, typename WeightMap>
double temperature(typename graph_traits<network_t>::vertex_descriptor v,
		   const network_t&n, const WeightMap&w)
{ typename graph_traits<network_t>::in_edge_iterator first, last;
  double temp = 0; 
  for(tie(first, last) = in_edges(v,n); first != last; first++)
    temp += get(w,*first);
  return temp;
}

template<typename network_t>
double median_temperature(const network_t&n)
{ stochastic_edge_weight_map<network_t> wm(n);
  multiset<double> temps;
  typename graph_traits<network_t>::vertex_iterator vi,vend;
  for (tie(vi,vend) = vertices(n); vi != vend; ++vi)
    temps.insert(temperature(*vi,n,wm));
  return median(temps.begin(),temps.size());
}
#endif//0

// central point dominance measures how many paths through the graph
// go through the same central point.  it's provided by boost, so here
// it is.
template<typename network_t>
double central_point_dominance_ind(const network_t&n)
{ typedef
    map<typename graph_traits<network_t>::vertex_descriptor,double> map_t;
  map_t cmap;
  associative_property_map<map_t> centralitymap(cmap);
  brandes_betweenness_centrality(n,centralitymap);
  relative_betweenness_centrality(n,centralitymap);
  return central_point_dominance(n,centralitymap);
}

// some machinery used in finding diameter and mean shortest path length
//#include <boost/numeric/ublas/matrix.hpp>
#include <vnl/vnl_matrix.h>
#include <boost/graph/floyd_warshall_shortest.hpp>
template<typename network_t>
vnl_matrix<int> shortest_paths_matrix_raw(const network_t&n)
{ typename graph_traits<network_t>::vertices_size_type nv = num_vertices(n);
  typedef vnl_matrix<int> matrix_t;
  matrix_t distance_matrix(nv,nv);
  typedef constant_property_map<
    typename graph_traits<network_t>::edge_descriptor, int > wm_t;
  wm_t weight_map(1);
  floyd_warshall_all_pairs_shortest_paths(n,distance_matrix,
					  boost::weight_map(weight_map));
  return distance_matrix;
}

template<typename network_t>
vnl_matrix<int> shortest_paths_matrix(const network_t&n)
{ vnl_matrix<int> (*spmr)(const network_t&)
    = &shortest_paths_matrix_raw<network_t>;
  static deterministic_indicator_caching<network_t, vnl_matrix<int> >
    shortest_paths_caching(spmr);
  return shortest_paths_caching(n);
}

template<typename network_t>
double diameter_ind(const network_t&n)
{ vnl_matrix<int> distance_matrix = shortest_paths_matrix(n);
  return distance_matrix.array_inf_norm();
}

// `relative diameter' here is the diameter divided by the max value
// it could have, i.e. by (#vertices - 1)
template<typename network_t>
double relative_diameter_ind(const network_t&n)
{ double rd = diameter_ind(n) / (num_vertices(n) - 1);
  return rd < 1 ? rd : 1;
}

// mean [shortest] path length is average over all vertex pairs of the
// shortest distance between.
template<typename network_t>
double mean_path_length(const network_t&n)
{ vnl_matrix<int> distance_matrix = shortest_paths_matrix(n);
  double sum = 0;
  for (int i = 0; i < distance_matrix.rows(); ++i)
    for (int j = 0; j < distance_matrix.cols(); ++j)
      if (i != j)
	sum += distance_matrix(i,j);
  typename graph_traits<network_t>::vertices_size_type nv = num_vertices(n);
  double mpl = sum / (double)(nv*(nv-1));
  return mpl;
}

// relative mean path length is divided by the largest it could be, again
// (#vertices - 1)
template<typename network_t>
double relative_mean_path_length(const network_t&n)
{ double mpl = mean_path_length(n) / (num_vertices(n) - 1);
  return mpl < 1 ? mpl : 1;
}

// mutuality [or symmetry?] is proportion of the [directed] graph's
// edges that are matched by an edge in the opposite direction
template<typename network_t>
double mutuality_ind(const network_t&n)
{ typename graph_traits<network_t>::edge_iterator ei,eend;
  unsigned mut_count = 0; 
  for(tie(ei,eend) = edges(n); ei != eend; ++ei)
  { typename graph_traits<network_t>::edge_descriptor reverse;
    bool found;
    tie(reverse,found) = edge(target(*ei,n),source(*ei,n),n);
    if (found)
      ++mut_count;
  }
  return mut_count/(double)num_edges(n);
}

// transitivity is proportion, over all consecutive pairs of edges
// (a->b->c), of those whose triangle is closed by a third edge c->a.
template<typename network_t>
double transitivity_ind(const network_t&n)
{ typename graph_traits<network_t>::vertex_iterator   vi, vend;
  typename graph_traits<network_t>::out_edge_iterator e1i,e1end;
  typename graph_traits<network_t>::out_edge_iterator e2i,e2end;
  unsigned pair_count = 0, hypot_count = 0;
  for (tie(vi,vend) = vertices(n); vi != vend; ++vi)
    for (tie(e1i,e1end) = out_edges(*vi,n); e1i != e1end; ++e1i)
      for (tie(e2i,e2end) = out_edges(target(*e1i,n),n); e2i != e2end; ++e2i)
        if (target(*e2i,n) != *vi)
        { ++pair_count;
          typename graph_traits<network_t>::edge_descriptor hypot;
          bool found;
          tie(hypot,found) = edge(target(*e2i,n),*vi,n);
          if (found)
            ++hypot_count;
        }
  if (pair_count == 0)
    return 0;
  else
    return hypot_count / (double)pair_count;
}

// not sure what wavefront is, but it's provided by the library.
// I think it's used in balancing a graph so its matrix is easy to
// invert, or something.
template<typename network_t>
double relative_mean_wavefront_ind(const network_t&n)
{ return aver_wavefront(n) / (double)num_vertices(n);
}

// out density of a vertex v is its out degree, divided by the maximum
// possible out degree, which is (#vertices-1).  mean out density is
// always equal to the graph's density; so comparing median out
// density to graph density gives an indication of what the out degree
// distribution is like.  might be good to look at its variance, though.
template<typename network_t>
double median_out_density(const network_t&n)
{ typedef multiset<typename network_t::degree_size_type> distrib_t;
  distrib_t distrib; 
  typename network_t::vertices_size_type nv = num_vertices(n);
  typename graph_traits<network_t>::vertex_iterator vi,vend;
  for(tie(vi,vend) = vertices(n); vi != vend; ++vi)
    distrib.insert(out_degree(*vi,n));
//   unsigned sz = distrib.size();
//   typename distrib_t::iterator di = distrib.begin();
//   for (unsigned i = 0; 2*i < sz; ++i)
//     ++di;
//   double med = *di;
//   if (sz % 2 == 1)
//     med = (med + *++di)/2.0;
  return median(distrib.begin(),distrib.size()) / nv;
}

// analogous to median out density
template<typename network_t>
double median_in_density(const network_t&n)
{ typedef multiset<typename network_t::degree_size_type> distrib_t;
  distrib_t distrib; 
  typename network_t::vertices_size_type nv = num_vertices(n);
  typename graph_traits<network_t>::vertex_iterator vi,vend;
  for(tie(vi,vend) = vertices(n); vi != vend; ++vi)
    distrib.insert(in_degree(*vi,n));
  return median(distrib.begin(),distrib.size()) / nv;
}

// temperature of a vertex is sum_j w_ji, i.e. sum of all in-edges'
// temperature (1 over the edge's source vertex's out degree).
// the stochastic edge weight map should be provided to the function,
// to give it edge temperatures.
template<typename network_t, typename WeightMap>
double temperature(typename graph_traits<network_t>::vertex_descriptor v,
		   const network_t&n, const WeightMap&w)
{ typename graph_traits<network_t>::in_edge_iterator first, last;
  double temp = 0; 
  for(tie(first, last) = in_edges(v,n); first != last; first++)
    temp += get(w,*first);
  return temp;
}

// median temperature over all vertices.  mean temperature is always
// 1, so this gives a view into temperature distribution.  as with
// degree distributions, probably wise to look at variance instead.
template<typename network_t>
double median_temperature(const network_t&n)
{ stochastic_edge_weight_map<network_t> wm(n);
  multiset<double> temps;
  typename graph_traits<network_t>::vertex_iterator vi,vend;
  for (tie(vi,vend) = vertices(n); vi != vend; ++vi)
    temps.insert(temperature(*vi,n,wm));
  return median(temps.begin(),temps.size());
}


// this just gives the moran probability of the # vertices of the
// graph
class moran_indicator
{ double r;
public:
  moran_indicator(double _r) : r(_r) {}
  template<typename network_t>
  double operator()(const network_t&n)
  { typename network_t::vertices_size_type nv = num_vertices(n);
    return moran_probability(r,nv);
  }
};

// it's not stochastic
bool indicator_is_stochastic(moran_indicator&mind)
{ return false; }

// before doing fixation simulation call this to evaluate
//  whether there's a possibility of nonzero fixation probability
template<typename network_t>
bool is_fixation_candidate(const network_t &n)
{
  // try in 2 phases.  first quick look for isolated or source vertices.
  // then sophisticated look for isolated or source components.
  // when is each worth the time cost?
  typename graph_traits<network_t>::vertex_iterator vi,vend;
  int n_sources = 0;
  for (tie(vi,vend) = vertices(n); vi != vend; ++vi)
  { typename inv_adjacency_iterator_generator<network_t>::type ini, inend;
    tie(ini,inend) = inv_adjacent_vertices(*vi,n);
    // if there are no edges into *vi
    if (ini == inend)
    { typename graph_traits<network_t>::adjacency_iterator outi,outend;
      tie(outi,outend) = adjacent_vertices(*vi,n);
      // and none out either, it's an isolated vertex, no fixation
      if (outi == outend)
	return 0;
      // if out edges only, a source, if more than one source no fixation
      else if (n_sources > 0)
	return 0;
      else
	++n_sources;
    }
  }
  // if we get through that without returning 0, try the sophisticated
  // check

  // first get the strongly connected components of the network.
  typedef typename network_t::vertices_size_type comp_t;
  vector<comp_t> component_map;
  component_map.resize(num_vertices(n));
  comp_t ncomponents = strong_components(n,component_map.data());
  if (ncomponents == 1)
    return true;
  // if there's more than one component with no in-edges,
  // or even one isolated component, no fixation.
  vector<bool> has_in_edges(ncomponents), has_out_edges(ncomponents);
  typename graph_traits<network_t>::edge_iterator ei,eend;
  for (tie(ei,eend) = edges(n); ei != eend; ++ei)
  { comp_t source_comp = component_map[source(*ei,n)],
      target_comp = component_map[target(*ei,n)];
    if (source_comp != target_comp)
    { has_in_edges[target_comp] = true;
      has_out_edges[source_comp] = true;
    }
  }
  bool found_source_component = false;
  for(int i = 0; i < ncomponents; ++i)
    if (!has_in_edges[i])
    { if (!has_out_edges[i])
        return false; // isolated component - no fixation
      else if (found_source_component)
	return false; // multiple source components - no fixation
      else
	found_source_component = true;
    }
//   if (print_stuff)
//   {
//     cout << "has multiple components but passes fixation test...\n"
// 	 << "has_in_edges: ";
//     for(int i=0; i<ncomponents; ++i)
//       cout << has_in_edges[i]? 1:0;
//     cout << "\nhas_out_edges: ";
//     for(int i=0; i<ncomponents; ++i)
//       cout << has_out_edges[i]? 1:0;
//     cout << endl;
//   }
  return true;
}

// before doing fixation simulation use this machinery to evaluate
//  whether there's a possibility of fixation, in case the graph is
//  disconnected or otherwise incompetent.
typedef enum
  { STRONGLY_CONNECTED, DISCONNECTED,
    ONE_SOURCE_COMPONENT, MULTIPLE_SOURCE_COMPONENTS }
  ncsf_t;
class network_component_status
{
public:
  ncsf_t flag;
  int source_size; // used in case of one source component
  network_component_status(ncsf_t _f, int _s=0)
    : flag(_f), source_size(_s) {}
};

#include <boost/graph/graph_utility.hpp>
// some trouble in the boost version of this fn, so I have an altered
// version here.  no record of what is changed.
template <typename VertexListGraph, typename VertexColorMap>
inline bool is_connected(const VertexListGraph& g, VertexColorMap color)
{
  typedef typename property_traits<VertexColorMap>::value_type ColorValue;
  typedef color_traits<ColorValue> Color;
  typename graph_traits<VertexListGraph>::vertex_iterator 
    ui, ui_end, vi, vi_end, ci, ci_end;
  for (tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui)
    for (tie(vi, vi_end) = vertices(g); vi != vi_end; ++vi)
      if (*ui != *vi) {
	for (tie(ci, ci_end) = vertices(g); ci != ci_end; ++ci) 
	  put(color, *ci, Color::white());
	if (! is_reachable(*ui, *vi, g, color))
	  return false;
      }
  return true;
}

// very simple quick exponentiation for just the cases we need
double index_power(int i, int n)
{ switch(n) 
  { case 0: return 1;
    case 1: return i;
    case 2: return i*i;
    case -1: return 1.0 / i;
    case -2: return 1.0 / (i*i);
    default: return 0;
  }
}

// moments of degree distribution
// degree_moment(i,j) = sum_x in(x)^i out(x)^j / N
template<typename network_t>
double degree_moment(const network_t&n, int i, int j)
{ double accum = 0;
  typename graph_traits<network_t>::vertex_iterator 
    vi, vi_end;
  for (tie(vi, vi_end) = vertices(n); vi != vi_end; ++vi)
    accum += index_power(in_degree(*vi,n),i)*index_power(out_degree(*vi,n),j);
    //accum += pow(double(in_degree(*vi,n)),i)*pow(double(out_degree(*vi,n)),j);
  return accum / num_vertices(n);
}

// degree-weighted moments of network state
// weighted_moment(m,n)
//   = sum_x state(x) in(x)^m out(x)^n / sum_x in(x)^m out(x)^n
template<typename network_t, typename state_t>
double weighted_moment(const network_t&net, int m, int n, state_t &state)
{ double numer = 0, denom = 0;
  typename graph_traits<network_t>::vertex_iterator 
    vi, vi_end;
  for (tie(vi, vi_end) = vertices(net); vi != vi_end; ++vi)
  { double ixmoxn = 
      index_power(in_degree(*vi,net),m)*index_power(out_degree(*vi,net),n);
    numer += state[*vi] * ixmoxn;
    denom += ixmoxn;
  }
  return numer / denom;
}

// functions to evaluate whether a graph is capable of fixation.
// first, for undirected graphs, they just need to be connected.
template<typename A, typename B,  typename D,
         typename E, typename F, typename G>
network_component_status
test_fixation_candidate(
  const boost::adjacency_list<A,B,boost::undirectedS,D,E,F,G> &n)
{
  typedef boost::adjacency_list<A,B,boost::undirectedS,D,E,F,G> network_t;
  map<typename graph_traits<network_t>::vertex_descriptor,
      default_color_type>
    colors;
  if (!::is_connected(n,make_assoc_property_map(colors)))
    return DISCONNECTED;
  return STRONGLY_CONNECTED;
}

// now for directed graphs, more complex.  in order to be really
// interesting they need to be strongly connected.  if only weakly
// connected there is a component that can't be reached from anywhere
// else, and fixation can only happen if the mutation arises there and
// fixates there, and then the rest of the network will necessarily
// fixate eventually.  if there are multiple such components, or the
// graph isn't connected, fixation is impossible.
template<typename network_t>
network_component_status test_fixation_candidate_directed(const network_t &n)
{
  // try in 2 phases.  first quick look for isolated or source vertices,
  // then sophisticated look for isolated or source components.
  // when is each worth the time cost?
  typename graph_traits<network_t>::vertex_iterator vi,vend;
  int n_sources = 0;
  for (tie(vi,vend) = vertices(n); vi != vend; ++vi)
  { typename inv_adjacency_iterator_generator<network_t>::type ini, inend;
    tie(ini,inend) = inv_adjacent_vertices(*vi,n);
    // if there are no edges into *vi
    if (ini == inend)
    { typename graph_traits<network_t>::adjacency_iterator outi,outend;
      tie(outi,outend) = adjacent_vertices(*vi,n);
      // and none out either, it's an isolated vertex, no fixation
      if (outi == outend)
	return DISCONNECTED;
      // if out edges only, a source, if more than one source no fixation
      else if (n_sources > 0)
	return MULTIPLE_SOURCE_COMPONENTS;
      else
	++n_sources;
    }
  }
  // if we get through that without returning failure, try the
  // more sophisticated check

  // first get the strongly connected components of the network.
  typedef typename network_t::vertices_size_type comp_t;
  comp_t component_map[num_vertices(n)];
  comp_t ncomponents = strong_components(n,&component_map[0]);
  if (ncomponents == 1)
    return STRONGLY_CONNECTED;
  // if there's more than one component with no in-edges,
  // or even one isolated component, no fixation.
  vector<bool> has_in_edges(ncomponents), has_out_edges(ncomponents);
  typename graph_traits<network_t>::edge_iterator ei,eend;
  for (tie(ei,eend) = edges(n); ei != eend; ++ei)
  { comp_t source_comp = component_map[source(*ei,n)],
      target_comp = component_map[target(*ei,n)];
    if (source_comp != target_comp)
    { has_in_edges[target_comp] = true;
      has_out_edges[source_comp] = true;
    }
  }
  bool found_source_component = false;
  int size;
  for(int i = 0; i < ncomponents; ++i)
    if (!has_in_edges[i])
    { if (!has_out_edges[i])
        return DISCONNECTED; // isolated component - no fixation
      else if (found_source_component)
         // multiple source components - no fixation
	return MULTIPLE_SOURCE_COMPONENTS;
      else
      { found_source_component = true;
        size = count(component_map,component_map+num_vertices(n),i);
      }
    }
  // if there's a single source component, fixation is possible,
  // though it's less likely
  if (found_source_component)
    return network_component_status(ONE_SOURCE_COMPONENT,size);
//   if (print_stuff)
//   {
//     cout << "has multiple components but passes fixation test...\n"
// 	 << "has_in_edges: ";
//     for(int i=0; i<ncomponents; ++i)
//       cout << has_in_edges[i]? 1:0;
//     cout << "\nhas_out_edges: ";
//     for(int i=0; i<ncomponents; ++i)
//       cout << has_out_edges[i]? 1:0;
//     cout << endl;
//   }
  return STRONGLY_CONNECTED;
}

template<typename A, typename B,  typename D,
         typename E, typename F, typename G>
network_component_status
test_fixation_candidate(
  const boost::adjacency_list<A,B,boost::directedS,D,E,F,G> &n)
{ return test_fixation_candidate_directed(n); }

template<typename A, typename B,  typename D,
         typename E, typename F, typename G>
network_component_status
test_fixation_candidate(
  const boost::adjacency_list<A,B,boost::bidirectionalS,D,E,F,G> &n)
{ return test_fixation_candidate_directed(n); }

// this is generally useful to network_selection_indicator:
//  choose one of a sequence of things, with probability proportional to
//  some weighting
// requires the type of _Map[*_Iterator] to be able to cast to double
// a positive total_weights can be provided, to speed this up
// it's used to find a source vertex with probability proportional to
//  fitness, of course.
template<typename _Iterator, typename _Map, typename _RNG>
typename iterator_traits<_Iterator>::value_type
choose_weighted(_Iterator first, _Iterator end, _Map &weights,
		_RNG &rng, double total_weight=-1)
{ if (total_weight <= 0)
    total_weight = accumulate_mapped(first,end,0.0,weights);
  ur_t choose_01(rng,uniform_real<>());
  double cumul = choose_01() * total_weight;
  for (; first != end; ++first)
    if ((cumul -= weights[*first]) <= 0)
      return *first;
  // should never come to here
  return *end;
}
template<typename network_t, typename RNG_t, typename params_t>
class network_selection_indicator;

// calling the monte carlo code for different kinds of adjacency_list
// needs to be outside the class, it seems, to make the template code work
template<typename A,  typename C, typename D, typename E,
         typename F, typename G, typename RNG_t, typename params_t>
fixation_result 
call_try_fixation(network_selection_indicator<
                    boost::adjacency_list<A,boost::vecS,C,D,E,F,G>,
                    RNG_t,params_t>*ts,
                  const boost::adjacency_list<A,boost::vecS,C,D,E,F,G>&n)
{ vector<double> state(num_vertices(n));
  return ts->fixation_once(n,random_vertex(n,ts->rng),state);
}

template<typename A,  typename C, typename D, typename E,
         typename F, typename G, typename RNG_t, typename params_t>
fixation_result 
call_try_fixation(network_selection_indicator<
                    boost::adjacency_list<A,boost::setS,C,D,E,F,G>,
                    RNG_t,params_t>*ts,
                  const boost::adjacency_list<A,boost::setS,C,D,E,F,G>&n)
{ // not implemented
  // set storage in this simulation seems to be really inefficient
}

// specialized fixation_stats for network -- includes recording of
// fixation probability as a function of starting vertex
template<typename network_t>
class net_fixation_stats : public fixation_stats
{
public:
  typedef typename graph_traits<network_t>::vertex_descriptor vertex_descriptor;
  typedef map<vertex_descriptor, fixation_record> fixations_by_vertex_t;
  fixations_by_vertex_t fixations_by_vertex;    
};

#ifdef DISPLAY
#include "BoostDotDisplay.h"
#endif
// class for monte carlo fixation simulation for graphs.
// this specializes the general monte carlo fixation code in
// indicators-general.h, qv.
template<typename network_t, typename RNG_t, typename params_t>
class network_selection_indicator
  : public selection_indicator<network_t,RNG_t,params_t,
                               net_fixation_stats<network_t> >
{ typedef selection_indicator< network_t,RNG_t,params_t,
                               net_fixation_stats<network_t> > selind;
public:
  using ObjectWithParameters<params_t>::params;
  using selind::currently_processing;

  // cache parameters values for quickness in looping
  Parameters::update_rule_t update_rule;
  int maxsteps;
public:
  typedef net_fixation_stats<network_t> fixation_stats_t;

  network_selection_indicator(RNG_t _rng, params_t*_p=0)
    : selind(_rng,_p) {}

  string cache_key(const network_t&n)
  { const string *pvs = params.get("perturb_r_at_vertex");
    if (!pvs)
      return canonical(n);
    else
      return canonical(n) + '\n' + *pvs;
  }

  // evaluate network for possibility of fixation
  // one special case: when exactly one vertex is the unique source
  // component, fixation probability is 1/#vertices, since fixation is
  // guaranteed when the mutation starts there and impossible otherwise.
  bool test_fixation_candidate(const network_t&n, double *answer_p)
  { network_component_status status = ::test_fixation_candidate(n);
    switch(status.flag)
    {
    case STRONGLY_CONNECTED:
      //cout << "strongly connected\n";
      return false;
    case DISCONNECTED:
    case MULTIPLE_SOURCE_COMPONENTS:
      cout << "multiple source components or disconnected\n";
      *answer_p = 0;
      return true;
    case ONE_SOURCE_COMPONENT:
      cout << "single source component, size = "
           << status.source_size << '\n';
      if (status.source_size == 1)
      { *answer_p = 1.0 / num_vertices(n);
        return true;
      }
      //else
    default:
      if (params.require_strongly_connected())
      { cout << "Graph is not strongly connected.  Rejecting." << endl;
        *answer_p = -1;
        return true;
      }
      // else
      return false;
    }
  }

  // simulate to fixation/extinction once
  // call out to templated function, which calls back in to fixation_once
  fixation_result try_fixation(const network_t&n)
  { return ::call_try_fixation(this,n);
  }

  // machinery to choose an edge to reproduce along.
  // 3 cases, depending on the update_rule parameter:
  // egt (aka birth-death): choose a source vertex with probability
  //  proportional to fitness among all vertices; then choose among
  //  all its out edges equally
  // vm (aka death-birth): choose a target vertex with probability
  //  inversely proportional to fitness, then choose among all its in
  //  edges equally
  // skye (also a death-birth process): choose a target vertex among
  //  all vertices equally, then choose among all the target's
  //  in-edges, with probability proportional to the fitnesses of
  //  their source vertices.
  typedef pair<typename graph_traits<network_t>::vertex_descriptor,
	       typename graph_traits<network_t>::vertex_descriptor>
  parent_child_pair;

  // egt case
  // don't call if there are no edges, it'll search forever
  template<typename fitness_map_t>
  parent_child_pair
  choose_parent_child_egt(const network_t &n, double total_fitness,
                          fitness_map_t&state)
  {
    typename graph_traits<network_t>::vertex_descriptor parent,child;
    while(1) // loop until a successful match found
    {
      // select a parent vertex, with probability proportional to fitness
      typename graph_traits<network_t>::vertex_iterator vi,vend;
      tie(vi,vend) = vertices(n);
      parent = choose_weighted(vi,vend,state,this->rng,total_fitness);
      // now choose a child vertex, with prob. proportional to edge weight
      //  which are all assumed to be 1
      typename network_t::degree_size_type n_children = out_degree(parent,n);
      typename graph_traits<network_t>::adjacency_iterator ci,cend;
      if (n_children <= 0)
	continue; // no out edges, start over
      else if (n_children == 1)
      { tie(ci,cend) = adjacent_vertices(parent,n);
        child = *ci;
      }
      else
      { ui_t choose_a_child(this->rng,uniform_int<>(0,n_children-1));
        int child_selection = choose_a_child();
	for(tie(ci,cend) = adjacent_vertices(parent,n);
	    ci != cend; ++ci)
	  if (child_selection-- == 0)
	    break;
	if (child_selection >= 0)
	{ cerr << "Error finding a child!" << endl; continue; } 
	child = *ci;
      }
      // if we get to here, we found a parent and child.
      return make_pair(parent,child);
    }
  }
  
  // skye case
  // don't call if there are no edges, it'll search forever
  template<typename fitness_map_t>
  parent_child_pair
  choose_parent_child_skye(const network_t &n, fitness_map_t&state)
  {
    typename network_t::vertices_size_type nv = num_vertices(n);
    typename network_t::vertices_size_type n_infected = 0;
    typename graph_traits<network_t>::vertex_descriptor parent, child;
    while(1) // loop until a successful match found
    {
      // select a child at random
      child = random_vertex(n,this->rng);
      // select parent via in edges
      typename inv_adjacency_iterator_generator<network_t>::type pi, pend;
      tie(pi,pend) = inv_adjacent_vertices(child,n);
      if (pi == pend) // if child has no in edges, try again
        continue;
      // weighted by fitness
      parent = choose_weighted(pi,pend,state,this->rng);
      // sanity check
      typename graph_traits<network_t>::edge_descriptor ed;
      bool ed_found;
      tie(ed,ed_found) = edge(parent,child,n);
      if (!ed_found)
      { cerr << "found invalid edge " << parent << ',' << child << endl;
        continue;
      }
      // if we get here, we're good
      return make_pair(parent,child);
    }
  }

  // vm case
  // don't call if there are no edges, it'll search forever
  template<typename rate_map_t>
  parent_child_pair
  choose_parent_child_vm(const network_t &n, double total_rate,
                         rate_map_t&rate)
  { typename graph_traits<network_t>::vertex_descriptor parent,child;
    while(1) // loop until a successful match found
    { // select a parent vertex, with probability proportional to fitness
      typename graph_traits<network_t>::vertex_iterator vi,vend;
      tie(vi,vend) = vertices(n);
      child = choose_weighted(vi,vend,rate,this->rng,total_rate);
      // now choose a parent vertex, with prob. proportional to edge weight
      //  which are all assumed to be 1
      typename network_t::degree_size_type n_parents = in_degree(child,n);
      typename inv_adjacency_iterator_generator<network_t>::type pi, pend;
      if (n_parents <= 0)
        continue; // no out edges, start over
      else if (n_parents == 1)
      { tie(pi,pend) = inv_adjacent_vertices(child,n);
        parent = *pi;
      }
      else
      { ui_t choose_a_parent(this->rng,uniform_int<>(0,n_parents-1));
        int parent_selection = choose_a_parent();
	for(tie(pi,pend) = inv_adjacent_vertices(child,n);
	    pi != pend; ++pi)
	  if (parent_selection-- == 0)
	    break;
	if (parent_selection >= 0)
	{ cerr << "Error finding a parent!" << endl; continue; } 
	parent = *pi;
      }
      // if we get to here, we found a parent and child.
      return make_pair(parent,child);
    }
  }
  
  // assign one mutant in a resident population, simulate
  //  whether it fixates (within a given threshold of patience).
  // which place to mutate is specified by the caller, so it can either
  //  mutate the same vertex over and over, or do a random one each time.

//   fixation_result
//   fixation_once(const network_t &n,
//                 typename graph_traits<network_t>::vertex_descriptor mutation)
//   { statevector.resize(num_vertices(n));
//     return fixation_once(n,mutation,statevector);
//   }

  // yes, it would be better to do these tests outside the loop,
  // leave me alone
  template<typename state_t>
  fixation_result
  fixation_once(const network_t &n,
                typename graph_traits<network_t>::vertex_descriptor mutation,
                state_t &state)
  { typedef typename graph_traits<network_t>::vertex_descriptor vertex_t;
    constant_indicator rfi(params.residentFitness());

    // is there a fitness-perturbed vertex?
    const string *pvs = params.get("perturb_r_at_vertex");
    if (pvs) // yes, there is one
    { vertex_t pv = string_to_unsigned(*pvs);
      typedef one_of_these_things_indicator<vertex_t,double> mfi_t;
      mfi_t mfi(params.mutantFitness(),
                pv, params.perturbed_mutant_fitness());
      fixation_result res = fixation_once(n,mutation,state,rfi,mfi);
      return res;
    }
    else // no
    { constant_indicator mfi(params.mutantFitness());
      return fixation_once(n,mutation,state,rfi,mfi);
    }
  }

  template<typename state_t,
           typename res_fitness_ind_t, typename mut_fitness_ind_t> 
  fixation_result
  fixation_once(const network_t &n,
                typename graph_traits<network_t>::vertex_descriptor mutation,
                state_t &state, res_fitness_ind_t res_fitness_ind,
                mut_fitness_ind_t mut_fitness_ind)
  { typename network_t::vertices_size_type nv = num_vertices(n);
    typename network_t::edges_size_type ne = num_edges(n);
    if (nv <= 0 || ne < 0)
      return fixation_result(ERROR,0);

    // store the parameters, to speed up the inner loops
    update_rule = params.update_rule();
    int maxsteps = params.maxStepsToFixation();

    fill(state.begin(),state.end(),0);

    // is there a special initial arrangement of mutants?
    const char *special_init = params.get("initial_configuration");
    if ( !special_init ) // no, initial state is all residents but one
      state[mutation] = 1;
    else 
    { unsigned thresh = params.initial_configuration_threshold();
      if (special_init == "low-in")
      { typename graph_traits<network_t>::vertex_iterator vi,vend;
        for (tie(vi,vend) = vertices(n); vi != vend; ++vi)
          if (in_degree(*vi,n) < thresh)
            state[*vi] = 1;
      } else if (special_init =="high-in")
      { typename graph_traits<network_t>::vertex_iterator vi,vend;
        for (tie(vi,vend) = vertices(n); vi != vend; ++vi)
          if (in_degree(*vi,n) >= thresh)
            state[*vi] = 1;
      } else if (special_init =="low-out")
      { typename graph_traits<network_t>::vertex_iterator vi,vend;
        for (tie(vi,vend) = vertices(n); vi != vend; ++vi)
          if (out_degree(*vi,n) < thresh)
            state[*vi] = 1;
      } else if (special_init =="high-out")
      { typename graph_traits<network_t>::vertex_iterator vi,vend;
        for (tie(vi,vend) = vertices(n); vi != vend; ++vi)
          if (out_degree(*vi,n) >= thresh)
            state[*vi] = 1;
      }
    }

    // below, we use this to map vertex to fitness
    typedef map_to_indicator<state_t> state_indicator_t;
    state_indicator_t stind(state);
    typedef switching_indicator_t<res_fitness_ind_t,
      mut_fitness_ind_t, state_indicator_t> fitness_ind_t;
    fitness_ind_t fitness_ind(res_fitness_ind, stind, 1,
                              mut_fitness_ind);

    // counter of number of infected vertices
    typename network_t::vertices_size_type n_infected = 1;
    // counter of reproduction events
    int n_steps = 0;
    // clock runs slower than n_steps, providing a time axis that
    // matches the analytical calculation of absorption time: fitness
    // equates to expected number of reproduction events per vertex
    // per time unit, so average number of events per time unit is the
    // sum of all the fitnesses.
    double clock = 0;
    // keep track of biggest number infected
    int peak = 0;

#ifdef DISPLAY
    // devices for animating the fixation process
    bool animate = params.animate_fixation();
    // running count of infected vertices
    IndicatorsDisplayController<network_t,params_t> fixationDisplay;
    fixationDisplay.inheritParametersFrom(params);
    fixationDisplay.setoutfilenamebase("fixation");
    fixationDisplay.params.setdisplayToScreen(animate);
    // moving histogram of probability of fixation given reaching a
    // certain number infected
    HistogramDisplayController<int,params_t> hist_display;
    hist_display.inheritParametersFrom(params);
    hist_display.params.setdisplayToScreen(animate);
    if(animate)
    { fixationDisplay.params.setrecordEvery(0);
      fixationDisplay.installIndicator(variableWatcher(n_infected),"mutants");
      hist_display.params.setrecordEvery(0);
      hist_display.setTitle(0,"conditional prob of fixation");
    }

    // update raw data for histogram at each step
    //typename selind::fixation_stats &nc = selind::cache[cache_key(n)];
    vector<int> history;
    if(animate)
    {//  nc.histogram_success.resize(nv+1);
//       nc.histogram_failure.resize(nv+1);
      history.push_back(n_infected);
    }
#endif
    static bool recorded = false;
    bool record_dynamics = !recorded && params.record_dynamics();
    recorded = true;
    CSVDisplay dynamics_csv;
    if (record_dynamics)
    { dynamics_csv.openFile(params.outputDirectory() + "/dynamics.csv");
      dynamics_csv << "t";
      for (int i = -2; i <= 2; ++i)
        for (int j = -2; j <= 2; ++j)
          dynamics_csv << stringf("omega.%d.%d", i, j);
      dynamics_csv.newRow();
      dynamics_csv << clock;
      for (int i = -2; i <= 2; ++i)
        for (int j = -2; j <= 2; ++j)
          dynamics_csv << weighted_moment(n,i,j,state);
      dynamics_csv.newRow();
    }

    // do infection process until fixation, extinction or timeout
    while (0 < n_infected && n_infected < nv &&
           n_steps++ < maxsteps)
    { //cout << n_infected << ' ';

//       double total_fitness = nv*residentFitness +
//         n_infected*(mutantFitness-residentFitness);
      typename graph_traits<network_t>::vertex_iterator vi, vend;
      tie(vi,vend) = vertices(n);
      double total_rate;

      // choose edge for propagation
      typename graph_traits<network_t>::vertex_descriptor parent, child;
      typedef relative_indicator_t<constant_indicator, fitness_ind_t>
        inv_fitness_ind_t;
      switch(update_rule)
      {
      case Parameters::EGT:
        { indicator_to_map<fitness_ind_t> fitness_map(fitness_ind);
          total_rate = accumulate_transformed(vi,vend,0.0,fitness_ind);
          tie(parent,child) = choose_parent_child_egt(n,total_rate,
                                                      fitness_map);
        }
        break;
      case Parameters::VM:
        { constant_indicator one_ind(1.0);
          inv_fitness_ind_t inv_fitness_ind(one_ind,fitness_ind);
          indicator_to_map<inv_fitness_ind_t> inv_fitness_map(inv_fitness_ind);
          total_rate = accumulate_transformed(vi,vend,0.0,inv_fitness_ind);
          tie(parent,child) = choose_parent_child_vm(n,total_rate,
                                                     inv_fitness_map);
        }
        break;
      case Parameters::SKYE:
        { indicator_to_map<fitness_ind_t> fitness_map(fitness_ind);
          // ?? is this really the total rate ??
          total_rate = accumulate_transformed(vi,vend,0.0,fitness_ind);
          tie(parent,child) = choose_parent_child_skye(n,fitness_map);
        }
        break;
      }
      //cout << parent << ',' << child << ' ';

      // now copy the parent's state to the child.
      if (state[child] == 1)
      { if (state[parent] == 0)
        { state[child] = 0;
          --n_infected;
        }
      }
      else
      { if (state[parent] == 1)
        { state[child] = 1;
          ++n_infected;
        }
      }

      // update counters
      if (n_infected > peak)
        peak = n_infected;
      clock += 1/total_rate;
      // this may be obsolete: if there's a single source, we can get
      // stuck with everything else infected.  extinction is
      // theoretically inevitable, but it can take cosmic amounts of
      // time.  detect this condition and quit.
      if (n_infected == nv - 1)
      { stochastic_edge_weight_map<network_t> wm(n);
        typename graph_traits<network_t>::vertex_iterator it,iend;
        for (tie(it,iend) = vertices(n); it != iend; ++it)
          if (state[*it] == 1)
            if ( temperature(*it,n,wm) == 0 )
            { cout << " temperature of lone holdout is 0\n";
              return fixation_result(TIMEOUT, clock, peak);
            }
      }

#ifdef DISPLAY
      // update display machinery
      if(animate)
      { history.push_back(n_infected);
        fixationDisplay.update(clock,n);
      }
#endif
      if (record_dynamics)
      { dynamics_csv << clock;
        for (int i = -2; i <= 2; ++i)
          for (int j = -2; j <= 2; ++j)
            dynamics_csv << weighted_moment(n,i,j,state);
        dynamics_csv.newRow();
      }
    }// while n_infected

    // we're out of the loop, update everything and return
    //cout << endl;
#ifdef DISPLAY
    if (animate)
      fixationDisplay.updateDisplay(clock);
#endif
    fixation_result
      result(((n_infected == nv) ?
              FIXATION : ((n_infected == 0) ?
                          EXTINCTION : TIMEOUT)),
             clock, peak );
//     if (reallyDisplay
//         && (result.result == FIXATION
//             || result.result == EXTINCTION))
//     { // record history to histogram
//       using namespace boost::lambda;
//       vector<int>&histogram =
//         (result.result == FIXATION
//          ? nc.histogram_success : nc.histogram_failure);
//       for_each(history.begin(), history.end(), ++var(histogram)[_1]);
//       // push updated stats to display
//       transform(histogram.begin(),histogram.end(),
//                 hist_display.global_x_inserter(),
//                 &_1 - &histogram[0]);
//       transform(nc.histogram_success.begin(),
//                 nc.histogram_success.end(),
//                 nc.histogram_failure.begin(),
//                 hist_display.y_inserter(0),
//                 _1 / ret<double>(_1 + _2));
//       hist_display.updateDisplay(0);
//     }

    fixation_stats_t &nc = selind::cache[cache_key(n)];
    fixation_record &vstats = nc.fixations_by_vertex[mutation];
    switch(result.result) {
    case FIXATION:
      ++vstats.n_fixations;
      break;
    case EXTINCTION:
      ++vstats.n_extinctions;
      break;
    }
    ++vstats.n_trials;

    return result;
  }
};

// monte carlo fixation indicator is stochastic
template<typename network_t, typename RNG_t, typename params_t>
bool indicator_is_stochastic(network_selection_indicator<network_t,RNG_t,
                                                         params_t>&i)
{ return true; }

// variant of monte carlo network selection, beginning from the same
// vertex every time.
template<typename network_t, typename RNG_t, typename params_t>
class selection_indicator_fixed_origin
  : public selection_indicator<network_t, RNG_t, params_t>
{
public:
  selection_indicator_fixed_origin
    (typename graph_traits<network_t>::vertex_descriptor v, RNG_t _rng)
      : origin(v), selection_indicator<network_t,RNG_t,params_t>(_rng) {}
  bool try_fixation(const network_t&n)
  { return fixation_once(n,origin);
  }
protected:
  typename graph_traits<network_t>::vertex_descriptor origin;
};

// analytic fixation probability indicators are in here
#include "matrix-math.h"

// given fixation probability of a given graph; assume it's actually
// a point from a function of the form p = (1 - 1/r^K)/(1 - 1/r^KN) --
// if true, K expresses how much the graph amplifies selection since
// the Moran process corresponds to K=1.  So given p, r, and N, we solve
// for K.
// unfortunately, this calculation isn't very reliable, because of
// multiple polynomial roots.  better not to use it.
#include <vnl/algo/vnl_rpoly_roots.h>
#include <vnl/vnl_real_polynomial.h>
#include <boost/lambda/lambda.hpp>
template<typename pt>
class probability_to_amplification
{
  typedef ObjectWithParameters<pt> prob_indicator_t;
  prob_indicator_t *prob_indicator;
public:
  probability_to_amplification(prob_indicator_t*_i)
    : prob_indicator(_i)
  {}
  probability_to_amplification() : prob_indicator(new prob_indicator_t)
  {} 
  // given fitness r, and prob. of fixation p; assuming
  // the form p = (1 - 1/r^K) / (1 - 1/r^(K*N)),
  // what is K?
  // it becomes (a = r^K)
  // (p - 1) a^N + a^(N-1) - p = 0,
  // so solve for a...
  template<typename network_t>
  double operator()(const network_t&n)
  { using namespace boost::lambda;
    typename graph_traits<network_t>::vertices_size_type N
      = num_vertices(n);
    double p = (*prob_indicator)(n);
    double r = prob_indicator->params.mutantFitness();
    double moran = //(1 - 1/r) / (1 - 1/pow(r,N));
      moran_probability(r,N);
    // reverse order: c[0] x^N + ... + c[N]
    vnl_vector<double> poly_coeffs(N+1,0);
    poly_coeffs[0] = p-1;
    poly_coeffs[1] = 1;
    poly_coeffs[N] = -p;
    vnl_real_polynomial P(poly_coeffs);
    vnl_rpoly_roots roots(P);
    vnl_vector<double> realroots(roots.realroots());
    set<double> answers;
    for(int i = 0; i < realroots.size(); ++i)
    { double a = realroots[i];
//       if (print_stuff)
//         cout << "P(" << a << ") = "
// 	     << P.evaluate(a) << endl
// 	     << "(1 - 1/a)/(1 - 1/a^" << N << ") = "
// 	     << (1 - 1/a)/(1 - 1/pow(a,(int)N)) << endl;
      if(a >= 0 && abs(a - 1) > 1e-5)
	answers.insert(a);
    }
    switch (answers.size())
    {
    case 0:
      answers.insert(1);
      break;
    case 1:
      break; // all well
    default:
      cerr << answers.size() << " roots!";
      for (set<double>::iterator ai = answers.begin(); ai != answers.end();
	   ++ai)
	cerr << ' ' << *ai;
      cerr << endl;
      break;
    }
    double a = *answers.rbegin();
    return log(a) / log(r);
  }
};

// indicator for selecting bad graph structures --
// generally just inverts the selection probability but also selects
// against graphs that can't fixate at all for most starting nodes.
template<typename network_t,typename pc>
class bad_at_selection_but_strongly_connected_indicator :
  public analytic_selection_indicator<network_t,pc>
{
public:
  bad_at_selection_but_strongly_connected_indicator(pc*_p)
    : analytic_selection_indicator<network_t,pc>(_p) {}

  double operator()(const network_t&network)
  { map<typename graph_traits<network_t>::vertex_descriptor,
         default_color_type>
      colors;
    // this criterion maybe should be different?
    if (!::is_connected(network,make_assoc_property_map(colors)))
      return 0;
    else
      return 1 - analytic_selection_indicator<network_t,pc>::operator()(network);
  }
};

#if 0
// indicator for selecting bad graph structures --
// generally just inverts the selection probability but also selects
// against graphs that can't fixate at all for most starting nodes.
template<typename network_t>
class bad_at_selection_but_strongly_connected_indicator :
  public analytic_selection_indicator<network_t>
{
public:
  double operator()(const network_t&network)
  {  map<typename graph_traits<network_t>::vertex_descriptor,
       default_color_type>
       colors;
    if (!::is_connected(network,make_assoc_property_map(colors)))
      return 0;
    else
      return 1 - analytic_selection_indicator<network_t>::operator()(network);
  }
};
#endif//0

#endif //__indicators_network_h__
