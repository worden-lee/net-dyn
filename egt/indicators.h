/*
 * indicator classes for network optimization
 *
 * an indicator is anything that you can do this to:
 *   double value = my_indicator(my_network)
 * in other words it needs to have an operator() that takes a network
 *  and returns a double.
 *
 * you give an indicator to improveNetwork(), and it finds a network
 *  that optimizes the indicator.
 *
 * indicators are also used to gather statistics on the networks,
 *  since each statistical measurement is an indicator.
*/
#ifndef __indicators_h__
#define __indicators_h__
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/random.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/betweenness_centrality.hpp>
#include <boost/vector_property_map.hpp>
#include <numeric>
#include <sstream>
#include "network-optimize.h"
using namespace boost;
using namespace std;

// represent a graph as a bit string, which is really the incidence
// matrix strung out in a line.  used for caching results.
template<typename network_t>
vector<bool> canonical(const network_t&n)
{ typename network_t::vertices_size_type nv = num_vertices(n);
  vector<bool> it(nv*nv,false);
  typename graph_traits<network_t>::edge_iterator ei,ee;
  for(tie(ei,ee) = edges(n); ei!=ee; ++ei)
    it[nv*source(*ei,n) + target(*ei,n)] = true;
  return it;
}

// write out a vector of bits as a string of 0's and 1's.
// useful for debugging caching.
string vbstring(vector<bool>&vb)
{ string vs(vb.size(),'0');
  string::iterator si = vs.begin();
  for(vector<bool>::iterator vi = vb.begin(); vi != vb.end(); ++vi,++si)
    if (*vi)
      *si = '1';
  return vs;
}

template<typename subject_t, typename ret_t = double,
	 typename indicator_t = ret_t(*)(const subject_t&)>
class deterministic_indicator_caching
{
  indicator_t &indicator;
  typedef map<vector<bool>,ret_t> cache_t;
  cache_t cache;
public:
  // memory leak here
  deterministic_indicator_caching() : indicator(*new indicator_t())
  {}
  deterministic_indicator_caching(indicator_t&_ind) : indicator(_ind)
  {}
  ret_t operator()(const subject_t&subject)
  { vector<bool> key = canonical(subject);
    typename cache_t::iterator past = cache.find(key);
    if (past == cache.end())
    { ret_t &answer = cache[key] = indicator(subject);
    //cout << answer << '\n';
      return answer;
    }
    else
    { //cout << past->second << '\n';
      return past->second;
    }
  }
};

template<typename network_t>
double zero_indicator(const network_t&n)
{ return 0;
}

// density of a graph is not found in the boost headers
template<typename network_t>
double density(const network_t &n)
{ typename network_t::vertices_size_type nv = num_vertices(n);
  return double(num_edges(n)) / (nv*(nv-1));
}

template<typename ind_t, typename ref_t>
class relative_indicator_t
{ ind_t &ind;
  ref_t &ref;
public:
  relative_indicator_t(ind_t&_ind,ref_t&_ref) : ind(_ind), ref(_ref)
  {}
  template<typename network_t>
  double operator()(const network_t&n)
  { return ind(n) / (double)ref(n);
  }
};

template<typename ind_t, typename ref_t>
relative_indicator_t<ind_t,ref_t>
relative_indicator(ind_t&_i,ref_t&_r)
{ return relative_indicator_t<ind_t,ref_t>(_i,_r);
}

class constant_indicator
{ double c;
public:
  constant_indicator(double d) : c(d)
  {}
  template<typename network_t>
  double operator()(const network_t&n)
  { return c; } 
};
    
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

//this is modeled after identity_property_map
template<typename _key_t, typename _value_t>
struct constant_property_map
  : public put_get_helper<_value_t, constant_property_map<_key_t,_value_t> >
{ typedef _key_t key_type;
  typedef _value_t value_type;
  typedef _key_t reference;
  typedef boost::readable_property_map_tag category;
  value_type value;
  constant_property_map(value_type _v) : value(_v) {}
  inline value_type operator[](const key_type& _t) const
  { return value; }
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

// this gives the median of a range of sorted values.
// it makes unstated assumptions about the iterator type.
template<typename iterator_t>
double median(iterator_t begin, unsigned n)
{ for (unsigned i = 0; 2*i < n; ++i)
    ++begin;
  double med = *begin;
  if (n % 2 == 1)
    med = (med + *++begin)/2.0;
  return med;
}

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

// an extension of the std::accumulate idea
// accumulates a transform of the values in [first,last) using operator+,
// starting with the value init.
// i.e. returns init + f(*first) + f(*(first+1)) + ...
template<typename _Iterator, typename _Rtp, typename _Fn>
_Rtp accumulate_transformed(_Iterator first, _Iterator end, _Rtp init, _Fn f)
{
  for (; first != end; ++first)
    init = init + f(*first);
  return init;
}

// this is like accumulate_transformed but uses [] instead of ()
// i.e. returns init + map[*first] + map[*(first+1)] + ...
template<typename _Iterator, typename _Rtp, typename _Map>
_Rtp accumulate_mapped(_Iterator first, _Iterator end, _Rtp init, _Map map)
{
  for (; first != end; ++first)
    init = init + map[*first];
  return init;
}

// this is generally useful to selection_indicator:
//  choose one of a sequence of things, with probability proportional to
//  some weighting
// requires the type of _Map[*_Iterator] to be able to cast to double
// a positive total_weights can be provided, to speed this up
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

template<typename network_t, typename RNG_t>
class selection_indicator
{
protected:
  class fixation_stats 
  { public:
    int fixations;
    int trials;
    fixation_stats() : fixations(0), trials(0) {}
  };
  typedef map< vector<bool>, fixation_stats > cache_t;
  cache_t cache;

  double cache_trials(const network_t &n, int fixations, int trials)
  { vector<bool> ncan = canonical(n);
    fixation_stats&ncache = cache[ncan];
    ncache.fixations += fixations;
    ncache.trials += trials;
//     cout << vbstring(ncan) << ' '
// 	 << ncache.fixations << ':' << ncache.trials << endl;
    return double(ncache.fixations)/ncache.trials;
  }
  RNG_t rng;
  // reusable state vector for fixation simulation
  vector<double>state;
public:
  selection_indicator(RNG_t _rng): rng(_rng) {}
  
  virtual bool do_one_fixation(const network_t&n)
  { return fixation_once(n,random_vertex(n,rng));
  }

  double operator()(const network_t &n)
  { // first, run analytic test: if the graph isn't connected, or
    // there is more than one source component of the graph, fixation
    // is impossible.
    if (!is_fixation_candidate(n))
      return 0;

    int n_fixations = 0;
    //cout << "testing for fixation: ";
    for (int i = 0; i < timesToTestFixation; ++i)
    { bool fixed = do_one_fixation(n);
    //cout << (fixed ? 'y':'n') << flush;
      if (fixed)
	++n_fixations;
    }
    //cout << '\n';
    return cache_trials(n,n_fixations,timesToTestFixation);
  }

  typedef pair<typename graph_traits<network_t>::vertex_descriptor,
	       typename graph_traits<network_t>::vertex_descriptor>
  parent_child_pair;

  parent_child_pair
  choose_parent_child_by_child(const network_t &n)
  {
    typename network_t::vertices_size_type nv = num_vertices(n);
    typename network_t::vertices_size_type n_infected = 0;
    typename graph_traits<network_t>::vertex_descriptor parent, child;
    while(1) // loop until a successful match found
    {
      // select a child at random
      child = random_vertex(n,rng);
      // select parent via in edges
      typename inv_adjacency_iterator_generator<network_t>::type pi, pend;
      tie(pi,pend) = inv_adjacent_vertices(child,n);
      if (pi == pend) // if child has no in edges, move on
	continue;
      // weighted by fitness
      parent = choose_weighted(pi,pend,state,rng);
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
      
  parent_child_pair
  choose_parent_child_by_parent(const network_t &n, double total_fitness)
  {
    typename graph_traits<network_t>::vertex_descriptor parent,child;
    while(1) // loop until a successful match found
    {
      // select a parent vertex, with probability proportional to fitness
      typename graph_traits<network_t>::vertex_iterator vi,vend;
      tie(vi,vend) = vertices(n);
      parent = choose_weighted(vi,vend,state,rng,total_fitness);
      // now choose a child vertex, with prob. proportional to edge weight
      //  which are all assumed to be 1
      typename network_t::degree_size_type n_children
	= out_degree(parent,n);
      typename graph_traits<network_t>::adjacency_iterator ci,cend;
      if (n_children <= 0)
	continue; // no out edges, start over
      else if (n_children == 1)
      { tie(ci,cend) = adjacent_vertices(parent,n);
        child = *ci;
      }
      else
      { ui_t choose_a_child(rng,uniform_int<>(0,n_children-1));
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
  
  // assign one mutant in a resident population, simulate
  //  whether it fixates (within a given threshold of patience).
  bool fixation_once(const network_t &n,
		     typename graph_traits<network_t>::vertex_descriptor mutation)
  { typename network_t::vertices_size_type nv = num_vertices(n);
    typename network_t::edges_size_type ne = num_edges(n);
    if (nv <= 0 || ne <= 0)
      return false;
    // initial state is all residents but one
    state.resize(nv);
    state.assign(nv,residentFitness);

    state[mutation] = mutantFitness;

    typename network_t::vertices_size_type n_infected = 1;
    int n_steps = 0;
    // do infection process until fixation reached
    while (0 < n_infected && n_infected < nv && n_steps++ < maxStepsToFixation)
    { //cout << n_infected << ' ';
      typename graph_traits<network_t>::vertex_descriptor parent, child;
      if (choose_parent_first)
	tie(parent,child) = choose_parent_child_by_parent(n,
			       nv*residentFitness +
			       n_infected*(mutantFitness-residentFitness));
      else
	tie(parent,child) = choose_parent_child_by_child(n);
      //cout << parent << ',' << child << ' ';
      // now copy the parent's state to the child.
      if (state[child] == mutantFitness)
      { if (state[parent] == residentFitness)
	{ state[child] = state[parent];
	  --n_infected;
        }
      }
      else
      { if (state[parent] == mutantFitness)
        { state[child] = state[parent];
	  ++n_infected;
	}
      }
    }// while n_infected
    //cout << endl;
    return (n_infected == nv);
  }
};

template<typename network_t, typename RNG_t>
class selection_indicator_fixed_origin
  : public selection_indicator<network_t, RNG_t>
{
public:
  selection_indicator_fixed_origin
    (typename graph_traits<network_t>::vertex_descriptor v, RNG_t _rng)
    : origin(v), selection_indicator<network_t,RNG_t>(_rng) {}
  bool do_one_fixation(const network_t&n)
  { return fixation_once(n,origin);
  }
protected:
  typename graph_traits<network_t>::vertex_descriptor origin;
};

// analytic fixation probability indicators are in here
#include "matrix-math.h"

// template<typename network_t, typename RNG_t>
// class comparing_indicator
// {
//   selection_indicator<network_t,RNG_t> si;
//   analytic_selection_indicator<network_t> ai;
// public:
//   comparing_indicator(RNG_t&rng)
//     : si(rng), ai() {}
  
//   double operator()(const network_t&network)
//   { double sv = si(network);
//     double av = ai(network);
//     if (print_stuff)
//       cout << " [" << av << ' ' << sv << "]\n";
// //  if (abs(sv - av) > 0.15)
// //    av = ai(network);
//     return av;
//   }
// };

template<typename i0_t, typename i1_t>
class comparing_indicator
{
  i0_t &i0;
  i1_t &i1;
public:
  comparing_indicator(i0_t &_i0, i1_t &_i1)
    : i0(_i0), i1(_i1) {}

  template<typename network_t>
  double operator()(const network_t&network)
  { double v0 = i0(network);
    double v1 = i1(network);
    if (print_stuff)
      cout << " [" << v0 << ' ' << v1 << "]\n";
//  if (abs(sv - av) > 0.15)
//    av = ai(network);
    return v0;
  }
};

// template<typename network_t, typename RNG_t>
// class selection_comparing_indicator
//   : public comparing_indicator< analytic_selection_indicator<network_t>,
// 				selection_indicator<network_t,RNG_t> >
// {
// public:
//   selection_comparing_indicator(RNG_t&rng)
//     : comparing_indicator< analytic_selection_indicator<network_t>,
// 			   selection_indicator<network_t,RNG_t> >
//        ( analytic_selection_indicator<network_t>(),
// 	 selection_indicator<network_t,RNG_t>(rng) )
//   {}
// };

// template<typename network_t, typename RNG_t>
// class selection_comparing_indicator_fixed_origin
//   : public comparing_indicator<
//       analytic_selection_indicator_fixed_origin<network_t>,
//       selection_indicator_fixed_origin<network_t,RNG_t> >
// {
// public:
//   selection_comparing_indicator_fixed_origin(
//     typename graph_traits<network_t>::vertex_descriptor v, RNG_t&rng)
//     : comparing_indicator<
//         analytic_selection_indicator_fixed_origin<network_t>,
//         selection_indicator_fixed_origin<network_t,RNG_t> >
//           ( analytic_selection_indicator_fixed_origin<network_t>(v),
// 	    selection_indicator<network_t,RNG_t>(v,rng) )
//   {}
// };

#include <boost/graph/graph_utility.hpp>
// some trouble in this fn
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

// given fixation probability of a given graph; assume it's actually
// a point from a function of the form p = (1 - 1/r^K)/(1 - 1/r^KN) --
// if true, K expresses how much the graph amplifies selection since
// the Moran process corresponds to K=1.  So given p, r, and N, we solve
// for K.
#include <vnl/algo/vnl_rpoly_roots.h>
#include <vnl/vnl_real_polynomial.h>
#include <boost/lambda/lambda.hpp>
template<typename prob_indicator_t>
class probability_to_amplification
{
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
    double r = mutantFitness;
    double moran = (1 - 1/r) / (1 - 1/pow(r,12));
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
    return log(a) / log(mutantFitness);
  }
};

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

#endif //__indicators_h__
