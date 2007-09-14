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
*/
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/random.hpp>
#include <boost/graph/strong_components.hpp>
#include <numeric>
#include <sstream>
#include "network-optimize.h"
using namespace std;

// represent a graph as a bit string, which is really the incidence
// matrix strung out in a line.  used by selection_indicator.
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
string vbstring(vector<bool>&vb)
{ string vs(vb.size(),'0');
  string::iterator si = vs.begin();
  for(vector<bool>::iterator vi = vb.begin(); vi != vb.end(); ++vi,++si)
    if (*vi)
      *si = '1';
  return vs;
}

// density of a graph is not found in the boost headers
template<typename network_t>
double density(const network_t &n)
{ typename network_t::vertices_size_type nv = num_vertices(n);
  return double(num_edges(n)) / (nv*nv);
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

#if 0
// 2 fns used by is_fixation_candidate
template<typename network_t, typename visited_t>
void dfs(typename graph_traits<network_t>::vertex_descriptor v,
		 visited_t &visited, const network_t&n)
{ visited[v] = true;
  typename graph_traits<network_t>::adjacency_iterator ni, nend;
  for (tie(ni,nend) = adjacent_vertices; ni != nend; ++ni)
    if (!visited[*ni])
      reverse_dfs(*ni,visited,n);
}

template<typename network_t, typename visited_t>
void reverse_dfs(typename graph_traits<network_t>::vertex_descriptor v,
		 visited_t &visited, const network_t&n)
{ visited[v] = true;
  typename inv_adjacency_iterator_generator<network_t>::type ni, nend;
  for (tie(ni,nend) = inv_adjacent_vertices; ni != nend; ++ni)
    if (!visited[*ni])
      reverse_dfs(*ni,visited,n);
}

// before doing fixation simulation call this to evaluate
//  whether there's a possibility of nonzero fixation probability
template<typename network_t>
bool is_fixation_candidate(const network_t &n)
{
  typedef typename graph_traits<network_t>::vertex_descriptor vertex;
  set<vertex> found_sources;
  typename graph_traits<network_t>::vertex_iterator vi,vend;
  // for each vertex, is it a part of an isolated or source component?
  typedef vector<bool> visited_t;
  visited_t visited(num_vertices(n));
  for (tie(vi,vend) = vertices(n); vi != vend; ++vi)
  { visited.assign(false);
    // mark all vertices that you can reach *vi FROM
    reverse_dfs(*vi,visited,n);
    for (visited_t::iterator vsi = visited.begin();
	 vsi != visited.end(); ++vsi)
      if (!*vsi)
	;
    ///???
  }
  
}
#endif//0

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
      // if out edges only, a source, if more than one, no fixation
      else if (n_sources > 0)
	return 0;
      else
	++n_sources;
    }
  }
  // if we get through that without returning 0, try the sophisticated
  // check

  // first get the strongly connected components of the network.
  typedef typename network_t::vertices_size_type valtype;
  //map<typename graph_traits<network_t>::vertex_descriptor,valtype>
  vector<valtype>
    component_map;
  component_map.resize(num_vertices(n));
  valtype ncomponents = strong_components(n,component_map.data());
  if (ncomponents == 1)
    return true;
  // if there's more than one component with no in-edges,
  // or even one isolated component, no fixation.
  vector<bool> has_in_edges(ncomponents), has_out_edges(ncomponents);
  typename graph_traits<network_t>::edge_iterator ei,eend;
  for (tie(ei,eend) = edges(n); ei != eend; ++ei)
  { has_in_edges[component_map[target(*ei,n)]] = true;
    has_out_edges[component_map[source(*ei,n)]] = true;
  }
  bool found_source_component = false;
  for(int i = 0; i < ncomponents; ++i)
    if (!has_in_edges[i])
    { if (!has_out_edges[i])
        return false;
      else if (found_source_component)
	return false;
      else
	found_source_component = true;
    }
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

// this is supposed to give the same as selection_indicator, but by
// calculating the fixation probability, not by monte carlo sampling.
//
// the theory: let s be a state of the network, i.e. a set of vertices
// that are infected by the mutant type, and P(s) be probability of
// reaching fixation of the mutant type from state s.  Then
//   P(s) = sum_t p(transition from s to t) P(t);
// or P(s) = sum_t m_st P(t),
// with P(empty state) = 0, P(full state) = 1; rewrite as
// sum_t n_st = sum_t (delta_st - m_st) P(t) = 0 for all s not empty or
// full;  then n P = v, with v a vector of zeros except 1 in the full
// state's row.  Solve.
// Fixation probability given a single random starting vertex is done by
// averaging those states' P.
//#include <boost/numeric/ublas/matrix_sparse.hpp>
#include "vnl/vnl_sparse_matrix_linear_system.h"
#include "vnl/algo/vnl_lsqr.h"

// used to go recursively through all subsets of the vertex set,
// i.e. through all possible graph states, and apply the visitor
// operator to each state.
template<typename element_t, typename visit_f>
void traverse_sets(set<element_t> &given, set<element_t> &draw_from,
		   visit_f &visitor)
{//  cout << "traverse: ";
//   copy(given.begin(), given.end(), ostream_iterator<element_t>(cout));
//   cout << " | ";
//   copy(draw_from.begin(), draw_from.end(), ostream_iterator<element_t>(cout));
//   cout << endl;
  if (draw_from.empty())
    visitor(given);
  else
  { typename set<element_t>::iterator ii = draw_from.begin();
    element_t i = *ii;
    draw_from.erase(ii);
    ii = given.insert(given.end(),i);
    traverse_sets(given,draw_from,visitor);
    given.erase(ii);
    traverse_sets(given,draw_from,visitor);
    draw_from.insert(draw_from.begin(),i);
  }
}

template<typename int_t>
unsigned long two_to_the(const int_t &i)
{ return 1ul << i; } 

template<typename T>
string asString(vnl_sparse_matrix<T> &m)
{
  ostringstream os;
  m.reset();
  while (m.next())
  { os << '(' << m.getrow() << ',' << m.getcolumn() << ": "
       << m.value() << ")\n";
  }
  return os.str();
}

template<typename network_t>
class analytic_selection_indicator
{
protected:

  typedef typename graph_traits<network_t>::vertex_descriptor vertex_t;

  // this is a visitor operation for traverse_sets (see above): for
  // each state s, put all the terms m_st into the matrix.
  template<typename matrix_t>
  class construct_n_visitor
  {
    const network_t &network;
    matrix_t &n;
  public:
    construct_n_visitor(const network_t &_network, matrix_t &_n)
      : network(_network), n(_n) {}
    
    void operator()(set<vertex_t> &state)
    { // don't do the first + last rows
      if (state.size() == 0 || state.size() == num_vertices(network))
        return;
      // this mapping of state to unique index only works if vertex_t
      // is an integral type
      unsigned long
	state_index = accumulate_transformed(state.begin(),state.end(),
					     0,two_to_the<vertex_t>);
//       if (print_stuff)
// 	cout << '{' << state_index << "}\n";
      double R = num_vertices(network)*residentFitness
  	           + state.size()*(mutantFitness-residentFitness);
      // every edge j -> i contributes to either n_s,s-i, n_ss or
      // n_s,s+i
      typename graph_traits<network_t>::edge_iterator ei,eend;
      for(tie(ei,eend) = edges(network); ei != eend; ++ei)
      { // could maybe optimize this if edges sorted by source?
        vertex_t j = source(*ei,network), i = target(*ei,network);
	typename network_t::degree_size_type 
	  outdeg_j = out_degree(j,network);
        if (state.count(j))
	{ if (state.count(i))
	    n(state_index,state_index) -= mutantFitness / (R*outdeg_j);
	  else
	    n(state_index,state_index + two_to_the(i))
	      -= mutantFitness / (R*outdeg_j);
	}
	else
	{ if (state.count(i))
	    n(state_index,state_index - two_to_the(i))
	      -= residentFitness / (R*outdeg_j);
  	  else
	    n(state_index,state_index) -= residentFitness / (R*outdeg_j);
	}
      }
      // correction for isolated vertices: when parent has no out
      // edges, no change of state
      typename graph_traits<network_t>::vertex_iterator vi,vend;
      for(tie(vi,vend) = vertices(network); vi != vend; ++vi)
	if (out_degree(*vi,network) == 0)
	{ double fitness = state.count(*vi) ? mutantFitness:residentFitness;
	  n(state_index,state_index) -= fitness / R;
	}
      // the sum of those contributions from all edges + isolated
      // nodes makes 1 (per state)
      // so the row sums of n should all be 0 (except 1 for the first
      // and last)
    }
  };

  // store the matrix and vector used in the calculation
  // use sparse matrix?? not sure
  // for now
//   typedef boost::numeric::ublas::mapped_matrix<double> matrix_t;
//   matrix_t n;
//   boost::numeric::ublas::mapped_vector<double> v;
  // trouble with boost (solver is missing)
  // using vnl for now
  typedef vnl_sparse_matrix<double> matrix_t;
  matrix_t n;
  vnl_vector<double> v, P;
  
public:
  void calculate(const network_t &network)
  { typename network_t::vertices_size_type nv = num_vertices(network);
    // caution with sizeable networks!  big matrix!!
    int n_states = two_to_the(nv);
    n = matrix_t(n_states,n_states); // resize and clear
    v = vnl_vector<double>(n_states,0);

    // v is easy: just one entry at the end.
    v[n_states-1] = 1;

    // n is like I said above, I - m.

    // don't know if just assigning to an identity matrix would mess
    // up the sparseness
    for (int i = 0; i < n_states; ++i)
      n(i,i) = 1;
    
    //  m_st = / if t = s+{i}  sum_{j in s, j -> i} r / (R outdegree(j))
    //         | if t = s      sum_{j -> i, i,j same state} f(j)/(R outdeg(j))
    //         \ if t = s-{i}  sum_{j in \s, j -> i} 1 / (R outdeg(j))
    //   with f(j) = r if j infected, 1 if not;
    //   R = sum_{vertices i} f(i)
    // so for each s we need to account for every edge once.
    // this is done by visiting all the states recursively.
    // the visitor iterates through the edges to assemble the above terms.
    construct_n_visitor<matrix_t> visitor(network,n);
    typename graph_traits<network_t>::vertex_iterator vi,vend;
    tie(vi,vend) = vertices(network);
    set<vertex_t> addto, takefrom(vi,vend);
    traverse_sets(addto,takefrom,visitor);
    //traverse_sets(set<vertex_t>(),set<vertex_t>(vi,vend),visitor);
    
    // now n and v are assembled, find the answer!
    //solve(n, v, unknown_storage_tag);
    vnl_sparse_matrix_linear_system<double> system(n, v);
    P = vnl_vector<double>(n_states,0);
    // doesn't seem to do the right thing
    //solver.multiply(v,P);
    vnl_lsqr solver(system);
    solver.minimize(P);
//     if (print_stuff)
//     { print_network(network);
//       cout << asString(n) << '\n' << v << '\n' << P << endl;
//     }
    // if (print_stuff)
    //  solver.diagnose_outcome(cout);
  }
  double operator()(const network_t &network)
  { // do simple graph criteria, in case of zero
    if (!is_fixation_candidate(network))
      return 0;
    // if not, do the full matrix operation
    calculate(network);
    // return mean prob of fixation across singleton states
    typename network_t::vertices_size_type nv = num_vertices(network);
    typename graph_traits<network_t>::vertex_iterator vi,vend;
    double prob = 0;
    for(tie(vi,vend) = vertices(network); vi != vend; ++vi)
    { // the index of a given state is sum_{i in s} 2^i
      unsigned long index = two_to_the(*vi);
      prob += P[index];
    }
    prob /= nv;
    return prob;
  }
};

template<typename network_t, typename RNG_t>
class comparing_indicator
{
  selection_indicator<network_t,RNG_t> si;
  analytic_selection_indicator<network_t> ai;
public:
  comparing_indicator(RNG_t&rng)
    : si(rng), ai() {}
  
  double operator()(const network_t&network)
  { double sv = si(network);
    double av = ai(network);
    if (print_stuff)
      cout << " [" << av << ' ' << sv << "]\n";
    return av;
  }
};

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

