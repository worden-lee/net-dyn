#ifndef MATRIX_MATH_H
#define MATRIX_MATH_H
//#include <boost/numeric/ublas/matrix_sparse.hpp>
#include "vnl/vnl_sparse_matrix_linear_system.h"
#include "vnl/algo/vnl_sparse_lu.h"
#include "vnl/algo/vnl_lsqr.h"

// analytic calculations of fixation probability

// up front, general utility functions

// all through the file, a state, i.e. a set of graph vertices, is
// mapped to an index by the function index = sum_{v in s} 2^v, since
// vertices are represented by ints in [0 .. n-1].
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

// function to go recursively through all subsets of the vertex set,
// i.e. through all possible graph states, and apply the visitor
// operation to each state.
template<typename element_t, typename visit_f>
void traverse_sets(set<element_t> &given, set<element_t> &draw_from,
		   visit_f &visitor)
{
// cout << "traverse: ";
// copy(given.begin(), given.end(), ostream_iterator<element_t>(cout));
// cout << " | ";
// copy(draw_from.begin(), draw_from.end(), ostream_iterator<element_t>(cout));
// cout << endl;
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

// ==== first discrete-time markov process calculations ====
// skip ahead for continuous time version which is probably better

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

// this is a visitor operation for traverse_sets (see above): for
// any given state s, put all the terms m_st into the matrix.
//
// first the one for when we choose an edge by parent
// 
//  m_st = / if t = s+{i}  sum_{j in s, j -> i} r / (R outdegree(j))
//         | if t = s      sum_{j -> i, i,j same state} f(j)/(R outdeg(j))
//         \ if t = s-{i}  sum_{j in \s, j -> i} 1 / (R outdeg(j))
//   with f(j) = r if j infected, 1 if not;
//   R = sum_{vertices i} f(i)
// so for each s we need to account for every edge once.
// this is done by visiting all the states recursively.
// the matrix to be constructed is n = I - m.
// first the matrix is set to an identity matrix; then
// the visitor iterates through the edges to sutract the above terms.
template<typename network_t, typename matrix_t>
class construct_n_visitor_by_parent
{
  const network_t &network;
  matrix_t &n;
 public:
  construct_n_visitor_by_parent(const network_t &_network, matrix_t &_n)
    : network(_network), n(_n) {}

  typedef typename graph_traits<network_t>::vertex_descriptor vertex_t;
  
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

// now the visitor for when we choose edges by child
//
//  m_st = / if t = s+{i}  sum_{j in s, j -> i} r / (n R(i))
//         | if t = s      sum_{j -> i, i,j same state} f(j)/(n R(j))
//         \ if t = s-{i}  sum_{j in \s, j -> i} 1 / (n R(j))
//   with f(j) = r if j infected, 1 if not;
//   R(i) = sum_{j->i} f(j)
// so for each s we need to account for every edge once; in this case,
//  for each vertex, account for each of its in-edges.
// the matrix to be constructed is n = 1 - m.
// first the matrix is set to an identity matrix; then
// the visitor iterates through the edges to sutract the above terms.
unsigned long watch_state = -1;
template<typename network_t, typename matrix_t>
class construct_n_visitor_by_child
{
  const network_t &network;
  matrix_t &n;
public:
  construct_n_visitor_by_child(const network_t &_network, matrix_t &_n)
    : network(_network), n(_n) {}
    
  typedef typename graph_traits<network_t>::vertex_descriptor vertex_t;

  void operator()(set<vertex_t> &state)
  { typename network_t::vertices_size_type nv = num_vertices(network);
    // don't do the first + last rows
    if (state.size() == 0 || state.size() == nv)
      return;
    // this mapping of state to unique index only works if vertex_t
    // is an integral type
    unsigned long
      state_index = accumulate_transformed(state.begin(),state.end(),
					   0,two_to_the<vertex_t>);

#define DO(nn,s0,s1) \
if (print_stuff && state_index == watch_state) \
{ cout << #nn "(" << s0 << ',' << s1 \
       << ") = " << nn(s0,s1) << '\n'; \
}
    DO(n,state_index,state_index)

    typename graph_traits<network_t>::vertex_iterator ii,iend;
    for(tie(ii,iend) = vertices(network); ii != iend; ++ii)
    { typename inv_adjacency_iterator_generator<network_t>::type ji, jend;
      tie(ji,jend) = inv_adjacent_vertices(*ii,network);
      if (ji == jend)
      { DO(n,state_index,state_index);
        n(state_index,state_index) -= 1.0/nv;
	DO(n,state_index,state_index);
      }
      else
      { double Ri = 0;
        for (; ji != jend; ++ji)
	{ Ri += ((state.count(*ji) > 0) ? mutantFitness : residentFitness); }
        for (tie(ji,jend) = inv_adjacent_vertices(*ii,network);
	     ji != jend; ++ji)
	{ if (state.count(*ji))
  	  { if (state.count(*ii))
	    { DO(n,state_index,state_index);
	      n(state_index,state_index) -= mutantFitness / (nv*Ri);
	      DO(n,state_index,state_index);
	    }
	    else
	    { DO(n,state_index,state_index + two_to_the(*ii));
	      n(state_index,state_index + two_to_the(*ii))
		-= mutantFitness / (nv*Ri);
	      DO(n,state_index,state_index + two_to_the(*ii));
	    }
	  }
	  else
	  { if (state.count(*ii))
	    { DO(n,state_index,state_index - two_to_the(*ii));
	      n(state_index,state_index - two_to_the(*ii))
		-= residentFitness / (nv*Ri);
	      DO(n,state_index,state_index - two_to_the(*ii));
	    }
	    else
	    { DO(n,state_index,state_index);
	      n(state_index,state_index) -= residentFitness / (nv*Ri);
	      DO(n,state_index,state_index);
	    }
	  }
	}
      }
    }
    // the sum of those contributions from all edges + isolated
    // nodes makes 1 (per state)
    // so the row sums of n should all be 0 (except 1 for the first
    // and last)
  }
};

// ==== now the calculations in continuous time using generator matrix ====

template<typename network_t, typename matrix_t>
class construct_generator_visitor
{
protected:
  const network_t &network;
  matrix_t &G;

  typedef typename graph_traits<network_t>::vertex_descriptor vertex_t;
  
public:
  construct_generator_visitor(const network_t &_network, matrix_t &_G)
    : network(_network), G(_G) {}
    
  // () called once for each state: construct that state's row in the
  // matrix, all rates and the diagonal term.
  void operator()(set<vertex_t> &state)
  { unsigned long state_index
      = accumulate_transformed(state.begin(), state.end(),
			       0, two_to_the<vertex_t>);
    double v = 0.0;
    typename graph_traits<network_t>::edge_iterator ei, eend;
    for (tie(ei,eend) = edges(network); ei != eend; ++ei)
    { bool source_infected = (state.count(source(*ei,network)) > 0);
      bool target_infected = (state.count(target(*ei,network)) > 0);
      if (source_infected && !target_infected)
      { double rij = r(*ei,state);
        if (rij != 0.0)
	{ G(state_index, state_index + two_to_the(target(*ei,network)))
	    += rij;
  	  v += rij;
	}
      }
      else if (!source_infected && target_infected)
      { double rij = r(*ei,state);
        if (rij != 0.0)
	{ G(state_index, state_index - two_to_the(target(*ei,network)))
	    += rij;
  	  v += rij;
	}
      }
    }
    if (v != 0.0)
      G(state_index,state_index) = -v;
  }
  double r(typename graph_traits<network_t>::edge_descriptor e,
	   set<vertex_t> &s)
  { if (choose_parent_first)
      return by_parent_edge_rate(e,network,s);
    else
      return by_child_edge_rate(e,network,s);
  }
};

// a state is a set recording which nodes have mutant fitness
template<typename node_t, typename state_t>
double node_fitness(node_t i, state_t &state)
{ return (state.count(i) > 0) ? mutantFitness : residentFitness;
}

// this is the rate of the 'exponential clock' attached to edge e
// for parent-first process: rate(i->j) = f(i)/outdegree(i)
template<typename network_t, typename state_t>
double by_parent_edge_rate(typename graph_traits<network_t>::edge_descriptor e,
			   network_t &n, state_t &s)
{ double fi = node_fitness(source(e,n),s);
  return fi / out_degree(source(e,n),n);
}

// for child-first process: rate(i->j) = f(i)/sum_{k->j}f(k)
template<typename network_t, typename state_t>
double by_child_edge_rate(typename graph_traits<network_t>::edge_descriptor e,
			  const network_t &n, state_t &s)
{ double fi = node_fitness(source(e,n),s);
  double total_in_fitness = 0;
  typename inv_adjacency_iterator_generator<network_t>::type ki, kend;
  for(tie(ki,kend) = inv_adjacent_vertices(target(e,n),n);
      ki != kend; ++ki)
    total_in_fitness += node_fitness(*ki,s);
  return fi / total_in_fitness;
}

// ==== before we get to the indicator classes, we need a function
//      that vnl_sparse_matrix, for some reason, doesn't provide :-) ====

// have to use a subclass to get access to the protected members.
template<typename T>
class vnl_sparse_matrix_utils : public vnl_sparse_matrix<T>
{
public:
//   static void extract_submatrix(const vnl_sparse_matrix<T> &M,
// 				vnl_sparse_matrix<T> &N,
// 				unsigned long row0, unsigned long col0,
// 				unsigned long rows, unsigned long cols)
//   { N.set_size(rows, cols);
//     M.reset();
//     // step through the nonzero matrix entries using the internal iterator,
//     // see vnl_sparse_matrix.h
//     while (M.next())
//     { int r = M.getrow();
//       int c = M.getcolumn();
//       T v = M.value();
//       if (r >= row0 && r < row0 + rows &&
// 	  c >= col0 && c < col0 + cols &&
// 	  v != 0)
// 	N(r - row0, c - col0) = v;
//     }
//   }

//   static void extract_subcolumn(const vnl_sparse_matrix<T> &M,
// 				vnl_vector<T> &N, unsigned long col,
// 				unsigned long row0, unsigned long rows)

  // do them together, way more efficient
  static void extract_submatrix_and_subcolumn(const vnl_sparse_matrix<T> &M,
					      vnl_sparse_matrix<T> &N,
					      unsigned long N_row0,
					      unsigned long N_col0,
					      unsigned long N_rows,
					      unsigned long N_cols,
					      vnl_vector<T> &V,
					      unsigned long V_col,
					      unsigned long V_row0,
					      unsigned long V_rows)
  { N.set_size(N_rows, N_cols);
    V.set_size(V_rows);
    unsigned int v_ctr = 0;
    // step through the nonzero matrix entries using the internal iterator,
    // see vnl_sparse_matrix.h, and take the ones we want.
    M.reset();
    while (M.next())
    { int r = M.getrow();
      int c = M.getcolumn();
      T v = M.value();
      if (r >= N_row0 && r < N_row0 + N_rows &&
	  c >= N_col0 && c < N_col0 + N_cols &&
	  v != 0)
	N(r - N_row0, c - N_col0) = v;
      if (c == V_col && r >= V_row0 && r < V_row0 + V_rows)
      { while (v_ctr < r - V_row0)
	  V[v_ctr++] = 0;
        V[r - V_row0] = v;
	v_ctr = r - V_row0 + 1;
      }
    }
    while(v_ctr < V_rows)
      V[v_ctr++] = 0;
  }
};

// ==== now the indicator classes that use the above ====

template<typename network_t>
class analytic_selection_indicator_discrete
{
protected:

  typedef typename graph_traits<network_t>::vertex_descriptor vertex_t;

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
    unsigned long n_states = two_to_the(nv);
    n = matrix_t(n_states,n_states); // resize and clear
    v = vnl_vector<double>(n_states,0);

    // v is easy: just one entry at the end.
    v[n_states-1] = 1;

    // n is like I said above, I - m.

    // don't know if just assigning to an identity matrix would mess
    // up the sparseness
    for (unsigned long i = 0; i < n_states; ++i)
      n(i,i) = 1;
    
    typename graph_traits<network_t>::vertex_iterator vi,vend;
    tie(vi,vend) = vertices(network);
    set<vertex_t> addto, takefrom(vi,vend);
    if (choose_parent_first)
    { construct_n_visitor_by_parent<network_t,matrix_t> visitor(network,n);
      traverse_sets(addto,takefrom,visitor);
    }
    else
    { construct_n_visitor_by_child<network_t,matrix_t> visitor(network,n);
      traverse_sets(addto,takefrom,visitor);
    }
    
    // now n and v are assembled, find the answer!
#if 1
    //solve(n, v, unknown_storage_tag);
    vnl_sparse_matrix_linear_system<double> system(n, v);
    P = vnl_vector<double>(n_states,0);
    // doesn't seem to do the right thing
    //solver.multiply(v,P);
    vnl_lsqr solver(system);
    solver.minimize(P);
//     if (print_stuff)
//     { cout << asString(n) << v << '\n' << P << endl;
//       solver.diagnose_outcome(cout);
//     }
#else
    vnl_sparse_lu solver(n,vnl_sparse_lu::estimate_condition);
    solver.solve(v,&P);
    if (print_stuff)
    { cout << asString(n) << v << '\n' << P << endl;
      cout << "det " << solver.determinant()
	   << " rcond " << solver.rcond()
	   << " upbnd " << solver.max_error_bound() << endl;
    }
#endif
  }
  double operator()(const network_t &network)
  { // if (print_stuff)
//     { print_network(network); cout << flush; } 
    // do simple graph criteria, in case of zero
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

template<typename network_t>
class analytic_selection_indicator
{
protected:
  // store the generator
  typedef vnl_sparse_matrix<double> matrix_t;
  matrix_t G;
  // and intermediate calculations
  matrix_t Q;
  vnl_vector<double> B;
  // and the answer
  // p[state_index-1] is the negative of prob of fixation from that state
  // not including the empty and full states
  vnl_vector<double> p;

public:
  void construct_matrix(const network_t &network)
  { typename network_t::vertices_size_type nv = num_vertices(network);
    // caution with sizeable networks!  big matrix!!
    unsigned long n_states = two_to_the(nv);
    G = matrix_t(n_states,n_states); // resize and clear

    typename graph_traits<network_t>::vertex_iterator vi,vend;
    tie(vi,vend) = vertices(network);
    set<typename graph_traits<network_t>::vertex_descriptor>
      addto, takefrom(vi,vend);
    construct_generator_visitor<network_t,matrix_t> visitor(network,G);
    traverse_sets(addto,takefrom,visitor);
  }

  double compute_p(const network_t &network)  
  {
    construct_matrix(network);

    // now G is constructed.
    // see the wiki for derivation: the matrix decomposes as
    //     [ 0 0 0 ]
    // G = [ A Q B ] and the fixation probability is 
    //     [ 0 0 0 ]
    // p_N(infinity) = p_N(0) - p_c(0) Q^-1 B
    // where p(t) = [ p_0(t) p_c(t) p_N(t) ]
    // and p(0) is the initial distribution

    // here we take p_N(0) = 0 and p_c(0) has 1/n for each singleton state
    // generate Q^-1 B and do the last bit informally
    vnl_sparse_matrix_utils<double>::
      extract_submatrix_and_subcolumn(G,Q,1,1,G.rows()-2,G.cols()-2,
				      B,G.cols()-1,1,G.cols()-2);
    p.set_size(B.size());
#if 1
    vnl_sparse_matrix_linear_system<double> system(Q,B);
    vnl_lsqr solver(system);
    solver.minimize(p);
#else
    vnl_sparse_lu solver(Q,vnl_sparse_lu::estimate_condition);
    solver.solve(B,&p);
    if (0 && print_stuff)
    { //cout << asString(Q);
      cout << "B: " << B << "\np: " << p << endl;
      cout << "det " << solver.determinant()
	   << " rcond " << solver.rcond()
	   << " upbnd " << solver.max_error_bound() << endl;
    }
    if (solver.determinant() == 0.0)
      cerr << "Q is singular!!" << endl;
#endif
  }

  virtual double operator()(const network_t &network)
  {
    if (!is_fixation_candidate(network))
      return 0;

    compute_p(network);

    double pfix = 0;
    typename network_t::vertices_size_type nv = num_vertices(network);
    for(unsigned long singleton = 1; singleton < G.rows(); singleton <<= 1)
    { pfix -= p[singleton - 1] / nv;
    }

    return pfix;
  }
};

template<typename network_t>
class analytic_selection_indicator_fixed_origin
  : public analytic_selection_indicator<network_t>
{
public:
  unsigned long start_state;

  analytic_selection_indicator_fixed_origin(
    typename graph_traits<network_t>::vertex_descriptor start_vertex)
    : start_state(two_to_the(start_vertex))
  {} 

  virtual double operator()(const network_t &network)
  { if (!is_fixation_candidate(network))
      return 0;

    compute_p(network);

//     print_network(network);

//     cout << "Q:\n" << asString(analytic_selection_indicator<network_t>::Q)
// 	 << "B:\n" << analytic_selection_indicator<network_t>::B << endl
// 	 << "p: " << analytic_selection_indicator<network_t>::p << endl;

    return - analytic_selection_indicator<network_t>::p[start_state-1];
  }
};

#endif
