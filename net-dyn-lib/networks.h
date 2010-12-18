/*  -*- C++ -*-
** networks.h
** 
** Made by <wonder>
** 
** Started on  Sun Dec 30 19:26:53 2007 
** Last update Sun Dec 30 19:26:53 2007 
*/

#ifndef   	NETWORKS_H_
# define   	NETWORKS_H_

#include <iostream>
#include <vector>
#include <deque>
#include <string>
using namespace std;
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/random.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/graph/plod_generator.hpp>
//#include "powerlaw.h"
#include "molloy-reed.h"
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include "network-mutation.h"
#include "util.h"

// represent a graph as a bit string, which is really the incidence
// matrix strung out in a line.  used for caching results.
template<typename network_t>
vector<bool> canonical_vb(const network_t&n)
{ typename network_t::vertices_size_type nv = num_vertices(n);
  vector<bool> it(nv*nv,false);
  //typename graph_traits<network_t>::edge_iterator ei,ee;
  //  for(tie(ei,ee) = edges(n); ei!=ee; ++ei)
  //    it[nv*source(*ei,n) + target(*ei,n)] = true;
  typename graph_traits<network_t>::vertex_iterator si,send; 
  for (tie(si,send) = vertices(n); si != send; ++si)
  { typename graph_traits<network_t>::out_edge_iterator ei,eend;
    for (tie(ei,eend) = out_edges(*si,n); ei != eend; ++ei)
      it[nv*(*si) + target(*ei,n)] = true;
  }
  return it;
}

template<typename A,  typename C, typename D,
         typename E, typename F, typename G>
string canonical(const boost::adjacency_list<A,boost::vecS,C,D,E,F,G>&n)
{ return vbstring(canonical_vb(n));
}
  
template<typename A,  typename C, typename D,
         typename E, typename F, typename G>
string canonical(const boost::adjacency_list<A,boost::setS,C,D,E,F,G>&n)
{ ostringstream os;
  typedef boost::adjacency_list<A,boost::setS,C,D,E,F,G> graph_t;
  typename graph_traits<graph_t>::edge_iterator ei,eend;
  for (tie(ei,eend) = edges(n); ei != eend; ++ei)
    os << *ei << '\n';
  return os.str();
}

template<typename A, typename B, typename C, typename D,
         typename E, typename F, typename G,
         typename params_t>
void write_snapshot(boost::adjacency_list<A,B,C,D,E,F,G>&n,
                    double time, params_t&params)
{ string basename
    = params.outputDirectory()+"/snapshot-"+double_to_string(time);
  string graphname = basename + ".adjacency";
  ofstream os((basename+".settings").c_str());
  os << "initial_graph_type CUSTOM\n";
  os << "n_vertices " << num_vertices(n) << "\n";
  os << "custom_graph_file " << graphname << "\n";
  os.close();
  os.open(graphname.c_str());
  print_object(n,os);
}

// this edge_inserter is modeled on stl's front_inserter, etc.
// it's good for copying an iterator that gives edges into a graph
// instead of giving that iterator to the graph's constructor
template<typename graph_t>
class edge_insert_iterator
{
  graph_t&graph;
public:
  typedef typename graph_traits<graph_t>::vertices_size_type
    vertices_size_type;

  typedef typename graph_traits<graph_t>::edge_descriptor
                                                      edge_type;
  typedef pair<vertices_size_type,vertices_size_type> other_edge_type;

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

template<typename network_t>
struct custom_types
{
};

// little fn to turn an edge on/off
template<typename graph_t,typename rng_t, typename params_t>
void toggle_edge
(typename network_mutation_generator<graph_t,rng_t,params_t>::possible_edge ft,
   graph_t &n)
{ typename graph_traits<graph_t>::edge_descriptor ed;
  bool ed_found;
  tie(ed,ed_found) = edge(ft.first,ft.second,n);
  if (ed_found)
    remove_edge(ed,n);
  else
    add_edge(ft.first,ft.second,n);
}

// fill graph object with edges.
// assumes num_vertices(n) is already as it should be
// reads from parameters some of:
//  initial_graph_type
//  superstar_m
//  superstar_l
//  graph_density
//  pl_in_exp
//  pl_out_exp
//  pl_alpha
//  pl_beta
//  custom_graph_file
//  reverse_initial_graph
//  flip_one_initial_edge
template<typename network_t, typename params_t, typename RNG_t>
void construct_network(network_t&n, params_t&parameters, RNG_t&rng)
{
  string igt = parameters.initial_graph_type();
  if (igt == "COMPLETE")
  { // complete graph
    typename graph_traits<network_t>::vertex_iterator si,send,ti,tend;
    for (tie(si,send) = vertices(n); si != send; ++si)
      for (tie(ti,tend) = vertices(n); ti != tend; ++ti)
        if (*si != *ti)
          add_edge(*si,*ti,n);
  }
  else if (igt == "STAR")
  { typename graph_traits<network_t>::vertex_iterator ti,tend;
    for (tie(ti,tend) = vertices(n); ti != tend; ++ti)
      if (*ti > 0)
      { add_edge(0,*ti,n);
        add_edge(*ti,0,n);
      }
  }
  else if (igt == "SUPERSTAR")
  { typename graph_traits<network_t>::vertex_iterator vi,vend;
    typedef typename graph_traits<network_t>::vertex_descriptor vertex;
    deque<vertex> q;
    tie(vi,vend) = vertices(n);
    copy(vi,vend,back_inserter(q));
    vertex center = q.front(); q.pop_front();
    // m is # of lobes, and # of leaves per lobe
    // (though there may be a few missing at the end)
    // l is length of chain from leaves to center
    unsigned m = parameters.superstar_m();
    unsigned l = parameters.superstar_l();
    // total size is m(m+l)+1 = m^2 + ml + 1
    if (m == 0)
    { if (l == 0)
      { m = static_cast<unsigned>(sqrt((num_vertices(n)-1)/2));
        l = m;
      }
      else
        for (m = 1; m*(m+l)+1<num_vertices(n); ++m)
          ;
    }
    else if (l == 0)
      l = (num_vertices(n) - 1) / m - m;
    while(!q.empty())
    { vertex collector = q.front(); q.pop_front();
      add_edge(collector,center,n);
      int i;
      for (i = 1; i < l && !q.empty(); ++i)
      { vertex nc = q.front(); q.pop_front();
        add_edge(nc,collector,n);
        collector = nc;
      }
      for (i = 0; i < m && !q.empty(); ++i)
      { vertex leaf = q.front(); q.pop_front();
        add_edge(leaf,collector,n);
        add_edge(center,leaf,n);
      }
      if (i == 0)
        add_edge(center,collector,n);
    }
  }
  else if (igt == "RANDOM")
  { // Erdos-Renyi random graph with given size and density
    typedef boost::erdos_renyi_iterator<boost::minstd_rand, network_t> ERgen;
    double density = parameters.graph_density();
    if (density > 0)
      copy(ERgen(rng, num_vertices(n), density), ERgen(), edge_inserter(n));
  }
  else if (igt == "POWER")
  { // scale-free random graph with specified exponents for
    // in and out degree distributions
#if 0
    typedef powerlaw_iterator_fixed_density<rng_t,network_t> PLGen;
    copy(PLGen(n, rng, num_vertices(n),
               parameters.graph_density(), parameters.pl_in_exp(),
	             parameters.pl_out_exp()),
         PLGen::end(), edge_inserter(n));
#else
    typedef power_law_molloy_reed_iterator<rng_t,network_t> PLGen;
    copy(PLGen(rng, num_vertices(n), parameters.graph_density(),
          parameters.pl_in_exp(), parameters.pl_out_exp()),
         PLGen::end(), edge_inserter(n));
#endif
  }
  else if (igt == "PLOD")
  { // power law out-degree random graph with given size and density
    // (poisson in-degree)
    typedef plod_iterator<rng_t,network_t> PLGen;
    double pl_beta = parameters.pl_beta();
    if (pl_beta <= 0)
    { pl_beta = pow(num_vertices(n),parameters.pl_alpha()+0.5);
      if (parameters.print_stuff())
        cout << "pl_beta: " << pl_beta << '\n';
    }
    copy(PLGen(rng, num_vertices(n),
               parameters.pl_alpha(), pl_beta),
	       PLGen(), edge_inserter(n));
  }
  else if (igt == "CUSTOM")
  { const string *graphfile = parameters.get("custom_graph_file");
    // we assume the num_vertices is set correctly
    if (graphfile && *graphfile != "")
    { ifstream is(graphfile->c_str());
      string line;
      int source = 0;
      while( std::getline((istream&)is,line) )
      { vector<string> entries;
        split(entries, line, boost::algorithm::is_space());
        for (int target = 0; target < entries.size();  ++target)
          if (entries[target] == "1")
            add_edge(source,target,n);
        ++source;
      }
    }       
  }
  else // leave graph empty, I guess    
  { cerr << "construct_network: don't recognize initial_graph_type "
         << igt << endl;
  }

  if (parameters.reverse_initial_graph()) // reverse all the edges
  { typename graph_traits<network_t>::edge_iterator ei, eend;
    tie(ei,eend) = edges(n);
    typedef deque< typename graph_traits<network_t>::edge_descriptor > eq_t;
    eq_t eq;
    copy(ei,eend,back_inserter(eq));
    n.clear();
    while(!eq.empty())
    { typename graph_traits<network_t>::edge_descriptor
        ed = eq.front(); eq.pop_front();
      add_edge(target(ed,n),source(ed,n),n);
    }
    // does this work? [yes]
//     cout << "reversed edges:\n";
//     print_object(n,cout);
  }
  if (parameters.flip_one_initial_edge()) // perturb one edge
  { int v1 = random_vertex(n,rng), v2;
    do v2 = random_vertex(n,rng); while (v2 == v1);
    toggle_edge_mutation<network_t>(make_pair(v1,v2))(n);
  }
}

template<typename network_t>
void print_object(network_t &n, ostream&os)
{
  typename graph_traits<network_t>::vertex_iterator ii,iend,ji,jend;
  for (tie(ii,iend) = vertices(n); ii != iend; ++ii)
  { for (tie(ji,jend) = vertices(n); ji != jend; ++ji)
    { typename graph_traits<network_t>::edge_descriptor ed;
      bool ed_found;
      tie(ed,ed_found) = edge(*ii,*ji,n);
      os << (ed_found ? '1':'0') << ' ';
    }
    os << endl;
  }
}

#endif 	    /* !NETWORKS_H_ */
