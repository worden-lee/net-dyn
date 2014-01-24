/*  -*- c++ -*-
** 
** network-mutation.h
** 
** Made by 
** Login   <wonder@inducement>
** 
** Started on  Mon Jan 28 14:46:39 2008 
** Last update Mon Jan 28 14:46:39 2008 
*/

#ifndef   	NETWORK_MUTATION_H_
# define   	NETWORK_MUTATION_H_

#include "util.h"
#include <boost/tuple/tuple.hpp>
using namespace boost;

// class hierarchy of network mutations
// these are for internal use by network_mutation_iterator
// that class has a holder class that represents these to the world
template<typename network_t>
class network_mutation
{public:
  typedef pair< typename graph_traits<network_t>::vertex_descriptor,
		typename graph_traits<network_t>::vertex_descriptor >
    possible_edge;
  network_mutation()
  {}
  // 3 virt funs reqd of descendants
  // note this is different from the api reqd of mutation objects
  virtual void operator()(network_t&n) = 0;
  // this has to return a new object that can be delete'd
  virtual network_mutation *undo() = 0;
  virtual string asString() = 0;
};

template<typename network_t>
class toggle_edge_mutation : public network_mutation<network_t>
{public:
  typedef typename network_mutation<network_t>::possible_edge possible_edge;
  possible_edge e;
  bool had; // true if the network has the edge when the operation is done
            // only used in printing
            // this design means asString only works right after the
            // mutation is applied
  toggle_edge_mutation(const possible_edge &_e)
    : e(_e)
  { //typename graph_traits<network_t>::edge_descriptor ed;
    //tie(ed,has) = edge(e.first,e.second,n);
  }
  void operator()(network_t&n)
  { typename graph_traits<network_t>::edge_descriptor ed;
    tie(ed,had) = edge(e.first,e.second,n);
    if (had)
      remove_edge(e.first,e.second,n);
    else
      add_edge(e.first,e.second,n);
  }
  network_mutation<network_t> *undo()
  { return new toggle_edge_mutation(e);
  }
  string asString()
  { return (had ? '-':'+')
      + stringf("(%d->%d)", e.first, e.second);
  }
};

template<typename network_t>
class double_toggle_mutation : public toggle_edge_mutation<network_t>
{ typedef typename toggle_edge_mutation<network_t>::possible_edge possible_edge;
  using toggle_edge_mutation<network_t>::e;
  using toggle_edge_mutation<network_t>::had;
public:
  bool had_backwards;
  double_toggle_mutation(const possible_edge &_e)
    : toggle_edge_mutation<network_t>(_e) {}
  void operator()(network_t&n)
  { swap(e.first,e.second);
    toggle_edge_mutation<network_t>::operator()(n);
    had_backwards = had;
    swap(e.first,e.second);
    toggle_edge_mutation<network_t>::operator()(n);
    //cout << "double toggle?\n";
  }
  network_mutation<network_t> *undo()
  { return new double_toggle_mutation(e);
  }
  string asString()
  { return (had ? (had_backwards ?
                   stringf("-(%d--%d)", e.first, e.second)
                   : stringf("reverse (%d->%d)", e.first, e.second))
            : (had_backwards ?
               stringf("reverse (%d->%d)", e.second, e.first)
               : stringf("+(%d--%d)", e.first, e.second)));
  }
};

// tricky: how to store undo information when removing a vertex
// since it invalidates vertex descriptors and iterators?
// inefficient but for now I'll just take a copy of the network

// this class not to be used directly, but used as an undo operation
// for complex mutations
template<typename network_t>
class restore_network_mutation : public network_mutation<network_t>
{ network_t&before;
public:
  restore_network_mutation(network_t&_n) : before(_n) {}
  void operator()(network_t&n)
  { n = before; }
  network_mutation<network_t> *undo() // this class doesn't undo
  { cerr << "restore_network_mutation::undo shouldn't be called!\n";
    return this;
  }
  string asString()
  { return "(restore network to stored state)"; }
};

// an abstract superclass
template<typename network_t>
class mutation_with_brute_force_undo : public network_mutation<network_t>
{public:
  bool _store_before; // if not set, we won't keep undo information
                      // default is set
  network_t before;
  mutation_with_brute_force_undo(bool _store) : _store_before(_store) {}
  // children should implement this, not ()
  virtual void do_mutation(network_t&n) = 0;
  void operator()(network_t&n)
  { if (_store_before)
      before = n;
    do_mutation(n);
  }
  network_mutation<network_t> *undo()
  { return new restore_network_mutation<network_t>(before);
  }
};

template<typename network_t>
class remove_vertex_mutation
  : public mutation_with_brute_force_undo<network_t>
{public:
  typedef typename graph_traits<network_t>::vertex_descriptor vertex;
  vertex v;
  remove_vertex_mutation(vertex _v, bool _store=true)
    : v(_v), mutation_with_brute_force_undo<network_t>(_store)
  {}
  void do_mutation(network_t&n)
  { clear_vertex(v,n);
    remove_vertex(v,n);
  }
  string asString()
  { return stringf("remove %d", v);
  }
};

template<typename network_t>
class merge_vertices_mutation
  : public mutation_with_brute_force_undo<network_t>
{public:
  typedef typename graph_traits<network_t>::vertex_descriptor vertex;
  vertex v, w;
  merge_vertices_mutation(vertex _v, vertex _w, bool _store=true)
    : v(_v), w(_w), mutation_with_brute_force_undo<network_t>(_store)
  {}
  void do_mutation(network_t&n)
  {
    // all edges in/out of w, duplicate them at v
    typename graph_traits<network_t>::adjacency_iterator ti,tend;
    for(tie(ti,tend) = adjacent_vertices(w,n); ti != tend; ++ti)
      if (*ti != v)
        add_edge(v,*ti,n); 

    typename inv_adjacency_iterator_generator<network_t>::type si,send;
    for(tie(si,send) = inv_adjacent_vertices(w,n); si != send; ++si)
      if (*si != v)
        add_edge(*si,v,n);
    // now get rid of w
    clear_vertex(w,n);
    remove_vertex(w,n);
    //cout << "merge?\n";
  }
  string asString()
  { return stringf("merge %d, %d", v, w);
  }
};

template<typename network_t>
class split_edge_mutation : public network_mutation<network_t>
{public:
  network_t&n;
  typename network_mutation<network_t>::possible_edge e;
  typename graph_traits<network_t>::vertex_descriptor new_vertex;
  split_edge_mutation(typename graph_traits<network_t>::edge_descriptor _e,
                      network_t&_n)
    : e(make_pair(source(_e,_n),target(_e,_n))), n(_n)
  {}
  void operator()(network_t&_n)
  { new_vertex = add_vertex(_n);
    pair<typename graph_traits<network_t>::edge_descriptor,bool>
      ee = edge(e.first,e.second,_n);
    if (ee.second) // should always be true
    { add_edge(e.first,new_vertex,_n);
      add_edge(new_vertex,e.second,_n);
      remove_edge(ee.first,_n);
    }
    ee = edge(e.second,e.first,_n);
    if (ee.second) // not always
    { add_edge(e.second,new_vertex,_n);
      add_edge(new_vertex,e.first,_n);
      remove_edge(ee.first,_n);
    }
    //cout << "split?\n";
  }
  network_mutation<network_t> *undo()
  { return new merge_vertices_mutation<network_t>(e.first,new_vertex,false);
  }
  string asString()
  { return stringf("split %d-%d-%d", e.first, new_vertex, e.second);
  }
};

template<typename network_t>
class attach_new_vertex_mutation : public network_mutation<network_t>
{public:
  typedef typename graph_traits<network_t>::vertex_descriptor vertex;
  vertex new_vertex, attach_to;
  attach_new_vertex_mutation(vertex _a) : attach_to(_a) {}
  void operator()(network_t&_n)
  { new_vertex = add_vertex(_n);
    add_edge(attach_to,new_vertex,_n);
    add_edge(new_vertex,attach_to,_n);
  }
  network_mutation<network_t> *undo()
  { return new remove_vertex_mutation<network_t>(new_vertex,false);
  }
  string asString()
  { return stringf("attach vertex %d at %d", new_vertex, attach_to);
  }
};

template<typename network_t>
class reattach_vertex_mutation
  : public mutation_with_brute_force_undo<network_t>
{ // this one removes an individual from the group and reattaches it
  // to a random vertex
public:
  typedef typename graph_traits<network_t>::vertex_descriptor vertex;
  vertex victim, // this one gets moved
    beneficiary, // this one inherits all its connections
    new_friend;  // it gets reattached to this one
  reattach_vertex_mutation(vertex _v, vertex _b, vertex _f, bool _r=true)
    : mutation_with_brute_force_undo<network_t>(_r),
      victim(_v), beneficiary(_b), new_friend(_f) {}
  void do_mutation(network_t&_n)
  { // all edges in/out of victim, duplicate them at beneficiary
    typename graph_traits<network_t>::adjacency_iterator ti,tend;
    for(tie(ti,tend) = adjacent_vertices(victim,_n); ti != tend; ++ti)
      if (*ti != beneficiary)
        add_edge(beneficiary,*ti,_n); 

    typename inv_adjacency_iterator_generator<network_t>::type si,send;
    for(tie(si,send) = inv_adjacent_vertices(victim,_n); si != send; ++si)
      if (*si != beneficiary)
        add_edge(*si,beneficiary,_n);
    // relocate the victim
    clear_vertex(victim,_n);
    add_edge(victim,new_friend,_n);
    add_edge(new_friend,victim,_n);
    //cout << "merge/reattach?\n";
  }
  string asString()
  { return stringf("merge %d to %d and reattach to %d",
                   victim, beneficiary, new_friend);
  }
};

template<typename network_t>
class merge_and_clone_mutation
  : public mutation_with_brute_force_undo<network_t>
{ // this one takes vertex 'ind' from the group, gives its links to
  // 'recipient' and makes it into a clone into 'donor'
public:
  typedef typename graph_traits<network_t>::vertex_descriptor vertex;
  vertex ind,  // this one gets moved
    recipient, // this one inherits all its connections
    donor;     // it gets copies of all this vertex's connections
  merge_and_clone_mutation(vertex _v, vertex _b, vertex _f, bool _r=true)
    : mutation_with_brute_force_undo<network_t>(_r),
      ind(_v), recipient(_b), donor(_f) {}
  void do_mutation(network_t&_n)
  { // all edges in/out of victim, duplicate them at beneficiary
    typename graph_traits<network_t>::adjacency_iterator ti,tend;
    for(tie(ti,tend) = adjacent_vertices(ind,_n); ti != tend; ++ti)
      if (*ti != recipient)
        add_edge(recipient,*ti,_n); 

    typename inv_adjacency_iterator_generator<network_t>::type si,send;
    for(tie(si,send) = inv_adjacent_vertices(ind,_n); si != send; ++si)
      if (*si != recipient)
        add_edge(*si,recipient,_n);
    // relocate the victim
    clear_vertex(ind,_n);
    // all edges in/out of donor, duplicate them at ind
    //typename graph_traits<network_t>::adjacency_iterator ti,tend;
    for(tie(ti,tend) = adjacent_vertices(donor,_n); ti != tend; ++ti)
      add_edge(ind,*ti,_n); 

    //typename inv_adjacency_iterator_generator<network_t>::type si,send;
    for(tie(si,send) = inv_adjacent_vertices(donor,_n); si != send; ++si)
      add_edge(*si,ind,_n);
    //cout << "merge/clone?\n";

    // could include this edge or not, between the donor and clone
    add_edge(donor,ind,_n);
    add_edge(ind,donor,_n);
  }
  string asString()
  { return stringf("merge %d to %d and clone %d",
                   ind, recipient, donor);
  }
};

// declared for real in Parameters.h
template<typename params_t>
class ObjectWithParameters;


template<typename network_t, typename rng_t, typename params_t>
class network_mutation_generator : public ObjectWithParameters<params_t>
{
public:
  // mutation generators have to have a mutation type
  // also that class-hierarchy has to have a holder - here we go
  class mutation
  { network_mutation<network_t>*mutation_object;
  public:
    mutation(network_mutation<network_t>*_m) : mutation_object(_m)
    {} // we are owner of the object, will deallocate it
    void operator()(network_t&n)
    { (*mutation_object)(n); }
    static mutation undo(mutation&m)
    { return mutation(m.mutation_object->undo()); }
    string asString()
    { return mutation_object->asString(); }
    ~mutation()
    { delete mutation_object; }
  };
  string cached;
  typedef typename network_mutation<network_t>::possible_edge possible_edge;
  typedef vector<possible_edge> edges_t;
  edges_t edges;

  rng_t &rng;

  network_mutation_generator(rng_t&_) : rng(_) {}

  // this iterator is used to enumerate all the possible mutations
  class mutation_iterator : public edges_t::iterator
  {public:
    typedef enum { TOGGLE, MERGE, SPLIT, END_OF_POSSIBILITIES } mut_op;
    mut_op next_op;
    network_t*n;
    mutation_iterator(typename edges_t::iterator sit, network_t&_n)
      : edges_t::iterator(sit), n(&_n)
    {}
    mutation_iterator(typename edges_t::iterator sit)
      : edges_t::iterator(sit), n(0)
    {}
    mutation operator*()
    { possible_edge &p_edge = this->edges_t::iterator::operator*();
      switch(next_op)
      {
      default:
        cerr << "Bad operation in network mutation iterator!\n";
        // fall through
      case TOGGLE:
        return mutation(new toggle_edge_mutation<network_t>(p_edge));
      case MERGE:
        return mutation(new merge_vertices_mutation<network_t>
                        (p_edge.first,p_edge.second,n));
      case SPLIT:
        return mutation(new split_edge_mutation<network_t>(p_edge,n));
      }
    }
    mutation_iterator&operator++()
    { while(1)
        switch(++next_op)
        { // advance through possibilities until a valid one
        default:
        case END_OF_POSSIBILITIES: // skip ahead to next possible_edge
          next_op = TOGGLE;
          edges_t::iterator::operator++();
          break;
        case SPLIT:  // split only if the edge is actually there
          { possible_edge &e = this->edges_t::iterator::operator*();
            typename graph_traits<network_t>::edge_descriptor ed;
            bool has;
            tie(ed,has) = edge(e.first,e.second,n);
            if (has)
              return *this;
            else
              break;
          }
        case MERGE: // merge any 2 vertices
        case TOGGLE: // toggle any pair
          return *this; // this includes the end() case
        }
    }
  };

  mutation_iterator begin(network_t&n)
  { string can = canonical(n);
    if (can != cached)
    { typename network_t::vertices_size_type nv = num_vertices(n);
      if (edges.size() != nv)
      { edges.reserve(nv*nv);
        edges.clear();
        typename graph_traits<network_t>::vertex_iterator fb,fe,tb,te;
        for (tie(fb,fe) = vertices(n); fb!=fe; ++fb)
          for (tie(tb,te) = vertices(n); tb!=te; ++tb)
            if (*tb!=*fb)
              edges.push_back(make_pair(*fb,*tb));
      }
      cached = can;
    }
    // shuffle the candidate edges each time.  this matters because
    //  the first one found that's an improvement will be used, or if
    //  none, the first one found that's neutral.
    random_shuffle(edges.begin(),edges.end());
    return mutation_iterator(edges.begin(),n);
  }
  mutation_iterator end()
  { return mutation_iterator(edges.end());
  }
  mutation random_mutation(network_t&n)
  { // where should this really go?
		typedef boost::variate_generator<rng_t&, boost::uniform_real<> > ur_t;
    ur_t choose_01(rng,uniform_real<>());
    double choice = choose_01();
    double toggle_edge_prob = this->toggle_edge_mutation_probability();
    choice -= toggle_edge_prob;
    if (choice < 0)
    { typename graph_traits<network_t>::vertex_descriptor
        s = random_vertex(n,rng),
        t;
      do
        t = random_vertex(n,rng);
      while(t == s);
      typename graph_traits<network_t>::edge_descriptor ed;
      return mutation(new toggle_edge_mutation<network_t>(possible_edge(s,t)));
    }
    double double_toggle_prob = this->double_toggle_mutation_probability();
    choice -= double_toggle_prob;
    if (choice < 0)
    { typename graph_traits<network_t>::vertex_descriptor
        s = random_vertex(n,rng),
        t;
      do
        t = random_vertex(n,rng);
      while(t == s);
      typename graph_traits<network_t>::edge_descriptor ed;
      return
        mutation(new double_toggle_mutation<network_t>(possible_edge(s,t)));
    }
    double remove_prob = this->remove_vertex_mutation_probability();
    choice -= remove_prob;
    if (choice < 0)
    { typename graph_traits<network_t>::vertex_descriptor
        v = random_vertex(n,rng);
      return mutation(new remove_vertex_mutation<network_t>(v));
    }
    double reattach_prob = this->reattach_vertex_mutation_probability();
    choice -= reattach_prob;
    if (choice < 0)
    { typename graph_traits<network_t>::vertex_descriptor
        v = random_vertex(n,rng), b, f;
      do // b can't be v
        b = random_vertex(n,rng);
      while(b == v);
      do // f can't be v but it can be b
        f = random_vertex(n,rng);
      while(f == v);
      return mutation(new reattach_vertex_mutation<network_t>(v,b,f));
    }
    double clone_prob = this->merge_and_clone_mutation_probability();
    choice -= clone_prob;
    if (choice < 0)
    { typename graph_traits<network_t>::vertex_descriptor
        v = random_vertex(n,rng), b, f;
      do // b can't be v
        b = random_vertex(n,rng);
      while(b == v);
      do // f can't be v but it can be b
        f = random_vertex(n,rng);
      while(f == v);
      return mutation(new merge_and_clone_mutation<network_t>(v,b,f));
    }
    double attach_prob = this->attach_vertex_mutation_probability();
    choice -= attach_prob;
    if (choice < 0)
    { return mutation(new attach_new_vertex_mutation<network_t>
                        (random_vertex(n,rng)));
    }
    double merge_prob = this->merge_vertices_mutation_probability();
    //double split_prob;
    choice -= merge_prob;
    typename graph_traits<network_t>::edge_descriptor e = random_edge(n,rng);
    if (choice < 0)
      return mutation(new merge_vertices_mutation<network_t>(source(e,n),
                                                             target(e,n)));
    else
      return mutation(new split_edge_mutation<network_t>(e,n));
  }
};

#endif 	    /* !NETWORK_MUTATION_H_ */
