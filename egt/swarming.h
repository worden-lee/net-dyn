//  -*- C++ -*-
//
// types for swarming
#include "SwarmingExperimentParameters.h"
#include "indicators-general.h"
#include "indicators-network.h"
#include "networks.h"
#include "util.h"
#include <vnl_matrix.h>
#include <vnl_vector.h>
#include <iostream>
#include <sstream>
#include <boost/multi_array.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/tuple/tuple.hpp>
#include <set>
using namespace std;
using namespace boost;

class swarming_pattern
  : public ObjectWithParameters<SwarmingExperimentParameters>
{
public:
  // unique string for each possible pattern
  virtual string canonical() const = 0;
};

// this swarming pattern is just parameters; but it's useful for it to
// be a class.
class dynamic_lattice_swarming_pattern
  : public swarming_pattern
{ string canonical() const
  { ostringstream os;
    unsigned maxdim = lattice_dimensions();
    os << maxdim << ':';
    for (int i = 0; i < maxdim; ++i)
    { const string *ds = get(string("lattice_dim")+fstring("%u",i));
      if (!ds)
        ds = get("lattice_dim_generic");
      os << ' ' << *ds;
    }
    os << '\n' << n_individuals() << '\n' << lattice_move_probability() << '\n';
    //return "<dynamic lattice swarming pattern>\n";
    return os.str();
  }
};

// implementation of a collections of individuals on a rectangular
// lattice. if they're neighbors, they communicate.
// this object stores a list of `individual' objects, plus the lattice
// with some sites occupied, and also a graph structure recording the
// resultant connections.
//  redesign 2/3/08 lw: no grid structure, just use the locations attached
//  to the individuals
template <std::size_t Dimension>
class lattice_population
  : public ObjectWithParameters<SwarmingExperimentParameters>
{
public:
  class lattice_location : public boost::array<int,Dimension> 
  { typedef boost::array<int,Dimension> base;
  public:
    lattice_location() {}
    lattice_location(const base &_x) :
      base(_x)
    {}
    string asString() const
    { string s;
      for (int i = 0; i < Dimension; ++i)
        s += fstring((i == 0 ? "(%u" : ",%u"), base::at(i));
      return s+")";
    }
    lattice_location operator+(const lattice_location &_y) const
    { lattice_location sum = *this;
      for (int i = 0; i < Dimension; ++i)
        sum[i] += _y[i];
      return sum;
    }
  };
  typedef lattice_location lattice_location_diff;
  typedef lattice_location lattice_size;

  class individual;
  // the individuals are stored in this list
  // it's a list so its iterators can be used as persistent pointers
  typedef list<individual> population_t;
  population_t population;

  typedef typename population_t::iterator individual_indicator;

  // each vertex of the graph has attached to it
  // a pointer to the corresponding individual
  typedef adjacency_list<setS,vecS,undirectedS,individual_indicator> graph_t;

  class individual
  {public:
    double fitness;
    lattice_location location;
    typename graph_traits<graph_t>::vertex_descriptor vertex;
  };

  // this is the dimensions of the lattice
  lattice_size shape;

  // this is true if the lattice is a torus, i.e. if the coordinates
  // wrap around edges
  bool wrap;

//   class lattice_site
//   {public:
//     bool occupied;
//     individual_indicator occupant;
//     lattice_site() : occupied(false) {}
//   };

  // this is a two dimensional lattice structure, occupied sites point
  // to the individuals in the vector
//   typedef multi_array<lattice_site,Dimension> lattice_t;
//   lattice_t lattice;

  // this graph is synchronized with the above, records connections
  // between neighbors
  graph_t graph;

  // this is the list of occupied locations, sorted, for the 'canonical'
  // representation
  typedef std::set<lattice_location> canonical_t;
  canonical_t canonical_list;
  // this is set when the list is out of date (because of a move,
  // insertion, or deletion)
  bool canonical_dirty;

  //default constructor seems okay but explicit = is needed
  lattice_population<Dimension>&operator=(const lattice_population<Dimension>&l)
  { population = l.population;
    // whatta pain
//     array<size_t, Dimension> shape;
//     copy(l.lattice.shape(),l.lattice.shape()+Dimension,shape.begin());
//     lattice.resize(shape);
//     lattice    = l.lattice;
    graph      = l.graph;
    wrap       = l.wrap;
    return *this;
  }

  unsigned distance(lattice_location&x, lattice_location&y)
  { // max dist (infinity norm)
    unsigned md = 0;
    for (int i = 0; i < Dimension; ++i)
    { unsigned d = abs(x[i]-y[i]);
      if (wrap) // modular version of abs
      { if (d > shape[i]/2)
          d = shape[i] - d;
      }
      if (d > md)
        md = d;
    }
    return md;
  }

  bool occupied(lattice_location &lx) const
  { for (typename population_t::const_iterator pi = population.begin();
         pi != population.end(); ++pi)
      if (lx == pi->location)
        return true;
    return false;
  }

  // this function - in fact a lot of this design - requires the
  // population_t be a list (or possibly a set) so that iterators
  // don't go invalid when individuals come or go or move
  void remove_individual(individual_indicator ii)
  { clear_vertex(ii->vertex,graph);
    remove_vertex(ii->vertex,graph);
//     lattice(ii->location).occupied = false;
    population.erase(ii);
    canonical_dirty = true;
  }

  // this doesn't check that the location is unoccupied
  individual_indicator add_individual(lattice_location loc)
  { individual_indicator
      ii = population.insert(population.end(), individual());
//     lattice(loc).occupied = true;
//     lattice(loc).occupant = ii;
    ii->location = loc;
    ii->vertex = add_vertex(graph);
    graph[ii->vertex] = ii;
    record_neighbors(ii,true);
    canonical_dirty = true;
    return ii;
  }

  // this either
  void move_individual(individual_indicator ii, lattice_location loc)
  {//  if (lattice(loc).occupied)
//     { cout << "ERROR asked to move " << ii->location.asString()
//            << " to an occupied location " << loc.asString() << "\n";
//       return;
//     }
//     lattice(ii->location).occupied = false;
    ii->location = loc;
//     lattice(loc).occupied = true;
//     lattice(loc).occupant = ii;
    clear_vertex(ii->vertex,graph);
    record_neighbors(ii,true);
    canonical_dirty = true;
  }

  // find the neighbors of an individual, insert edges where
  // appropriate into the graph.  if `all' is false, it will only look
  // at those to the right of (or directly below) the given location,
  // to avoid redundancy.
  //  that is, in general, at sites for which the first coordinate that
  //  differs from the focal site is bigger.  (actually, we do a few extra.)
  // assumes the `individual' object is completely and correctly initialized
  inline void record_neighbors(const individual_indicator who, bool all)
  { int horizon = this->lattice_neighbor_distance();

    for (individual_indicator ii = population.begin();
         ii != population.end(); ++ii)
      if (ii != who)
      { unsigned dist = distance(ii->location, who->location);
        if (dist == 0)
          cerr << "collision at " << ii->location.asString()
               << " in record_neighbors!\n";
        if (dist <= horizon)
          add_edge(who->vertex, ii->vertex, graph);
      }
    //     const typename lattice_t::size_type *shape = lattice.shape();
    //cout << "record_neighbors " << x << ' ' << y << '\n';

//     // locate neighbors by searching all nearby lattice locations
//     lattice_location initial_neighbor;
//     initial_neighbor[0] = all ? -horizon : 0;
//     for (int i = 1; i < Dimension; ++i)
//       initial_neighbor[i] = -horizon;

//     //for (int dy = -horizon; dy <= horizon; ++dy)
//       //if (0 <= y+dy && y+dy < shape[1])
//       //for (int dx = dx_start_offset; dx <= horizon; ++dx)
//           //if ((dx != 0 || dy != 0)
//     for(lattice_neighborhood_iterator it
//           = lattice_neighborhood_iterator::begin(who->location,*this,horizon,
//                                                  false,initial_neighbor);
//         it != lattice_neighborhood_iterator::end(); ++it)
//     {
//           //{  if (!all && dy <= 0 && dx <= 0)
//                //continue;
//             //cout << "try " << x+dx << ' ' << y+dy;
//             // lattice_site &neighbor = lattice[x+dx][y+dy];
//       bool invalid = false, different = false;
//       for (int i = 0; i < Dimension; ++i)
//         if ((*it)[i] != who->location[i])
//         { different = true;
//           if ((*it)[i] < 0 || (*it)[i] >= shape[i])
//             invalid = true;
//         }
//       if (different && !invalid)
//       { lattice_site &neighbor = lattice(*it);
//         if (neighbor.occupied)
//         { //cout << " *";
//           add_edge(who->vertex, neighbor.occupant->vertex, graph);
//         }
//       //cout << '\n';
//       }
//     }
  }

  void set_wrap(bool _)
  { wrap = _; }

  void set_size(lattice_size sz)
  {//  lattice_population::lattice_t::extent_gen extents;
//     lattice.resize(extents[wid][hgt]);
//     lattice.resize(sz);
//     if (lattice.num_elements() == 0)
//       cerr << "lattice didn't allocate!! danger!!\n";
    shape = sz;
   }

  void clear()
  { population.clear();
//     for (typename lattice_t::element *ep = lattice.data();
//          ep < lattice.data()+lattice.num_elements(); ++ep)
//       ep->occupied = false;
    graph.clear();
    canonical_dirty = true;
  }

  template<typename RNG_t>
  void place_n_random_individuals(unsigned n_to_fill, RNG_t &rng)
  {//  const typename lattice_t::size_type *shape = lattice.shape();
//     ui_t choose_x(rng,uniform_int<>(0,shape[0]-1));
//     ui_t choose_y(rng,uniform_int<>(0,shape[1]-1));
    while(n_to_fill > 0)
    {// int x = choose_x(), y = choose_y();
      lattice_location loc;
      for (int i = 0; i < Dimension; ++i)
      { ui_t choose_xi(rng,uniform_int<>(0,shape[i]-1));
        loc[i] = choose_xi();
      }
//       if (lattice(loc).occupied)
//         continue;
      --n_to_fill;
      individual_indicator ii
        = population.insert(population.end(), individual());
      ii->location = loc;
      ii->vertex = add_vertex(graph);
      graph[ii->vertex] = ii;
//       lattice(loc).occupied = true;
//       lattice(loc).occupant = ii;

      // currently (with just a list and a graph) this iterative
      // record_neighbors is sufficient
      record_neighbors(ii,false);
    }
    //    print_object(graph,cout);
  }

  string canonical() const
  { ostringstream os;
    // it's not really const, sorry
    canonical_t &_cl = const_cast<canonical_t&>(canonical_list);
    if (canonical_dirty || _cl.size() != population.size())
    { _cl.clear();
      for (typename population_t::const_iterator ii = population.begin();
           ii != population.end(); ++ii)
      { _cl.insert(ii->location);
      }
    }
    for (typename canonical_t::iterator ci = _cl.begin();
         ci != _cl.end(); ++ci)
    { os << ci->asString() << ' '; }
    os << '\n';

//     if (Dimension == 2)
//     {//  const typename lattice_t::size_type *shape
//         = lattice.shape();
//       lattice_location loc;
//       for (loc[1] = shape[1] - 1; loc[1] >= 0; --loc[1])
//       { for (loc[0] = 0; loc[0] < shape[0]; ++loc[0])
//         //os << /*(y == 0 ? "" : " ") << */(lattice[x][y].occupied ? 1:0);
//           os << (lattice(loc).occupied ? '*':'.');
//         os << '\n';
//       }
//     }
//     else if (Dimension == 3)
//     { const typename lattice_t::size_type *shape
//         = lattice.shape();
//       lattice_location loc;
//       for (loc[2] = shape[2] - 1; loc[2] >= 0; --loc[2])
//         for (loc[1] = shape[1] - 1; loc[1] >= 0; --loc[1])
//         { for (int i = 0; i < shape[1] - 1 - loc[1]; ++i)
//             os << ' ';
//           for (loc[0] = 0; loc[0] < shape[0]; ++loc[0])
//             //os << /*(y == 0 ? "" : " ") << */(lattice[x][y].occupied ? 1:0);
//             os << (lattice(loc).occupied ? '*':'.');
//           os << '\n';
//         }
//     }
//     else
//     { for (const typename lattice_t::element *ep = lattice.data();
//            ep < lattice.data()+lattice.num_elements(); ++ep)
//         os << (ep->occupied ? '*':'.');
//       os << '\n';
//     }
    return os.str();
  }
};

class static_lattice_swarming_pattern
  : public swarming_pattern
{
public:
  static const unsigned Dimension = 10;
  lattice_population<Dimension> lattice;
  typedef lattice_population<Dimension>::graph_t graph_t;

  static_lattice_swarming_pattern()
  { //lattice.inheritParametersFrom(this);
  }

  template<typename RNG_t>
  void initialize(RNG_t&rng)
  { typename lattice_population<Dimension>::lattice_size size;
    lattice.inheritParametersFrom(this);
    unsigned maxdim = lattice_dimensions();
    unsigned i;
    for(i = 0; i < maxdim; ++i)
    { const string *ds = get(string("lattice_dim")+fstring("%u",i));
      if (!ds)
        ds = get("lattice_dim_generic");
      size[i] = ds ? string_to_unsigned(*ds) : 0;
    }
    for (; i < Dimension; ++i)
      size[i] = 1;
    lattice.set_size(size);
    lattice.clear();
    lattice.set_wrap(lattice_is_torus());
    lattice.place_n_random_individuals(n_individuals(),rng);
    cout << canonical();
  }
  string canonical() const
  { return lattice.canonical();
  }
};

template<typename params_t>
void write_snapshot(static_lattice_swarming_pattern&swp,
                    double time, params_t&params)
{ // deal with later
}

void print_object(static_lattice_swarming_pattern&swp,ostream&os)
{ string can = swp.canonical(); 
  if (can.length() < 1600) // 20 packed lines
    os << can;
}
template<std::size_t Dimension>
class lattice_mutation
{public:
  typedef typename lattice_population<Dimension>::individual_indicator
    individual_iterator;
  typedef typename lattice_population<Dimension>::lattice_location location;
  typedef enum { REMOVE, ADD, MOVE } lattice_operation;
  lattice_population<Dimension> *lat;
  lattice_operation operation;
  individual_iterator ind;
  location loc;
  location _old_location;
  lattice_mutation()
  {}
  lattice_mutation(lattice_population<Dimension> *p,
                   lattice_operation op, individual_iterator i)
    : lat(p), operation(op), ind(i)
  {}
  lattice_mutation(lattice_population<Dimension> *p,
                   lattice_operation op, location l)
    : lat(p), operation(op), loc(l)
  {}
  lattice_mutation(lattice_population<Dimension> *p,
                   lattice_operation op, individual_iterator i, location l)
    : lat(p), operation(op), ind(i), loc(l), _old_location(i->location)
  {}
  void operator()(lattice_population<Dimension>&p)
  { if (&p != lat)
      cout << "lattice_mutation: wrong lattice_population!\n";
    switch (operation)
    {
    case REMOVE:
      p.remove_individual(ind);
      break;
    case ADD:
      ind = p.add_individual(loc);
      break;
    case MOVE:
      p.move_individual(ind,loc);
      break;
    }
  }
  void operator()(static_lattice_swarming_pattern &swp)
  { return operator()(swp.lattice); }
  // this has to be called after the mutation has been applied and
  // before it's undone
  static lattice_mutation undo(const lattice_mutation&m)
  { switch(m.operation)
    {
    case REMOVE:
      return lattice_mutation(m.lat, ADD,    m.ind->location);
    case ADD:
      return lattice_mutation(m.lat, REMOVE, m.ind);
    case MOVE:
      return lattice_mutation(m.lat, MOVE,   m.ind, m._old_location);
    }
  }
  string asString()
  { ostringstream oss;
    switch(operation)
    {
    case ADD:    oss << '+' << loc.asString(); break;
    case REMOVE: oss << '-' << ind->location.asString(); break;
    case MOVE:   oss << _old_location.asString() << "->" << loc.asString();
                 break;
    }
    return oss.str();
  }
};

// lattice_neighborhood_iterator enumerates all the neighbors of a given
// lattice location.
template<unsigned Dimension>
class lattice_neighborhood_iterator
{ typedef lattice_population<Dimension> lattice_population_t;
  typedef typename lattice_population_t::lattice_location lattice_location;
  typedef typename
    lattice_population_t::lattice_location_diff lattice_location_diff;

  lattice_location focus; // neighbors of this location
  lattice_population_t *lattice; // reference to the lattice_population
  int horizon;            // extent of neighborhood
  lattice_location_diff dx;   // current offset
  lattice_location location;  // current value of iterator
  // (invariant: location == focus + dx)
  // uninitialized iterator. for use only by end().
  lattice_neighborhood_iterator()
    : lattice(0) {}
public:
  lattice_neighborhood_iterator(lattice_location&_focus,
                                lattice_population_t*_lattice,
                                int _horizon,
                                bool _wrap,
                                lattice_location_diff&_dx)
    : focus(_focus), lattice(_lattice), horizon(_horizon), dx(_dx)
  { search_from_dx();
  }
  // given dx, which may or may not point to a valid place, advance
  // it as necessary to a valid place, or end(), return iterator with
  // location guaranteed to match, or end().
  // this function takes care of wrapping.  when wrapping occurs dx is
  // not adjusted, only location.
  void search_from_dx()
  { const typename lattice_population_t::lattice_size &shape = lattice->shape;
    for (int i = Dimension - 1; i >= 0; --i)
    { // if dx too small, advance
      int mn = lattice->wrap ? -horizon : std::max(-horizon, -focus[i]);
      int mx = lattice->wrap ? horizon :
        std::min(horizon, int(shape[i]) - 1 - focus[i]);
      if (lattice->wrap && shape[i] == 1) // dimension not in use
        mn = mx = 0;
      if (dx[i] < mn)
        dx[i] = mn;
      // if dx too large, carry to next entry
      if (dx[i] > mx)
      { dx[i] = mn;
        if (i > 0) // carry
          ++dx[i-1];
        else       // no more carry, we've reached the end
          horizon = -1;
      }
      // assign location
      location[i] = focus[i] + dx[i];
      if (lattice->wrap) // wrap location if needed
      { if (location[i] < 0 || location[i] >= shape[i])
          location[i] = mod_correct(location[i],(int)shape[i]);
      }
    }
    // special case: skip dx == 0
    if (location == focus)
    { ++dx[Dimension-1];
      search_from_dx();
    }
  }
  static lattice_neighborhood_iterator begin(lattice_location&focus,
                                             lattice_population_t*lattice,
                                             int horizon,
                                             bool wrap,
                                             lattice_location_diff&init_dx)
  { return lattice_neighborhood_iterator(focus,lattice,horizon,wrap,init_dx);
  }
  static lattice_neighborhood_iterator begin(lattice_location&focus,
                                             lattice_population_t*lattice,
                                             int horizon,
                                             bool wrap)
  { lattice_location_diff dx;
    for (int i = 0; i < Dimension; ++i)
      dx[i] = -horizon;
    return lattice_neighborhood_iterator(focus,lattice,horizon,wrap,dx);
  }
  static lattice_neighborhood_iterator end()
  { lattice_neighborhood_iterator et;
    et.horizon = -1; // this flags an invalid iterator
    return et;
  }
  lattice_neighborhood_iterator &operator++()
  { ++dx[Dimension-1];
    search_from_dx();
    return *this;
  }
  const lattice_location &operator*() const
  { return location; }
  // invalid iterators are equal 
  bool operator==(const lattice_neighborhood_iterator &other)
  { return (horizon < 0 && other.horizon < 0)
      || (horizon == other.horizon && focus == other.focus
          && location == other.location);
  }
  bool operator!=(const lattice_neighborhood_iterator &other)
  { return !operator==(other);
  }
};

template<size_t Dimension, typename RNG_t=void>
class lattice_moves_generator
{public:
  typedef typename lattice_population<Dimension>::lattice_location
    lattice_location;
  typedef lattice_neighborhood_iterator<Dimension>
    neighbor_iterator;
  typedef lattice_mutation<Dimension>                 mutation;
  typedef vector< lattice_mutation<Dimension> >       mutations_t;
  typedef typename mutations_t::iterator              mutation_iterator;

  mutations_t mutations;
  RNG_t &rng;

  // you can use this constructor when begin() .. end() interface is
  // wanted
  lattice_moves_generator() : rng(*(RNG_t*)0) {}
  // but use this constructor when random_mutation() interface is
  // wanted
  lattice_moves_generator(RNG_t&rng) : rng(rng) {}

  mutation_iterator begin(lattice_population<Dimension>&lat)
  { mutations.reserve(8*lat.population.size());
    mutations.clear();
    const typename lattice_population<Dimension>::lattice_size
      &shape = lat.shape;
    for (typename lattice_population<Dimension>::individual_indicator
           ii = lat.population.begin(); ii != lat.population.end(); ++ii)
    { lattice_location &x = ii->location;
      for (neighbor_iterator ni = neighbor_iterator::begin(x,lat,1,false);
           ni != neighbor_iterator::end(); ++ni)
      { bool invalid = false, different = false;
        for (int i = 0; i < Dimension; ++i)
        { if ((*ni)[i] < 0 || (*ni)[i] >= shape[i])
            invalid = true;
          if ((*ni)[i] != x[i])
            different = true;
        }
        if (different && !invalid)
        { if (!lat.lattice(*ni).occupied)
            mutations.push_back(mutation(&lat, mutation::MOVE, ii, *ni));
        }
      }
    }
    random_shuffle(mutations.begin(),mutations.end());
    return mutations.begin();
  }
  mutation_iterator begin(static_lattice_swarming_pattern&swp)
  { return begin(swp.lattice); }
  mutation_iterator end()
  { return mutations.end(); }
  mutation random_mutation(static_lattice_swarming_pattern&swp)
  { return random_move(swp.lattice,rng); }
};

template<std::size_t Dimension, typename RNG_t>
lattice_mutation<Dimension> random_move(lattice_population<Dimension> &lattice,
                                        RNG_t &rng)
{ ui_t choose_dx(rng,uniform_int<>(-1,1));
  const typename lattice_population<Dimension>::lattice_size &shape
    = lattice.shape;
  while(1)
  { // pick an individual and a direction
    typename lattice_population<Dimension>::graph_t::vertex_descriptor rvx
      = random_vertex(lattice.graph,rng);
    typename lattice_population<Dimension>::individual_indicator ind
      = lattice.graph[rvx];
    typename lattice_population<Dimension>::lattice_location_diff dx;
    for (int i = 0; i < Dimension; ++i)
      dx[i] = choose_dx();
    bool nonzero = false;
    for (int i = 0; i < Dimension; ++i)
    { if (dx[i] != 0)
        nonzero = true;
    }
    if (nonzero)
    { typename lattice_population<Dimension>::lattice_location tx
        = ind->location + dx;
      bool invalid = false;
      for (int i = 0; i < Dimension; ++i)
        if (tx[i] < 0 || tx[i] >= shape[i])
        { if (lattice.wrap)
            tx[i] = mod_correct(tx[i], (int)shape[i]);
          else
            invalid = true;
        }
      if (!invalid && !lattice.occupied(tx))
        return lattice_mutation<Dimension>(&lattice,
                                           lattice_mutation<Dimension>::MOVE,
                                           ind, tx);
    }
  }
}
  
// these are of use with BoostDotGraphDisplay
template<std::size_t Dimension, typename network_t>
class vertex_location_map
{ public:
  const network_t&n;
  double scale;
  unsigned nDimensions; // how many to print
  vertex_location_map(const network_t&_n, double _scale, unsigned _nd)
    : n(_n), scale(_scale), nDimensions(_nd) {}
  string operator[](typename graph_traits<network_t>::vertex_descriptor const&v)
  { return '"' + operator()(v,n) + '"'; }
  string operator()(typename graph_traits<network_t>::vertex_descriptor const&v,
                    const network_t&c)
  { typename lattice_population<Dimension>::lattice_location &loc
      = c[v]->location;
    string str;
    for (int i = 0; i < nDimensions; ++i)
      str += (i?",":"") + double_to_string(scale*loc[i]);
    return str;
  }
};

template<std::size_t Dimension, typename network_t>
string get(vertex_location_map<Dimension,network_t> &m,
           typename graph_traits<network_t>::vertex_descriptor const&v)
{ return m[v]; }

//template<std::size_t Dimension>
vertex_location_map< static_lattice_swarming_pattern::Dimension,
                     static_lattice_swarming_pattern::graph_t >
make_vertex_id(const static_lattice_swarming_pattern::graph_t&n,
               double _scale=1,
               unsigned _nDim=static_lattice_swarming_pattern::Dimension)
{ return
    vertex_location_map< static_lattice_swarming_pattern::Dimension,
      static_lattice_swarming_pattern::graph_t >(n,_scale,_nDim);
}

template<std::size_t Dimension, typename network_t>
class place_vertex
{
public:
  const network_t&n;
  double scale;
  vector< vector<double> > offsets;
  place_vertex(const network_t&_n, double _scale)
    : n(_n), scale(_scale), offsets(Dimension)
  { offsets[0].resize(2);
    offsets[0][0] = scale;
    offsets[0][1] = 0;
    offsets[1].resize(2);
    offsets[1][0] = 0;
    offsets[1][1] = scale;
    offsets[2].resize(2);
    offsets[2][0] = scale/8;
    offsets[2][1] = scale/12;
    offsets[3].resize(2);
    offsets[3][0] = - scale/8;
    offsets[3][1] =   scale/8;
    // ?? for 'other' dimensions
    for (int i = 4; i < Dimension; ++i)
    { offsets[i].resize(2);
      offsets[i][0] = - sin(0.1*i)*scale/4;
      offsets[i][1] = - cos(0.1*i)*scale/4;
    }
  }
  string operator()(typename graph_traits<network_t>::vertex_descriptor const&v,
                    const network_t&c)
  { typename lattice_population<Dimension>::lattice_location &loc
      = c[v]->location;
    vector<double> xy(2);
    for (int d = 0; d < Dimension; ++d)
    { xy[0] += loc[d]*offsets[d][0];
      xy[1] += loc[d]*offsets[d][1];
    }
    ostringstream str;
    str << xy[0] << "," << xy[1];
    return str.str();
  }
};

// the rooms thing is seeming unpromising
// in fact I think I can prove it's always isothermal
class rooms_swarming_pattern
  : public swarming_pattern
{
public:
  // number fixed in each room
  // should add up to n_individuals - n_mobile
  vnl_vector<unsigned int> fixed;
  // implementation dependent:
  // prob of switching from which room to which
  virtual double p_switch(int from, int to) const = 0;
};

class complete_swarming_pattern
  : public rooms_swarming_pattern
{
public:
  // prob a mobile individual in room i will switch to room j
  // (given a room switch event, see room_switch_rate)
  // each row should add to 1
  vnl_matrix<double> _p_switch;
  double p_switch(int from, int to) const
  { return _p_switch[from][to]; }
  string canonical() const
  { ostringstream os;
    os << mutantFitness() << ' ' << residentFitness() << '\n'
       << n_rooms() << ' '
       << n_individuals() << ' '
       << n_mobile() << ' ' << fixed << '\n'
       << room_switch_rate() << ' ' << _p_switch;
    return os.str();
  }
};

#include "networks.h"

class graph_swarming_pattern
  : public rooms_swarming_pattern
{
public:
  typedef adjacency_list<setS,vecS,bidirectionalS> network_t;
  network_t n;

  double p_switch(int from, int to) const
  { graph_traits<network_t>::edge_descriptor ed;
    bool ed_found;
    tie(ed,ed_found) = edge(from,to,n);
    if (!ed_found) return 0;
    //else
    static int cache_v = -1; static double cache_temp;
    if(cache_v == from)
      return cache_temp;
    //else
    graph_traits<network_t>::degree_size_type
      outdeg = out_degree(from,n);
    cache_v = from;
    return (cache_temp = 1.0 / outdeg);
  }

  string canonical() const
  { ostringstream os;
    os << mutantFitness() << ' ' << residentFitness() << '\n'
       << n_rooms() << ' '
       << n_individuals() << ' '
       << n_mobile() << ' ' << fixed << '\n'
       << room_switch_rate() << '\n';
    print_object(n, os);
    return os.str();
  }
};

string canonical(const swarming_pattern &sw)
{ return sw.canonical();
}


// this object attaches to a rooms_swarming_pattern, and records
// the momentary state of the swarming process
class swarming_state
{
public:
  // pattern to use
  const rooms_swarming_pattern &pattern;
  // how many mobile ind in each room
  vnl_vector<unsigned int> mobiles;
  // how many infected mobile ind in each room
  vnl_vector<unsigned int> infected_mobiles;
  // how many infected fixed ind in each room
  vnl_vector<unsigned int> infected_fixed;
  // total infected
  unsigned int total_infected;
  // total infected mobile
  unsigned int n_infected_mobile;

  // --- values to cache for quick looping
  int n_mobile, n_individuals, n_rooms;
  double residentFitness, mutantFitness, mob_factor, room_switch_rate;
  bool print_stuff, choose_parent_first;

  void _fill_cache()
  { n_mobile = pattern.n_mobile();
    n_individuals = pattern.n_individuals();
    n_rooms = pattern.n_rooms();
    residentFitness = pattern.residentFitness();
    mutantFitness = pattern.mutantFitness();
    mob_factor =
      pattern.mobilityIncreasesWithFitness() ?
        mutantFitness / residentFitness : 1;
    room_switch_rate = pattern.room_switch_rate();
    print_stuff = pattern.print_stuff();
    choose_parent_first = pattern.choose_parent_first();
  }

  template<typename RNG_t>
  swarming_state(const rooms_swarming_pattern&_pattern, RNG_t &rng)
    : pattern(_pattern), mobiles(_pattern.n_rooms(),0),
      infected_mobiles(_pattern.n_rooms(),0),
      infected_fixed(_pattern.n_rooms(),0), total_infected(0),
      n_infected_mobile(0)
  { _fill_cache();
    // all vectors are created full of zeros, mobiles need to be placed.
    randomize_mobiles(rng);
  }

  template<typename RNG_t>
  void randomize_mobiles(RNG_t&rng)
  { // what distribution to use?
    ui_t choose_room(rng,uniform_int<>(0,n_rooms-1));
    mobiles.fill(0);
    for(unsigned int i = n_mobile; i > 0; --i)
      ++mobiles[choose_room()];
  }

  class individual_descriptor
  {
  public:
    unsigned int room;
    bool fixed;
    bool infected;
    bool exists; // set to false if no individual is represented
    individual_descriptor(int _r, bool _f, bool _i, bool _e=true)
      : room(_r), fixed(_f), infected(_i), exists(_e) {}
    individual_descriptor() {}
  };
  
  template<typename RNG_t>
  individual_descriptor choose_any_individual(RNG_t&rng)
  { ui_t chooser(rng,uniform_int<>(0,n_individuals-1));
    int choice = chooser();
    int seen = 0;
    if (choice < n_mobile)
      for (int i = 0; i < n_rooms; ++i)
      { if (choice < seen + infected_mobiles[i])
          return individual_descriptor(i,false,true);
        else if (choice < seen + mobiles[i])
	  return individual_descriptor(i,false,false);
        else
          seen += mobiles[i];
      }
    else
    { seen += n_mobile;
      for (int i = 0; i < n_rooms; ++i)
	if (choice < seen + infected_fixed[i])
          return individual_descriptor(i,true,true);
        else if (choice < seen + pattern.fixed[i])
	  return individual_descriptor(i,true,false);
        else
          seen += pattern.fixed[i];
    }
    cerr << "error in choose_any_individual" << endl;
  }

  template<typename RNG_t>
  individual_descriptor choose_mobile_individual(RNG_t&rng)
  { ur_t chooser(rng,uniform_real<>(0, n_mobile
                                        + (mob_factor-1)*n_infected_mobile));
    double choice = chooser();
    double seen = 0;
    for (int i = 0; i < n_rooms; ++i)
    { if (choice < seen + mob_factor*infected_mobiles[i])
        return individual_descriptor(i,false,true);
      else if (choice < seen + mobiles[i] + (mob_factor-1)*infected_mobiles[i])
        return individual_descriptor(i,false,false);
      seen += mobiles[i] + (mob_factor-1)*infected_mobiles[i];
    }
    cerr << "error in choose_mobile_individual" << endl;
  }

  template<typename RNG_t>
  individual_descriptor choose_source_individual(RNG_t&rng)
  { double r = mutantFitness / residentFitness;
    ur_t chooser(rng,uniform_real<>(0, n_individuals +(r-1)
                                        *total_infected));
    double choice = chooser();
    double seen = 0;
    // first consider all the mobile individuals
    for (int i = 0; i < n_rooms; ++i)
    { double residentFitnesses
        = mobiles[i] - infected_mobiles[i];
      double mutantFitnesses
        = r*infected_mobiles[i];
      if (choice < seen + mutantFitnesses)
        return individual_descriptor(i,false,true);
      else if (choice < seen + mutantFitnesses + residentFitnesses)
        return individual_descriptor(i,false,false);
      seen += mutantFitnesses + residentFitnesses;              
    }
    // now all the fixed ones
    for (int i = 0; i < n_rooms; ++i)
    { double residentFitnesses
        = pattern.fixed[i] - infected_fixed[i];
      double mutantFitnesses = r*infected_fixed[i];
      if (choice < seen + mutantFitnesses)
        return individual_descriptor(i,true,true);          
      else if (choice < seen + mutantFitnesses + residentFitnesses)
        return individual_descriptor(i,true,false);
      seen += mutantFitnesses + residentFitnesses;
    }
    cerr << "error in choose_source_individual" << endl;
  }

  template<typename RNG_t>
  individual_descriptor choose_in_same_room(individual_descriptor as,
                                            bool weighted,RNG_t&rng)
  { double _r = weighted? mutantFitness / residentFitness : 1;
    unsigned int i = as.room;
    double total_r = (pattern.fixed[i]+mobiles[i]
                      +(_r-1)*(infected_fixed[i]
                               +infected_mobiles[i]));
    if (as.infected) total_r -= _r;
    else             total_r -= 1;
    if (total_r <= 0)
      return individual_descriptor(i,false,false,false);
    ur_t choose_target(rng,uniform_real<>(0,total_r));
    double target = choose_target();
    double seen = 0;
    // is target an infected mobile?
    seen += infected_mobiles[i];
    if (as.infected && !as.fixed)
      --seen;
    if (target < seen)
      return individual_descriptor(i,false,true);
    // is target an uninfected mobile?
    seen += mobiles[i] - infected_mobiles[i];
    if (!as.infected && !as.fixed)
      --seen;
    if (target < seen)
      return individual_descriptor(i,false,false);
    // is target an infected fixed?
    seen += infected_fixed[i];
    if (as.infected && as.fixed)
      --seen;
    if (target < seen)
      return individual_descriptor(i,true,true);
    // is target an uninfected fixed?
    seen += pattern.fixed[i] - infected_fixed[i];
    if (!as.infected && as.fixed)
      --seen;
    if (target < seen)
      return individual_descriptor(i,true,false);
    cerr << "error in choose_in_same_room" << endl;
  }
  
  // this won't work right unless no one is infected yet
  template<typename RNG_t>
  void infect_random_individual(RNG_t&rng)
  { individual_descriptor ind = choose_any_individual(rng);
    if (print_stuff)
      cout << "infect a " << (ind.fixed ? "fixed":"mobile")
           << " individual in room " << ind.room << endl;
    if (ind.fixed)
      ++infected_fixed[ind.room];
    else
    { ++infected_mobiles[ind.room];
      ++n_infected_mobile;
    }
    ++total_infected;
  }
  
  template<typename RNG_t>
  void step(RNG_t&rng)
  { ur_t choose_01(rng,uniform_real<>());
    // room switches happen at rate room_switch_rate (or r times that)
    //  per mobile individual
    // infections happen at rate 1 or r per source, or 1 per target
    double total_switch_rate =
      room_switch_rate * (n_mobile + (mob_factor-1)*n_infected_mobile);
    double total_inf_rate =
      (choose_parent_first
       ? (residentFitness * n_individuals
          + (mutantFitness-residentFitness)*total_infected)
       : n_individuals);
    double p_switch
      = total_switch_rate / (total_switch_rate + total_inf_rate);
    // either try a room switch
    if (choose_01() < p_switch)
    { individual_descriptor ind = choose_mobile_individual(rng);
      double move_choice = choose_01();
      double seen = 0.0;
      for (int j = 0; j < n_rooms; ++j)
      { seen += pattern.p_switch(ind.room,j);
        if (move_choice < seen)
	{ if (ind.infected)
          { --infected_mobiles[ind.room];
	    ++infected_mobiles[j];
	  }
          --mobiles[ind.room];
	  ++mobiles[j];
          if (print_stuff)
            cout << "switch " << ind.room << "->" << j
                 << (ind.infected?" (infected)":" (uninfected)") << endl;
	  return;
	}
      }
      if(print_stuff)
         cout << "-";
//         cout << "no switch\n";
    }
    // or try an infection
    else
    { individual_descriptor source, target;
      if (choose_parent_first)
      { source = choose_source_individual(rng);
        // choose a target in same room (excluding the source)
        target = choose_in_same_room(source,false,rng);
      }
      else
      { target = choose_any_individual(rng);
        source = choose_in_same_room(target,true,rng);
      }
      if (source.exists && target.exists
          && source.infected && !target.infected)
      { if (print_stuff)
          cout << "infect " << (source.fixed ? 'f':'m')
               << "->" << (target.fixed ?'f':'m') << " in room "
               << source.room << endl;
        if (target.fixed)
          ++infected_fixed[target.room];
        else
        { ++infected_mobiles[target.room];
          ++n_infected_mobile;
        }
        ++total_infected;
      }
      else if (source.exists && target.exists
               && !source.infected && target.infected)
      { if (print_stuff)
          cout << "uninfect " << (source.fixed ? 'f':'m')
               << "->" << (target.fixed ?'f':'m') << " in room "
               << source.room << endl;
        if (target.fixed)
          --infected_fixed[target.room];
        else
        { --infected_mobiles[target.room];
          --n_infected_mobile;
        }
        --total_infected;
      }
      else
        if(print_stuff)
           cout << '.';
//           cout << "no infection\n";
    }
  }
};
