//  -*- c++ -*-
//
#ifndef BOOST_DOT_DISPLAY_H
#define BOOST_DOT_DISPLAY_H

#include "DotDisplay.h"
#include <boost/graph/graphviz.hpp>
#include "util.h"
#include "networks.h"
#include "indicators-network.h"
using namespace boost;

template<typename network_t>
class color_edges : public put_get_helper< string,
					   color_edges<network_t> >
{
protected:
  network_t &n;
public:
  typedef typename graph_traits<network_t>::edge_descriptor key_type;
  typedef string value_type;
  color_edges(network_t&_n) : n(_n) {}

  virtual value_type operator[](const key_type &e) const = 0;
};

template<typename network_t>
bool belongs_to_a_2_loop(typename graph_traits<network_t>::edge_descriptor e,
			 const network_t &n)
{ typename graph_traits<network_t>::edge_descriptor f;
  bool found;
  tie(f,found) = edge(target(e,n),source(e,n),n);
  return found;
}

template<typename network_t>
bool belongs_to_a_3_loop(typename graph_traits<network_t>::edge_descriptor e,
			 const network_t &n)
{ typename graph_traits<network_t>::adjacency_iterator ai,aend;
  for (tie(ai,aend) = adjacent_vertices(target(e,n),n); ai!=aend; ++ai)
  { typename graph_traits<network_t>::edge_descriptor f;
    bool connect;
    tie(f,connect) = edge(*ai,source(e,n),n);
    if (connect)
      return true;
  }
  return false;
}

template<typename network_t>
class color_2_loops : public color_edges<network_t>
{
public:
  typedef typename color_edges<network_t>::key_type key_type;
  typedef typename color_edges<network_t>::value_type value_type;
  color_2_loops(network_t&n) : color_edges<network_t>(n) {}
  using color_edges<network_t>::n;
  value_type operator[](const key_type &e) const
  { return belongs_to_a_2_loop(e,n)? "red":"black";
  }
};

template<typename network_t>
class color_3_loops : public color_edges<network_t>
{
public:
  typedef typename color_edges<network_t>::key_type key_type;
  typedef typename color_edges<network_t>::value_type value_type;
  color_3_loops(network_t&n) : color_edges<network_t>(n) {}
  using color_edges<network_t>::n;
  value_type operator[](const key_type &e) const
  { return belongs_to_a_3_loop(e,n)? "blue":"black";
  }
};

template<typename network_t>
class color_loops_23 : public color_edges<network_t>
{
public:
  typedef typename color_edges<network_t>::key_type key_type;
  typedef typename color_edges<network_t>::value_type value_type;
  color_loops_23(network_t&n) : color_edges<network_t>(n) {}
  using color_edges<network_t>::n;
  value_type operator[](const key_type &e) const
  { return belongs_to_a_2_loop(e,n)? "red":
            belongs_to_a_3_loop(e,n)? "blue":"black";
  }
};

template<typename network_t>
color_2_loops<network_t> make_color_2_loops(network_t&n)
{ return color_2_loops<network_t>(n);
}
  
template<typename network_t>
color_3_loops<network_t> make_color_3_loops(network_t&n)
{ return color_3_loops<network_t>(n);
}

template<typename network_t>
color_loops_23<network_t> make_color_loops_23(network_t&n)
{ return color_loops_23<network_t>(n);
}
  
template<typename network_t>
class color_vertices : public put_get_helper< string,
					      color_vertices<network_t> >
{
protected:
  network_t &n;
public:
  typedef typename graph_traits<network_t>::vertex_descriptor key_type;
  typedef string value_type;
  color_vertices(network_t&_n) : n(_n) {}

  virtual value_type operator[](const key_type &e) const = 0;
};

template<typename network_t>
class color_temperature : public color_vertices<network_t>
{
//   typedef constant_property_map<
//     typename graph_traits<network_t>::edge_descriptor, int > wm_t;
public:
  typedef typename color_vertices<network_t>::key_type key_type;
  typedef typename color_vertices<network_t>::value_type value_type;
  using color_vertices<network_t>::n;

  color_temperature(network_t&n) : color_vertices<network_t>(n)/*, wm(n)*/ {}
  value_type operator[](const key_type &e) const
  { return (*this)(e,n); }

  value_type operator()(const key_type &e, const network_t&c) const
  { typedef stochastic_edge_weight_map<network_t> wm_t;
    wm_t wm(c);
    typename graph_traits<network_t>::vertices_size_type nv = num_vertices(c);
    double temp = temperature(e, c, wm);
    // if all temps are 1, it's an isothermic graph, in which case
    // fixation process is equivalent to the Moran process

    // this scales 0 -> 1, 1 -> ~1/2, nv-1 -> 1
    double scaled = temp / (1 + ((double)nv-2)/(nv-1)*temp);
    // this makes the middle brighter
    scaled = sqrt(scaled);
    unsigned int byte =
      static_cast<unsigned int>(255*scaled);
    return stringf("#%02x%02x%02x",byte,0,255-byte); // rgb
  }
};
  
template<typename network_t>
color_temperature<network_t> make_color_temperature(network_t&n)
{ return color_temperature<network_t>(n);
}

template<typename network_t, typename indicator_t, typename params_t>
class color_influence : public color_vertices<network_t>
{
//   typedef constant_property_map<
//     typename graph_traits<network_t>::edge_descriptor, int > wm_t;
public:
  typedef typename color_vertices<network_t>::key_type key_type;
  typedef typename color_vertices<network_t>::value_type value_type;
  using color_vertices<network_t>::n;
  indicator_t &ind;
  params_t &params;

  // note _p has to be in the chain of parameters that _i uses
  color_influence(network_t&_n, indicator_t&_i, params_t&_p)
    : color_vertices<network_t>(_n),
      ind(_i), params(_p) {}
  value_type operator[](const key_type &e) const
  { return (*this)(e,n); }

  value_type operator()(const key_type &e, const network_t&c) const
  { double dr
      = params.perturbed_mutant_fitness() - params.mutantFitness();
    double p0 = ind(c);
    params.setperturb_r_at_vertex(e);
    double p1 = ind(c);
    params.unset("perturb_r_at_vertex");
    double inf = (p1 - p0) / dr;
//    cout << "<< " << inf << " >>";
    if (inf > 1)
      inf = 1;
    double scaled = inf;
    // this makes the middle brighter
    scaled = sqrt(scaled);
    unsigned int byte =
      static_cast<unsigned int>(255*scaled);
    // color from black to yellow
    return stringf("#%02x%02x%02x",byte,byte,0); // rgb
  }
};
  
template<typename network_t, typename indicator_t, typename params_t>
color_influence<network_t,indicator_t,params_t>
  make_color_influence(network_t&n,indicator_t&i,params_t&p)
{ return color_influence<network_t,indicator_t,params_t>(n,i,p);
}
  
template<typename property_t>
class string_property_map
{
public:
  string name;
  const property_t &property;
  typedef typename property_t::key_type key_type;
  typedef vector< pair<string,typename property_t::value_type> > value_type;
  string_property_map(string _nm, const property_t&_pr)
    : name(_nm), property(_pr) {}
  value_type operator[](const key_type &e) const
  { value_type v;
    v.push_back(make_pair(name,property[e]));
    return v;
  }
};

using namespace boost::lambda;

template<typename arg_t, typename context_t>
class multiple_string_property_map
{
protected:
  class prop_shim_base
  {public:
    //virtual string operator[](const arg_t&a) = 0;
    virtual string operator()(const arg_t&a,const context_t&c) = 0;
  };
  template<typename prop_t>
  class prop_shim : public prop_shim_base
  {public:
    prop_t prop;
    prop_shim(prop_t _p) : prop(_p) {}
//     string operator[](const arg_t&a)
//     { return prop[a]; }
    string operator()(const arg_t&a,const context_t&c)
    { return prop(a,c); }
  };
public:
  typedef vector< pair<string,string> > value_type;

  typedef vector< pair<string,prop_shim_base*> >props_t;
  props_t props;
  template<typename prop_t>
  void push_property(string name, prop_t _p)
  { props.push_back(make_pair(name, new prop_shim<prop_t>(_p))); }

  const context_t*context;
  void set_context(const context_t&c)
  { context = &c; }
  value_type operator[](const arg_t &a) const
  { value_type v;
    for (int i = 0; i < props.size(); ++i)
      v.push_back(make_pair(props[i].first,(*props[i].second)(a,*context)));
    return v;
  }
};

template<typename bool_indicator, typename context_t>
class bool_string_property
{
public:
  string if_yes, if_no;
  bool_indicator &indicator;
  context_t &context;
  //typedef key_t key_type;
  //typedef vector< pair<string,string> > value_type;
  bool_string_property(string _y, string _n,
                       bool_indicator &_i, context_t &_c)
    : if_yes(_y), if_no(_n), indicator(_i), context(_c)
  {}
  template<typename key_t>
  string operator[](const key_t&k) const
  { return indicator(k,context) ? if_yes : if_no;
  }
  template<typename key_t>
  string operator()(const key_t&k, const context_t&c) const
  { return indicator(k,c) ? if_yes : if_no; }
};

template<typename bool_indicator, typename subject_t>
bool_string_property<bool_indicator, subject_t>
make_bool_string_property(string _y, string _n, bool_indicator&_i,
                          subject_t _context)
{ return bool_string_property<bool_indicator,subject_t>
    (_y, _n, _i,_context);
}

template<typename bool_indicator, typename context_t>
class bool_string_property_map
{
public:
  string name;
  string if_yes, if_no;
  bool_indicator &indicator;
  context_t &context;
  //typedef key_t key_type;
  typedef vector< pair<string,string> > value_type;
  bool_string_property_map(string _name, string _y, string _n,
			   bool_indicator &_i, context_t &_c)
    : name(_name), if_yes(_y), if_no(_n), indicator(_i), context(_c)
  {}
  template<typename key_t>
  value_type operator[](const key_t&k) const
  { value_type v;
    v.push_back(make_pair(name,indicator(k,context) ? if_yes : if_no));
    return v;
  }
};

template<typename bool_indicator, typename subject_t>
bool_string_property_map<bool_indicator, subject_t>
make_bool_string_map(string _name, string _y, string _n, bool_indicator&_i,
		     subject_t _context)
{ return bool_string_property_map<bool_indicator,subject_t>
    (_name, _y, _n, _i,_context);
}

template<typename colorer>
class colormap : public string_property_map<colorer>
{
public:
  typedef string_property_map<colorer> super;
  typedef typename super::key_type   key_type;
  typedef typename super::value_type value_type;
  colormap(const colorer&_ce) : string_property_map<colorer>("color",_ce) {}
};

template<typename colorer>
colormap<colorer> make_colormap(const colorer&co)
{ return colormap<colorer>(co);
}

template<typename attrs>
attributes_writer<attrs> make_attributes_writer(const attrs&a)
{ return attributes_writer<attrs>(a);
}

class identity_map
{ public:
  //template<typename T>
  //T &operator[](T&t) { return t; }
};
template<typename T>
T get(identity_map &i, const T &t) { return t; }

// ================================================================
// how to label vertices: in the case of vector storage,
//  vertex_descriptor makes a fine label
// ================================================================
template<typename A,  typename C, typename D,
         typename E, typename F, typename G>
identity_map make_vertex_id(const boost::adjacency_list<A,boost::vecS,C,D,E,F,G>&n)
{ return identity_map(); }

// this is more serious
// and won't work in general. hack.
template<typename network_t>
class vertex_pointer_map
{ public:
  const network_t&n;
  vertex_pointer_map(const network_t&_n) : n(_n) {}
//   void *operator[](typename graph_traits<network_t>::vertex_descriptor &v)
//   { return &n[v]; }
};

template<typename network_t>
string get(vertex_pointer_map<network_t> &m,
           typename graph_traits<network_t>::vertex_descriptor const&v)
{ const void *pointer = &m.n[v];
  ostringstream os;
  os << pointer;
  // string now begins with '0x' which trips up graphviz
  return '"' + os.str() + '"';
}

template<typename A,  typename C, typename D,
         typename E, typename F, typename G>
vertex_pointer_map< boost::adjacency_list<A,boost::setS,C,D,E,F,G> >
make_vertex_id(const boost::adjacency_list<A,boost::setS,C,D,E,F,G>&n)
{ return
    vertex_pointer_map< boost::adjacency_list<A,boost::setS,C,D,E,F,G> >(n);
}

template<typename network_t, typename PC, bool Lattice_p=false>
class BoostDotGraphController : public DotGraphController<PC>
{
public:
  using ObjectWithParameters<PC>::params;

  BoostDotGraphController() : n(0) {}

  // this may be a departure from the boost idiom, but here it is:
  // each `property' inserted into this thing has to be able to do:
  //  string val = prop(v,n)
  multiple_string_property_map<typename network_t::vertex_descriptor, network_t>
    vertex_properties;
  //  string val = prop(e,n)
  multiple_string_property_map<typename network_t::edge_descriptor, network_t>
    edge_properties;

  template<typename prop_t>
  void add_vertex_property(string name, prop_t _p)
  { vertex_properties.push_property(name,_p); }
  template<typename prop_t>
  void add_edge_property(string name, prop_t _p)
  { edge_properties.push_property(name,_p); }

  void createDisplay();
  string epsfile(void)
  { return this->dotdir()+"/network.eps"; }
  string filenamebase(void)
  { return filebase; }
  void update(double t, const network_t &_n)
  { n = &_n;
    bool update;
    if (params.update_only_on_network_change())
    { string new_n = canonical(_n);
      update = (new_n != last_n);
      last_n = new_n;
    }
    else
      update = true;
    if (update)
    { filebase = stringf("%g",t);
      DisplayController<DotDisplay<PC>,PC>::update(t);
    }
  }
  bool writes_positions()
  { return Lattice_p; }

  void writeGraph(ostream &os)
  {// write_graphviz(os, *n);
    // assign colors to vertices
    //default_writer vw;
    // assign colors to edges
    default_writer gw;
    vertex_properties.set_context(*n);
    edge_properties.set_context(*n);
    write_graphviz(os,*n,
                   make_attributes_writer(vertex_properties),
                   make_attributes_writer(edge_properties),
                   gw,make_vertex_id(*n));
  }
protected:
  string filebase;
  const network_t *n;
  string last_n;
};

template<typename network_t, typename PC>
class BoostDotDisplay : public DotDisplay<PC>
{
public:
  BoostDotDisplay(DotGraphController<PC>*dc) : DotDisplay<PC>(dc) {}
  void writeGraph(ostream&ost)
  { static_cast<BoostDotGraphController<network_t,PC>*>(DotDisplay<PC>::xs)
      ->writeGraph(ost);
  }
};

template<typename network_t, typename PC, bool L>
void BoostDotGraphController<network_t,PC,L>::createDisplay()
{
  DotGraphController<PC>::display
    = new BoostDotDisplay<network_t,PC>(this);
}

#endif//BOOST_DOT_DISPLAY_H
