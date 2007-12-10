#include "DotDisplay.h"
#include <boost/graph/graphviz.hpp>
#include "util.h"
#include "network-optimize.h"
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
  typedef stochastic_edge_weight_map<network_t> wm_t;
  wm_t wm;
public:
  typedef typename color_vertices<network_t>::key_type key_type;
  typedef typename color_vertices<network_t>::value_type value_type;
  using color_vertices<network_t>::n;

  color_temperature(network_t&n) : color_vertices<network_t>(n), wm(n) {}
  value_type operator[](const key_type &e) const
  { typename graph_traits<network_t>::vertices_size_type nv = num_vertices(n);
    double temp = temperature(e, n, wm);
    // this scales 0 -> 1, 1 -> ~1/2, nv-1 -> 1
    // if all temps are 1, it's an isothermic graph, which is
    // equivalent to the Moran process
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

template<class network_t>
class BoostDotGraphController : public DotGraphController
{
public:
  BoostDotGraphController(double p = 1, double fp = 1) :
    DotGraphController(p,fp), n(0)
  {} 
  void createDisplay();
  string outdir()
  {
    return "out";
  }
  string epsfile(void)
  {
    return "network.eps";
  }
  string graphfile(void)
  {
    return filebase + ".dot";
  }
  void update(double t, const network_t &_n)
  {
    n = &_n;
    filebase = stringf("%g",t);
    DisplayController<DotDisplay>::update(t);
  }
  void writeGraph(ostream &os)
  {// write_graphviz(os, *n);
    // assign colors to vertices
    //default_writer vw;
    // assign colors to edges
    write_graphviz(os,*n,
      make_attributes_writer(make_colormap(make_color_temperature(*n))),
      //make_attributes_writer(make_colormap(make_color_2_loops(*n))));
      make_attributes_writer(make_bool_string_map("style","bold","",
				belongs_to_a_2_loop<network_t>,*n)));
  }
protected:
  string filebase;
  const network_t *n;
};

template<class network_t>
class BoostDotDisplay : public DotDisplay
{
public:
  BoostDotDisplay(DotGraphController*dc) : DotDisplay(dc) {}
  void writeGraph(ostream&ost)
  {
    dynamic_cast<BoostDotGraphController<network_t>*>(xs)->writeGraph(ost);
  }
};

template<class network_t>
void BoostDotGraphController<network_t>::createDisplay()
{
  display = new BoostDotDisplay<network_t>(this);
}
