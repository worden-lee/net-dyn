#include "DotDisplay.h"
#include <boost/graph/graphviz.hpp>
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
  
template<typename colorer>
class colormap
{
public:
  const colorer &ce;
  typedef typename colorer::key_type key_type;
  typedef vector< pair<string,typename colorer::value_type> > value_type;
  colormap(const colorer&_ce) :ce(_ce) {}
  value_type operator[](const key_type &e) const
  {
    value_type v;
    v.push_back(make_pair("color",ce[e]));
    return v;
  }
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
    return ".";
  }
  string filebasename(void)
  {
    return filebase;
  }
  void update(double t, network_t &_n, string _filebase)
  {
    n = &_n;
    filebase = _filebase;
    //_update();
    DisplayController<DotDisplay>::update(t);
  }
  void writeGraph(ostream &os)
  {// write_graphviz(os, *n);
    // write with colors?
    default_writer vw;
    write_graphviz(os,*n, vw,
      make_attributes_writer(make_colormap(make_color_2_loops(*n))));
  }
  
protected:
  string filebase;
  network_t *n;
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
