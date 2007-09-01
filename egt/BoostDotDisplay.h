#include "DotDisplay.h"
#include <boost/graph/graphviz.hpp>

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
  void update(network_t &_n, string _filebase)
  {
    n = &_n;
    filebase = _filebase;
    _update();
  }
  void writeGraph(ostream &os)
  { write_graphviz(os, *n);
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
