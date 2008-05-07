//  -*- c++ -*-
//
#ifndef SENSITIVITY_DISPLAY_H
#define SENSITIVITY_DISPLAY_H

#include "BoostDotDisplay.h"
using namespace boost;

template<typename network_t, typename RNG_t, typename params_t>
class color_influence_stochastic : public color_vertices<const network_t>
{
public:
  typedef typename color_vertices<const network_t>::key_type key_type;
  typedef typename color_vertices<const network_t>::value_type value_type;
  using color_vertices<const network_t>::n;
  typedef
    typename sensitivity_simulation<network_t,RNG_t,params_t>::track_case_t
    track_case_t;
  track_case_t &track_baseline;
  vector<track_case_t>&track_perturbations;
  double r, dr;

  color_influence_stochastic(track_case_t&_tb, vector<track_case_t>&_tp,
                             double _r, double _dr, const network_t&_n)
    : color_vertices<const network_t>(_n),
      track_baseline(_tb), track_perturbations(_tp),
      r(_r), dr(_dr) {}
  value_type operator[](const key_type &e) const
  { return (*this)(e,n); }

  value_type operator()(const key_type &e, const network_t&c) const
  { double p0 = track_baseline.fixation_probability();
    double p1 = track_perturbations[e].fixation_probability();
    double inf = (p1 - p0) / dr;
    if (inf < 0)
      inf = 0;
//    cout << "<< " << inf << " >>";
    // an estimate of the characteristic size of these numbers
    double scaled = inf * (r*r*num_vertices(c)/2);
    if (scaled > 1)
      scaled = 1;
    else // this makes the middle brighter
      scaled = sqrt(scaled);
    unsigned int byte =
      static_cast<unsigned int>(255*scaled);
    // color from black to red
    return stringf("#%02x%02x%02x",byte,0,0); // rgb
  }
};
  
template<typename network_t, typename RNG_t, typename params_t>
color_influence_stochastic<network_t,RNG_t,params_t>
make_color_influence_stochastic(
  typename sensitivity_simulation<network_t,RNG_t,params_t>::track_case_t&tb,
  vector<typename sensitivity_simulation<network_t,RNG_t,params_t>::track_case_t>&tp,
  double r, double dr,const network_t&n,RNG_t&rng,params_t&p)
{ return color_influence_stochastic<network_t,RNG_t,params_t>(tb,tp,r,dr,n);
}

#endif//SENSITIVITY_DISPLAY_H
