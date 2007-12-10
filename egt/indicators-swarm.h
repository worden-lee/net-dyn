// -*- C++ -*-
//
// indicators to do with patterns of swarming individuals
// 
#include "indicators-general.h"

template <typename sw_t, typename RNG_t>
class swarm_selection_indicator : public selection_indicator<sw_t, RNG_t>
{
  using selection_indicator<sw_t,RNG_t>::rng;
  
public:
  swarm_selection_indicator(RNG_t&_rng)
    : selection_indicator<sw_t,RNG_t>(_rng) {}

  bool try_fixation(const sw_t&sw)
  {
    rooms_swarming_state sws(sw);
    //sws.randomize_mobiles(rng);
    sws.infect_random_individual(rng);
    int tries = 0;
    while (0 < sws.total_infected && sws.total_infected < sw.n
           /*&& tries < maxStepsToFixation*/)
    { sws.step(rng);
      ++tries;
    }
    return (sws.total_infected == sw.n);
  }
};
