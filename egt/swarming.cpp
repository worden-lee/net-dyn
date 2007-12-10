// try the swarm-fixation code
#include "network-optimize.h"
#include "swarming.h"
#include "indicators-swarm.h"
#include "IndicatorsDisplay.h"

int main()
{
  rng_t rng;
  rng.seed((unsigned)time(0));
  // rand() used in random_shuffle()
  srand((unsigned)time(0));

  print_stuff = false;//true;
  timesToTestFixation = 10000;

  rooms_swarming_pattern swp;
  switch(1)
  { // set up the network and transition rates
  case 0: // complete dyad
    swp.r = 1.2;
    swp.n_rooms = 2;
    swp.n = 11;
    swp.n_mobile = 1;
    swp.fixed.set_size(swp.n_rooms);
    swp.fixed[0] = swp.fixed[1] = 5;
    swp.room_switch_rate = 0.01;
    swp.p_switch.set_size(swp.n_rooms,swp.n_rooms);
    swp.p_switch(0,1) = 0.01;
    swp.p_switch(1,0) = 1;
    break;
  case 1: // 5-vertex star graph
          // egt fixation probability is 1.09322 * moran prob
    swp.r = 1.2;
    swp.n_rooms = 5;
    swp.n = swp.n_mobile = 10;
    swp.fixed.set_size(swp.n_rooms);
    for (int i = 0; i < swp.n_rooms; ++i)
    { swp.n += swp.fixed[i] = 0;//5;
    }
    swp.room_switch_rate = 0.01;
    swp.p_switch.set_size(swp.n_rooms,swp.n_rooms);
    swp.p_switch(0,1) = 0.25;
    swp.p_switch(0,2) = 0.25;
    swp.p_switch(0,3) = 0.25;
    swp.p_switch(0,4) = 0.25;
    swp.p_switch(1,0) = 1;
    swp.p_switch(2,0) = 1;
    swp.p_switch(3,0) = 1;
    swp.p_switch(4,0) = 1;
    break;
  }    
  swarm_selection_indicator<rooms_swarming_pattern,rng_t> sind(rng);
  double moranprob =
      (1 - 1/mutantFitness) /
      (1 - 1/pow(mutantFitness,(double)swp.n));
  constant_indicator mind(moranprob);
  relative_indicator_t<
    swarm_selection_indicator<rooms_swarming_pattern,rng_t>,
    constant_indicator > rind(sind,mind);

  IndicatorsDisplayController<rooms_swarming_pattern>
    ind_display(animate_optimization?0:-1);
  ind_display.installIndicator(rind,"fixation/Moran prob");

//   cout << "fixation / Moran probability is " << rind(swp) << endl;
  for (swp.room_switch_rate = 0.001; swp.room_switch_rate <= 1000;
       swp.room_switch_rate *= 1.01)
//   for (swp.room_switch_rate = 1000; swp.room_switch_rate >= 0.001;
//        swp.room_switch_rate /= 1.01)
  { double p_switch
      = swp.room_switch_rate / (1 + swp.room_switch_rate);
    ind_display.update(p_switch,swp);
  }
}
