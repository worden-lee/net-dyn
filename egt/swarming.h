//  -*- C++ -*-
//
// types for swarming
#include "network-optimize.h"
#include <vnl_matrix.h>
#include <vnl_vector.h>
#include <iostream>
#include <sstream>
using namespace std;
using namespace boost;

class rooms_swarming_pattern
{
public:
  // mutant fitness
  double r;
  // number of rooms
  unsigned int n_rooms;
  // number of total individuals (the below should add up to this)
  unsigned int n;
  // number of mobile individuals
  unsigned int n_mobile;
  // number fixed in each room
  vnl_vector<unsigned int> fixed;
  // this controls the rate of room switching relative to
  // the rate of infection.
  double room_switch_rate;
  // prob a mobile individual in room i will switch to room j
  vnl_matrix<double> p_switch;
};

string canonical(const rooms_swarming_pattern &sw)
{ ostringstream os;
  os << sw.r << ' ' << sw.n_rooms << ' ' << sw.n << ' '
     << sw.n_mobile << ' ' << sw.fixed << ' '
     << sw.room_switch_rate << ' ' << sw.p_switch;
  return os.str();
}

class rooms_swarming_state
{
public:
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
  // pattern to use
  const rooms_swarming_pattern &pattern;
  
  rooms_swarming_state(const rooms_swarming_pattern&_pattern)
    : mobiles(_pattern.n_rooms,0), infected_mobiles(_pattern.n_rooms,0),
      infected_fixed(_pattern.n_rooms,0), total_infected(0),
      n_infected_mobile(0), pattern(_pattern)
  { // all vectors are created full of zeros, mobiles need to be placed.
    // just put them all in room 0, I guess.
    mobiles[0] = pattern.n_mobile;
  }

  template<typename RNG_t>
  void randomize_mobiles(RNG_t&rng)
  { // what distribution to use?
    //for(unsigned int i = 0; i < pattern.n_mobile; ++i)
    cerr << "error: randomize_mobiles not implemented!" << endl;
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
  { ui_t chooser(rng,uniform_int<>(0,pattern.n-1));
    int choice = chooser();
    int seen = 0;
    if (choice < pattern.n_mobile)
      for (int i = 0; i < pattern.n_rooms; ++i)
      { if (choice < seen + infected_mobiles[i])
          return individual_descriptor(i,false,true);
        else if (choice < seen + mobiles[i])
	  return individual_descriptor(i,false,false);
        else
          seen += mobiles[i];
      }
    else
    { seen += pattern.n_mobile;
      for (int i = 0; i < pattern.n_rooms; ++i)
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
  { double factor = (mobilityIncreasesWithFitness ?
                     pattern.r : residentFitness);
    ur_t chooser(rng,uniform_real<>(0, pattern.n_mobile
                                        + (factor-1)*n_infected_mobile));
    double choice = chooser();
    double seen = 0;
    for (int i = 0; i < pattern.n_rooms; ++i)
    { if (choice < seen + factor*infected_mobiles[i])
        return individual_descriptor(i,false,true);
      else if (choice < seen + mobiles[i] + (factor-1)*infected_mobiles[i])
        return individual_descriptor(i,false,false);
      seen += mobiles[i] + (factor-1)*infected_mobiles[i];
    }
    cerr << "error in choose_mobile_individual" << endl;
  }

  template<typename RNG_t>
  individual_descriptor choose_source_individual(RNG_t&rng)
  { ur_t chooser(rng,uniform_real<>(0,pattern.n+(pattern.r-1)*total_infected));
    double choice = chooser();
    double seen = 0;
    for (int i = 0; i < pattern.n_rooms; ++i)
    { if (choice < seen + pattern.r*infected_mobiles[i])
        return individual_descriptor(i,false,true);
      else if (choice < (seen + mobiles[i] - infected_mobiles[i]
                         + pattern.r*infected_mobiles[i]))
        return individual_descriptor(i,false,false);
      seen += (mobiles[i] - infected_mobiles[i]
               + pattern.r*infected_mobiles[i]);
    }
    for (int i = 0; i < pattern.n_rooms; ++i)
    { if (choice < seen + pattern.r*infected_fixed[i])
        return individual_descriptor(i,true,true);          
      else if (choice < (seen + pattern.fixed[i] - infected_fixed[i]
                         + pattern.r*infected_fixed[i]))
        return individual_descriptor(i,true,false);
      seen += (pattern.fixed[i] - infected_fixed[i]
               + pattern.r*infected_fixed[i]);
    }
    cerr << "error in choose_source_individual" << endl;
  }

  template<typename RNG_t>
  individual_descriptor choose_in_same_room(individual_descriptor as,
                                            bool weighted,RNG_t&rng)
  { double r = weighted? pattern.r : 1;
    unsigned int i = as.room;
    double total_r = (pattern.fixed[i]+mobiles[i]
                      +(r-1)*(infected_fixed[i]
                              +infected_mobiles[i]));
    if (as.infected) total_r -= r;
    else             total_r -= 1;
    if (total_r == 0)
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
    double mobfactor = (mobilityIncreasesWithFitness ?
                        pattern.r : residentFitness);
    double total_switch_rate =
      pattern.room_switch_rate * (pattern.n_mobile
                                  + (mobfactor-1)*n_infected_mobile);
    double total_inf_rate =
      (choose_parent_first ?
       pattern.n + (pattern.r-1)*total_infected : pattern.n);
    double p_switch
      = total_switch_rate / (total_switch_rate + total_inf_rate);
    // either try a room switch
    if (choose_01() < p_switch)
    { individual_descriptor ind = choose_mobile_individual(rng);
      double move_choice = choose_01();
      double seen = 0.0;
      for (int j = 0; j < pattern.n_rooms; ++j)
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
