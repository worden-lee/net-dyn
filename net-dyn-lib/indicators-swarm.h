// -*- C++ -*-
//
// indicators to do with patterns of swarming individuals
// 
#include "indicators-general.h"
#include "indicators-network.h"
#include "BoostDotDisplay.h"

// implementation of monte carlo fixation sampling for network defined
// by neighbor relations on an n-dimensional lattice
// see indicators-general for the superclass
template<typename RNG_t, typename pt>
class static_lattice_swarm_selection_indicator
  : public selection_indicator<static_lattice_swarming_pattern, RNG_t, pt>
{ typedef static_lattice_swarming_pattern sw_t;
  typedef selection_indicator<sw_t,RNG_t,pt> sel_ind_t;
  using sel_ind_t::rng;
  using sel_ind_t::currently_processing;

  // this object does its work by handing the derived network to
  // network-style indicators; either analytic or monte carlo
  // (written before switching_indicator, which could be used here
  // now).
  network_selection_indicator<sw_t::graph_t,RNG_t,pt> nsind;
  typedef analytic_selection_indicator<sw_t::graph_t,pt> naind_t;
  naind_t naind;
  deterministic_indicator_caching<sw_t::graph_t,double,naind_t> dnaind;

public:
  static_lattice_swarm_selection_indicator(RNG_t&_rng)
    : sel_ind_t(_rng), nsind(_rng,&this->params),
      naind(&this->params), dnaind(naind)
  {}

  // this was used when caching statistics here, now we use the one
  // provided by the network indicator directly.
//   bool test_fixation_candidate(const sw_t&a, double*answer_p)
//   { return nsind.test_fixation_candidate(a.lattice.graph,answer_p); }

  // try_fixation isn't used, same as test_fixation_candidate
//   { //state_property_map mp(sw.graph);
//     vector<double> state(num_vertices(sw.lattice.graph));
//     return nsind.fixation_once(sw.lattice.graph,
//                                random_vertex(sw.lattice.graph,rng),state);
//   }

  // evaluate the lattice
  indicator_value operator()(const sw_t&sw)
  { this->inheritParametersFrom(sw.params);

    // what we do to get fixation probability on the lattice is, we
    // just hand the graph derived from the lattice configuration to
    // an indicator for graphs.  if it's small enough to do
    // analytically, do that; otherwise, do a series of trials using
    // the network indicator class for each single trial, and keep the
    // // statistics in this object, not that one.
    // that's dumb, better to keep them there, that way we get more
    // precision by amalgamating samples from equivalent graphs.
    // note they cache using the canonical string for the
    // graph, not the one for the lattice, because many lattice
    // configurations map to the same graph
    if (num_vertices(sw.lattice.graph) <= this->max_n_vertices_for_matrix())
      return dnaind(sw.lattice.graph);
    else
      //return sel_ind_t::operator()(sw);
      return nsind(sw.lattice.graph);
  }

  // use this when vertex_descriptors aren't an integer type
  // but don't because it seems to run really slow, just use vecS for
  // your graph type and the default implementation is fine
  class state_property_map
  { public:
    typedef static_lattice_swarming_pattern::graph_t graph_t;
    graph_t&graph;
    graph_traits<graph_t>::vertex_iterator _begin,_end;
    state_property_map(const graph_t&_g) : graph(const_cast<graph_t&>(_g))
    { tie(_begin,_end) = vertices(graph); }
    class state_iterator : public graph_traits<graph_t>::vertex_iterator
    { public:
      graph_t&graph;
      state_iterator(graph_traits<graph_t>::vertex_iterator _i,
                     graph_t&_g) : graph_traits<graph_t>::vertex_iterator(_i),
                                   graph(_g) {}
      double &operator*()
      { return graph[graph_traits<graph_t>::vertex_iterator::operator*()]->fitness;
      }
    };
    state_iterator begin()
    { return state_iterator(_begin,graph); }
    state_iterator end()
    { return state_iterator(_end,graph); }
    double &operator[](graph_traits<graph_t>::vertex_descriptor v)
    { return graph[v]->fitness; }
  };
};

// -- fixation probability on dynamic lattice --

// this selection indicator doesn't call out to the network one, it
// does it all itself, because movement on the lattice changes the
// network structure during the fixation process.
template <typename RNG_t, typename pt>
class dynamic_lattice_swarm_selection_indicator
  : public selection_indicator<dynamic_lattice_swarming_pattern, RNG_t, pt>
{ typedef selection_indicator<dynamic_lattice_swarming_pattern,RNG_t,pt>
    sel_ind_t;
  using sel_ind_t::rng;
  using sel_ind_t::currently_processing;
  typedef dynamic_lattice_swarming_pattern sw_t;

  static const unsigned Dimension = 2;
public:
  dynamic_lattice_swarm_selection_indicator(RNG_t&_rng)
    : sel_ind_t(_rng)
  {}

  // unlike in the general case (so far), this object comes with
  // parameters attached, which we should use in preference to the
  // global parameters.
  double operator()(const sw_t&sw)
  { this->inheritParametersFrom(sw.params);
    double answer = sel_ind_t::operator()(sw);
    return answer;
  }
  // sample the fixation process.  infect something, do lattice moves
  // and reproductions, until extinction, fixation, or timeout.
  fixation_result try_fixation(const sw_t&sw)
  { // the sw object doesn't provide a lattice, only specifications.
    // create a lattice and populate it before beginning.
    lattice_population<Dimension> lattice;
    lattice.inheritParametersFrom(this);
    typename lattice_population<Dimension>::lattice_size size;
    for(unsigned i = 0; i < Dimension; ++i)
    { const string *ds = sw.get(string("lattice_dim")+fstring("%u",i));
      size[i] = ds ? string_to_unsigned(*ds) : 0;
    }
    lattice.set_size(size);
    //lattice.set_size(sw.lattice_width(), sw.lattice_height());
    lattice.clear();
    lattice.place_n_random_individuals(sw.n_individuals(),rng);

    if (/*sw.n_rooms() < 1 || */sw.n_individuals() < 1)
      return fixation_result(ERROR,0);

    int total_infected = 1;

#ifdef DISPLAY
    // set up fixation animation machinery
    bool reallyDisplay = sw.animate_fixation();
    IndicatorsDisplayController<sw_t,ParamsClass> fixationDisplay;
    fixationDisplay.params.setdisplayToScreen(reallyDisplay);
    fixationDisplay.params.setdisplayEvery(100);
    if(reallyDisplay)
    { fixationDisplay.inheritParametersFrom(sw.params);
      fixationDisplay.params.setrecordEvery(0);
      fixationDisplay.installIndicator(variableWatcher(total_infected),
                                       "infected");
//      fixationDisplay.installIndicator(variableWatcher(sws.n_infected_mobile),
//                                       "mobile infected");
    }
#endif

    // uninfect the general population
    double resident_fitness = sw.residentFitness();
    double mutant_fitness = sw.mutantFitness();

    for (lattice_population<Dimension>::individual_indicator pi
           = lattice.population.begin(); pi != lattice.population.end(); ++pi)
      pi->fitness = resident_fitness;

    // infect a random individual
    graph_traits<lattice_population<Dimension>::graph_t>::vertex_descriptor rv
      = random_vertex(lattice.graph,rng);
    lattice.graph[rv]->fitness = mutant_fitness;

    // record parameters
    int    total_n   = sw.n_individuals(); 
    int    maxsteps  = sw.maxStepsToFixation();
    bool   cpf       = sw.choose_parent_first();
    double move_prob = sw.lattice_move_probability();
    // set up counters
    int tries = 0;
    double total_weight
      = total_n * resident_fitness
         + total_infected * (mutant_fitness - resident_fitness);
    // this object is used to generate random numbers in (0,1)
    ur_t choose_01(rng,uniform_real<>());

    // loop, doing moves and infections
    while (0 < total_infected && total_infected < total_n
           && tries < maxsteps)
    { int last_inf = total_infected;

      if (choose_01() < move_prob)
      { // choose to do a move
        lattice_mutation<Dimension> mut = random_move(lattice,rng);
        mut(lattice);
        //cout << "move " << mut.asString() << '\n' << lattice.canonical();
      }
      else // choose to do an infection
      { lattice_population<Dimension>::individual_indicator source, target;
        if (cpf)
        { // choose source first, weighted by fitness
          double choice = choose_01() * total_weight;
          for (lattice_population<Dimension>::individual_indicator pi
                 = lattice.population.begin();
               pi != lattice.population.end(); ++pi)
            if ((choice -= pi->fitness) <= 0)
            { source = pi;
              break;
            }
          // choose target second, equally among neighbors
          typename lattice_population<Dimension>::graph_t::degree_size_type
            n_out = out_degree(source->vertex,lattice.graph);
          typename graph_traits<lattice_population<Dimension>::graph_t>
            ::adjacency_iterator ci,cend;
          if (n_out <= 0)
            continue; // no neighbors, start over
          else if (n_out == 1)
          { tie(ci,cend) = adjacent_vertices(source->vertex,lattice.graph);
            target = lattice.graph[*ci];
          }
          else
          { ui_t choose_target(rng,uniform_int<>(0,n_out-1));
            int target_selection = choose_target();
            for(tie(ci,cend)
                  = adjacent_vertices(source->vertex, lattice.graph);
                ci != cend; ++ci)
              if (target_selection-- == 0)
                break;
            if (target_selection >= 0)
            { cerr << "Error finding a target!" << endl; continue; } 
            target = lattice.graph[*ci];
          }
        }
        else
        { // choose target first, equally among all individuals
          graph_traits<lattice_population<Dimension>::graph_t>::vertex_descriptor
            tv = random_vertex(lattice.graph,rng);
          target = lattice.graph[tv];
          // choose source second, weighted by fitness among neighbors
          typename inv_adjacency_iterator_generator<lattice_population<Dimension>::graph_t>::type
            pi, pend;
          tie(pi,pend) = inv_adjacent_vertices(target->vertex, lattice.graph);
          if (pi == pend) // if target has no in edges, move on
            continue;
          double total_weight = 0;
          for (; pi != pend; ++pi)
            total_weight += lattice.graph[*pi]->fitness;
          double choice = choose_01() * total_weight;
          for (tie(pi,pend)
                 = inv_adjacent_vertices(target->vertex,lattice.graph);
               pi != pend; ++pi)
            if ((choice -= lattice.graph[*pi]->fitness) <= 0)
              break;
          source = lattice.graph[*pi];
        }
        // now source and target are set
        // sanity check
        graph_traits<lattice_population<Dimension>::graph_t>::edge_descriptor
          ed;
        bool ed_found;
        tie(ed,ed_found) = edge(source->vertex,target->vertex,lattice.graph);
        if (!ed_found)
        { cerr << "found invalid edge " << source->location.asString()
               << "--" << target->location.asString() << endl;
          continue;
        }
        // if we get here, we're ready to do the event
        if (target->fitness != source->fitness)
        { if (target->fitness == resident_fitness)
            ++total_infected;
          else
            --total_infected;
          target->fitness = source->fitness;
          total_weight += source->fitness - target->fitness;
        }
      }// end of if () { move } else { infection }
      // now we've done the event
      ++tries;
#ifdef DISPLAY
      if (reallyDisplay)
        fixationDisplay.update(tries,sw);
#endif
    }// end of for loop
    // now we're done doing events, because of extinction, fixation,
    // or timeout
#ifdef DISPLAY
      if (reallyDisplay)
        fixationDisplay.updateDisplay(tries);
#endif
    return fixation_result(((total_infected == total_n) ?
                            FIXATION :
                            ((total_infected == 0) ?
                             EXTINCTION : TIMEOUT)),
                           tries, /*peak*/0 );
   }
};

// -- fixation probability for rooms experiment

// unfortunately, it turns out always to be equal to moran
// probability, so don't waste too much time on this code.
//
// it models a network of 'rooms' with multiple individuals per room,
// infecting each other and moving from room to room.  similar to the
// dynamic lattice model, but 'neighbors' are individuals in the same
// place, not neighboring places, and space is defined by a graph
// object instead of a lattice object.
template <typename sw_t, typename RNG_t, typename pt>
class swarm_selection_indicator
  : public selection_indicator<sw_t, RNG_t, pt>
{ typedef selection_indicator<sw_t,RNG_t,pt> sel_ind_t;
  using sel_ind_t::rng;
  using sel_ind_t::currently_processing;

public:
  swarm_selection_indicator(RNG_t&_rng)
    : sel_ind_t(_rng)
  {}

  double operator()(const sw_t&sw)
  { inheritParametersFrom(sw.params);
    double answer = sel_ind_t::operator()(sw);
    return answer;
  }
  fixation_result try_fixation(const sw_t&sw)
  {
    swarming_state sws(sw,rng);
    if (sw.n_rooms() < 1 || sw.n_individuals() < 1)
      return fixation_result(ERROR,0);
#ifdef DISPLAY
    bool reallyDisplay = sw.animate_fixation();
    IndicatorsDisplayController<sw_t,ParamsClass> fixationDisplay;
    fixationDisplay.params.setdisplayToScreen(reallyDisplay);
    fixationDisplay.params.setdisplayEvery(100);
    if(reallyDisplay)
    { fixationDisplay.inheritParametersFrom(sw.params);
      fixationDisplay.params.setrecordEvery(0);
      fixationDisplay.installIndicator(variableWatcher(sws.total_infected),
                                       "infected");
      fixationDisplay.installIndicator(variableWatcher(sws.n_infected_mobile),
                                       "mobile infected");
    }
#endif

    //sws.randomize_mobiles(rng);
    sws.infect_random_individual(rng);
    int tries = 0;
    while (0 < sws.total_infected && sws.total_infected < sw.n_individuals()
           /*&& tries < maxStepsToFixation*/)
    { sws.step(rng);
      ++tries;
#ifdef DISPLAY
      if (reallyDisplay)
        fixationDisplay.update(tries,sw);
#endif
    }
#ifdef DISPLAY
      if (reallyDisplay)
        fixationDisplay.updateDisplay(tries);
#endif
    return fixation_result(((sws.total_infected == sw.n_individuals()) ?
                            FIXATION :
                            ((sws.total_infected == 0) ?
                             EXTINCTION : TIMEOUT)),
                           tries, /*peak*/0 );
   }
};

template <typename sw_t, typename RNG_t, typename PC>
bool indicator_is_stochastic(swarm_selection_indicator<sw_t,RNG_t,PC>&i)
{ return true; }
