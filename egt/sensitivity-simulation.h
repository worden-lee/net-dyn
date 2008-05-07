/*  -*- C++ -*-
** 
** sensitivity-simulation.h
** 
** Made by 
** Login   <wonder@stripper>
** 
** Started on  Wed Mar 26 17:10:34 2008 
** Last update Wed Mar 26 17:10:34 2008 
*/

#ifndef   	SENSITIVITY_SIMULATION_H_
# define   	SENSITIVITY_SIMULATION_H_

#include <set>
#include <boost/lambda/lambda.hpp>
#include "indicators-general.h"

// this object does an extended simulation, keeping track of both what
// happens with homogeneous mutant fitnesses on all nodes, and with
// any one of the nodes' mutant fitness perturbed.
// see the wiki for details.
// It doesn't use operator(), it only descends from
// selection_indicator in order to share some of its machinery.
// It requires these parameters:
//  mutant_fitness
//  resident_fitness
//  perturbed_mutant_fitness
//  desired_accuracy

template<typename network_t, typename RNG_t, typename params_t>
class sensitivity_simulation
  : public selection_indicator<network_t,RNG_t,params_t>
{
  typedef selection_indicator<network_t, RNG_t, params_t> selind_t;
  using selind_t::params;

  typedef typename graph_traits<network_t>::vertex_descriptor
    vertex_descriptor;

  // these represent the "logical clauses" discussed in the wiki.
  // store one of these per node, instead of just a boolean for
  // resident or mutant.
  class node_state
  {
  public:
    bool is_T, is_F;
    //typedef set< vertex_descriptor > nodelist_t;
    typedef vector<bool> nodelist_t;
    nodelist_t nodelist;

    node_state(bool b=false)
      : is_T(b), is_F(!b)
    {}

    node_state&operator=(bool b)
    { is_F = !(is_T = b);
      // don't shrink the nodelist, just 0 it
      fill(nodelist.begin(),nodelist.end(),false);
    }
    node_state&operator=(const node_state&that)
    { is_T = that.is_T;
      is_F = that.is_F;
      if (!is_T && !is_F) // don't copy if not used
        nodelist = that.nodelist;
    }

    void nv_is(unsigned nv)
    { nodelist.resize(nv); }

    // evaluate the clause to T or F by specifying the r's
    // either baseline -- all unperturbed, or one r perturbed to r+
    bool value(vertex_descriptor perturbed=-1)
    { if (is_T) return true;
      else if (is_F) return false;
      else if (perturbed == -1) // == because it's probably unsigned
        return false;
      else
        //return (nodelist.count(perturbed) > 0);
        return nodelist[perturbed];
    }

    // replace *this with (*this or that)
    void or_with_other_clause(const node_state&that)
    { // (x v T) is T
      // (x v F) is x
      if (is_T || that.is_F)
        return;
      else if (that.is_T || is_F)
      { *this = that;
        return;
      }
      else // both are list clauses
      { // replace this->nodelist with the union of both lists
//         nodelist_t tmp = nodelist;
//         nodelist.clear();
//         set_union(tmp.begin(), tmp.end(), that.nodelist.begin(),
//                   that.nodelist.end(),
//         inserter(nodelist,nodelist.end()));
        nodelist_t::iterator nb = nodelist.begin(), ne = nodelist.end();
        nodelist_t::const_iterator tnb = that.nodelist.begin(),
          tne = that.nodelist.end();
        while(nb != ne)
        { *nb = *nb | *tnb;
          ++nb; ++tnb;
        }
      }
    }
    // replace *this with (*this or (that and r_x>=r+))
    void or_with_other_clause_and_r_term(const node_state&that,
                                         vertex_descriptor x)
    { if (is_T || that.is_F)
        return;
      // at this point, this is either F or a list,
      // that is either T or a list
      { // if that is (r_x >= r+ or ...),
        //  (that ^ r_x >= r+) is just r_x >= r+
        // if that is T, the same.
        // if that is a disjunction not including r_x,
        // (that ^ r_x>=r+) is F.
        // so the 'and' result is either r_x>=r+ or F,
        // so we either add the x clause or do nothing.
        //bool add_x = (that.is_T || that.nodelist.count(x) > 0);
        bool add_x = (that.is_T || //that.nodelist.count(x) > 0);
                      that.nodelist[x]);
        if (add_x && is_F) // list might be out of date
          nodelist.assign(nodelist.size(),false);
        if (add_x)
        {// nodelist.insert(x);
          nodelist[x] = true;
          is_F = false;
        }
      }
    }
    string asString()
    { if (is_T) return "T";
      if (is_F) return "F";
      ostringstream ss;
      for (int i = 0; i < nodelist.size(); ++i)
      { //cout << *ni << ' ';
        if (nodelist[i]) ss << i << ' ';
      }
      return ss.str();
    }
  };

public:
  //using selind_t::fixation_stats;

  // make one of these for each case (baseline, some node perturbed)
  class track_case_t : public selind_t::fixation_stats
  { // use the superclass to keep overall track of how many fixated/extinct 
  public:
    bool finished; // short-term: has this case come to fix/ext this
                   // time thru
    bool current_state; // even shorter term:
    bool mixed;         // used while checking for fix/ext
    track_case_t() : finished(false) {}
  };
  // this does colors for the display, reads the track_case's directly
  template<typename _network_t, typename _RNG_t, typename _params_t>
  friend class color_influence_stochastic;

  sensitivity_simulation(RNG_t&_rng, params_t*_params)
    : selind_t(_rng,_params)
  { //inheritParametersFrom(_params);
  }

  void do_simulation(const network_t &n);
};

#include "sensitivity-display.h"

template<typename network_t, typename RNG_t, typename params_t>
void sensitivity_simulation<network_t, RNG_t, params_t>
  ::do_simulation(const network_t &n)
{ network_component_status status = ::test_fixation_candidate(n);
  switch(status.flag)
  {
  case DISCONNECTED:
  case MULTIPLE_SOURCE_COMPONENTS:
    cout << "unconnected network\n";
    return;
  case ONE_SOURCE_COMPONENT:
    cout << "single source component, size = "
         << status.source_size << '\n';
    //       if (status.source_size == 1)
    //       { *answer_p = 1.0 / num_vertices(n);
    //         return true;
    //       }
    //       else
    //         break;
    return;
  default:// otherwise: go for it.
    break;
  }

  unsigned nv = num_vertices(n);
  // state is a vector of clause objects.
  vector<node_state> state(nv);

  for (int i = 0; i < nv; ++i)
    state[i].nv_is(nv);

  // there is one baseline case and nv perturbed cases.
  // we need to track all of them watching for fixation/extinction.
  // here the long-term fixation record, and a bool for whether it's
  // done in the current run.
  track_case_t track_baseline;
  // have to change this collection if vertex_descriptor isn't
  // unsigned int.
  vector<track_case_t> track_perturbations(nv);
 
  // get the parameters.
  double u = params.residentFitness();
  double r = params.mutantFitness();
  double rp = params.perturbed_mutant_fitness();
  double desired_accuracy = params.monte_carlo_selection_accuracy();

#if DISPLAY
  BoostDotGraphController<network_t,params_t> influence_display;
  influence_display.inheritParametersFrom(params);
  influence_display.add_vertex_property("color",
    make_color_influence_stochastic(track_baseline,track_perturbations,
                                    r,(rp-r),n,rng,params));
  influence_display.add_edge_property("style",
    make_bool_string_property("bold","",belongs_to_a_2_loop<network_t>,n));
  //influence_display.setdisplayEvery(1000);
  influence_display.setupdate_only_on_network_change(false);
#endif
  // random variable for reuse
  ur_t choose_xi(rng,uniform_real<>(0,rp));

  // do runs until accuracy is guaranteed small enough.
  unsigned max_trials = Bernoulli_estimate_n(desired_accuracy);
  cout << max_trials << " trials\n";

  while (track_baseline.n_trials < max_trials)
  { unsigned n_steps = 0;
    unsigned n_cases_to_finish = 1 + nv;

    // initialize to all F except one T.
    fill(state.begin(),state.end(),false);
    state[random_vertex(n,rng)] = true;

    track_baseline.finished = false;
    for (unsigned i = 0; i < track_perturbations.size(); ++i)
      track_perturbations[i].finished = false;

    // do events until all cases are finished.
    while (n_cases_to_finish > 0)
    { // TODO -- death-birth case
      // choose source
      vertex_descriptor source = random_vertex(n,rng);
      // choose target
      vertex_descriptor target;
      typename network_t::degree_size_type n_children
        = out_degree(source,n);
      typename graph_traits<network_t>::adjacency_iterator ci,cend;
      if (n_children <= 0)
        continue; // no out edges, start over
      else if (n_children == 1)
      { tie(ci,cend) = adjacent_vertices(source,n);
        target = *ci;
      }
      else
      { ui_t choose_a_child(this->rng,uniform_int<>(0,n_children-1));
        int child_selection = choose_a_child();
        for(tie(ci,cend) = adjacent_vertices(source,n);
            ci != cend; ++ci)
          if (child_selection-- == 0)
            break;
        if (child_selection >= 0)
        { cerr << "Error finding a child!" << endl; continue; } 
        target = *ci;
      }
      // have chosen source and target.  now do, or do not.
      double xi = choose_xi();
      // 3 cases depending on random variable xi.
      // case 1: uninfection/infection event.  copy s->t.
      if (xi <= u)
      { state[target] = state[source];
      }
      // case 2: infection only event.  t is infected if s is or t was.
      else if (xi <= r)
      { state[target].or_with_other_clause(state[source]);
      }
      // case 3: conditional infection, only in the case the source
      // is perturbed
      else
      {//  cout << source << " -> " << target << ": "
        //                << state[source].asString() << " -> "
        //                << state[target].asString() <<": ";
        state[target].or_with_other_clause_and_r_term(state[source],
                                                      source);
        //           cout << state[target].asString() << endl;
      }
      // now the event is done.
      // now check for fixations.
      typename graph_traits<network_t>::vertex_iterator vi,vend;
      tie(vi,vend) = vertices(n);
      if (!track_baseline.finished)
      { track_baseline.current_state = state[*vi].value();
        track_baseline.mixed = false;
      }
      for (unsigned i = 0; i < state.size(); ++i)
        if (!track_perturbations[i].finished)
        { track_perturbations[i].current_state = state[*vi].value(i);
          track_perturbations[i].mixed = false;
        }
      for (++vi; vi != vend; ++vi)
      { if (!track_baseline.finished && !track_baseline.mixed)
          if (state[*vi].value() != track_baseline.current_state)
            track_baseline.mixed = true;
        for (unsigned i = 0; i < state.size(); ++i)
          if (!track_perturbations[i].finished
              && !track_perturbations[i].mixed)
            if (state[*vi].value(i) != track_perturbations[i].current_state)
              track_perturbations[i].mixed = true;
      }
      if (!track_baseline.finished && !track_baseline.mixed)
      { if (track_baseline.current_state)
          ++track_baseline.n_fixations;
        else
          ++track_baseline.n_extinctions;
        ++track_baseline.n_trials;
        track_baseline.finished = true;
        --n_cases_to_finish;
      }
      for (unsigned i = 0; i < state.size(); ++i)
        if (!track_perturbations[i].finished
            && !track_perturbations[i].mixed)
        { if (track_perturbations[i].current_state)
            ++track_perturbations[i].n_fixations;
          else
            ++track_perturbations[i].n_extinctions;
          ++track_perturbations[i].n_trials;
          track_perturbations[i].finished = true;
          --n_cases_to_finish;
        }

      ++n_steps;
      if (0)//n_steps % 100 == 0 || n_cases_to_finish == 0)
      { cout << "after " << n_steps << " steps\nstate:\n";
        for (int i = 0; i < state.size(); ++i)
          cout << state[i].asString() << '\n';
        cout << "\nfixations:\nbaseline "
             << (track_baseline.finished ?
                 (track_baseline.current_state?'T':'F') : '-')
             << '\n';
        for (int i = 0; i < track_perturbations.size(); ++i)
        { cout //<< i << ' '
            << (track_perturbations[i].finished ?
                (track_perturbations[i].current_state?'T':'F') : '-')
            << ' ';
        }
        cout << '\n' << endl;
      }
    }
    cout << track_baseline.n_trials << ' '
         << track_baseline.fixation_probability()
         << "±" << track_baseline.fixation_accuracy() << endl;
#if DISPLAY
    if (track_baseline.n_trials % 1000 == 0 ||
        track_baseline.n_trials >= max_trials)
    { cout << "=== " << track_baseline.n_trials << " ===\n";
      influence_display.update(n_steps,n);
    }
#endif
  }
  //accuracy achieved.
  cout << "\ndone. fixation probabilities:\nbaseline: "
       << track_baseline.fixation_probability()
       << "±" << track_baseline.fixation_accuracy() << '\n';
  for(int i = 0; i < track_perturbations.size(); ++i)
    cout << track_perturbations[i].fixation_probability()
         << "±" << track_perturbations[i].fixation_accuracy() << ' ';
  cout << "\n\nchanged outcomes:\n";
  for(int i = 0; i < track_perturbations.size(); ++i)
    cout << track_perturbations[i].n_fixations - track_baseline.n_fixations
         << ' ';
  cout << "\n\nsensitivities:\n";
  vector<double> dpdri(track_perturbations.size());
  for(int i = 0; i < track_perturbations.size(); ++i)
  { dpdri[i] = (track_perturbations[i].fixation_probability()
                - track_baseline.fixation_probability()) / (rp - r);
    if (abs(dpdri[i]) < desired_accuracy/1000)
      dpdri[i] = 0;
    cout << dpdri[i] << ' ';
  }
  cout << endl;

  // these "displays" are for creating a data base to mine for correlations
  CSVDisplay network_csv("out/network.csv"), nodes_csv("out/nodes.csv");

  network_csv << "p" << "accuracy" << "density" << "mutuality"
              << "transitivity" << "median in density" << "median out density";
  network_csv.newRow();
  network_csv << track_baseline.fixation_probability()
              << track_baseline.fixation_accuracy()
              << density(n) << mutuality_ind(n) << transitivity_ind(n)
              << median_in_density(n) << median_out_density(n);
  network_csv.newRow();

  nodes_csv << "dp/dr[i]" << "n_changed" << "in degree" << "out degree"
            << "temperature";
  nodes_csv.newRow();
  stochastic_edge_weight_map<network_t> wm(n);
  for(int i = 0; i < track_perturbations.size(); ++i)
  { nodes_csv << dpdri[i]
              << (track_perturbations[i].n_fixations
                  - track_baseline.n_fixations)
              << in_degree(i,n) << out_degree(i,n)
              << temperature(i,n,wm);
    nodes_csv.newRow();
  }
}

#endif 	    /* !SENSITIVITY-SIMULATION_H_ */
