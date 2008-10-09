/*  -*- C++ -*-
** 
** search.h
** 
** Made by lw
** Login   <wonder@Jamison>
** 
** Started on  Sat Jan 12 17:26:14 2008 
** Last update Sat Jan 12 17:26:14 2008 
*/

#ifndef   	SEARCH_H_
# define   	SEARCH_H_

// search algoriths, pretty generically defined.
// in all these cases:
// 
// args:
//   n is the object to be improved (in place)
//   indicator is something that can be used like this:
//     double fitness = indicator(n);
//   mutagen is something that generates mutation objects, which
//     are things (`mut') that can do
//       mut(n);
//       ...
//       cout << mut.asString();
//       undo(mut)(n);
//     how mutagen is used depends on the algorithm
//   parameters and df can be null
// returns: void, but modifies n.

// hill_climb(): improve the object until the indicator can't get any
// bigger (within one mutation).
//
// here, mutagen is used like this:
//     for(typename mutation_generator_t::mutation_iterator
//         mi = mutagen.begin(n); mi != mutagen.end(); ++mi)
//     { typename mutation_generator_t::mutation mut = *mi;
//       (...)
//     }
//
// requires the following to be defined:
//   write_snapshot(n,(double)t,*parameters)
//
// uses parameters:
//   trialsForLocalOptimum
//   maxNeutralSteps
//   print_stuff
//
template<typename object_t, typename indicator_t,
         typename mutation_generator_t, typename params_t,
         typename displayfarm_t>
void hill_climb(object_t &n, indicator_t &indicator,
                mutation_generator_t &mutagen,
                params_t*parameters,
                displayfarm_t*df=0)
{ int nSteps = 0; 
  int consecutiveTimesFoundLocalOptimum = 0;
  int trialsForLocalOptimum
    = (indicator_is_stochastic(indicator)
       ? parameters->trialsForLocalOptimum() : 1);
  double currentIndicator;

  currentIndicator = indicator(n);
  if (parameters->print_stuff())
    cout << "indicator: " << currentIndicator << '\n';

  write_snapshot(n,nSteps,*parameters);
  
  while (consecutiveTimesFoundLocalOptimum < trialsForLocalOptimum)
  {
    bool foundImprovement = true;
    bool foundNeutralNeighbor = false;
    bool nadd;
    typename mutation_generator_t::mutation neutral_neighbor;
    int neutralSteps = 0;
    while(foundImprovement ||
	  (foundNeutralNeighbor
           && neutralSteps < parameters->maxNeutralSteps()))
    {
      // stochastic indicator? call it again each time through, for
      // precision
      if (indicator_is_stochastic(indicator))
      { currentIndicator = indicator(n);
        if (parameters->print_stuff())
          cout << "indicator: " << currentIndicator << '\n';
      }

//     cout << "current network\n";
      print_object(n,cout);
#ifdef DISPLAY
      if(df)
        df->update(nSteps,n);
#endif//DISPLAY
      
      foundImprovement = foundNeutralNeighbor = false;
      // assess each candidate for change
      for (typename mutation_generator_t::mutation_iterator
             mi = mutagen.begin(n);
	   !foundImprovement && mi != mutagen.end(); ++mi)
      { typename mutation_generator_t::mutation mut = *mi;
        // do the mutation
        mut(n);
        // see if the indicator is better
        double newIndicator = indicator(n);
	// for stochastic indicators: if it's big, call it again
	//  and make sure it's really big
	// (with the caching I've got in there, the more you call it
	//  the more accurate it is)
	if (indicator_is_stochastic(indicator) &&
	    newIndicator >= currentIndicator)
	{ newIndicator = indicator(n);
	} 
	if (newIndicator >= currentIndicator &&
	    // need <= in case of infinity
	    newIndicator <= currentIndicator + 1e-7)
	{ if (!foundNeutralNeighbor)
  	  { foundNeutralNeighbor = true;
	    neutral_neighbor = mut;
 	  }
	}
	else if (newIndicator > currentIndicator)
	{ currentIndicator = newIndicator;
          foundImprovement = true;
	  if (parameters->print_stuff())
	  { cout << "found improvement " << mut.asString() << ": "
		 << currentIndicator << '\n';
//             if (nv <= 40)
//            print_object(n,cout);
            write_snapshot(n,nSteps,*parameters);
	  }
	  break;
	}

	// if we arrive here, we didn't find an improvement, put it back
        mutation_generator_t::mutation::undo(mut)(n);
      }
      if (foundImprovement)
      { ++nSteps;
	neutralSteps = 0;
	consecutiveTimesFoundLocalOptimum = 0;
      }
      else if (foundNeutralNeighbor
               && neutralSteps++ < parameters->maxNeutralSteps())
      { neutral_neighbor(n);
        if (parameters->print_stuff())
  	{ cout << "taking neutral step " << neutral_neighbor.asString()
               << ": " << currentIndicator << '\n';
//           if (nv <= 40)
//          print_object(n,cout);
          write_snapshot(n,nSteps,*parameters);
	}
        ++nSteps;
	consecutiveTimesFoundLocalOptimum = 0;
      }
    }
    if (neutralSteps >= parameters->maxNeutralSteps())
    { if (parameters->print_stuff())
	cout << "reached max neutral steps\n";
      break;
    }
    else
    { if (parameters->print_stuff())
        cout << "at local optimum\n";
      ++consecutiveTimesFoundLocalOptimum;
    }
#ifdef DISPLAY
    if (df)
      df->update(nSteps,n);
#endif//DISPLAY
  }
  if (parameters->print_stuff())
    cout << "indicator: " << currentIndicator << '\n';
}

#ifndef SKIP_SIM_ANNEAL

// simulated_annealing(): take random walk biased toward steps that
// improve the indicator; temperature gradually decreases (strength of
// bias increases) until a final state is reached.
//
// here, mutagen is used like this:
//   typename mutation_generator_t::mutation mut =
//     mutagen.random_mutation(n,rng);
//
// requires the following to be defined:
//   write_snapshot(n,(double)t,*parameters)
//
// uses parameters:
//   annealing_initial_temperature     starting temp
//   annealing_cooling_rate            amt to decrease temp per time unit
//   rate_factor                       # of mutation events per time unit
//   print_stuff                       whether to print messages
//
template<typename object_t, typename indicator_t,
         typename mutation_generator_t, typename params_t,
         typename rng_t, typename displayfarm_t>
void simulated_annealing(object_t &n, indicator_t &indicator,
                         mutation_generator_t &mutagen,
                         params_t*parameters, rng_t&rng,
                         displayfarm_t*df=0)
{ double temperature;
  double initial_temperature = parameters->annealing_initial_temperature();
  double coolingrate = parameters->annealing_cooling_rate();
  double ratefactor  = parameters->annealing_rate_factor();
  coolingrate /= ratefactor;
  double maxSteps = initial_temperature / coolingrate;
  int nSteps = 0;
  double multiplier = exp( - log(maxSteps) / maxSteps );
  // this figures it will take maxSteps steps to get down to coolingrate
  bool print = parameters->print_stuff();
#define progress() (nSteps/maxSteps)

  indicator_value currentIndicator;

  object_t best;
  indicator_value bestIndicator = -HUGE;
  double best_t;

  currentIndicator = indicator(n);

  write_snapshot(n,progress(),*parameters);

  ur_t choose_01(rng,uniform_real<>());

  temperature = initial_temperature;
 
  if (print)
    cout << temperature << ": " << currentIndicator << '\n';

  // consider making timesToTestFixation() rise as temp decreases?
  for (nSteps = 0; 
       temperature > 0; temperature -= coolingrate/*, ++nSteps*/)
    //temperature >= coolingrate; temperature *= multiplier)
  {
// want accuracy of both new and old indicators to be small enough
// so have to sample enough times
// if f' = f - T the prob drops from 1 to 1/e
// so try requiring accuracy = T
// but; if we use (f' - f)/f we just need accuracy = Tf
    while(currentIndicator.accuracy() > temperature * currentIndicator)
    {// cout << '/';
      currentIndicator = indicator(n);
    }
//     if (print)
//       cout << temperature << ": " << currentIndicator << '\n';

    //long t_bef = time(0);

    typename mutation_generator_t::mutation mut = mutagen.random_mutation(n);
    // do the mutation
    mut(n);

    // see if the new indicator is better
    indicator_value newIndicator = indicator(n);

    while(newIndicator.accuracy() > temperature * currentIndicator)
    {// cout << '\\';
      newIndicator = indicator(n);
    }
    //long t_aft = time(0);
    //cout << (t_aft - t_bef) << " seconds to get an estimate" << endl;

    double prob;
    if (newIndicator >= currentIndicator)
      prob = 1;
    else
      prob = exp(((newIndicator - currentIndicator)
                  /(currentIndicator.value()?:1.0))/temperature);
//     if (currentIndicator == 0)
//       prob = 1;
//     else
//       prob = 1 - (1 - newIndicator/currentIndicator)/temperature;
    if (prob == 1 || choose_01() < prob)
    { currentIndicator = newIndicator;
      ++nSteps;
      if (print)
      {//  if (prob < 1)
//           cout << "drop - prob = " << prob << '\n';
        cout << temperature << ": " << mut.asString() << ' '
             << currentIndicator << '\n';
        write_snapshot(n,progress(),*parameters);
      }
      if (currentIndicator > bestIndicator)
      { best = n;
        bestIndicator = currentIndicator;
        best_t = progress();
      }
    }
    else
    { //if (print)
      //  cout << '(' << newIndicator << ")\n";
      mutation_generator_t::mutation::undo(mut)(n);
      ++nSteps;
    }
#ifdef DISPLAY
    if (df)
      df->update(progress(),n);
#endif//DISPLAY
  }
  // done
  if (print)
    cout << "best at " << best_t << ": " << bestIndicator << '\n'
         << canonical(best) << '\n';
#ifdef DISPLAY
  if (df)
    df->update(progress(),n);
#endif//DISPLAY
}

#endif // SKIP_SIM_ANNEAL

#endif 	    /* !SEARCH_H_ */
