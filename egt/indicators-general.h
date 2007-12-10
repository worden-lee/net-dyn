/* -*- C++ -*-
 * 
 * general purpose indicator code
 *
 * an indicator is anything that you can do this to:
 *   double value = my_indicator(my_object)
 * in other words it needs to have an operator() that takes an appropriate
 *  thing and returns a double.
 *
 * they are used, for one thing, to provide a fitness function for
 *  optimization, and for another thing, to embody various measures
 *  that change over time.
 */
#ifndef __indicators_general_h__
#define __indicators_general_h__
#include <string>
#include <vector>
#include <map>
#include "network-optimize.h"
using namespace std;

// write out a vector of bits as a string of 0's and 1's.
// useful for debugging caching.
string vbstring(const vector<bool>&vb)
{ string vs(vb.size(),'0');
  string::iterator si = vs.begin();
  for(vector<bool>::const_iterator vi = vb.begin(); vi != vb.end(); ++vi,++si)
    if (*vi)
      *si = '1';
  return vs;
}

// provides a cached version of another indicator, so it
// won't be calculated more than once for the same argument.
// the argument type must have a 'canonical' function
// that maps it to the key type, which defaults to vector<bool>
template<typename subject_t, typename ret_t = double,
	 typename indicator_t = ret_t(*)(const subject_t&),
	 typename key_t = string >
class deterministic_indicator_caching
{
  indicator_t &indicator;
  typedef map<key_t,ret_t> cache_t;
  cache_t cache;
public:
  // memory leak here
  deterministic_indicator_caching() : indicator(*new indicator_t())
  {}
  deterministic_indicator_caching(indicator_t&_ind) : indicator(_ind)
  {}
  ret_t operator()(const subject_t&subject)
  { key_t key = canonical(subject);
    typename cache_t::iterator past = cache.find(key);
    if (past == cache.end())
    { ret_t &answer = cache[key] = indicator(subject);
    //cout << answer << '\n';
      return answer;
    }
    else
    { //cout << past->second << '\n';
      return past->second;
    }
  }
};

class constant_indicator
{ double c;
public:
  constant_indicator(double d) : c(d)
  {}
  template<typename network_t>
  double operator()(const network_t&n)
  { return c; } 
};

template<typename subject_t>
double zero_indicator(const subject_t&n)
{ return 0;
}

template<typename ind_t, typename ref_t>
class relative_indicator_t
{ ind_t &ind;
  ref_t &ref;
public:
  relative_indicator_t(ind_t&_ind,ref_t&_ref) : ind(_ind), ref(_ref)
  {}
  template<typename network_t>
  double operator()(const network_t&n)
  { return ind(n) / (double)ref(n);
  }
};

template<typename ind_t, typename ref_t>
relative_indicator_t<ind_t,ref_t>
relative_indicator(ind_t&_i,ref_t&_r)
{ return relative_indicator_t<ind_t,ref_t>(_i,_r);
}

// this uses one indicator as its value, but first compares
// it to the value of another one, so you can see if the two
// are giving the same result
template<typename i0_t, typename i1_t>
class comparing_indicator
{
  i0_t &i0;
  i1_t &i1;
public:
  comparing_indicator(i0_t &_i0, i1_t &_i1)
    : i0(_i0), i1(_i1) {}

  template<typename network_t>
  double operator()(const network_t&network)
  { double v0 = i0(network);
    double v1 = i1(network);
    if (print_stuff)
      cout << " [" << v0 << ' ' << v1 << "]\n";
//  if (abs(sv - av) > 0.15)
//    av = ai(network);
    return v0;
  }
};

// given an argument type that is capable of 'fixation' this
// class samples the probability of fixation, and caches its
// results so that multiple calls improve the estimate.
template<typename argument_t, typename RNG_t>
class selection_indicator
{
protected:
  class fixation_stats 
  { public:
    int fixations;
    int trials;
    int sum_fixtime;
    fixation_stats() : fixations(0), trials(0), sum_fixtime(0) {}
  };
  typedef map< string, fixation_stats > cache_t;
  cache_t cache;

  double cache_trials(const argument_t &a, int fixations, int trials,
                      int fixtime=0)
  { string ncan = canonical(a);
    fixation_stats&ncache = cache[ncan];
    ncache.fixations += fixations;
    ncache.trials += trials;
    ncache.sum_fixtime += fixtime;
//     cout << vbstring(ncan) << ' '
// 	 << ncache.fixations << ':' << ncache.trials << endl;
    return double(ncache.fixations)/ncache.trials;
  }
  RNG_t rng;

public:
  selection_indicator(RNG_t _rng): rng(_rng) {}

  virtual bool is_fixation_candidate(const argument_t&a)
  { return true; }

  virtual bool try_fixation(const argument_t&) = 0;

  // do a number of fixation samples, return sample probability.
  double operator()(const argument_t &a)
  { // first, run analytic test: if the graph isn't connected, or 
    // whatever, fixation is impossible.
    if (!is_fixation_candidate(a))
      return 0;

    int n_fixations = 0;
    //cout << "testing for fixation: ";
    for (int i = 0; i < timesToTestFixation; ++i)
    { bool fixed = try_fixation(a);
//       if (print_stuff)
//         cout << (fixed ? "":"no ") << "fixation" << endl;
      if (fixed)
	++n_fixations;
    }
    //cout << '\n';
    return cache_trials(a,n_fixations,timesToTestFixation);
  }
  double report_mean_fixation_time(const argument_t&a)
  { fixation_stats&ncache = cache[canonical(a)];
    return double(ncache.sum_fixtime)/ncache.trials;
  }
};

// the selection_indicator classes are capable of collecting
// mean fixation time as they go, but they don't report it to
// callers.  this class attaches to that one giving you an indicator
// that accesses the recorded mean fixation time.
template<typename sel_ind_class>
class mean_fixation_time_indicator
{
protected:
  sel_ind_class&sel_ind;
public:
  mean_fixation_time_indicator(sel_ind_class&_si) : sel_ind(_si) {}
  template<typename arg_t>
  double operator()(const arg_t&a)
  { return sel_ind.report_mean_fixation_time(a);
  }
};

#endif //__indicators_general_h__
