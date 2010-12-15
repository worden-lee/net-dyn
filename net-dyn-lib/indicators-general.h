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
#include "IndicatorsDisplay.h"
#include "HistogramDisplay.h"
//#include "network-optimize.h"
using namespace std;
#include <boost/lambda/lambda.hpp>
using namespace boost;
using namespace boost::lambda;
#include <iostream>
#include <limits>

// it is sometimes useful for an indicator to return this type.
// it allows you to query how accurate the value is (since sometimes
// it's a monte carlo sample).
class indicator_value
{ double _value;
  double _accuracy;
public:
  indicator_value(double v=0, double a=0) : _value(v), _accuracy(a) {}

  double value()
  { return _value; }
  double accuracy()
  { return _accuracy; }
  operator double()
  { return _value; }
};

// indicator_values print out like doubles, but with ± part after the value
ostream&operator<<(ostream&os,indicator_value&val)
{ os << val.value();
  if (val.accuracy() != 0)
    os << "±" << val.accuracy();
  return os;
}

// useful to know whether you'll gain information by calling the
// indicator more than once.  call indicator_is_stochastic(indicator)
// to find out.  though the indicator_value machinery may make this
// redundant, we'll see.
// was having trouble with default template, now trying saying default
// is false for all pure functions
template<typename arg_t>
bool indicator_is_stochastic(double(*f)(const arg_t&))
{ return false; }

// canonical (see below) doesn't have a default implementation.  it
// was making things hard.

// template<typename xt>
// string canonical(const xt&x)
// { return "<this type needs a canonical function!>\n"; }

// boost functions sometimes use a map api rather than an operator() api
// so it's useful to convert back or forth
template<typename ind_t, typename ret_t=double>
class indicator_to_map
{ ind_t&ind;
public:
  indicator_to_map(ind_t&_i) : ind(_i) {}
  template<typename arg_t>
  ret_t operator[](const arg_t&arg)
  { return ind(arg); }
};

template<typename map_t, typename ret_t=double>
class map_to_indicator
{ map_t&map;
public:
  map_to_indicator(map_t&_m) : map(_m) {}
  template<typename arg_t>
  ret_t operator()(const arg_t&arg)
  { return map[arg]; }
};


// provides a cached version of another indicator, so it won't be
// calculated more than once for the same argument.  the argument type
// must have a 'canonical' function that maps it uniquely to the key
// type, which is string in all cases so far
template<typename subject_t, typename ret_t = double,
	 typename indicator_t = ret_t(*)(const subject_t&),
	 typename key_t = string >
class deterministic_indicator_caching
{
  indicator_t &indicator;
  typedef map<key_t,ret_t> cache_t;
  cache_t cache;
public:
  // memory leak here: this indicator is never freed
  deterministic_indicator_caching() : indicator(*new indicator_t())
  {}
  deterministic_indicator_caching(indicator_t&_ind) : indicator(_ind)
  {}
  ret_t operator()(const subject_t&subject)
  { key_t key = canonical(subject);
    // if the value isn't cached, calculate it and cache it
    typename cache_t::iterator past = cache.find(key);
    if (past != cache.end()) // if the value is cached, just return it
    { //cout << past->second << '\n';
      return past->second;
    }
    else // if not, calculate it, store it in the cache, return it
    { ret_t &answer = cache[key] = indicator(subject);
    //cout << answer << '\n';
      return answer;
    }
  }
};

template<typename S, typename R, typename I, typename K>
bool indicator_is_stochastic(deterministic_indicator_caching<S,R,I,K>&d)
{ return false; }

// an indicator that has the same value no matter what you apply it to.
// useful as the denominator of a quotient of two indicators, for instance.
// (maybe lambda_indicator makes this obsolete?)
class constant_indicator
{ double c;
public:
  constant_indicator(double d) : c(d)
  {}
  template<typename arg_t>
  double operator()(const arg_t&n)
  { return c; } 
};

// constant 0 indicator
template<typename subject_t>
double zero_indicator(const subject_t&n)
{ return 0;
}

// given two indicators, to be used under different circumstances (for
// instance, use monte carlo sampling when graph is large, but
// calculate fixation probability directly when it's small), this
// indicator calls one or the other of the two, depending on whether
// the 'threshold indicator' is below the threshold.
template<typename lo_ind_t, typename hi_ind_t, typename th_ind_t>
class switching_indicator_t
{
public:
  lo_ind_t&lo_ind;
  hi_ind_t&hi_ind;
  th_ind_t&th_ind;
  double threshold;
  switching_indicator_t(lo_ind_t&_li, th_ind_t&_ti, double _th, hi_ind_t&_hi)
    : lo_ind(_li), hi_ind(_hi), th_ind(_ti), threshold(_th)
  {}
  template<typename arg_t>
  indicator_value operator()(const arg_t&a)
  { if (th_ind(a) < threshold)
      return lo_ind(a);
    else
      return hi_ind(a);
  }
};

// function to create one of these without having to say all the
// template types
template<typename lo_ind_t, typename hi_ind_t, typename th_ind_t>
switching_indicator_t<lo_ind_t, hi_ind_t, th_ind_t>
switching_indicator(lo_ind_t&_li,th_ind_t&_ti,double _th,hi_ind_t&_hi)
{ return switching_indicator_t<lo_ind_t, hi_ind_t, th_ind_t>(_li,_ti,_th,_hi);
}

// switching indicator is stochastic if any of the things it calls is
template<typename lo_ind_t, typename hi_ind_t, typename th_ind_t>
bool indicator_is_stochastic(switching_indicator_t<lo_ind_t,hi_ind_t,th_ind_t>&ind)
{ return indicator_is_stochastic(ind.lo_ind)
    || indicator_is_stochastic(ind.hi_ind)
    || indicator_is_stochastic(ind.th_ind);
}

// return v0 for all arguments, except one special one, which is
// assigned v1
// arg_t has to have == operator
template<typename arg_t, typename ret_t>
class one_of_these_things_indicator
{ arg_t &special;
  ret_t v0, v1;
public:
  one_of_these_things_indicator(ret_t _0, arg_t &_s, ret_t _1)
    : special(_s), v0(_0), v1(_1)
  {}
  ret_t operator()(const arg_t&arg)
  { return arg == special ? v1 : v0; }
};

template<typename arg_t, typename ret_t>
bool indicator_is_stochastic(one_of_these_things_indicator<arg_t,ret_t>&ind)
{ return false;
}

// this indicator is attached to a variable in memory, and every time
// you call it it returns that variable's value.  so you can plot it
// easily in an IndicatorsDisplay, etc.
template<typename var_t>
class variableWatchIndicator
{ const var_t &v;
public:
  variableWatchIndicator(const var_t &_v) : v(_v) {}
  template<typename arg_t>
  double operator()(const arg_t& a)
  { return double(v); }
};

template<typename var_t>
variableWatchIndicator<var_t> variableWatcher(var_t&v)
{ return variableWatchIndicator<var_t>(v); }

// indicator returning the quotient of two other indicators.
// for instance, fixation probability divided by moran probability.
template<typename ind_t, typename ref_t>
class relative_indicator_t
{
public:
  ind_t &ind;
  ref_t &ref;
  relative_indicator_t(ind_t&_ind,ref_t&_ref) : ind(_ind), ref(_ref)
  {}
  template<typename network_t>
  // todo: how to figure accuracy when denominator is stochastic?
  indicator_value operator()(const network_t&n)
  { indicator_value num = ind(n);
    indicator_value den = ref(n);
    return indicator_value(num.value()/den.value(),
                           den.accuracy()==0 ? num.accuracy()/den.value() : 0);
  }
};

template<typename ind_t, typename ref_t>
relative_indicator_t<ind_t,ref_t>
relative_indicator(ind_t&_i,ref_t&_r)
{ return relative_indicator_t<ind_t,ref_t>(_i,_r);
}

// quotient is stochastic if numerator or denominator is.
template<typename ind_t, typename ref_t>
bool indicator_is_stochastic(relative_indicator_t<ind_t,ref_t>&rind)
{ return indicator_is_stochastic(rind.ind)
    || indicator_is_stochastic(rind.ref); 
}

// see boost/lambda/lambda.hpp
// Lambda objects are anonymous functions.  this indicator apparatus
// could be very powerful.
template<typename l_t, typename x0_t, typename x1_t>
class LambdaIndicator_t
{
protected:
  l_t lambda;
  const x0_t &x0;
  const x1_t &x1;
public:
  LambdaIndicator_t(l_t &_l, const x0_t &_x0, const x1_t &_x1)
    : lambda(_l), x0(_x0), x1(_x1)
  {}
  template<typename arg_t>
  indicator_value operator()(const arg_t&a)
  { return lambda(x0,x1); }
};

template<typename l_t, typename x0_t, typename x1_t>
LambdaIndicator_t<l_t,x0_t,x1_t>
LambdaIndicator(l_t&_l,const x0_t &_x0, const x1_t &_x1)
{ return LambdaIndicator_t<l_t,x0_t,x1_t>(_l,_x0,_x1); }

// generic_indicator is a shell for some other indicator
// allowing you to put them together in a collection object, on either
// side of a ?: operator, what not.  IndicatorsDisplay should use
// this, but (at the time of writing this comment) uses its own,
// nearly identical class.

// generic_indicator uses this hierarchy of template classes to house
// the different indicator types
// pure virtual base class
template<typename arg_t>
class generic_indicator_ref_base
{
public:
  virtual indicator_value operator()(const arg_t&arg) = 0;
  // indicator_is_stochastic has to go through this virtual member
  // function to get to the actual indicator type
  virtual bool _indicator_is_stochastic() = 0;
};

// a derived class for each indicator type here
template<typename ind_t,typename arg_t>
class generic_indicator_ref : public generic_indicator_ref_base<arg_t>
{ 
protected:
  ind_t&ref;
public:
  generic_indicator_ref(ind_t&ind_ref)
    : ref(ind_ref) {}
  indicator_value operator()(const arg_t&arg)
  { return ref(arg); }
  bool _indicator_is_stochastic()
  { return indicator_is_stochastic(ref);}
};

// the generic_indicator, which is a single type that uses the above
// to handle the differences in type
template<typename arg_t>
class generic_indicator
{
//protected:
public:
  generic_indicator_ref_base<arg_t>*ref;
public:
  generic_indicator() : ref(0) {}
  template<typename ind_t>
  generic_indicator(ind_t&ind)
    : ref(new generic_indicator_ref<ind_t,arg_t>(ind))
  {}
  indicator_value operator()(const arg_t&a)
  { return (*ref)(a); }
};

template<typename arg_t>
bool indicator_is_stochastic(generic_indicator<arg_t>&i)
{ return i.ref->_indicator_is_stochastic(); }

// this uses one indicator as its value, but first compares it to the
// value of another one [by printing both], so you can see if the two
// are giving the same result.  helpful with debugging.
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
    cout << " [" << v0 << ' ' << v1 << "]\n";
//  if (abs(sv - av) > 0.15)
//    av = ai(network);
    return v0;
  }
};

// moran fixation probability for well mixed populations.
// not strictly an indicator, but generally useful.
double moran_probability(double r, double n)
{ return (1 - 1/r) / (1 - pow(r,-n));
}

#if 0  // now defined by boost library
#include "boost/graph/property_maps/constant_property_map.hpp"
#else
//this is modeled after boost's identity_property_map
// it can, for instance, assign weight 1 to all edges,
// which is necessary in finding shortest paths
#include "boost/property_map.hpp"
template<class _key_t, class _value_t>
struct constant_property_map
  : public boost::put_get_helper<_value_t, 
      constant_property_map<_key_t,_value_t> >
{ typedef _key_t key_type;
  typedef _value_t value_type;
  typedef _key_t reference;
  typedef boost::readable_property_map_tag category;
  value_type value;
  constant_property_map(value_type _v) : value(_v) {}
  inline value_type operator[](const key_type& _t) const
  { return value; }
};
#endif

// for monte carlo fixation probability simulations.  try a single
// sample, see if it goes to fixation or extinction, or some other
// condition arises.  also record how long it took, and what was the
// maximum number infected (in case of extinction).  record the
// answers in this object and return it.
typedef enum { FIXATION, EXTINCTION, TIMEOUT, ERROR } fixation_result_type;
class fixation_result
{
public:
  fixation_result_type result;
  double time;
  int peak;
  fixation_result(fixation_result_type _r, double _t, int _p=0)
    : result(_r), time(_t), peak(_p) {}
};

// superclass for fixation_stats
class fixation_record
{ public:
  fixation_record() : n_fixations(0), n_extinctions(0), n_trials(0) {}
  unsigned n_fixations;
  unsigned n_extinctions;
  unsigned n_trials;
};

// in selection_indicator, below, an object of this type is stored for
// each argument.  if the same argument is received again, this object
// will be updated to reflect all results seen so far for that
// argument.  in this way, estimates improve until the desired
// accuracy is reached.
class fixation_stats : public fixation_record
{ public:
  double sum_fixationtime;
  double sum_extinctiontime;
  unsigned sum_peakbeforeextinction;
  unsigned max_peakbeforeextinction;
  // for figuring prob. of fixation given population of x
  // record how many times x is hit by all trajectories leading to
  // fixation and to extinction.
  //     vector<int> histogram_success;
  //     vector<int> histogram_failure;
  fixation_stats() : sum_fixationtime(0), sum_extinctiontime(0),
                     sum_peakbeforeextinction(0),
                     max_peakbeforeextinction(0)
  {}
  unsigned fixation_N()
  { return n_fixations + n_extinctions; }
  double fixation_probability()
  { return double(n_fixations)/fixation_N(); }
  double fixation_accuracy()
  { return Bernoulli_accuracy(fixation_probability(), fixation_N());
  }
};

// given an argument type that is capable of 'fixation' this
// class samples the probability of fixation, and caches its
// results so that multiple calls improve the estimate.
//
// for a specific scenario, construct a subclass that provides
//   try_fixation() and optionally test_fixation_candidate().
template<typename argument_t, typename RNG_t, typename params_t,
         typename fst=fixation_stats>
class selection_indicator : public ObjectWithParameters<params_t>
{
public:
  typedef fst fixation_stats_t;
  using ObjectWithParameters<params_t>::params;
protected:
  typedef map< string, fixation_stats_t > cache_t;
  cache_t cache;

public:
  RNG_t rng;

  // this is sometimes useful, points to the cache for the last thing
  // passed in for evaluation
  fixation_stats_t *currently_processing;

  selection_indicator(RNG_t _rng, params_t*_p=0)
    : rng(_rng), currently_processing(0)
  { inheritParametersFrom(_p); }

  // return true if there's an easy answer for fixation probability
  // (and provide it if so)
  // derived classes that can make this function perform should do so
  virtual bool test_fixation_candidate(const argument_t&a, double*answer_p)
  { return false; }

  // do a single realization, see whether it comes to fixation,
  // extinction, or some kind of trouble; return the result.
  // this is the main thing for derived classes to provide.
  // it would be pure virtual except sometimes it's convenient to
  // override operator() instead and not implement this function.
  virtual fixation_result try_fixation(const argument_t&)
  {}

  virtual string cache_key(const argument_t&a)
  { return canonical(a); }

  // do a number of fixation samples, return sample probability.
  indicator_value operator()(const argument_t &a)
  { // first, run analytic test: if the graph isn't connected, or 
    // whatever, fixation is impossible.
    { double answer;
      if (test_fixation_candidate(a,&answer))
        // try and verify
        //cout << "predict fixation probability: " << answer << '\n';
        return answer;
    }

    currently_processing = &cache[cache_key(a)];
#ifdef DISPLAY
    IndicatorsDisplayController<argument_t,params_t> trialsDisplay;
    if (params.animate_estimates())
    { trialsDisplay.inheritParametersFrom(this->params);
      trialsDisplay.params.setdisplayEvery(params.timesToTestFixation()/20);
      trialsDisplay.setoutfilenamebase("estimates");
      trialsDisplay.installIndicator(
                       LambdaIndicator(_1/ret<double>(_2),
                                       currently_processing->n_fixations,
                                       currently_processing->n_trials),
                       "sample fixation probability");
    }
#endif
    //cout << "testing for fixation: ";
    int times = params.timesToTestFixation();
    double desired_accuracy = params.monte_carlo_selection_accuracy();
    double p;
    double accuracy = HUGE;
    int tests = 0;
    while ((desired_accuracy > 0 && accuracy > desired_accuracy)
           || (desired_accuracy <= 0 && tests < times))
    { fixation_result trial_result = try_fixation(a);
      switch (trial_result.result)
      {
      case FIXATION:
        ++currently_processing->n_fixations;
        currently_processing->sum_fixationtime += trial_result.time;
        break;
      case EXTINCTION:
        ++currently_processing->n_extinctions;
        currently_processing->sum_extinctiontime += trial_result.time;
        currently_processing->sum_peakbeforeextinction += trial_result.peak;
        if (trial_result.peak > currently_processing->max_peakbeforeextinction)
          currently_processing->max_peakbeforeextinction = trial_result.peak;
        break;
      case TIMEOUT:
      default:
        break;
      }
      ++currently_processing->n_trials;
#ifdef DISPLAY
      if (params.animate_estimates())
      { trialsDisplay.update(currently_processing->n_trials,a);
      }
#endif
      ++tests;
      double N =
        currently_processing->n_fixations + currently_processing->n_extinctions;
      p = currently_processing->n_fixations/N;
      // N = var / (accuracy * p)^2
      // or accuracy = sqrt(var/N) / p
// 
// for example, if p = 0.25, accuracy = 0.1, we find N = 300
// (Bernoulli distribution, mean = p, var = p(1-p))
      if (p == 0 || p == 1)
        accuracy = 1/N; // this is a problem case, but is this the solution?
      else
        accuracy = sqrt(p*(1-p)/N)/p;
    }
    //cout << '\n';

    return indicator_value(p,accuracy);
  }

  fixation_stats_t&provide_stats(const argument_t&a)
  { return cache[cache_key(a)]; }
};

// this may not be needed since it's also provided for subclasses
template<typename argument_t, typename RNG_t, typename params_t>
bool indicator_is_stochastic(selection_indicator<argument_t,RNG_t,params_t>&i)
{ return true; }

// the selection_indicator classes are capable of collecting mean
// fixation time and related stats as they go, but they don't report
// it to callers.  these indicators attach to that one giving you an
// indicator that accesses the recorded mean
// fixation/absorption/extinction time.  But maybe this is a waste
// of typing...
template<typename sel_ind_class>
class mean_fixation_time_indicator
{
protected:
  sel_ind_class&sel_ind;
public:
  mean_fixation_time_indicator(sel_ind_class&_si) : sel_ind(_si) {}
  template<typename arg_t>
  double operator()(const arg_t&a)
  { typename sel_ind_class::fixation_stats_t&ncache = sel_ind.provide_stats(a);
    return ncache.sum_fixationtime/ncache.n_fixations;
  }
};

template<typename sel_ind_class>
class mean_extinction_time_indicator
{
protected:
  sel_ind_class&sel_ind;
public:
  mean_extinction_time_indicator(sel_ind_class&_si) : sel_ind(_si) {}
  template<typename arg_t>
  double operator()(const arg_t&a)
  { typename sel_ind_class::fixation_stats_t&ncache = sel_ind.provide_stats(a);
    return ncache.sum_extinctiontime/ncache.n_extinctions;
  }
};

template<typename sel_ind_class>
class mean_absorption_time_indicator
{
protected:
  sel_ind_class&sel_ind;
public:
  mean_absorption_time_indicator(sel_ind_class&_si) : sel_ind(_si) {}
  template<typename arg_t>
  double operator()(const arg_t&a)
  { typename sel_ind_class::fixation_stats_t&ncache = sel_ind.provide_stats(a);
    return (ncache.sum_extinctiontime + ncache.sum_fixationtime)
      /(ncache.n_extinctions + ncache.n_fixations);
  }
};

template<typename sel_ind_class>
class mean_peak_before_extinction_indicator
{
protected:
  sel_ind_class&sel_ind;
public:
  mean_peak_before_extinction_indicator(sel_ind_class&_si) : sel_ind(_si) {}
  template<typename arg_t>
  double operator()(const arg_t&a)
  { typename sel_ind_class::fixation_stats_t&ncache = sel_ind.provide_stats(a);
    return ncache.sum_peakbeforeextinction
      / double(ncache.n_extinctions);
  }
};

// write out a vector of bits as a string of 0's and 1's.
// useful for debugging caching in some cases.
string vbstring(const vector<bool>&vb)
{ string vs(vb.size(),'0');
  string::iterator si = vs.begin();
  for(vector<bool>::const_iterator vi = vb.begin(); vi != vb.end(); ++vi,++si)
    if (*vi)
      *si = '1';
  return vs;
}

#endif //__indicators_general_h__
