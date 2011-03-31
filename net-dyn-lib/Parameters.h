/* -*- C++ -*- 
 */
#ifndef PARAMETERS_H
#define PARAMETERS_H
#include "util.h"
#include <fstream>
#include <map>
using std::map;
#include <time.h>
#include <math.h>

#ifdef __unix__
#define UNIX
#define DIR_SEP '/'
#else
#define WINDOWS
#define DIR_SEP '\\'
#endif

typedef map<string, string> dictionary_t;

class Parameters : protected dictionary_t
{
//protected:
public:
  Parameters *inheritFrom;

  inline void set(string key, string val)
  { operator[](key) = val; }
  inline void unset(string key)
  { erase(key); }
  inline const string *get(string key) const
  { const_iterator it=find(key); \
    if (it != end())      return &it->second; \
    else if (inheritFrom) return inheritFrom->get(key); \
    else                  return 0; \
  }
  // this is not needed by the program but useful in debugger
  inline const string *get(const char *key) const
  { return get(string(key));
  }

// Declaring parameters at compile time:
// they are all represented internally by associating a
// value string with a key string in the dictionary,
// accessed via set() and get().
// But we also provide typed set/get functions for
// the rest of the program to use.
// So for instance,
// 
// DECLARE_PARAM(bool,abortOnFullMoon)
// 
// would expand to
// 
// inline void setabortOnFullMoon(bool val)
// { set("abortOnFullMoon",bool_to_string(val)); }
// inline bool abortOnFullMoon()
// { const string *str = get("abortOnFullMoon");
//   return str? string_to_bool(*str) : bool();
// }
//
// the various string_to_* and *_to_string functions are
// declared in util.h.  If some other type is needed those
// may have to be added to.
  
// RTYPE is helpful in case of enum types, apparently
#define DECLARE_PARAM_3(RTYPE,TYPE,NAME) \
  inline void set##NAME(TYPE val) \
  { set(#NAME, TYPE##_to_string(val)); } \
  inline RTYPE NAME() const \
  { const string *str = get(#NAME);\
    return str? string_to_##RTYPE(*str) : RTYPE(); \
  }

#define DECLARE_PARAM(TYPE,NAME) \
  DECLARE_PARAM_3(TYPE,TYPE,NAME)

// here's some nasty stuff for making the get/set
// business work with enum types.  We should
// probably just disallow enums and use string
// values instead.

// brought over from the climate project, 9/13/2006

// macros here to generate string_to_* and *_to_string
// functions for enum types
//
// say you have typedef enum { CONTINUE, PANIC } todo_t;
// then ENUM_DECLARATIONS(todo_t) gives you
//  map<string,todo_t> s_todo_t_map;
//  map<todo_t,string> todo_t_s_map;
//  inline todo_t string_to_todo_t(const string &x)
//  { return s_todo_t_map[x]; }
//  inline string todo_t_to_string(en_t x)
//  { return todo_t_s_map[x]; }
//
// and
//  DECLARE_ENUM_VALUE(todo_t,CONTINUE)
//  DECLARE_ENUM_VALUE(todo_t,PANIC)
// gives you
//  s_todo_t_map["CONTINUE"] = CONTINUE;
//  todo_t_s_map[CONTINUE] = "CONTINUE";
//  s_todo_t_map["PANIC"] = PANIC;
//  todo_t_s_map[PANIC] = "PANIC";

// this one has to go in the class declaration,
// once for each enum type
// notice at() in place of operator[]: necessary because of
// const, but have to be confident it won't raise a 'not found'
// exception.
#define ENUM_DECLARATIONS(en_t) \
  map<string,en_t> s_##en_t##_map; \
  map<en_t,string> en_t##_s_map; \
  inline en_t string_to_##en_t(const string &x) const \
  { return s_##en_t##_map.at(x); }                    \
  inline string en_t##_to_string(en_t x) const \
  { return en_t##_s_map.at(x); }

// and this one has to go in the Parameters class's constructor,
// once for each possible value of each enum type.
#define DECLARE_ENUM_VALUE(en_t,val) \
  s_##en_t##_map[#val] = val; \
  en_t##_s_map[val] = #val;

public:
  Parameters() : inheritFrom(0)
  { DECLARE_ENUM_VALUE(update_rule_t,EGT);
    DECLARE_ENUM_VALUE(update_rule_t,VM);
    DECLARE_ENUM_VALUE(update_rule_t,SKYE);
  }

  //  --- global stuff

  // settings/ directory is presumed to be where the executable is
  DECLARE_PARAM(string, dirnameForSettings)

  // whether to print information to the terminal
  DECLARE_PARAM(bool,print_stuff)

  // if so, including the big matrix?
  DECLARE_PARAM(bool,print_adjacency_matrix)

  // what experiment to do
  //typedef enum { OPTIMIZE, SAMPLE } experiment_t;
  DECLARE_PARAM(string,experiment)
  // output
  DECLARE_PARAM(string, outputDirectory)
  //DECLARE_PARAM(double, reopenOutputFilesInterval)
  // name of file to record all parameter values to
  DECLARE_PARAM(string, paramsDumpFile)
  // not used - whether a displayController should put its timeseries
  // into a file
  DECLARE_PARAM(bool, writeToFiles)
  // whether a displayController should display animation
  DECLARE_PARAM(bool, displayToScreen)
  // how often displayController should sample data
  DECLARE_PARAM(double, recordEvery)
  // how often displayController should update its animation
  DECLARE_PARAM(double, displayEvery)
  // how much time is visible at most at once
  DECLARE_PARAM(double, scrollingWindow)

  DECLARE_PARAM(string, gnuplotCommand)
  DECLARE_PARAM(string, ghostviewCommand)
  // if this is "signal" we cause the eps viewer to update by sending it
  // a signal (i.e. gv); if not we do it by invoking the viewer program again
  // (i.e. evince).
  DECLARE_PARAM(string, ghostviewUpdateStrategy)
  // if this is set for any given display object, it overrides
  // the display class's default filenames
  DECLARE_PARAM(string, outfilenamebase)
  // set this to repeat something.  if it's unset, it'll use the system
  // clock for a unique seed.
  DECLARE_PARAM(long, randSeed)

  //  --- for animations

  DECLARE_PARAM(bool, update_only_on_network_change)

  //  --- for the main function

  // start with a random graph with this many vertices and edges
  DECLARE_PARAM(unsigned, n_vertices)

  // kind of graph to start with:
  // COMPLETE, RANDOM, STAR, SUPERSTAR, POWER, PLOD
  DECLARE_PARAM(string, initial_graph_type)

  // density (used with only some graph types)
  DECLARE_PARAM(double, graph_density)

  // exponent of indegree distribution (for POWER graph type)
  DECLARE_PARAM(double, pl_in_exp)

  // exponent of outdegree distribution (for POWER graph type)
  DECLARE_PARAM(double, pl_out_exp)

  // exponent for PLOD random graph (or is it?)
  DECLARE_PARAM(double, pl_alpha)

  // y intercept for PLOD random graph (or is it?)
  // if this is set to 0 it's replaced by a default value
  DECLARE_PARAM(double, pl_beta)

  // # of lobes for superstar graph
  // if left to zero, program will choose based on # of vertices
  DECLARE_PARAM(int, superstar_m)
  // # of leaves per lobe for superstar graph
  // if left to zero, program will choose based on # of vertices
  DECLARE_PARAM(int, superstar_l)

  // if true, turn all the edges backward after creating initial graph
  DECLARE_PARAM(bool, reverse_initial_graph)

  // if true, add or remove one edge of initial graph before beginning
  DECLARE_PARAM(bool, flip_one_initial_edge)

  //  --- fixation experiments

  // if true, record the time series of weighted moments of the
  // network state for a representative realization
  DECLARE_PARAM(bool, record_dynamics)

  // if true, display the evolutionary dynamics as the population
  // walks to fixation or extinction
  DECLARE_PARAM(bool, animate_fixation)

  // if true, display partial estimates of fixation probability
  // as individual trials are sampled
  DECLARE_PARAM(bool, animate_estimates)

  // for deciding whether to find fixation probability
  // analytically (by solving a matrix equation)
  // or by monte carlo simulation
  DECLARE_PARAM(unsigned, max_n_vertices_for_matrix)

  //  --- for the fixation process

  DECLARE_PARAM(bool, skipFixation)

  // if true, the contact process is:
  //   choose a parent weighted by fitness
  //   choose a child among the parent's out-neighbors, weighed by edge weight
  //   (this is the evolutionary graph theory process)
  // if false, instead it's:
  //   choose a child at random
  //   choose a parent among the child's in-neighbors, weighted by
  //     fitness and edge weight
  // (though so far, edge weights are assumed to be 1)
  //DECLARE_PARAM(bool, choose_parent_first)
  typedef enum { EGT, VM, SKYE } update_rule_t;
  ENUM_DECLARATIONS(update_rule_t)
  DECLARE_PARAM(update_rule_t, update_rule)

  // fitness at the 'normal' vertices
  DECLARE_PARAM(double, residentFitness)

  // fitness of the mutant type, which might go to fixation
  DECLARE_PARAM(double, mutantFitness)

  // if so, don't evaluate fixation prob unless the graph is
  // strongly connected
  DECLARE_PARAM(bool, require_strongly_connected)

  // if this is assigned, one of the vertices has a slightly different
  // mutant fitness
  DECLARE_PARAM(unsigned, perturb_r_at_vertex)

  // and then this is that mutant fitness
  DECLARE_PARAM(double, perturbed_mutant_fitness)

  // if this is assigned, don't mutate one vertex, do all the vertices
  // above/below a certain degree.  values are "low-in", "high-out", etc.
  DECLARE_PARAM(string, initial_configuration)

  // and if you set that, set this too
  DECLARE_PARAM(unsigned, initial_configuration_threshold)

  //  --- monte carlo fixation

  // simulate introduction and possible fixation this many times at a
  //  go, to estimate probability of fixation
  // but note, these results are stored and more trials are added in at
  //  later calls, so this is only the minimum sample size.
  DECLARE_PARAM(int, timesToTestFixation)

  // alternatively, if this is set to a positive number, instead of
  // doing a fixed number of trials, keep doing trials until this
  // level of accuracy is reached.
  DECLARE_PARAM(double, monte_carlo_selection_accuracy)

  // if no fixation or extinction after this many reproduction events,
  // quit trying.  NOTE this might be measured in other time units, such
  // as reproduction events PER VERTEX.
  DECLARE_PARAM(int, maxStepsToFixation)

  //  --- for the optimizer

  // to control whether it displays moving graphics of the optimization
  // process
  DECLARE_PARAM(bool,animate_optimization)
#define DISPLAY 1

  // probability distribution of different kinds of network mutations
  DECLARE_PARAM(double, toggle_edge_mutation_probability)
  DECLARE_PARAM(double, double_toggle_mutation_probability)
  DECLARE_PARAM(double, remove_vertex_mutation_probability)
  DECLARE_PARAM(double, merge_vertices_mutation_probability)
  DECLARE_PARAM(double, attach_vertex_mutation_probability)
  DECLARE_PARAM(double, reattach_vertex_mutation_probability)
  DECLARE_PARAM(double, merge_and_clone_mutation_probability)
  // and split_edge_mutation_probability is the remainder

  //  --- hill climbing

  // if no changes improve the network, but some are neutral, take this
  // many neutral steps before quitting
  DECLARE_PARAM(int, maxNeutralSteps)

  // when we think we have an optimum, we try this many times before
  //  we're sure.  >1 is only useful when the indicator is stochastic.
  DECLARE_PARAM(int, trialsForLocalOptimum)

  //  --- simulated annealing

  // temperature starts here
  DECLARE_PARAM(double, annealing_initial_temperature)

  // and declines by this much at each step, until zero
  DECLARE_PARAM(double, annealing_cooling_rate)

  // this many mutations per time step
  // this will be supplied by the program, not the user
  DECLARE_PARAM(double, annealing_rate_factor)

  virtual void inheritValuesFrom(Parameters*upstream)
  { inheritFrom = upstream; }

  // get some settings from a file
  virtual void loadDefaults(string file)
  { parseSettingsFile(file);
    // if not set explicitly, rand seed is time
    if (!get("randSeed"))
      setrandSeed(time(0)*getpid());
  }
  void parseSettingsFile(string filename);
  void parseSettings(istream &);
  virtual string cleanLine(string);
  virtual void parseLine(string);

  // write them back out to a file
  virtual void writeAllSettings(ostream&);

  // read settings and file inclusions from command line invocation
  virtual void handleArgs(int argc, char **argv, const char *init_settings=0);

  // this is called after all settings are read in
  virtual void afterSetting(void)
  { string dumpfile = 
      outputDirectory()+"/"+paramsDumpFile();
    int mp = mkpath(dumpfile.c_str());
    if (mp)
    { string perr = 
          "Creating output directory and settings file " + dumpfile;
      perror(perr.c_str());
    }
    ofstream dumpto(dumpfile.c_str());
    writeAllSettings(dumpto);
  }

  // this is called after all the objects are set up and initialized
  virtual void finishInitialize(void)
  { }

  virtual ~Parameters() {}
};

template<typename ParametersClass>
class ObjectWithParameters : public ParametersClass
{
public:
  ParametersClass &params;
  ObjectWithParameters() : params(static_cast<ParametersClass&>(*this)) {}
  // assignment operator needs to be defined explicitly, so it calls
  // the default constructor and gets 'params' right
  ObjectWithParameters<ParametersClass>
    &operator=(const ObjectWithParameters<ParametersClass>&other)
  {}

  void inheritParametersFrom(Parameters*p)
  { params.inheritValuesFrom(p); }
  void inheritParametersFrom(Parameters&p)
  { inheritParametersFrom(&p); }
  template<typename otherParametersClass>
  void inheritParametersFrom(ObjectWithParameters<otherParametersClass>*op)
  { params.inheritValuesFrom(op ? &op->params : 0); }
};

#endif //PARAMETERS_H
