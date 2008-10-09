//  -*- C++ -*-
#ifndef _INDICATORS_DISPLAY_H_
#define _INDICATORS_DISPLAY_H_

#include "TimeSeriesDisplay.h"
#include "CSVController.h"
//#include "indicators.h"

// an indicator is anything that has a double operator()
// could be a class or a function pointer
// this thing provides a common class to hold them in
template<typename subject_t>
class indicatorholder_base
{
public:
  virtual double operator()(const subject_t&subject) = 0;
};

template<typename subject_t, typename indicator_t>
class indicatorholder : public indicatorholder_base<subject_t>
{ indicator_t *inside;
public:
  indicatorholder(indicator_t&indicator) : inside(&indicator)
  {}
  double operator()(const subject_t&subject)
  { return (*inside)(subject); } 
};

static unsigned idc_serial_counter = 0;

template<typename subject_t,typename ParamsClass>
class IndicatorsDisplayController :
  public TimeSeriesController<int,ParamsClass>
{
  vector< indicatorholder_base<subject_t>* > indicators;
  vector<string> names;
  CSVDisplay *csvdisplay;
  unsigned serial;
  using ObjectWithParameters<ParamsClass>::params;
  using TimeSeriesController<int,ParamsClass>::display;
  using DisplayController< TimeSeriesDisplay<int,ParamsClass>,
                           ParamsClass >::recordCounter;
public:
  IndicatorsDisplayController()
    : csvdisplay(0), serial(idc_serial_counter++)
  {}

//   void inheritValuesFrom(Parameters*upstream)
//   { Parameters::inheritValuesFrom(upstream);
//     // take recordEvery to local storage
//     recordEvery = params.recordEvery();
//     // and get rid of the one the superclass uses
//     params.setrecordEvery(-1);
//   }
  template<typename indic_t>  
  void _installIndicator(indic_t*i,string name)
  { if (!display)
      createDisplay();
    display->setTitle(indicators.size(),name);
    indicators.push_back(new indicatorholder<subject_t,indic_t>(*i));
    names.push_back(name);
  }
  template<typename indic_t>  
  void installIndicator(indic_t*i,string name)
  { _installIndicator(i,name);
  }
  template<typename indic_t>
  void installIndicator(indic_t&i,string name)
  { _installIndicator(&i,name);
  }
  // install a copy if it's a temporary object, best to avoid
  template<typename indic_t>
  void installIndicator(const indic_t&i,string name)
  { _installIndicator(new indic_t(i),name);
  }
  template<typename PC>
  void installIndicator(ObjectWithParameters<PC>*o)
  { _installIndicator(o);
    o->inheritParametersFrom(this);
  }
  template<typename PC>
  void installIndicator(ObjectWithParameters<PC>&o)
  { installIndicator(&o);
  }
  
  string filenamebase(void)
  { string pfb = params.outfilenamebase();
    if (pfb != "") return pfb;
    // else
    ostringstream oss;
    oss << "indicators-" << serial;
    return oss.str();
  }
  int color(int k)
  { return k+1;
  }

  // this is called by the superclass' update() but does nothing
  void recordFile(double t)
  {}
  //{ cerr << "wrong IndicatorsDisplay::recordFile" << endl; }

  vector<double> lasttime;
  // this is called by our update(), does something
  // returns true if some value is changed since last time it was
  // called
  bool recordFile(double t, const subject_t&subject)
  { if (!display) createDisplay();
    lasttime.resize(indicators.size());
    bool worthdisplaying = false;
    *csvdisplay << t;
    for (int i = 0; i < indicators.size(); ++i)
    { double ind_i = (*indicators[i])(subject);
      display->record(i, t, ind_i);
      if (ind_i != lasttime[i])
      { lasttime[i] = ind_i;
        worthdisplaying = true;
      }
      *csvdisplay << ind_i;
    }
    csvdisplay->newRow();
    return worthdisplaying;
  }

  void update(double t, const subject_t&subject)
  { if (!display)
      createDisplay();
    if (t == 0)
    { for (int i = 0; i < names.size(); ++i)
        *csvdisplay << (i==0? "#"+names[i] : names[i]);
      csvdisplay->newRow();
    }
    bool worthredisplaying = false;
    if (recordCounter + params.recordEvery() <= t)
    { worthredisplaying = recordFile(t,subject);
      recordCounter = t;
    }
    if (worthredisplaying)
      DisplayController< TimeSeriesDisplay<int,ParamsClass>, ParamsClass >
        ::update(t);
  }
  void createDisplay()
  { if (!display)
    { display = new TimeSeriesDisplay<int,ParamsClass>(this);
      display->inheritParametersFrom(params);
      display->initialize();
    }
    if (!csvdisplay)
    { csvdisplay = new CSVDisplay(params.outputDirectory()
                                  +'/'+filenamebase()+".csv");
    }
  }
  ~IndicatorsDisplayController()
  { if (display)    delete display;
    if (csvdisplay) delete csvdisplay;
  }
};

#endif//_INDICATORS_DISPLAY_H_
