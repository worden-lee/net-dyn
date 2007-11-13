#include "TimeSeriesDisplay.h"
#include "CSVController.h"
#include "indicators.h"

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

template<typename subject_t>
class IndicatorsDisplayController : public TimeSeriesController<int>
{
  vector< indicatorholder_base<subject_t>* > indicators;
  vector<string> names;
  CSVDisplay *csvdisplay;
public:
  IndicatorsDisplayController(double displayPeriod, double wind=HUGE)
    : TimeSeriesController<int>(displayPeriod,0,wind)
  {} 
  template<typename indic_t>  
  void installIndicator(indic_t*i,string name)
  { if (!display)
      createDisplay();
    display->setTitle(indicators.size(),name);
    indicators.push_back(new indicatorholder<subject_t,indic_t>(*i));
    names.push_back(name);
  }
  template<typename indic_t>
  void installIndicator(indic_t&i,string name)
  { installIndicator(&i,name);
  }
  // install a copy if its a temporary object, best to avoid
  template<typename indic_t>
  void installIndicator(indic_t i,string name)
  { installIndicator(new indic_t(i),name);
  }
  
  string outdir()
  { return "out";
  }
  string filenamebase(void)
  { return "indicators";
  }
  int color(int k)
  { return k+1;
  }
  void recordFile(double t) // shouldn't be called
  { cerr << "wrong IndicatorsDisplay::recordFile" << endl; }
  void recordFile(double t, const subject_t&subject)
  { for (int i = 0; i < indicators.size(); ++i)
    { double ind_i = (*indicators[i])(subject);
      display->record(i, t, ind_i);
      *csvdisplay << ind_i;
    }
    csvdisplay->newRow();
  }
  void update(double t, const subject_t&subject)
  { if (!display)
      createDisplay();
    if (t == 0)
    { for (int i = 0; i < names.size(); ++i)
        *csvdisplay << (i==0? "#"+names[i] : names[i]);
      csvdisplay->newRow();
    }
    recordFile(t,subject);
    DisplayController< TimeSeriesDisplay<int> >::update(t);
  }
  void createDisplay()
  { if (!display)
    { display = new TimeSeriesDisplay<int>(this);
      display->initialize();
    }
    if (!csvdisplay)
    { csvdisplay = new CSVDisplay(outdir()+'/'+filenamebase()+".csv");
    }
  }
};
