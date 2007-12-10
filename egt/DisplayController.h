/* -*- C++ -*- */
#ifndef DISPLAYCONTROLLER_H
#define DISPLAYCONTROLLER_H
#include <math.h>
#include <string>
using std::string;
#include <vector>
using std::vector;

template<class DISPLAY>
class DisplayController// : public OutputController
{
protected:
  DISPLAY *display;
  double recordEvery, recordCounter;
  double displayEvery, displayCounter;
  
public:
  // your subclass needs to assign the site pointer, either
  // in its constructor or thereafter.
  DisplayController(double displayPeriod, double filePeriod)
    : display(0),
      recordEvery(filePeriod), recordCounter(-HUGE),
      displayEvery(displayPeriod), displayCounter(-HUGE)
      // if ...Every < 0 never do it, if ...Every == 0 always do it
    // set ...Counter=0 to skip plotting at time=0
  {}
  virtual ~DisplayController() 
  {
  }

  // your subclass constructor should provide this, to allocate
  //  and initialize an appropriate Display object,
  //  and assign it to display.
  // this gets called from update() so be warned if you don't use that.
  virtual void createDisplay() = 0;
  
  void period(int period)
  { displayEvery = recordEvery = period; }
  
  virtual string outdir()
  { return "."; }

  virtual void flush(void)
  { if (display) display->flush();
  }

  virtual void finish(void) {}

  virtual void updateDisplay(double t)
  { if (!display)
      createDisplay();
    display->updateDisplay(t);
  }
  virtual void recordFile(double t)
  { if (!display)
      createDisplay();
    display->recordFile(t);
  }

  virtual void update(double t)
  {
    if (!display)
      createDisplay();
    if (recordEvery == 0 ||
        (recordEvery > 0 && recordCounter+recordEvery <= t))
    {
      recordFile(t);
      recordCounter = t;
    }
    if (displayEvery == 0 ||
        (displayEvery > 0 && displayCounter+displayEvery <= t))
    {
      updateDisplay(t);
      displayCounter = t;
    }
  }
};

template<typename argument_t>
class DisplayFarm
{
protected:
  class DisplayControllerShimRoot
  {
  public:
    //virtual void update(double t) = 0;
    virtual void update(double t, const argument_t&a) = 0;
  };
  template<typename DC_t>
  class DisplayControllerShim : public DisplayControllerShimRoot
  {
  protected:
    DC_t*dc;
  public:
    DisplayControllerShim(DC_t*_dc) : dc(_dc){}
//     virtual void update(double t)
//     { dc->update(t); }
    virtual void update(double t, const argument_t&a)
    { dc->update(t,a); }
  };
  typedef vector<DisplayControllerShimRoot*> dcs_t;
  dcs_t displayControllers;

public:
  DisplayFarm(){}
  template<typename DC_t>
  void installController(DC_t*dc)
  { displayControllers.push_back(new DisplayControllerShim<DC_t>(dc)); }
  template<typename DC_t>
  void installController(DC_t&dc)
  { displayControllers.push_back(new DisplayControllerShim<DC_t>(&dc)); }
//   void update(double t)
//   { for (typename dcs_t::iterator dci = displayControllers.begin();
//          dci != displayControllers.end(); ++dci)
//       (*dci)->update(t);
//   }
  void update(double t, const argument_t&a)
  { for (typename dcs_t::iterator dci = displayControllers.begin();
         dci != displayControllers.end(); ++dci)
      (*dci)->update(t,a);
  }
};

#endif // DISPLAYCONTROLLER_H
