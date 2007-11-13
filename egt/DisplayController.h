/* -*- C++ -*- */
#ifndef DISPLAYCONTROLLER_H
#define DISPLAYCONTROLLER_H
#include <math.h>
#include <string>
using std::string;

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
    if (recordEvery > 0 && recordCounter+recordEvery <= t)
    {
      recordFile(t);
      recordCounter = t;
    }
    if (displayEvery > 0 && displayCounter+displayEvery <= t)
    {
      updateDisplay(t);
      displayCounter = t;
    }
  }
};

#endif // DISPLAYCONTROLLER_H