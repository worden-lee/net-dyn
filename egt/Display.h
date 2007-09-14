/* -*- C++ -*- */
#ifndef DISPLAY_H
#define DISPLAY_H

class Display
{
public:
  virtual ~Display()
  {// cout << "~Display()" << endl;
  }
  virtual void initialize(void) {}
  virtual void flush() {}
  // override one or both of these
  virtual void updateDisplay(void) {}
  virtual void recordFile(void) {}
};

#endif // DISPLAY_H
