/* -*- C++ -*- */
// created by Jesse Edelstein 2007
// 9/2007 tweaked and adapted for network-optimize project by LW
#ifndef _CSVDISPLAY_H_
#define _CSVDISPLAY_H_
#include <string>
#include <fstream>
#include <iostream>

class CSVDisplay : public Display
{
protected:
  ofstream logfile;
  bool onBlankLine;
  //  string myOutfile;

public:
  CSVDisplay () {
    // included for compatibility, but a file is not open!
    // use openFile if this must be used
  }

  // preferred constructor form - pass in the filename to use
  CSVDisplay(string filename) {
    logfile.open(filename.c_str());
    onBlankLine = true;
    //    myOutfile = filename;
  }

//   ~CSVDisplay () {
//     logfile.close();
//   }

  // open a file, if there's not one loaded already
  void openFile(string filename) {
    if (logfile.is_open()) {
      cerr << "CSVDisplay::openFile called, but a file is already open!\n";
    } else {
      logfile.open(filename.c_str());
      onBlankLine = true;
      //      myOutfile = filename;
    }
  }

  // write a whole line (i.e., header data)
  void writeLine(const string line) {
    if (!logfile) {
       cerr << "CSVDisplay::writeLine called, but no logfile is open!\n";
    }
    logfile << line;
    // check if the line included a newline, and write as needed
    string::const_reverse_iterator lastchar = line.rbegin();
    if (*lastchar != '\n') logfile << endl;
    onBlankLine = true;
  }

  // Write a new datum -- should do comma insertion automagically
  // this should automatically handle other numericals
  CSVDisplay& operator<<(const double& data) { 
    if (&data == 0 || data == HUGE) {
      // HUGE: the international code for "it's broken!"
      *this << string("NA");
      return *this;
    }
    if (!logfile) {
      cerr << "CSVDisplay::<< called, but no logfile is open!\n";
    }
    if (onBlankLine) {
      onBlankLine = false;
    } else {
      logfile << ",";
    }
    logfile << data;
    return *this;
  }

  CSVDisplay& operator<<(const string data) {
    if (!logfile) {
      cerr << "CSVDisplay::<< called, but no logfile is open!\n";
    }
    if (onBlankLine) {
      onBlankLine = false;
    } else {
      logfile << ",";
    }
    logfile << data;
    return *this;
  }

  // start a new row of data
  // writes a new line only if we're not on a fresh line
  void newRow() {
    logfile << endl;
    onBlankLine = true;
  }

public:
  // included for compatibility, but I don't think these should be
  // called directly...
  void updateDisplay(double t) {} 
  void recordCommunity(double t) {}

  void flush() {
    logfile.flush();
  }
};

template<typename ParamsClass>
class CSVController : public DisplayController<CSVDisplay,ParamsClass>
{
  typedef DisplayController<CSVDisplay,ParamsClass> superclass;
  using superclass::display;
public:
  // We're stretching the DisplayController spec a little, because
  // this class will control many CSVDisplay classes: one that is
  // continually updated with general data as the simulation runs, and
  // many that contain snapshots of data for each live species.
  // 
  // So *display only points to the community logfile CSVDisplay. We
  // initialize the snapshot writers as needed.
  //
  // Of course there is no displayPeriod needed, because the CSVs are
  // not displayed, only written to disk.
  CSVController()
  { // create the community logfile
    createDisplay();
    superclass::recordCounter = this->time();
  }

  void createDisplay()
  { mkdir(outdir().c_str(), S_IRWXU|S_IRWXG|S_IRWXO);
    string communityOutfile = outdir() + "/community-log.csv";
    display = new CSVDisplay(communityOutfile);
    display->writeLine("temp,eco_t,evo_t,n_equil,n_species,n_resources,dT/dmaxT");
  }

  string outdir(){
    return ParamsClass::outputDirectory()+"/csv";
  }
  
  virtual double time() = 0;

  void step(double t) {}

  // more junk for completeness
  void flush(void) {
    display->flush();
  }
  void finish(void) {
    delete display;
  }
  void recordCommunity(void) {
    update(this->time());
  }
  void updateDisplay(void) {}

  void update(double t)
  {
//     if (!display)
//       createDisplay();
    double recordEvery = ParamsClass::recordEvery();
    if (recordEvery > 0 && superclass::recordCounter+recordEvery <= t && t != 0)
    {
      // (record stuff...)
      superclass::recordCounter = t;
    }
  }
};

#endif//_CSVDISPLAY_H_
