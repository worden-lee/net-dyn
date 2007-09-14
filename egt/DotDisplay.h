/* -*- C++ -*- */
#ifndef DOTDISPLAY_H
#define DOTDISPLAY_H

#include "Display.h"
#include "DisplayController.h"
#include <sys/types.h> // mkdir
#include <sys/stat.h>  // mkdir
#include <signal.h> // kill
#include <sys/wait.h> // waitpid
#include <errno.h>
#include <string>
#include <fstream>
#include <iostream>
using namespace std;

class DotDisplay;

class DotGraphController : public DisplayController<DotDisplay>
{
public:
  DotGraphController(double period, double fperiod);
  virtual ~DotGraphController();
  
  void createDisplay();
  
  virtual bool approve(void) 
  { return true; }
  virtual bool allowDisplaying(void)
  { return approve(); } 
  virtual bool allowRecording(void)
  { return true; }

  virtual string title(void)
  { return "Community Graph"; }
  virtual string outdir() = 0;
  virtual string dotdir()
  { return outdir()+"/dot"; }
  virtual string filebasename(void)
  { return dotdir()+"/graph"; }
  virtual string graphfile(void)
  { return filebasename()+".dot"; }
  virtual string epsfile(void)
  { return filebasename()+".eps"; }
  
  //void _update(void);
  void recordFile();
  void updateDisplay(void);
};

class DotDisplay: public Display
{
protected:
  DotGraphController *xs;
  ofstream os;
  string graphfile;
  bool startedgv;
  int gvpid;
  bool wrote;
  
public:
  DotDisplay(DotGraphController *a)
    : xs(a), startedgv(false), gvpid(-1), wrote(false)
  {
    (void)mkdir(xs->dotdir().c_str(),S_IRWXU|S_IRWXG|S_IRWXO);
    //if(displayEvery>0)
      updateDisplay();
  }

  void recordFile()
  {
    graphfile = xs->dotdir()+'/'+xs->graphfile();
    os.open(graphfile.c_str());
    writeGraph(os);
    os.close();
    wrote = true;
  }
  virtual void writeGraph(ostream &ost)
  {
  }
  static string ghostviewCommand()
  {
    return "/usr/bin/gv";
  }
  
  virtual void ensureGV(void)
  { // start/restart child ghostview if needed

    // gvpid is -1 if we haven't forked yet
    // if we have, make sure it's still running
    if (gvpid > 0)
    {
      int checkpid = waitpid(gvpid,0,WNOHANG);
      switch (checkpid)
      {
      case  0: // it's still running
	break;
      case -1: // error
	cerr << "error in DotDisplay::ensureGV() : "
	     << strerror(errno) << endl;
	break;
      default: // if return value is gvpid, it died
	if (checkpid != gvpid)
	  cerr << "DotDisplay: why does waitpid return "
	       << checkpid << " instead of " << gvpid
	       << "?" << endl;
// 	cerr << "DotDisplay " << this << ": waitpid(" << gvpid << ")"
// 	     << " returned nonzero : restarting gv" << endl;
	gvpid = -1; // need to start a new one
	break;
      }
    }

    if (gvpid <= 0) // if we need to start a gv process
    {
      int pid = fork();
      if (pid == 0)
      { // child process in here
	const string &gvcomm = ghostviewCommand();
	string eps = xs->dotdir() + "/" + xs->epsfile();
	//cout << "DotDisplay child process: execlp "
	//     << gvcomm << ' ' << buf << endl;
	execlp(gvcomm.c_str(), gvcomm.c_str(), eps.c_str(), NULL);
	// if we get here there's a problem
	cerr<<"DotDisplay had error starting "
	    << gvcomm << ' ' << eps << endl;
	exit(-1);
      }
      // parent process out here
      if (pid < 0)
	cerr << "DotDisplay" << this
	     << ": can't fork gv subprocess" << endl;
      else
      { gvpid = pid;
//         cerr << "DotDisplay " << this
// 	     << ": gv subprocess is " << gvpid << endl;
      }
    }
    else // if gv is already running
    { // tell it to refresh itself
//       cerr << "DotDisplay " << this
// 	   << ": send HUP to " << gvpid << endl;
      kill(gvpid,SIGHUP);
    }
  }

  void updateDisplay(void)
  {
    if (wrote)
    {
      string dotcomm =
	//string("/usr/bin/dot -Tps -o")
	string("/usr/bin/neato -Goverlap=scale -Gsplines=true -Tps -o")
	+ xs->dotdir() + "/" + xs->epsfile()
	+ " " + graphfile;
      
      //cout << dotcomm << endl;
      if(system(dotcomm.c_str()))
	cerr << "DotDisplay had problem running " << dotcomm << endl;

      ensureGV();
    } 
  }

  ~DotDisplay()
  {
    if (gvpid > 0)
    {
//       cerr << "DotDisplay::~DotDisplay() : kill child process "
// 	   << gvpid << endl;

      //kill(gvpid,SIGINT);

      //cout << "Waiting for gv to terminate" << endl;
      //waitpid(gvpid,0,0);
//       cerr << "DotDisplay::~DotDisplay() : "
// 	   << gvpid << " terminated" << endl;
    }
  }
};

#endif // DOTDISPLAY_H
