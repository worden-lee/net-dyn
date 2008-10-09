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

template<typename PC>
class DotDisplay;

template<typename PC>
class DotGraphController
  : public DisplayController<DotDisplay<PC>,PC>
{
public:
  using ObjectWithParameters<PC>::params;

  void createDisplay();
  
  virtual bool approve(void) 
  { return true; }
  virtual bool allowDisplaying(void)
  { return params.displayToScreen(); } 
  virtual bool allowRecording(void)
  { return true; }

  virtual string title(void)
  { return "Community Graph"; }
  virtual string dotdir()
  { return params.outputDirectory()+"/dot"; }
  virtual string filenamebase(void)
  { string pfb = params.outfilenamebase();
    if (pfb == "") pfb = "graph";
    return pfb;
  }
  virtual string graphfile(void)
  { return dotdir()+'/'+filenamebase()+".dot"; }
  virtual string epsfile(void)
  { return dotdir()+'/'+filenamebase()+".eps"; }
  virtual bool writes_positions()
  { return false; }
  
  //void _update(void);
//   void recordFile(double t);
//   void updateDisplay(double t);

  virtual ~DotGraphController();
};

template<typename ParamsClass>
class DotDisplay: public Display
{
protected:
  DotGraphController<ParamsClass> *xs;
  ofstream os;
  string graphfile;
  bool startedgv;
  int gvpid;
  bool wrote;
  bool useHUP;
  
public:
  DotDisplay(DotGraphController<ParamsClass> *a)
    : xs(a), startedgv(false), gvpid(-1), wrote(false),
      useHUP(xs->ghostviewUpdateStrategy() == "signal")
  {
    (void)mkdir(xs->dotdir().c_str(),S_IRWXU|S_IRWXG|S_IRWXO);
    //if(displayEvery>0)
    //  updateDisplay(0);
  }

  void recordFile(double t)
  {
    graphfile = xs->graphfile();
    os.open(graphfile.c_str());
    writeGraph(os);
    os.close();
    wrote = true;
  }
  virtual void writeGraph(ostream &ost)
  {
  }
  string ghostviewCommand()
  {
    return xs->params.ghostviewCommand();
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

    if (gvpid <= 0 || !useHUP) // if we need to start a gv process
    {
      int pid = fork();
      if (pid == 0)
      { // child process in here
	const string &gvcomm = ghostviewCommand();
	string eps = xs->epsfile();
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

  void updateDisplay(double t)
  {
    if (wrote)
    {
      string dotcomm =
	//string("/usr/bin/dot -Tps -o")
	string("/usr/bin/neato ") + (xs->writes_positions()? "-n1 ":"")
        + "-Goverlap=scale -Gsplines=true -Tps -o"
	+ xs->epsfile() + " " + graphfile;
      
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

template<class PC>
void DotGraphController<PC>::createDisplay()
{ DisplayController<DotDisplay<PC>,PC>::display
    = new DotDisplay<PC>(this);
}

template<class PC>
DotGraphController<PC>::~DotGraphController()
{ delete DisplayController<DotDisplay<PC>,PC>::display; }

#endif // DOTDISPLAY_H
