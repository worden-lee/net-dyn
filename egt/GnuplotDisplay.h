/* -*- C++ -*- */
#ifndef GNUPLOT_DISPLAY_H
#define GNUPLOT_DISPLAY_H
#include "Display.h"
#include "exec-stream.h"
#include <fstream>
#include <iostream>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
using namespace std;

class GnuplotDisplay : public Display
{
public:
  ofstream gnuplot;        
protected:
  exec_stream_t estream;
  ostream &gpProcess;
  string openFilename;
  bool started;

  // gnuplot is an output stream for the logfile
  // -- subclasses should use it, as the instructions go to gnuplot
  // gpProcess is a wrapper for gnuplot's stdin

public:
  // don't open estream in constructor because of virtual
  //  gnuplotLogFile() -- do it in initialize()
  GnuplotDisplay() : gpProcess(estream.in()), started(false) {}
  
  virtual ~GnuplotDisplay()
  { try
    { estream.kill();
    } 
    catch( std::exception const & e )
    { std::cerr << "error: " << e.what() << "\n";
    }
  }
  
  virtual void initialize(void)
  { if (!started)
    { try
      { openFilename = gnuplotLogFile();
        string::size_type slash = openFilename.rfind('/');
	if (slash != string::npos)
	{ string ofdir = openFilename.substr(0,slash);
	  mkdir(ofdir.c_str(),S_IRWXU|S_IRWXG|S_IRWXO);
	}
        gnuplot.open(openFilename.c_str(), ios_base::binary);
        //if (!parameters.disableDisplaying())
	{ vector<string> inv = gnuplotInvocation();
	  estream.start(inv[0], inv.begin()+1, inv.end());
	}
      }
      catch( std::exception const & e )
      { std::cerr << "GnuplotDisplay caught exception: "
		  <<  e.what()  <<  "\n";
      }
      started = true;
    }
  }

  // appropriate filename for gnuplot commands
  // if this returns different values at different times,
  // multiple files will be created.
  virtual string gnuplotLogFile(void)
  {
    return "default.gp";
  }

  // Appropriate call for gnuplot (or to bit bucket if display is off)
  virtual vector<string> gnuplotInvocation(void)
  {
    vector<string> inv;
//     if (parameters.disableDisplaying())
//     { inv.push_back("sh");
//       inv.push_back("-c");
//       inv.push_back("cat >/dev/null");
//     }
//     else
    { inv.push_back("gnuplot");
      inv.push_back("-noraise");
    }
    return inv;
  }

  void flush()
  { gnuplot.flush(); }

protected:
  // rewind to the beginning of the file
  void rewindStream()
  { string filename = gnuplotLogFile();
    if (filename == openFilename)
      gnuplot.seekp(0, ios::beg);
    else
    { gnuplot.close();
      gnuplot.open(filename.c_str(), ios_base::binary);
      openFilename = filename;
    }
  }

  // flush data to file, and send it to gnuplot (if we're doing that)
  void plotBufferedData() {
    flush();
    streampos pos = gnuplot.tellp();
    if (truncate(openFilename.c_str(),pos))
      cerr << "Error truncating " << openFilename << ": "
	   << strerror(errno) << endl;
    //if (!parameters.disableDisplaying())
      gpProcess << "load \"" << openFilename << "\"" << endl;
  }

  // subclasses should use this pattern for updating the display:
  // rewindStream();
  // gnuplot << whatever << commands << for << Gnuplot;
  // plotBufferedData();

  // it is permitted to return a different value from gnuplotLogFile()
  // to cause rewindStream() to close the file and open another one

  // Subclasses should implement this if they want
  // to write instructions for Gnuplot, such as "set title ..."
  virtual void writeGnuplotHeaders() {}

  // Clear the Gnuplot graph
  void clear()
  {
    //if (!parameters.disableDisplaying())
      gpProcess << "clear" << endl;
  }
};

#endif // GNUPLOT_DISPLAY_H
