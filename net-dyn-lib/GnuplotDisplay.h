/* -*- C++ -*- */
#ifndef GNUPLOT_DISPLAY_H
#define GNUPLOT_DISPLAY_H
#include "Display.h"
#include "exec-stream.h"
#include "Parameters.h"
#include <fstream>
#include <iostream>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
using namespace std;

extern vector<exec_stream_t*> estream_boneyard;

template<typename params_t>
class GnuplotDisplay : public Display, public ObjectWithParameters<params_t>
{
public:
  using ObjectWithParameters<params_t>::params;
protected:
  bool started;
  // getting fancy: recycle estreams in a boneyard, so we don't waste
  // time killing and starting gnuplots
  exec_stream_t *estream;
  // gnuplot is an output stream for the logfile
  // -- subclasses should use it, as the instructions go to gnuplot
  // gpProcess is a wrapper for gnuplot's stdin
  ostream &gpProcess;
public:
  ofstream gnuplot;        
protected:
  string openFilename;

  static exec_stream_t *get_estream()
  { if (estream_boneyard.empty())
      return 0;
    //else
    exec_stream_t*bk = estream_boneyard.back();
    estream_boneyard.pop_back();
    return bk;
  }
  static void give_estream(exec_stream_t*es)
  { estream_boneyard.push_back(es); }
  
public:
  // don't open estream in constructor because of virtual
  //  gnuplotLogFile() -- do it in initialize()
  GnuplotDisplay()
    : started(true),
      estream(get_estream() ? : ((started=false), new exec_stream_t)),
      gpProcess(estream->in())
  {}
  
  virtual ~GnuplotDisplay()
  { try
    { give_estream(estream);
      //estream.kill();
      estream = 0;
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
        if (params.displayToScreen())
	{ vector<string> inv = gnuplotInvocation();
	  estream->start(inv[0], inv.begin()+1, inv.end());
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
    return params.outputDirectory() + "/default.gp";
  }

  // Appropriate call for gnuplot (or to bit bucket if display is off)
  virtual vector<string> gnuplotInvocation(void)
  {
    vector<string> inv;
    inv.push_back("gnuplot");
    inv.push_back("-noraise");
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
    if (params.displayToScreen())
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
    if (params.displayToScreen())
      gpProcess << "clear" << endl;
  }
};

#endif // GNUPLOT_DISPLAY_H
