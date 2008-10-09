/* -*- C++ -*- */
// This draws a periodically changing histogram in Gnuplot.
// The center of the design is how you provide the histogram data
// via two iterators, one for x data and one for y data
#ifndef HISTOGRAM_DISPLAY_H
#define HISTOGRAM_DISPLAY_H
#include "DisplayController.h"
#include "GnuplotDisplay.h"
#include <iostream>
#include <sstream>

template<class KC, class PC>
class HistogramDisplayController;

template<class key_t, class params_t=Parameters>
class HistogramDisplay : public GnuplotDisplay<params_t>
{
public:
  using GnuplotDisplay<params_t>::gnuplot;
  
protected:
  HistogramDisplayController<key_t,params_t> &xs;

public:
  class series_t
  {
  public:
    // display attributes
    bool vertical;
    bool lines;
    string title;
    // data
    vector<double> ydata;
    //  if xdata is empty, the global xdata will apply
    vector<double> xdata;
    series_t() : vertical(false), lines(true) {}
  };
  vector<double> global_xdata;

  typedef map<key_t,series_t> series_map_t;
  series_map_t series_map;

  HistogramDisplay(HistogramDisplayController<key_t,params_t>*con);

  virtual ~HistogramDisplay() {}

  virtual void writeGnuplotHeaders(void)
  { gnuplot << "set parametric\n";
  }

  void updateDisplay(double t)
  { this->rewindStream();
    writeGnuplotHeaders();

    bool first = true;
    ostringstream ss;
    for (typename series_map_t::iterator mit = series_map.begin();
         mit != series_map.end(); ++mit)
    { if (first)
      { gnuplot << "plot "; first = false; }
      else
      { gnuplot << ", "; }
      series_t &series = mit->second;
      gnuplot << "'-' ";
      if (series.title != "")
        gnuplot << "title \"" << series.title << "\" ";
      gnuplot << "w " << (series.lines?'l':'p');
      vector<double>&xdata
        = (series.xdata.size() > 0 ? series.xdata : global_xdata);
      double lastx = 0;
      ss << lastx << " 0\n";
      int size = xdata.size();
      if (size > series.ydata.size()) size = series.ydata.size();
      for (int i = 0; i < size; ++i)
      { if (isnan(xdata[i]) || isnan(series.ydata[i]))
          continue;
        ss << lastx    << ' ' << series.ydata[i] << '\n'
                << xdata[i] << ' ' << series.ydata[i] << '\n';
        lastx = xdata[i];
      }
      ss << lastx << " 0\n";
      ss << "e\n";
    }
    gnuplot << '\n';
    gnuplot << ss.str();
    this->plotBufferedData();
  }

  string gnuplotLogFile(void)
  {
    return this->outputDirectory() + "/histogram.gp";
  }

  series_t &series(key_t k)
  { return series_map[k]; }
};

template<class key_t, class params_t>
class HistogramDisplayController :
  public DisplayController< HistogramDisplay<key_t, params_t>, params_t >
{
  friend class HistogramDisplay<key_t,params_t>;
public:
  typedef DisplayController< HistogramDisplay<key_t,params_t>, params_t >
    superclass;
protected:
  // newer g++'s are very picky about these things
  using superclass::display;
  using ObjectWithParameters<params_t>::params;

public:

  virtual ~HistogramDisplayController()
  { if (display) delete display;
  }

  virtual void createDisplay()
  { if (!display)
    { display = new HistogramDisplay<key_t,params_t>(this);
      display->inheritParametersFrom(this);
      display->initialize();
    }
  }

  // subclass could provide this
  virtual void recordFile(double t) {}
//   // or...
//
//   // receive the histogram data as an x iterator and a y iterator
//   template<typename xit_t, typename yit_t>
//   void recordData(xit_t x0, xit_t xend, yit_t yi)
//   { xdata.resize(xend - x0);
//     ydata.resize(xend - x0);
//     for (xit_t xi = x0; xi != xend; ++xi, ++yi)
//     { xdata[xi-x0] = *xi;
//       ydata[yi-y0] = *yi;
//     }
//   }

  // better, use these output iterators to put the x and y data into
  back_insert_iterator< vector<double> > global_x_inserter()
  { if (!display) createDisplay();
    display->global_xdata.clear();
    return back_inserter(display->global_xdata);
  }
  back_insert_iterator< vector<double> > x_inserter(key_t k)
  { if (!display) createDisplay();
    display->series(k).xdata.clear();
    return back_inserter(display->series(k).xdata);
  }
  back_insert_iterator< vector<double> > y_inserter(key_t k)
  { if (!display) createDisplay();
    display->series(k).ydata.clear();
    return back_inserter(display->series(k).ydata);
  }

  virtual void updateDisplay(double t)
  { if (!display) createDisplay();
    display->updateDisplay(t);
  }
  virtual double displayThreshold()
  { return 0; }
  void setTitle(key_t k, string _t)
  { if (!display) createDisplay();
    display->series(k).title = _t;
  }
  virtual string filenamebase(void)
  { string pfb = params.outfilenamebase();
    if (pfb != "") return pfb;
    return "histogram";
  }
  virtual string gpfilename(void)
  { string fnb = filenamebase();
    if (fnb.length() == 0)
      return "";
    return params.outputDirectory() + '/' + fnb + ".gp";
  }
  virtual string datafilename(void)
  { string fnb = filenamebase();
    if (fnb.length() == 0)
      return "";
    return params.outputDirectory() + '/' + fnb + ".dat";
  }
};

// template function definitions

template <class KC, class PC>
HistogramDisplay<KC,PC>
::HistogramDisplay(HistogramDisplayController<KC,PC>*con)
  : xs(*con)
{}

#endif//HISTOGRAM_DISPLAY_H
