/* -*- C++ -*- */
// This is a parent class for Displays that plot time series
//  in Gnuplot.
// Either curves, or branching lineages,
// either horizontal time or vertical time
#ifndef TIMESERIESDISPLAY_H
#define TIMESERIESDISPLAY_H
#include "DisplayController.h"
#include "GnuplotDisplay.h"
#include <fstream>
#include <deque>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <string>
#include <sstream>
#include "util.h"
#include <unistd.h>

//template<class DISPLAY>
template<class K, class ParamsClass>
class TimeSeriesController;

class point_tx_t
{
public:
  double t, x;
  point_tx_t(double _t, double _x)
    : t(_t), x(_x) {}
};
double get_t(const point_tx_t&p);
double get_x(const point_tx_t&p);

typedef deque<point_tx_t> series_t;
typedef vector<double> stoptimes_t;
class timeseries_t
{
public:
  series_t series;
  stoptimes_t stoptimes;
  string basename;
  string title;
  int color;
  bool axis2;
  bool plot_ok;
  bool initialized;
  timeseries_t() : color(-1), axis2(false), plot_ok(false), initialized(false)
  {}
};

template<class K, class ParamsClass=Parameters, class lineage_t=K>
class TimeSeriesDisplay : public GnuplotDisplay<ParamsClass>
{
public:
  using GnuplotDisplay<ParamsClass>::gnuplot;
  typedef K key_t;
  
protected:
  TimeSeriesController<key_t,ParamsClass> &xs;
  // if datafile is used, it records the time series to filesystem
  ofstream datafile;
  // store the timeseries here
  typedef map<K,timeseries_t> series_map_t;
  series_map_t series;
  // aux structures for displaying timeseries sorted by lineage
  //  use of these is optional: if no lineages are ever set explicitly
  //  they won't be used.
  typedef set< lineage_t > lineages_t;
  lineages_t lineages;
  typedef map< lineage_t, set<K> > series_by_lineage_t;
  series_by_lineage_t series_by_lineage;
  typedef map< K, lineage_t > lineage_by_series_t;
  lineage_by_series_t lineage_by_series;
  // display attributes
  bool vertical;
  bool lines;
  bool y2rangeset;
  string plottitle;
  
public:
  TimeSeriesDisplay(TimeSeriesController<key_t,ParamsClass>*con);

  virtual ~TimeSeriesDisplay() {}

  void initialize(void);

  // this seems to be needed for compilation
  // but should not be used
  // should be overridden by subclasses
  void updateDisplay(double t)
  { cerr <<"TimeSeriesDisplay: UPDATEDISPLAY ERROR\n"; }

  // this will assume points before t0 are no longer to be stored
  void updateDisplay(double t0, double t1);
  virtual void writeGnuplotHeaders(void);

  virtual string gnuplotLogFile(void);

  virtual void descendsFrom(key_t daughter, key_t parent,
			    double birthTime) {}
  // descendsFrom can call this, if desired, but doesn't by default
  virtual void connectLineages(key_t daughter, key_t parent,
			       double birthTime);

  // this is called by connectLineages but must be called directly too
  //  if you want display to be ordered by lineage
  virtual void setLineage(key_t series, lineage_t lineage);
  
  // this returns false if k was not previously known
  virtual bool record(key_t k, double t, double x);
  // this t should be one also used in record()
  virtual void deactivate(key_t k, double t);

  void setVertical(bool v)
  { vertical = v; }
  void setLines(bool l)
  { lines = l; }
  void setColor(key_t k, int c)
  { series[k].color = c; }
  void setBasename(key_t k, string s)
  { series[k].basename = s; }
  void setTitle(key_t k, string s)
  { series[k].title = s; }
  void setAxis2(key_t k, bool a)
  { series[k].axis2 = a; }
  int color(K k)
  { return series[k].color; } 
protected:
  void announceSeries(ostream&os, timeseries_t &ts, double t);
  virtual string title(timeseries_t &ts)
  { return ts.title; }
  virtual lineage_t defaultLineage(K k);
  void removeKeyFromMemory(K i);
};

template<class K, class ParamsClass>
class TimeSeriesController :
  public DisplayController< TimeSeriesDisplay<K,ParamsClass>, ParamsClass >
{
  friend class TimeSeriesDisplay<K,ParamsClass>;
public:
  typedef K key_t;
  typedef DisplayController< TimeSeriesDisplay<K,ParamsClass>, ParamsClass >
    superclass;
protected:
  // newer g++'s are very picky about these things
  using superclass::display;
  using ObjectWithParameters<ParamsClass>::params;

public:

  virtual ~TimeSeriesController()
  { //if (display) delete display;
  }

  // subclass should provide this
  virtual void recordFile(double t) = 0;

  // and optionally, some of these
  virtual void updateDisplay(double t)
  {
    double window = params.scrollingWindow();
    double t0 = (window < t ? t - window : 0);
    display->updateDisplay(t0,t);
  }
  virtual bool vertical()
  { return false; }
  virtual double displayThreshold()
  { return 0; }
  virtual string title(void)
  { return ""; } 
  virtual string filenamebase(void)
  { return params.outfilenamebase(); }
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
  virtual int color(key_t k)
  { return -1; } 

  // things to do when a series is introduced
  virtual void newSeries(key_t k)
  { display->setColor(k, color(k)); }
};

// template function definitions

template <class K, class PC, class lineage_t>
TimeSeriesDisplay<K, PC, lineage_t>
::TimeSeriesDisplay(TimeSeriesController<K,PC>*con)
  : xs(*con),
    vertical(xs.vertical()), lines(true), y2rangeset(false)
{}

template<class K, class PC, class lineage_t>
void TimeSeriesDisplay<K, PC, lineage_t>::initialize(void)
{
  GnuplotDisplay<PC>::initialize();
  string _datafilename = xs.datafilename();
  if (_datafilename.length() > 0)
    datafile.open(_datafilename.c_str());
  plottitle = xs.title();
//   if (_datafilename.length() > 0)
//     gnuplot << "set title \"" << _datafilename << "\"\n";
}

template <class K, class PC, class lineage_t>
string TimeSeriesDisplay<K, PC, lineage_t>::gnuplotLogFile()
{
  return xs.gpfilename();
}

template <class K, class PC, class lineage_t>
void TimeSeriesDisplay<K, PC, lineage_t>::
  connectLineages(K daughter, K parent, double birth_t)
{ if (lineage_by_series.find(parent) != lineage_by_series.end())
    setLineage(daughter,lineage_by_series[parent]);
 
  // search the parent series for the last value no later than birth_t
  timeseries_t &pts = series[parent];
  typename series_t::iterator birth_pt = pts.series.end();
  for (typename series_t::iterator pi=pts.series.begin();
       pi != pts.series.end(); ++pi)
  {
    if (pi->t > birth_t) break;
    birth_pt = pi;
  }
  // that's the first point in the daughter series.
//  series[daughter].series.push_front(*birth_pt);
  if (birth_pt != pts.series.end())
    record(daughter,birth_pt->t,birth_pt->x);
}

template <class K, class PC, class lineage_t>
void TimeSeriesDisplay<K, PC, lineage_t>::
  setLineage(K series, lineage_t lineage)
{
  lineages.insert(lineage);
  series_by_lineage[lineage].insert(series);
  lineage_by_series[series] = lineage;
}

// if K isn't compatible with lineage_t you have to override this
template <class K, class PC, class lineage_t>
lineage_t TimeSeriesDisplay<K, PC, lineage_t>::defaultLineage(K k)
{ return lineage_t(k);
}

// copied over from util.h
template<class A, class B>
ostream &operator<<(ostream&o, const pair<A,B>&p)
{ return o << '(' << p.first << ',' << p.second << ')';
}

template <class K, class PC, class lineage_t>
bool TimeSeriesDisplay<K, PC, lineage_t>::record(K k, double t, double x)
{
//   typename series_map_t::iterator fk = series.find(k);
//   bool have = (fk != series.end()); // was it already there
//   timeseries_t &ts = (have ? fk->second : series[k]);
  timeseries_t &ts = series[k];
  bool have = ts.initialized;
  if (!have)
  { xs.newSeries(k); // call back the controller
    // if it doesn't have an explicitly set lineage, use the key
    if (lineage_by_series.find(k) == lineage_by_series.end())
      setLineage(k,defaultLineage(k));
    ts.initialized = true;
  }
  unsigned tss = ts.series.size();
  // if it's already there, fine
  if (tss >= 1 &&
      ts.series[tss-1].t == t &&
      ts.series[tss-1].x == x)
    return have;
  // optimize: instead of recording the same x three
  // times in a row, replace the second one 
  else if (tss >= 2 &&
      ts.series[tss-1].x == x &&
      ts.series[tss-2].x == x &&
      !(ts.stoptimes.size() > 0 &&
	ts.stoptimes.back() >= ts.series[tss-2].t))
    ts.series[tss-1].t = t;
  else
    ts.series.push_back(point_tx_t(t,x));
  if (datafile.is_open())
    if (ts.basename != "")
      datafile << ts.basename << ' ' << t << ' ' << x << endl;
    else
      datafile << k << ' ' << t << ' ' << x << endl;
  return have;
}

template <class K, class PC, class lineage_t>
void TimeSeriesDisplay<K, PC, lineage_t>::deactivate(K k, double t)
{
  typename series_map_t::iterator fk = series.find(k);
  if (fk!=series.end())
  { // deactivate if not already deactivated
    timeseries_t &ts = fk->second;
    if (!ts.series.empty() &&
	(ts.stoptimes.empty() ||
	 ts.series.back().t > ts.stoptimes.back()))
      ts.stoptimes.push_back(t);
  }
}

template <class K, class PC, class lineage_t>
void TimeSeriesDisplay<K, PC, lineage_t>::
  announceSeries(ostream& os, timeseries_t &ts, double t)
{
  os << " '-'";
  if (ts.axis2 && !vertical) // deal with vertical axis2 later
    os << " axes x1y2";
  string _title = title(ts); // I. Title.
  if (_title.length() > 0)
    os << " title '" << _title << '\'';
  else
    os << " notitle";
  os << " w l " << ts.color;
}

template <class K, class PC, class lineage_t>
void TimeSeriesDisplay<K, PC, lineage_t>::removeKeyFromMemory(K i)
{ lineage_t lin = lineage_by_series[i];
  lineage_by_series.erase(i);
  series_by_lineage[lin].erase(i);
  if (series_by_lineage[lin].empty())
  { series_by_lineage.erase(lin);
    lineages.erase(lin);
  }
}

// Write header info for Gnuplot
// Subclasses should overwrite with their own needs
// (e.g., y2tics, parametric, nokey)
template <class K, class PC, class lineage_t>
void TimeSeriesDisplay<K, PC, lineage_t>::writeGnuplotHeaders()
{
  gnuplot << "set title \"" << plottitle << "\"\n";
}

template <class K, class PC, class lineage_t>
void TimeSeriesDisplay<K, PC, lineage_t>::updateDisplay(double t0, double t1)
{
  double recordEvery = xs.recordEvery();

  GnuplotDisplay<PC>::rewindStream();
  writeGnuplotHeaders();

  // first, figure out what times to plot, to get all the 
  //  line segments that intersect [t0, t1]

  // throw away obsolete data
  if (t0 > 0)
  {
    // collect all times found in the time series into `tset'
    set<double> tset;
    for (typename series_map_t::iterator sri=series.begin();
	 sri != series.end(); ++sri)
    {
      series_t &srp = sri->second.series;
      transform(srp.begin(),srp.end(),inserter(tset,tset.begin()),&get_t);
    }
  
    // so far, all times encountered <= firstTime

    // firstTime is the latest time before the beginning of the window

    // firstTime is the earliest time that should be plotted

    // for each time series, weed out all points no longer being plotted
    //  and decide whether it's still worth keeping
    for (typename series_map_t::iterator sri=series.begin();
	 sri != series.end();)
    {
      timeseries_t &ts = sri->second;
      ts.plot_ok = false;

      for (unsigned int i = 0; i < ts.series.size(); )
      {
	double lastt = -HUGE;
	// erase if no longer needed
	if (ts.series[i].t < t0 &&
	    (i+1 == ts.series.size() || ts.series[i+1].t < t0))
	{
//  	  cout << "erase " << ts.series[i].t << " ..";
//  	  if (i+1 == ts.series.size())
//  	    cout << "(*)";
//  	  else cout << ts.series[i+1].t;
//  	  cout << " (" << t0 << "--" << t1 << ") ... ";
//  	  cout << i << "th of " << ts.series.size();
	  ts.series.erase(ts.series.begin()+i); // in place of ++i
//  	  cout << ".." << ts.series.size() << '\n';
	}
	// erase if it's too close to the last one
	//  (the spacing `recordEvery' grows as time passes)
 	else if (ts.series[i].t - lastt < recordEvery/2 /* /10 */)
	{
//  	  cout << "weed (" << lastt << "..) " << ts.series[i].t << "\n";
 	  ts.series.erase(ts.series.begin()+i);
	}
	// otherwise we keep it
	else
	{ // if it's big enough, we mark the whole record worth keeping
	  if (ts.series[i].x > xs.displayThreshold())
	    ts.plot_ok = true;
	  lastt = ts.series[i].t;
	  ++i;
	}
      }
//       cout << "series " << ts.basename << " (" << sri->first
// 	   << ") size is "
// 	   << ts.series.size() << "\n";

      // now ts->plot_ok is true iff there is a point in the series
      //  above the threshold for displaying
      typename series_map_t::iterator srx = sri++;
      // get rid of the time series if it's not worth plotting
      if (!ts.plot_ok)
      { removeKeyFromMemory(srx->first);
	series.erase(srx);
      }
    }
//     cout << '\n';
  }

  // now series contains only time series that contain a point above the
  // threshold and within the range of times to be plotted
  // special case: if just one series on data axis 2, see if it's at a
  // steady state and do a special range command to keep gnuplot from
  // complaining
  string y2range = "";
  {
    set<timeseries_t *> tss;
    for (typename series_map_t::iterator sri=series.begin();
	 sri != series.end(); ++sri)
      if (sri->second.axis2)
	tss.insert(&sri->second);
    // if just one series
    if (tss.size()==1)
    {
      timeseries_t *t2 = *tss.begin();
      set<double> x2set;
      series_t &srp = t2->series;

      // put all the x's in the set
      transform(srp.begin(),srp.end(),
		inserter(x2set,x2set.begin()),&get_x);
      set<double>::iterator
	maxx = max_element(x2set.begin(),x2set.end()),
	minx = min_element(x2set.begin(),x2set.end());
      if (maxx==x2set.end() || minx==x2set.end())
      { cerr<<"TimeSeriesDisplay: ERROR! No max or min!"<<endl;
        return;
      }
      if (*maxx-*minx<0.001)
	y2range = stringf(" [%g:%g]", *minx-0.001,*minx+0.001);
    }
  }

  // removed stuff to check for y2range having been previously set,
  // because now it will be a different file anyway
  //   that's true but it's still necessary to unset it in the running
  //   gnuplot program
  if (y2range != "")
  { gnuplot << "set y2range" << y2range << "\n";
    y2rangeset = true;
  }
  // this if() here seems to break the displaying
  else// if (y2rangeset)
  { gnuplot << "set y2range [*:*]\n";
  //cout << "set y2range [*:*]\n";
    y2rangeset = false;
  }

  // Now: loop through all timeseries, write a plot command and series of
  //  plot points.
  ostringstream datastream;
  
  gnuplot << "plot";
  if (t0 > 0)
    gnuplot << " [" << t0 << ":" << t1 << "]";

  bool comma=false;
  // do lineages in order, series within each lineage
  for (typename lineages_t::iterator li=lineages.begin();
       li != lineages.end(); ++li)
    for (typename set<K>::iterator sri=series_by_lineage[*li].begin();
	 sri != series_by_lineage[*li].end(); ++sri)
  {
    timeseries_t &ts = series[*sri];
    string title;
    if (ts.basename.length() > 0)
      title = ts.basename;
    else if (ts.title.length() > 0)
      title = ts.title;
    datastream << "# " << title << "\n";
    if (ts.series.size() == 0)
      continue;
  
    // for each time series to be plotted
    //  announce it in the plot command and buffer up its points

    stoptimes_t::iterator si = ts.stoptimes.begin(),
      se = ts.stoptimes.end();
    bool last_was_blank = false;
    double st = (si != se ? *si : HUGE);
    for (typename series_t::iterator pi=ts.series.begin();
	 pi != ts.series.end(); ++pi)
    {
      while (st < pi->t)
      { 
	datastream << "e\n";
	st = (++si != se ? *si : HUGE);

	gnuplot << (comma ? "," : ((comma=true),""));
	announceSeries(gnuplot,ts,t1);
      }
      // possible gnuplot bug workaround
      double pt = pi->t, px=pi->x;
      if (pt < t0)
      { px = px + (pi[1].x - px) * (t0 - pt) / (pi[1].t - pt);
        pt = t0;
      }
      if (isnan(pi->x) || isnan(pt))
      { if (!last_was_blank)
          datastream << '\n';
        last_was_blank = true;
      }
      else
      { if (vertical)
          datastream << pi->x << ' ' << pt << '\n';
        else
          datastream << pt << ' ' << pi->x << '\n';
        last_was_blank = false;
      }
    }
    while (st != HUGE)
    { 
      datastream << "e\n";
      st = (++si != se ? *si : HUGE);

      gnuplot << (comma ? "," : ((comma=true),""));
      announceSeries(gnuplot,ts,t1);
    }
    if (ts.stoptimes.empty()
	|| ts.stoptimes.back() < ts.series.back().t)
    { datastream << "e\n";

      gnuplot << (comma ? "," : ((comma=true),""));
      announceSeries(gnuplot,ts,t1);
    }
  }
  gnuplot << '\n';

  // plot command is done, follow it with all the points.
  gnuplot << datastream.str();
  GnuplotDisplay<PC>::plotBufferedData();
}

#endif//TIMESERIESDISPLAY_H
