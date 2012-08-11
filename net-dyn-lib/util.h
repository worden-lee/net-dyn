/*
 * util stuff, some from Numerical Recipes in C
 */
#ifndef __UTIL_H__
#define __UTIL_H__

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <errno.h>
#include <stdlib.h>
#include <limits>
//typedef _G_va_list va_list;
#include <stdarg.h>
#include <string.h>
#include <string>
using namespace std;

/*
 * a little of header file from Numerical Recipes
 */

inline double DMAX(double a, double b) 
{ return(a > b ? a : b);
}

inline double DMIN(double a, double b) 
{ return(a < b ? a : b);
}

inline long LMAX(long a, long b)
{ return(a > b ? a : b);
}

#define SIGN(x,y) ((y) >= 0.0 ? fabs(x) : -fabs(x))

//#define nrerror perror
inline void nrerror(const char *string)
{
  if (errno != 0)
    perror(string);
  else
  {
    cerr << string << endl; cerr.flush();
  }
}

inline int fsign(double x) {
  if (x==0) return 0;
  else if (x>0) return 1;
  else return -1;
}

// estimate accuracy of a Bernoulli sample
// N = var / (accuracy * p)^2
// or accuracy = sqrt(var/N) / p
// for example, if p = 0.25, accuracy = 0.1, we find N = 300
// (sample mean = p, var = p(1-p))
inline double Bernoulli_accuracy(double p, unsigned N)
{ if (N == 0)
    return std::numeric_limits<double>::infinity();
  if (p <= 0 || p >= 1) // this is a problem case, but is this the solution?
    return sqrt(1.0/N);
  else
    return sqrt(p*(1-p)/N)/p;
}

// conservative estimate of how many Bernoulli samples to take for the
// specified accuracy (p unknown)
// the maximum variance is 0.25, so we use that
inline unsigned Bernoulli_estimate_n(double accuracy)
{ return unsigned(1/(accuracy*accuracy));//0.25 / (accuracy * 0.5)^2
}

// the % operator returns negative when the moduland is negative.
// this always returns a number between zero and modulus, not
// including modulus
template<typename int_type, typename int_type_2>
inline int_type mod_correct(int_type moduland, int_type_2 modulus)
{ int_type pc = moduland % modulus;
  if (modulus > 0 && pc < 0)
    pc += modulus;
  else if (modulus < 0 && pc > 0)
    pc += modulus;
  return pc;
}

// in case of PRESERVE_SIGNS, 'sign' tells what sign of v should be
// otherwise 'sign' isn't used
// either way make sure its modulus is at least that of m
inline double clamp(double v, double m, int sign, bool preserveSigns=false) 
{
  m = fabs(m);
  if (preserveSigns)
    switch ( sign )
    {
    case 0:
      return 0;
    case -1:
      return DMIN(v, -m);
    case 1:
      return DMAX(v, m);
    }
  else 
    switch (fsign(v))
    {
    case 0:
      return m;
    case -1:
      return DMIN(v, -m);
    case 1:
      return DMAX(v, m);
    }
  cerr << "Error in clamp!\n";
  return -HUGE;
}

#include <glob.h>
const char *glob_string(const char *string, const char *sep);

// for printf style formatting
inline const char *vfstring(const char *fmt, va_list va)
{
  static char buffer[1024];
  vsnprintf(buffer, sizeof buffer, fmt, va);
  return buffer;
}

inline const char *fstring(const char*fmt,...)
     __attribute__ ((format (printf, 1, 2)));
inline const char *fstring(const char*fmt,...)
{
  va_list va;
  va_start(va,fmt);
  const char *s = vfstring(fmt, va);
  va_end(va);
  return s;
}

#include "string"
inline string stringf(const char*fmt,...)
{
  va_list va;
  va_start(va,fmt);
  const char *s = vfstring(fmt, va);
  va_end(va);
  return string(s);
}

inline int string_to_int(const string &x)
{ return atoi(x.c_str()); }
inline string int_to_string(int x)
{ return string(fstring("%d",x)); }

inline unsigned string_to_unsigned(const string &x)
{ return (unsigned)atol(x.c_str()); }
inline string unsigned_to_string(unsigned x)
{ return string(fstring("%u",x)); }

inline long string_to_long(const string &x)
{ return atol(x.c_str()); }
inline string long_to_string(long x)
{ return string(fstring("%ld",x)); }

inline double string_to_double(const string &x)
{ return atof(x.c_str()); }
inline string double_to_string(double x)
{ return string(fstring("%g",x)); }

inline bool string_to_bool(const string &x)
{ return x == "true" || x == "TRUE"; } 
inline string bool_to_string(bool x)
{ return string(x?"true":"false"); }

inline const string &string_to_string(const string &x)
{ return x; }

// an extension of the std::accumulate idea
// accumulates a transform of the values in [first,last) using operator+,
// starting with the value init.
// i.e. returns init + f(*first) + f(*(first+1)) + ...
template<typename _Iterator, typename _Rtp, typename _Fn>
_Rtp accumulate_transformed(_Iterator first, _Iterator end, _Rtp init, _Fn f)
{
  for (; first != end; ++first)
    init = init + f(*first);
  return init;
}

// this is like accumulate_transformed but uses [] instead of ()
// i.e. returns init + map[*first] + map[*(first+1)] + ...
template<typename _Iterator, typename _Rtp, typename _Map>
_Rtp accumulate_mapped(_Iterator first, _Iterator end, _Rtp init, _Map map)
{
  for (; first != end; ++first)
    init = init + map[*first];
  return init;
}

// template<class A, class B>
// ostream &operator<<(ostream&o, pair<A,B>&p)
// { return o << '(' << p.first << ',' << p.second << ')';
// }

// ensure path can be created, by creating any parent directories needed.
// path itself is not created.
#include <sys/stat.h>
inline int mkpath(const char *path)
{ char buf[strlen(path)+1];
  const char *it = path;
  char *bit = buf;
  int ret;
  while(*it)
  { if (*it == '/')
    { struct stat st;
      *bit=0;
      ret = 0;
      if (stat(buf, &st) && errno == ENOENT)
        ret = mkdir(buf,0755);
      if (ret)
        return ret;
    }
    *bit++ = *it++;
  }
  return 0;
}

// this thing is like a set<> that only needs to remember its min and max
template<typename T>
class range
{
public:
  T min;
  T max;
  //typedef void *iterator;
  range() : min(HUGE), max(-HUGE) {}
  void clear()
  { min = HUGE;
    max = -HUGE;
  }
  void insert(const T &_x)
  { if (_x < min) min = _x;
    if (_x > max) max = _x;
  }
};

// template<typename T>
// class std::insert_iterator< range<T> >
// {
// protected:
//   range<T> &r;
// public:
//   insert_iterator(range<T> &_r, typename range<T>::iterator it =0)
//     : r(_r) {} 
//   insert_iterator &operator*()
//   { return *this; }
//   insert_iterator &operator=(const T &_x)
//   { r.insert(_x);
//     return *this;
//   }
//   insert_iterator &operator++()
//   { return *this; }
//   insert_iterator &operator+(int i)
//   { return *this; }
// };
// template <typename T, class Iterator>
// std::insert_iterator< range<T> > inserter (range<T>&r, range<T>::iterator i)
// {
//   return std::insert_iterator< range<T> >(r);
// }

#endif //__UTIL_H__
