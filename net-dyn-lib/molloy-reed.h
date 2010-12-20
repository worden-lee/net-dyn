// This file is modified from 'plod_generator.hpp' provided by the
// Boost Graph Library, which is Copyright 2004 The Trustees of
// Indiana University.

// The Boost Graph License is attached at the end of this file.

// 2008 Lee Worden

#ifndef MOLLOY_REED_GRAPH_GENERATOR_HPP
#define MOLLOY_REED_GRAPH_GENERATOR_HPP

#include <iterator>
#include <utility>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/graph/graph_traits.hpp>
#include <vector>
#include <map>
#include <cmath>

// this template is for generating directed graphs whose in-degree and
// out-degree distributions are both specified.
template<typename RandomGenerator, typename Graph>
class molloy_reed_iterator
{
  typedef std::vector<unsigned> degree_dist_t;
  typedef std::vector<std::size_t> stubs_t;
  typedef typename
    boost::graph_traits<Graph>::directed_category directed_category;

public:
  typedef std::input_iterator_tag iterator_category;
  typedef std::pair<std::size_t, std::size_t> value_type;
  typedef const value_type& reference;
  typedef const value_type* pointer;
  typedef void difference_type;

  molloy_reed_iterator() 
    : in_stubs(), out_stubs(), current_index(-1),
      allow_self_loops(false) { }

  molloy_reed_iterator(RandomGenerator& gen,
                       const degree_dist_t &in_deg, 
                       const degree_dist_t &out_deg,
                       bool allow_self_loops = false)
    : in_stubs(new stubs_t), out_stubs(new stubs_t),
      current_index(0), allow_self_loops(allow_self_loops)
  { if (in_deg.size() != out_deg.size())
      return; // got to agree on # of vertices
    unsigned n = accumulate(in_deg.begin(), in_deg.end(), 0);
    if (n != accumulate(out_deg.begin(), out_deg.end(), 0))
      return; // got to agree on # of edges

    for (std::size_t i = 0; i < in_deg.size(); ++i) {
      //cout << in_deg[i] << ' '<< out_deg[i] << endl;
      in_stubs->insert(in_stubs->end(), in_deg[i], i);
      out_stubs->insert(out_stubs->end(), out_deg[i], i);
    }
    boost::uniform_int<std::size_t> ui;
    boost::variate_generator<RandomGenerator, boost::uniform_int<std::size_t> >
      vgen(gen, ui);
    random_shuffle(out_stubs->begin(), out_stubs->end(), vgen);
  }

  static molloy_reed_iterator &end(void)
  { static molloy_reed_iterator end_iterator;
    return end_iterator;
  }

  reference operator*() const
  { if (0 <= current_index && current_index < in_stubs->size())
    { //cout << "here's an edge: " << current_edge << endl;
      return current_edge;
    }
    cerr << "bad molloy-reed iterator dereference" << endl;
    return current_edge;
  }
    
  molloy_reed_iterator& operator++()
  { next();
    return *this;
  }

  virtual bool operator==(const molloy_reed_iterator& other) const
  { return current_index == other.current_index;
  }

  bool operator!=(const molloy_reed_iterator& other) const
  { return !(*this == other); }

private:
  void next()
  { do {
      ++current_index;
    }
    while (!allow_self_loops && 
        current_index < in_stubs->size() &&
        (*in_stubs)[current_index] == (*out_stubs)[current_index]);
    if (current_index < in_stubs->size())
      current_edge = 
        make_pair((*in_stubs)[current_index], (*out_stubs)[current_index]);
    else
      current_index = -1; // end
  }

  boost::shared_ptr<stubs_t> in_stubs;
  boost::shared_ptr<stubs_t> out_stubs;
  bool allow_self_loops;
  int current_index;
  value_type current_edge;
};

// specialize the above to produce directed graphs with independent power
// laws for in and out degree
template<typename RandomGenerator, typename Graph>
class power_law_molloy_reed_iterator : 
  public molloy_reed_iterator<RandomGenerator, Graph>
{
protected:
  static vector<unsigned> 
    make_degrees(std::size_t n, double density, double exp, 
      RandomGenerator &gen)
  { std::vector<unsigned> degrees(n);
    boost::uniform_real<double> ur01(0.0, 1.0);
    std::vector<double> rel_degrees(n);
    // this wouldn't serve everyone - I want strongly connected graphs,
    // so I don't allow any degrees to be 0.
    // these are x0^(1-exp) and x1^(1-exp) where the degree is in [x0,x1]
    // see http://mathworld.wolfram.com/RandomNumber.html
    double x0x = 1, x1x = pow(n,1-exp);
    for (std::size_t i = 0; i < n; ++i)
    { do
      { double y = ur01(gen);
        double px = pow(x0x + (x1x-x0x)*y, 1/(1-exp));
        rel_degrees[i] = px;
        //cout << "pow " << (x0x + (x1x-x0x)*y) << ' ' << (1/(1-exp)) 
        //    << " = " << px << endl;
      } while ( rel_degrees[i] == 0 );
    }
    //cout << "rel_degrees :";
    //for (std::size_t i = 0; i < n; ++i)
    //  cout << ' ' << rel_degrees[i];
    //cout << endl;
    // scale all degrees to produce correct # of edges
    unsigned target = n * (n-1) * density / 2;
    // binary search for a factor that will get us the right number
    // of edges...
    double low = 0;
    double high = 2 * (target - n) /
        (accumulate(rel_degrees.begin(), rel_degrees.end(), 0.0) - n);
    unsigned edges;
    do {
      double factor = 0.5 * (high + low);
      //cout << "factor: " << factor << endl;
      for (std::size_t i = 0; i < n; ++i)
      { degrees[i] = 1 + (rel_degrees[i] - 1) * factor;
      }
      edges = accumulate(degrees.begin(), degrees.end(), 0);
      //cout << "total " << edges
      //  << " of " << target << endl;
      if ( edges < target )
        low = factor;
      else if (edges > target)
        high = factor;
    }  while (edges != target);
    //cout << "degrees :";
    //for (std::size_t i = 0; i < n; ++i)
    //  cout << ' ' << degrees[i];
    //cout << endl;
    return degrees;
  }

public:
  power_law_molloy_reed_iterator(RandomGenerator &gen, std::size_t n,
      double density, double in_exp, double out_exp,
      bool allow_self_loops = false) :
    molloy_reed_iterator<RandomGenerator, Graph>(gen,
        make_degrees(n, density, in_exp, gen),
        make_degrees(n, density, out_exp, gen), allow_self_loops)
  {
#if 0
    // short lived hack here!
    extern NetworkExperimentParameters parameters;
    string fn = parameters.outputDirectory()+"/degree-dist.csv";
    cout << "Write to " << fn << endl;
    CSVDisplay nodes_csv(fn);
    nodes_csv << "vertex" << "in degree" << "out degree";
    nodes_csv.newRow();
    for (std::size_t i = 0; i != n; ++i) {
      nodes_csv << i << (*in_degree_dist)[i] << (*out_degree_dist)[i];
      nodes_csv.newRow();
    }
#endif
  }

  power_law_molloy_reed_iterator() :
    molloy_reed_iterator<RandomGenerator, Graph>()
  {}

  static power_law_molloy_reed_iterator &end(void)
  { static power_law_molloy_reed_iterator end_iterator;
    return end_iterator;
  }
};

/*

// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

//  Authors: Douglas Gregor
//           Andrew Lumsdaine

	Boost.Graph Licence
	-------------------

COPYRIGHT NOTICE:

Copyright 1997-2000, University of Notre Dame.
Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek

The Boost Graph Library "Artistic License"

Preamble

The intent of this document is to state the conditions under which a
Package may be copied, such that the Copyright Holder maintains some
semblance of artistic control over the development of the package,
while giving the users of the package the right to use and distribute
the Package in a more-or-less free fashion, plus the right to make
reasonable modifications.

Definitions

"Package" refers to the collection of files distributed by the
Copyright Holder, and derivatives of that collection of files created
through textual modification.

"Standard Version" refers to such a Package if it has not been
modified, or has been modified in accordance with the wishes of the
Copyright Holder as specified below.

"Copyright Holder" is whoever is named in the copyright or copyrights for the package. 

"You" is you, if you're thinking about copying or distributing this Package. 

"Reasonable copying fee" is whatever you can justify on the basis of
media cost, duplication charges, time of people involved, and so
on. (You will not be required to justify it to the Copyright Holder,
but only to the computing community at large as a market that must
bear the fee.)

"Freely Available" means that no fee is charged for the item itself,
though there may be fees involved in handling the item. It also means
that recipients of the item may redistribute it under the same
conditions they received it.

1. You may make and give away verbatim copies of the source form of
the Standard Version of this Package without restriction, provided
that you duplicate all of the original copyright notices and
associated disclaimers.

2. You may apply bug fixes, portability fixes and other modifications
derived from the Public Domain or from the Copyright Holder. A Package
modified in such a way shall still be considered the Standard Version. 

3. You may otherwise modify your copy of this Package in any way,
provided that you insert a prominent notice in each changed file
stating how and when you changed that file, and provided that you do
at least ONE of the following:

  a. place your modifications in the Public Domain or otherwise make
    them Freely Available, such as by posting said modifications to Usenet
    or an equivalent medium, or placing the modifications on a major
    archive site such as uunet.uu.net, or by allowing the Copyright Holder
    to include your modifications in the Standard Version of the Package.
  b. use the modified Package only within your corporation or organization. 
  c. rename any non-standard types and functions so the names do not
    conflict with Standard Vibrary, which must also be provided, and
    provide a separate documentation for each non-standard type of function
    that clearly documents how it differs from the Standard Version.
  d. make other distribution arrangements with the Copyright Holder. 

4. You may charge a reasonable copying fee for any distribution of this
Package. You may charge any fee you choose for support of this
Package. You may not charge a fee for this Package itself. However,
you may distribute this Package in aggregate with other (possibly
commercial) programs as part of a larger (possibly commercial)
software distribution provided that you do not advertise this Package
as a product of your own. 

5. The name of the Copyright Holder may not be used to endorse or
promote products derived from this software without specific prior
written permission.

DISCLAIMER:

LICENSOR MAKES NO REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPLIED.
By way of example, but not limitation, Licensor MAKES NO
REPRESENTATIONS OR WARRANTIES OF MERCHANTABILITY OR FITNESS FOR ANY
PARTICULAR PURPOSE OR THAT THE USE OF THE LICENSED SOFTWARE COMPONENTS
OR DOCUMENTATION WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS
OR OTHER RIGHTS.

The Authors and the University of Notre Dame du Lac shall not be held
liable for any liability nor for any direct, indirect or consequential
damages with respect to any claim by LICENSEE or any third party on
account of or arising from this Agreement or use of this software.

Any disputes arising out of this Agreement or LICENSEE'S use of the
software at any time shall be resolved by the courts of the state of
Indiana.  LICENSEE hereby consents to the jurisdiction of the Indiana
courts and waives the right to challenge the jurisdiction thereof in
any dispute arising out of this Agreement or Licensee's use of the
software.

*/

#endif // MOLLOY_REED_GRAPH_GENERATOR_HPP
