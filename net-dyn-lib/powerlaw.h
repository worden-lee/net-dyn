// This file is modified from 'plod_generator.hpp' provided by the
// Boost Graph Library, which is Copyright 2004 The Trustees of
// Indiana University.

// The Boost Graph License is attached at the end of this file.

// 2008 Lee Worden

#ifndef POWER_LAW_GRAPH_GENERATOR_HPP
#define POWER_LAW_GRAPH_GENERATOR_HPP

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
// out-degree distributions are both scale-free and with different
// exponents.
template<typename RandomGenerator, typename Graph>
class powerlaw_iterator
{
  typedef std::vector<double> degree_dist_t;
  typedef typename
    boost::graph_traits<Graph>::directed_category directed_category;

public:
  typedef std::input_iterator_tag iterator_category;
  typedef std::pair<std::size_t, std::size_t> value_type;
  typedef const value_type& reference;
  typedef const value_type* pointer;
  typedef void difference_type;

  powerlaw_iterator() 
    : gen(0), in_degree_dist(), out_degree_dist(), degrees_left(0), 
      allow_self_loops(false) { }

  powerlaw_iterator(RandomGenerator& gen, std::size_t n,  
                    double density, double in_exp, double out_exp,
                    bool allow_self_loops = false)
    : gen(&gen), n(n), in_degree_dist(new degree_dist_t),
      out_degree_dist(new degree_dist_t),
    degrees_left((size_t)(density*n*(n-1))), allow_self_loops(allow_self_loops)
  {
    using std::pow;

    if (n <= 1)
      return; // good luck with that

    boost::uniform_int<std::size_t> x(0, n-1);
    double in_acc, out_acc;
    for (std::size_t i = 0; i != n; ++i) {
      std::size_t xv;
      double rel_degree;
      //std::size_t degree = (xv == 0? 0 : std::size_t(beta * pow(xv, -alpha)));
      //if (degree != 0) {
      //  out_degrees->push_back(std::make_pair(i, degree));
      //}
      do { xv = x(gen); } while (xv == 0);
      in_acc += (rel_degree = pow(xv,-in_exp));
      in_degree_dist->push_back(rel_degree);
        
      do { xv = x(gen); } while (xv == 0);
      out_acc += (rel_degree = pow(xv,-out_exp));
      out_degree_dist->push_back(rel_degree);
    }

    for (std::size_t i = 0; i != n; ++i) {
      (*in_degree_dist)[i] /= in_acc;
      (*out_degree_dist)[i] /= out_acc;
    }
    next();
  }

  //pointer operator->() const { return &current; }
  reference operator*() const
  { return current;
  }
    
  powerlaw_iterator& operator++()
  { next();
    --degrees_left;
    return *this;
  }

  bool operator==(const powerlaw_iterator& other) const
  { 
    return degrees_left == other.degrees_left; 
  }

  bool operator!=(const powerlaw_iterator& other) const
  { return !(*this == other); }

private:
  void next()
  { boost::uniform_real<double> ur01(0.0,1.0);
    double in_p = ur01(*gen), out_p = ur01(*gen);
    for (int i = 0; i < n; ++i)
    { if (in_p > 0)
      { in_p -= (*in_degree_dist)[i];
        if (in_p <= 0)
          current.first = i;
      }
      if (out_p > 0)
      { out_p -= (*out_degree_dist)[i];
        if (out_p <= 0)
          current.second = i;
      }
      else if (in_p < 0)
        break;
    }
    // disqualify loop
    if (!allow_self_loops && current.second == current.first)
      next();
  }

  RandomGenerator* gen;
  std::size_t n;
  boost::shared_ptr<degree_dist_t> out_degree_dist;
  boost::shared_ptr<degree_dist_t> in_degree_dist;
  std::size_t degrees_left;
  bool allow_self_loops;
  value_type current;
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

#endif // POWER_LAW_GRAPH_GENERATOR_HPP
