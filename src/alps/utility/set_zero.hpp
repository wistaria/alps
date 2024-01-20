/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 1994-2010 by Matthias Troyer <troyer@comp-phys.org>,
*                            Beat Ammon <ammon@ginnan.issp.u-tokyo.ac.jp>,
*                            Andreas Laeuchli <laeuchli@comp-phys.org>,
*                            Synge Todo <wistaria@comp-phys.org>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

/* $Id: set_zero.hpp 7842 2017-07-26 14:27:58Z domischi $ */

#ifndef ALPS_UTILITY_SET_ZERO_HPP
#define ALPS_UTILITY_SET_ZERO_HPP

#include <alps/type_traits/is_sequence.hpp>
#include <alps/type_traits/element_type.hpp>

#include <boost/utility/enable_if.hpp>

#include <algorithm>

namespace alps {

template <class X> 
inline typename boost::disable_if<is_sequence<X>,void>::type
set_zero(X& x) { x=X();}

template <class X> 
inline typename boost::enable_if<is_sequence<X>,void>::type
set_zero(X& x) 
{
  std::fill(x.begin(),x.end(),typename element_type<X>::type());
}

//Enable set_zero for valarrays
template <class X> 
inline typename boost::enable_if<is_sequence<std::valarray<X> >,void>::type
set_zero(std::valarray<X> & x) 
{
  for(int i=0;i<x.size();++i) x[i]=X();
}


} // end namespace alps

#endif // ALPS_UTILITY_RESIZE_HPP
