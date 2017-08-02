#ifndef STANDARD_LIBRARY_EXTENSIONS_HPP
#define STANDARD_LIBRARY_EXTENSIONS_HPP

#include <numeric>
#include <algorithm>
#include <vector>

namespace mellotron {

/*!
 *  \file   StandardLibraryExtensions.hpp
 *  \author Joey Dumont      <joey.dumont@gmail.com>
 *  \since  2017-07-28
 *  \brief  Declares and defines extensions to the standard library functions.
 */

/// Finds the indices that that sorts a given 1D array. This should work for both
/// a std::vector and a boost::multi_array with one dimension.
/// Modified from
/// https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
template <class Array1D>
std::vector<size_t> sort_indices(Array1D &v)
{
  // Initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
            [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

/// Takes a 1D array and rearranges the elements according to a given
/// array of indices.
template <class Array1D>
void rearrange(Array1D &v, std::vector<size_t> &idx)
{
  // Deep copy of the array with the copy constructor.
  auto v_copy(v);

  // We rearrange the elements
  for (uint i=0; i<v.size(); i++)
  {
    v[i] = v_copy[idx[i]];
  }
}

} // namespace mellotron

#endif // STANDARD_LIBRARY_EXTENSIONS_HPP