
//    --------------------------------------------------------------------
//
//    This file is part of Luna.
//
//    LUNA is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    Luna is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with Luna. If not, see <http://www.gnu.org/licenses/>.
//
//    Please see LICENSE.txt for more details.
//
//    --------------------------------------------------------------------

// extension of the method of Knuth and Welford for computing standard
// deviation in one pass through the data; source: John D. Cook
// Consulting


#ifndef __LUNA_RUNNINGSTATS_H__
#define __LUNA_RUNNINGSTATS_H__

class running_stats_t
{

 public:
  running_stats_t();
  void clear();
  void push(double x);
  long long num_data_values() const;
  double mean() const;
  double variance() const;
  double standard_deviation() const;
  double skewness() const;
  double kurtosis() const;
  
  friend running_stats_t operator+(const running_stats_t a, const running_stats_t b);

  running_stats_t& operator+=(const running_stats_t &rhs);

private:
  long long n;
    double M1, M2, M3, M4;
};

#endif


