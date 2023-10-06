
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

#include "running-stats.h"
#include <cmath>
#include <vector>

running_stats_t::running_stats_t() 
{
  clear();
}

void running_stats_t::clear()
{
  n = 0;
  M1 = M2 = M3 = M4 = 0.0;
}

void running_stats_t::push(double x)
{
  double delta, delta_n, delta_n2, term1;
  
  long long n1 = n;
  n++;
  delta = x - M1;
  delta_n = delta / n;
  delta_n2 = delta_n * delta_n;
  term1 = delta * delta_n * n1;
  M1 += delta_n;
  M4 += term1 * delta_n2 * (n*n - 3*n + 3) + 6 * delta_n2 * M2 - 4 * delta_n * M3;
  M3 += term1 * delta_n * (n - 2) - 3 * delta_n * M2;
  M2 += term1;
}

long long running_stats_t::num_data_values() const
{
  return n;
}

double running_stats_t::mean() const
{
  return M1;
}

double running_stats_t::variance() const
{
  return M2/(n-1.0);
}

double running_stats_t::standard_deviation() const
{
  return sqrt( variance() );
}

double running_stats_t::skewness() const
{
  return sqrt(double(n)) * M3/ pow(M2, 1.5);
}

double running_stats_t::kurtosis() const
{
  return double(n)*M4 / (M2*M2) - 3.0;
}

running_stats_t operator+(const running_stats_t a, const running_stats_t b)
{
  running_stats_t combined;
  
  combined.n = a.n + b.n;
  
  double delta = b.M1 - a.M1;
  double delta2 = delta*delta;
  double delta3 = delta*delta2;
  double delta4 = delta2*delta2;
  
  combined.M1 = (a.n*a.M1 + b.n*b.M1) / combined.n;
  
  combined.M2 = a.M2 + b.M2 + 
    delta2 * a.n * b.n / combined.n;
  
  combined.M3 = a.M3 + b.M3 + 
    delta3 * a.n * b.n * (a.n - b.n)/(combined.n*combined.n);
  combined.M3 += 3.0*delta * (a.n*b.M2 - b.n*a.M2) / combined.n;
  
  combined.M4 = a.M4 + b.M4 + delta4*a.n*b.n * (a.n*a.n - a.n*b.n + b.n*b.n) / 
    (combined.n*combined.n*combined.n);
  combined.M4 += 6.0*delta2 * (a.n*a.n*b.M2 + b.n*b.n*a.M2)/(combined.n*combined.n) + 
    4.0*delta*(a.n*b.M3 - b.n*a.M3) / combined.n;
  
  return combined;
}

running_stats_t& running_stats_t::operator+=(const running_stats_t& rhs)
{ 
  running_stats_t combined = *this + rhs;
  *this = combined;
  return *this;
}
