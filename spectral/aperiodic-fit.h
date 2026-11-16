
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

#ifndef LUNA_APERIODIC_FIT_H
#define LUNA_APERIODIC_FIT_H

// Self-contained power-law fit; no EDF, writer, or global band dependencies.
#include <algorithm>
#include <cmath>
#include <vector>

namespace psd_aperiodic {

inline double median( std::vector<double> v )
{
  std::sort( v.begin(), v.end() );
  const size_t n = v.size();
  return n % 2 ? v[n/2] : (v[n/2-1] + v[n/2]) / 2;
}

struct fit_t {
  bool valid = false;
  bool fallback = false;
  const char * failure = "fit not evaluated";
  int input_bins = 0;
  int n = 0;
  double offset = 0, exponent = 0, rsq = 0;
  std::vector<double> background;
};

// Fill only invalid positions. No missing-value sentinel is required, even in
// temporary tracks. Interior gaps interpolate in output units; edges extend
// the nearest observed value. An entirely invalid track uses a finite constant.
inline void fill( std::vector<double> & data, const std::vector<double> & valid,
                  double empty_value )
{
  size_t first=0;
  while (first<data.size() && !valid[first]) ++first;
  if (first==data.size()) { std::fill(data.begin(),data.end(),empty_value); return; }
  std::fill(data.begin(),data.begin()+first,data[first]);
  size_t left=first;
  while (left+1<data.size())
    {
      size_t right=left+1;
      while (right<data.size() && !valid[right]) ++right;
      if (right==data.size())
        { std::fill(data.begin()+left+1,data.end(),data[left]); break; }
      for (size_t i=left+1; i<right; ++i)
        {
          const double weight=double(i-left)/(right-left);
          data[i]=(1-weight)*data[left]+weight*data[right];
        }
      left=right;
    }
}

inline fit_t complete_fit( const std::vector<double> & f, const std::vector<double> & x,
                          const std::vector<double> & y, const std::vector<double> & keep,
                          double a, double b, double lwr, double upr, double pivot = 1 )
{
  fit_t result;
  result.input_bins=x.size();
  double mean=0, ss=0, error=0;
  for (size_t i=0; i<x.size(); ++i)
    if (keep[i]) { ++result.n; mean+=y[i]; }
  mean/=result.n;
  for (size_t i=0; i<x.size(); ++i)
    if (keep[i])
      { ss+=(y[i]-mean)*(y[i]-mean); error+=std::pow(y[i]-a-b*x[i],2); }
  // Offset is the fitted log10 power at the pivot frequency (default 1 Hz,
  // i.e. the log10(f)=0 intercept); the fitted line itself is unchanged.
  result.offset=(a+b*std::log10(pivot))/10;
  result.exponent=-b/10;
  result.rsq=ss>1e-12 ? 1-error/ss : (error<1e-12 ? 1 : 0);
  result.background.assign(f.size(),0);
  for (size_t i=0; i<f.size(); ++i)
    if (f[i]>=lwr && f[i]<=upr)
      {
        result.background[i]=std::pow(10.0,(a+b*std::log10(f[i]))/10);
        if (!std::isfinite(result.background[i]) || result.background[i]<=0)
          { result.failure="non-positive or non-finite fitted background"; return result; }
      }
  result.valid=true;
  result.failure="";
  return result;
}

inline bool line( const std::vector<double> & x, const std::vector<double> & y,
                  const std::vector<double> & w, double & a, double & b )
{
  double sw = 0, sx = 0, sy = 0;
  for (size_t i=0; i<x.size(); ++i)
    { sw += w[i]; sx += w[i]*x[i]; sy += w[i]*y[i]; }
  if ( sw <= 0 ) return false;
  const double mx = sx/sw, my = sy/sw;
  double xx = 0, xy = 0;
  for (size_t i=0; i<x.size(); ++i)
    { xx += w[i]*(x[i]-mx)*(x[i]-mx); xy += w[i]*(x[i]-mx)*(y[i]-my); }
  if ( xx <= 1e-12 ) return false;
  b = xy/xx;
  a = my-b*mx;
  return std::isfinite(a) && std::isfinite(b);
}

// Huber IRLS limits influence of both peaks and unusually low bins. A second
// stage explicitly rejects positive peaks. The 1 dB rejection floor prevents
// numerical noise in an almost exact power law from being classified as peaks.
inline fit_t fit( const std::vector<double> & f, const std::vector<double> & p,
                  double lwr, double upr, double threshold = 3, int iterations = 5,
                  double pivot = 1 )
{
  fit_t result;
  if ( f.size() != p.size() || !std::isfinite(lwr) || !std::isfinite(upr)
       || lwr <= 0 || upr <= lwr || !std::isfinite(threshold)
       || threshold <= 0 || iterations < 1 || !std::isfinite(pivot) || pivot <= 0 )
    { result.failure = "invalid fitting arguments"; return result; }
  std::vector<double> x, y;
  for (size_t i=0; i<f.size(); ++i)
    if ( std::isfinite(f[i]) && f[i] >= lwr && f[i] <= upr
         && std::isfinite(p[i]) && p[i] > 0 )
      { x.push_back(std::log10(f[i])); y.push_back(10*std::log10(p[i])); }
  result.input_bins = x.size();
  if ( x.size() < 5 )
    { result.failure = "fewer than five positive finite frequency bins"; return result; }
  std::vector<double> keep(x.size(), 1), weights = keep, residual(x.size());
  fit_t initial;
  const auto failed = [&](const char * reason) {
    fit_t alternative=initial.valid ? initial : result;
    alternative.fallback=initial.valid;
    alternative.failure=reason;
    return alternative;
  };
  double a = 0, b = 0;
  for (int pass=0; pass<=iterations; ++pass)
    {
      weights = keep;
      if ( !line(x,y,weights,a,b) )
        return failed("degenerate regression on retained bins");
      double location = 0, scale = 0;
      for (int irls=0; irls<30; ++irls)
        {
          std::vector<double> selected;
          for (size_t i=0; i<x.size(); ++i)
            { residual[i] = y[i]-a-b*x[i]; if (keep[i]) selected.push_back(residual[i]); }
          location = median(selected);
          for (double & v : selected) v = std::fabs(v-location);
          scale = 1.4826*median(selected);
          const double huber = 1.345*std::max(scale, 1e-8);
          for (size_t i=0; i<x.size(); ++i)
            weights[i] = keep[i]*std::min(1.0, huber/std::max(std::fabs(residual[i]-location),1e-12));
          double next_a, next_b;
          if ( !line(x,y,weights,next_a,next_b) )
            return failed("degenerate robust weighted regression");
          const double change = std::fabs(next_a-a)+std::fabs(next_b-b);
          a = next_a; b = next_b;
          if ( change < 1e-8 ) break;
        }
      if (pass==0) initial=complete_fit(f,x,y,keep,a,b,lwr,upr,pivot);
      if ( pass == iterations ) break;
      std::vector<double> next = keep;
      int count = 0;
      for (size_t i=0; i<x.size(); ++i)
        {
          // Once rejected, a peak cannot pull subsequent fits upward.
          next[i] = keep[i] && y[i]-a-b*x[i] <= location+std::max(1.0,threshold*scale);
          count += next[i] != 0;
        }
      if ( count < 5 )
        return failed("peak rejection left fewer than five frequency bins");
      if ( next == keep ) break;
      keep.swap(next);
    }
  const auto final=complete_fit(f,x,y,keep,a,b,lwr,upr,pivot);
  return final.valid ? final : failed(final.failure);
}

} // namespace psd_aperiodic
#endif
