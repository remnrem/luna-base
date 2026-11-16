
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

#include "spectral/psd-aperiodic.h"
#include "param.h"
#include "edf/edf.h"
#include "fftw/bandaid.h"
#include "db/db.h"

extern writer_t writer;
extern logger_t logger;

psd_aperiodic_t::psd_aperiodic_t( param_t & param, double min, double max )
  : enabled(param.has("aperiodic")), db(param.has("dB")), silent(param.has("silent")),
    epoch_output(param.has("epoch") || param.has("epoch-spectrum")),
    spectrum_output(param.has("spectrum") || param.has("epoch-spectrum")),
    epoch_spectrum(param.has("epoch-spectrum")),
    adding(param.has("add") || param.has("add-spectrum")),
    adding_spectrum(param.has("add-spectrum")), lwr(min), upr(max), threshold(3),
    bin_fac(param.has("fac") ? param.requires_int("fac") : 1), iterations(5)
{
  if ( !enabled )
    {
      if (param.has("aperiodic-th") || param.has("aperiodic-iter") || param.has("aperiodic-pivot"))
        Helper::halt("PSD aperiodic-th/aperiodic-iter/aperiodic-pivot require aperiodic");
      return;
    }
  if ( !param.empty("aperiodic") )
    {
      const std::vector<double> range = param.dblvector("aperiodic");
      if (range.size()!=2) Helper::halt("expecting aperiodic=lwr,upr (or bare aperiodic)");
      lwr = range[0]; upr = range[1];
    }
  if (param.has("aperiodic-th")) threshold = param.requires_dbl("aperiodic-th");
  if (param.has("aperiodic-iter")) iterations = param.requires_int("aperiodic-iter");
  if (param.has("aperiodic-pivot")) pivot = param.requires_dbl("aperiodic-pivot");
  if (!std::isfinite(pivot) || pivot<=0)
    Helper::halt("invalid PSD aperiodic-pivot (expecting a positive frequency in Hz)");
  if (!std::isfinite(lwr) || !std::isfinite(upr) || lwr<=0 || upr<=lwr
      || !std::isfinite(threshold) || threshold<=0 || iterations<1 || iterations>100 || bin_fac<1)
    Helper::halt("invalid PSD aperiodic range, threshold, or iteration count (1-100)");
  if (adding)
    {
      prefix = param.empty("add") ? "" : param.value("add");
      const std::set<std::string> excluded = param.strset("skip-bands");
      skip.assign(excluded.begin(),excluded.end());
    }
  logger << "  PSD aperiodic power-law fit: " << lwr << " - " << upr
         << " Hz, positive-peak threshold " << threshold << " MAD SDs, offset pivot "
         << pivot << " Hz\n";
}

bool psd_aperiodic_t::band_in_range( int band ) const
{
  const freq_range_t & r = globals::freq_band[(frequency_band_t)band];
  return !frequencies.empty() && r.first>=lwr && r.second<=upr
    && r.first>=frequencies.front() && r.second<=frequencies.back();
}

void psd_aperiodic_t::values( double original, double background )
{
  if (background<=0 || !std::isfinite(background) || !std::isfinite(original) || original<0) return;
  if (!db || original>0)
    writer.value("APER_PSD",db ? 10*std::log10(original) : original);
  writer.value("APER",db ? 10*std::log10(background) : background);
  const double residual = original-background;
  writer.value("PER",residual); // Always signed, always linear.
  if (db && residual>0) writer.value("PER_DB",10*std::log10(residual));
  if (original>0) writer.value("PER_EXCESS_DB",10*(std::log10(original)-std::log10(background)));
}

void psd_aperiodic_t::report( const std::vector<double> & f, const std::vector<double> & p,
                             const psd_aperiodic::fit_t & fit, bool bands, bool spectrum )
{
  if (silent) return;
  writer.value("APER_VALID",fit.valid ? 1 : 0);
  writer.value("APER_FALLBACK",fit.valid && fit.fallback ? 1 : 0);
  if (!fit.valid) return;
  writer.value("APER_OFFSET",fit.offset);
  writer.value("APER_EXP",fit.exponent);
  writer.value("APER_RSQ",fit.rsq);
  writer.value("APER_N",fit.n);
  if (bands)
    {
      bandaid_t original, background;
      original.calc_bandpower(f,p);
      background.calc_bandpower(f,fit.background);
      for (auto band : original.bands)
        if (band_in_range(band))
          {
            writer.level(globals::band(band),globals::band_strat);
            values(original.fetch(band),background.fetch(band));
          }
      writer.unlevel(globals::band_strat);
    }
  if (spectrum)
    {
      // Match PSD frequency binning, but integrate/average linear components
      // before converting to dB. Never fit the display-binned spectrum.
      size_t start=0;
      while (start<f.size() && f[start]<lwr) ++start;
      const size_t width=(size_t)bin_fac;
      for (size_t i=start; i+width<=f.size(); i+=width)
        {
          if (f[i+width-1]>upr) break;
          double original=0, background=0;
          for (size_t j=i; j<i+width; ++j)
            { original+=p[j]; background+=fit.background[j]; }
          writer.level((f[i]+f[i+width-1])/2,globals::freq_strat);
          values(original/width,background/width);
        }
      writer.unlevel(globals::freq_strat);
    }
}

void psd_aperiodic_t::epoch( const std::vector<double> & f, const std::vector<double> & p, bool bands,
                            const std::vector<double> & data, const std::string & channel,
                            int epoch, double start, double stop )
{
  if (!enabled) return;
  channel_label=channel;
  if (frequencies.empty()) { frequencies=f; sum.assign(p.size(),0); }
  if (f!=frequencies || p.size()!=sum.size())
    Helper::halt("PSD aperiodic requires a consistent frequency grid across epochs");
  for (size_t i=0; i<p.size(); ++i) sum[i]+=p[i];
  ++count;
  // Exact constant input can leave minute rounding residuals after centering;
  // their window leakage must not be interpreted as an aperiodic background.
  bool varies=false, finite=true;
  for (double value : data)
    {
      if (!std::isfinite(value)) { finite=false; break; }
      if (value!=data.front()) varies=true;
    }
  any_variable = any_variable || (finite && varies);
  // Fit every window even for summary-only requests, so failure counts always
  // refer to the same analyzed windows, independently of output verbosity.
  auto fit = finite && varies ? psd_aperiodic::fit(f,p,lwr,upr,threshold,iterations,pivot)
                              : psd_aperiodic::fit_t();
  if (!finite) fit.failure = "non-finite input samples";
  else if (!varies) fit.failure = data.empty() ? "empty input window" : "constant input window";
  // Stage a complete window before appending any values. If one requested
  // component is undefined, replace the whole window and flag it consistently.
  bool usable=fit.valid;
  std::map<std::string,double> candidate;
  const auto stage = [&](const std::string & suffix, double op, double ap) {
    candidate["_APER"+suffix]=0;
    candidate["_PEX"+suffix]=0;
    if (!fit.valid) return;
    if (op<=0 || ap<=0 || !std::isfinite(op) || !std::isfinite(ap))
      { usable=false; return; }
    const double background=db ? 10*std::log10(ap) : ap;
    const double excess=10*(std::log10(op)-std::log10(ap));
    if (!std::isfinite(background) || !std::isfinite(excess))
      { usable=false; return; }
    candidate["_APER"+suffix]=background;
    candidate["_PEX"+suffix]=excess;
  };
  if (adding && bands)
    {
      bandaid_t original, background;
      if (fit.valid)
        { original.calc_bandpower(f,p); background.calc_bandpower(f,fit.background); }
      int index=0;
      for (auto band : original.bands)
        {
          if (std::find(skip.begin(),skip.end(),globals::band(band))!=skip.end()) continue;
          ++index; // Preserve the numbering of existing PSD _Bn channels.
          if (!band_in_range(band)) continue;
          const std::string suffix="_B"+Helper::int2str(index);
          stage(suffix,original.fetch(band),background.fetch(band));
        }
    }
  if (adding && adding_spectrum)
    for (size_t i=0; i<f.size(); ++i)
      if (f[i]>=lwr && f[i]<=upr)
        {
          const std::string suffix="_F"+Helper::int2str((int)i);
          stage(suffix,p[i],fit.valid ? fit.background[i] : 0);
        }
  if (fit.valid && !usable) fit.failure="undefined added component power";
  fit.valid=usable;
  if (!usable)
    {
      ++failed_count;
      if (finite && !varies) ++flat_count;
      max_gap=std::max(max_gap,++gap);
      if (failed_count==1)
        logger << "  PSD aperiodic: first invalid window for " << channel
               << ", epoch " << epoch << ", interval " << start << "-" << stop
               << " seconds: " << fit.failure << "\n";
    }
  else
    {
      gap=0;
      if (fit.fallback) ++fallback_count;
    }
  if (epoch_output) report(f,p,fit,bands,epoch_spectrum);
  if (adding)
    {
      // Fit parameters share component validity and the same 1 Hz alignment.
      // These stay in log10 power units regardless of the dB display option.
      candidate["_APER_SLOPE"]=-fit.exponent;
      candidate["_APER_OFFSET"]=fit.offset;
      candidate["_APER_RSQ"]=fit.rsq;
      valid_windows.push_back(usable ? 1 : 0);
      fallback_windows.push_back(usable && fit.fallback ? 1 : 0);
      for (const auto & value : candidate)
        tracks[value.first].push_back(usable ? value.second : 0);
    }
}

void psd_aperiodic_t::summary( bool bands )
{
  if (!enabled || !count) return;
  std::vector<double> mean=sum;
  for (double & v : mean) v/=count;
  const auto fit = any_variable ? psd_aperiodic::fit(frequencies,mean,lwr,upr,threshold,iterations,pivot)
                                : psd_aperiodic::fit_t();
  report(frequencies,mean,fit,bands,spectrum_output);
  if (!silent)
    {
      writer.value("APER_NWIN",count);
      writer.value("APER_NFAIL",failed_count);
      writer.value("APER_NFALLBACK",fallback_count);
      writer.value("APER_NFLAT",flat_count);
      writer.value("APER_MAXGAP",max_gap);
      writer.value("APER_NFILL",adding ? failed_count : 0);
    }
  logger << "  PSD aperiodic windows for " << channel_label << ": " << count
         << " total, " << failed_count << " failed (" << flat_count << " constant), "
         << fallback_count << " robust-fit fallbacks, longest invalid run " << max_gap
         << " windows" << (adding ? "; failed windows will be filled\n" : "\n");
}

void psd_aperiodic_t::add( edf_t & edf, const std::string & label,
                         int observed, int expected, int pad1, int pad2 )
{
  if (!enabled || !count) return;
  if (edf.timeline.num_epochs()!=edf.timeline.num_total_epochs())
    Helper::halt("PSD aperiodic add requires all epochs to be unmasked");
  for (const auto & track : tracks)
    {
      const bool parameter = track.first=="_APER_SLOPE"
        || track.first=="_APER_OFFSET" || track.first=="_APER_RSQ";
      std::vector<double> data=track.second;
      if ((int)data.size()!=observed || observed+pad1+pad2!=expected)
        Helper::halt("PSD aperiodic add requires all regularly spaced epochs; track length/padding mismatch");
      psd_aperiodic::fill(data,valid_windows,
                         !parameter && track.first.find("_APER")==0 ? (db ? -120 : 1e-12) : 0);
      for (double value : data)
        if (!std::isfinite(value)) Helper::halt("internal error: non-finite filled PSD aperiodic track");
      data.insert(data.begin(),pad1,0);
      data.insert(data.end(),pad2,0);
      if (failed_count==count && (parameter || track.first.find("_PEX")==0))
        edf.add_signal(prefix+label+track.first,1,data,0,1,0,32767);
      else if (failed_count==count && !db)
        // Anchor the physical minimum at the floor. Automatic EDF scaling of
        // an almost-constant tiny signal would otherwise round it to ~1e-5.
        edf.add_signal(prefix+label+track.first,1,data,1e-12,1,0,32767);
      else
        edf.add_signal(prefix+label+track.first,1,data);
    }
  // Quality channels use identical time alignment. Padding is invalid and
  // never counted as a failed analyzed window.
  for (const auto & quality : { std::make_pair("_APER_VALID",valid_windows),
                               std::make_pair("_APER_FALLBACK",fallback_windows) })
    {
      std::vector<double> data=quality.second;
      if ((int)data.size()!=observed || observed+pad1+pad2!=expected)
        Helper::halt("PSD aperiodic add: quality track length/padding mismatch");
      data.insert(data.begin(),pad1,0);
      data.insert(data.end(),pad2,0);
      edf.add_signal(prefix+label+quality.first,1,data,0,1,0,32767);
    }
  if (!tracks.empty())
    logger << "  added " << tracks.size() << " aperiodic/excess signals for " << label
           << " at 1 Hz, with " << pad1 << " leading and " << pad2 << " trailing padding samples\n";
}
