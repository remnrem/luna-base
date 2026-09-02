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

#include "stats/dpp-spec.h"

#include "helper/helper.h"
#include "helper/luna_io.h"

#include <cmath>
#include <fstream>

std::map<std::string, dpp_feature_t> dpp_specs_t::lab2ftr;
std::map<dpp_feature_t, std::string> dpp_specs_t::ftr2lab;

namespace {

// "lo-hi" -> (lwr,upr)
bool parse_range(const std::string &s, double *lwr, double *upr) {
  size_t p = s.find('-');
  if (p == std::string::npos || p == 0 || p == s.size() - 1)
    return false;
  if (!Helper::str2dbl(s.substr(0, p), lwr))
    return false;
  if (!Helper::str2dbl(s.substr(p + 1), upr))
    return false;
  return true;
}

std::vector<double> parse_window_list(const std::string &s) {
  std::vector<std::string> tok = Helper::parse(s, ',');
  std::vector<double> w;
  for (int i = 0; i < (int)tok.size(); i++) {
    double d;
    if (!Helper::str2dbl(tok[i], &d) || d <= 0)
      Helper::halt("bad windows= value: " + tok[i]);
    w.push_back(d);
  }
  return w;
}

} // namespace

bool dpp_spec_t::has(const std::string &key) const {
  return arg.find(key) != arg.end();
}

double dpp_spec_t::narg(const std::string &key) const {
  if (!has(key))
    return 0;
  double d;
  if (!Helper::str2dbl(arg.find(key)->second, &d))
    Helper::halt("problem converting string -> numeric: " + key);
  return d;
}

int dpp_spec_t::cols() const {
  switch (ftr) {
  case DPP_PSD:
    return 5; // canonical: delta, theta, alpha, sigma, beta
  case DPP_SLOPE:
    return 1;
  case DPP_HJORTH:
    return 3; // activity, mobility, complexity
  case DPP_SKEW:
    return 1;
  case DPP_KURTOSIS:
    return 1;
  case DPP_MSE:
    return 1;
  case DPP_ENVELOPE:
    return 3; // mean, sd, cv
  case DPP_PLV:
    return 2; // plv, mean_ipc
  case DPP_COH:
    return (band.size() == 1 && band[0] != "") ? 1 : 5;
  case DPP_PSI:
    return 1;
  case DPP_CATCH22:
    return (has("catch24") && Helper::yesno(arg.at("catch24"))) ? 24 : 22;
  case DPP_PAC:
    return 2; // normalized mean-vector-length, preferred phase (radians)
  }
  return 0;
}

std::string dpp_spec_t::label_root() const {
  const std::string &ftrlab = dpp_specs_t::ftr2lab[ftr];
  // "-", not "," (comma-separated channel pairs are the *input* grammar,
  // e.g. "PLV C3,C4" -- but a comma inside an output VAR level name is
  // fragile for downstream CSV-ish consumers)
  const std::string chlab = ch.size() == 2 ? (ch[0] + "-" + ch[1]) : ch[0];
  std::string s = ftrlab + "." + chlab;

  if (ftr == DPP_ENVELOPE || ftr == DPP_PLV || ftr == DPP_PAC) {
    if (band.size() == 2 && band[0] != band[1])
      s += "." + band[0] + "-" + band[1];
    else if (band.size() >= 1)
      s += "." + band[0];
  }
  // COH/PSI (optional/required single band=), and HJORTH/SKEW/KURTOSIS/MSE/
  // CATCH22's own optional band= (each expanded spec carries exactly one
  // band value, per read()'s band_expandable handling) -- same suffix rule
  else if (band.size() == 1 && band[0] != "")
    s += "." + band[0];

  return s;
}

void dpp_specs_t::init() {
  chs.clear();
  filters.clear();
  specs.clear();

  lab2ftr.clear();
  ftr2lab.clear();

  lab2ftr["PSD"] = DPP_PSD;
  lab2ftr["SLOPE"] = DPP_SLOPE;
  lab2ftr["HJORTH"] = DPP_HJORTH;
  lab2ftr["SKEW"] = DPP_SKEW;
  lab2ftr["KURTOSIS"] = DPP_KURTOSIS;
  lab2ftr["MSE"] = DPP_MSE;
  lab2ftr["ENVELOPE"] = DPP_ENVELOPE;
  lab2ftr["PLV"] = DPP_PLV;
  lab2ftr["COH"] = DPP_COH;
  lab2ftr["PSI"] = DPP_PSI;
  lab2ftr["CATCH22"] = DPP_CATCH22;
  lab2ftr["PAC"] = DPP_PAC;

  ftr2lab[DPP_PSD] = "PSD";
  ftr2lab[DPP_SLOPE] = "SLOPE";
  ftr2lab[DPP_HJORTH] = "HJORTH";
  ftr2lab[DPP_SKEW] = "SKEW";
  ftr2lab[DPP_KURTOSIS] = "KURTOSIS";
  ftr2lab[DPP_MSE] = "MSE";
  ftr2lab[DPP_ENVELOPE] = "ENVELOPE";
  ftr2lab[DPP_PLV] = "PLV";
  ftr2lab[DPP_COH] = "COH";
  ftr2lab[DPP_PSI] = "PSI";
  ftr2lab[DPP_CATCH22] = "CATCH22";
  ftr2lab[DPP_PAC] = "PAC";

  if (default_window_sec <= 0)
    default_window_sec = 30;
  if (default_step_sec <= 0)
    default_step_sec = 30;
}

void dpp_specs_t::read(const std::string &f) {
  if (specs.size())
    return;

  init();

  const std::string expanded = Helper::expand(f);
  if (!Helper::fileExists(expanded))
    Helper::halt("could not open " + f);

  std::ifstream IN1 = LunaIO::open_ifstream(expanded, std::ios::in);

  while (1) {
    std::string line;
    Helper::safe_getline(IN1, line);
    if (IN1.bad() || IN1.eof())
      break;
    if (line == "")
      continue;
    if (line[0] == '%')
      continue;

    std::vector<std::string> tok = Helper::parse(line, " \t");
    if (tok.size() == 0)
      continue;

    // CH <label> [prefilter=lo-hi]
    if (Helper::toupper(tok[0]) == "CH") {
      if (tok.size() < 2)
        Helper::halt("expecting: CH label [prefilter=lo-hi]");
      dpp_channel_t c;
      c.ch = tok[1];
      for (int i = 2; i < (int)tok.size(); i++) {
        std::vector<std::string> kv = Helper::parse(tok[i], '=');
        if (kv.size() != 2)
          Helper::halt("bad CH arg: " + tok[i]);
        if (Helper::toupper(kv[0]) == "PREFILTER") {
          double lwr, upr;
          if (!parse_range(kv[1], &lwr, &upr))
            Helper::halt("bad prefilter= range (expecting lo-hi): " + kv[1]);
          c.has_prefilter = true;
          c.prefilter_lwr = lwr;
          c.prefilter_upr = upr;
        } else
          Helper::halt("unrecognized CH arg: " + kv[0]);
      }
      chs[c.ch] = c;
      continue;
    }

    // FILTER <name> <lwr> <upr> [ripple=][tw=]
    if (Helper::toupper(tok[0]) == "FILTER") {
      if (tok.size() < 4)
        Helper::halt("expecting: FILTER name lwr upr [ripple=][tw=]");
      dpp_filter_t filt;
      filt.name = tok[1];
      if (!Helper::str2dbl(tok[2], &filt.lwr))
        Helper::halt("bad FILTER lwr: " + line);
      if (!Helper::str2dbl(tok[3], &filt.upr))
        Helper::halt("bad FILTER upr: " + line);
      filt.ripple = 0.02;
      filt.tw = 1;
      for (int i = 4; i < (int)tok.size(); i++) {
        std::vector<std::string> kv = Helper::parse(tok[i], '=');
        if (kv.size() != 2)
          Helper::halt("bad FILTER arg: " + tok[i]);
        if (Helper::toupper(kv[0]) == "RIPPLE") {
          if (!Helper::str2dbl(kv[1], &filt.ripple))
            Helper::halt("bad ripple=: " + kv[1]);
        } else if (Helper::toupper(kv[0]) == "TW") {
          if (!Helper::str2dbl(kv[1], &filt.tw))
            Helper::halt("bad tw=: " + kv[1]);
        } else
          Helper::halt("unrecognized FILTER arg: " + kv[0]);
      }
      filters[filt.name] = filt;
      continue;
    }

    // feature block line:  block: FEATURE ch[,ch2] [band=][windows=][args]
    const std::string &block1 = tok[0];
    if (block1.size() == 1 || block1[block1.size() - 1] != ':')
      Helper::halt("expecting colon after block name\n  block: <feature> "
                   "<channel(s)> [args]");
    const std::string block =
        Helper::toupper(block1.substr(0, block1.size() - 1));

    if (tok.size() < 3)
      Helper::halt("expecting: block: FEATURE channel(s) [args]");

    const std::string ftr_lab = Helper::toupper(tok[1]);
    if (lab2ftr.find(ftr_lab) == lab2ftr.end())
      Helper::halt("feature not recognized: " + tok[1]);
    const dpp_feature_t ftr = lab2ftr[ftr_lab];
    const bool two_channel =
        (ftr == DPP_PLV || ftr == DPP_COH || ftr == DPP_PSI || ftr == DPP_PAC);

    std::vector<std::string> chan_tokens;
    std::map<std::string, std::string> targs;

    for (int i = 2; i < (int)tok.size(); i++) {
      std::vector<std::string> kv = Helper::parse(tok[i], '=');
      if (kv.size() == 1)
        chan_tokens.push_back(tok[i]);
      else if (kv.size() == 2)
        targs[kv[0]] = kv[1];
      else
        Helper::halt("bad token: " + tok[i]);
    }

    if (chan_tokens.size() == 0)
      Helper::halt("no channel(s) specified: " + line);

    // PAC (phase-channel,amplitude-channel) is deliberately restricted to
    // one pair per line, unlike PLV/COH/PSI which allow several
    // space-separated pairs -- with a directional phase-first/amplitude-
    // second role, crossing multiple pairs on one line risks silently
    // pairing the wrong roles across pairs; one line, one pair keeps it
    // unambiguous
    if (ftr == DPP_PAC && chan_tokens.size() > 1)
      Helper::halt("PAC allows only one channel pair per line: " + line);

    // window list for this block: arg override, else global default
    std::vector<double> windows;
    if (targs.find("windows") != targs.end())
      windows = parse_window_list(targs["windows"]);
    else
      windows.push_back(default_window_sec);

    for (int c = 0; c < (int)chan_tokens.size(); c++) {
      std::vector<std::string> chlabs;

      if (two_channel) {
        std::vector<std::string> pair = Helper::parse(chan_tokens[c], ',');
        if (pair.size() != 2)
          Helper::halt(
              ftr_lab +
              " requires a comma-separated channel pair, e.g. C3,C4: " + line);
        chlabs = pair;
      } else
        chlabs.push_back(chan_tokens[c]);

      for (int k = 0; k < (int)chlabs.size(); k++)
        if (!has_channel(chlabs[k]))
          Helper::halt(chlabs[k] + " not specified via CH: " + line);

      // band= -- five different roles by feature:
      //  - PSD/SLOPE: never accepted (band= would silently do nothing at
      //    compute time otherwise -- reject loudly instead)
      //  - ENVELOPE: required, one shared value (always single-channel,
      //    so chlabs.size()==1 here -- the "one per channel" half of the
      //    old shared rule never actually applied to it)
      //  - PLV: required, and now *within-band only* -- a single value
      //    (shared by both channels) is the only accepted form; two
      //    distinct values (cross-band, same-channel "coupling") used to
      //    be allowed but is removed -- it measured phase-phase
      //    synchronization, not what SO-spindle-style coupling studies
      //    actually want (phase-*amplitude* coupling -- see PAC below)
      //  - PSI: required, exactly one value; COH: optional, exactly one
      //    value if given (else all 5 canonical bands) -- unchanged
      //  - PAC: required, exactly *two* values, band[0]=phase source,
      //    band[1]=amplitude source -- never a single shared value
      //    (phase/amplitude of the same band is a degenerate, not the
      //    intended, use of this feature)
      //  - HJORTH/SKEW/KURTOSIS/MSE/CATCH22: optional (default: computed
      //    on the raw/prefiltered trace, exactly as before band= support
      //    existed); if given, one *or more* comma-separated values --
      //    each expands into its own spec, exactly like windows= already
      //    crosses a block across multiple window lengths. Since these
      //    features are always single-channel (two_channel is false),
      //    chlabs.size()==1 here always, so this never collides with
      //    ENVELOPE/PLV's pair-vs-shared distinction above.
      const bool band_rejected = (ftr == DPP_PSD || ftr == DPP_SLOPE);
      const bool band_env_plv = (ftr == DPP_ENVELOPE || ftr == DPP_PLV);
      const bool band_one_only = (ftr == DPP_COH || ftr == DPP_PSI);
      const bool band_pac = (ftr == DPP_PAC);
      const bool band_expandable =
          (ftr == DPP_HJORTH || ftr == DPP_SKEW || ftr == DPP_KURTOSIS ||
           ftr == DPP_MSE || ftr == DPP_CATCH22);

      std::vector<std::string> bands;
      bool expand_bands = false;

      if (targs.find("band") != targs.end()) {
        if (band_rejected)
          Helper::halt(ftr_lab + " does not support band= : " + line);

        std::vector<std::string> blist = Helper::parse(targs["band"], ',');

        if (band_one_only) {
          if (blist.size() != 1)
            Helper::halt(ftr_lab + " band= takes exactly one value: " + line);
          bands = blist;
        } else if (band_pac) {
          if (blist.size() != 2)
            Helper::halt("PAC requires band=<phase-band>,<amplitude-band> "
                         "(exactly two values): " +
                         line);
          bands = blist;
        } else if (band_env_plv) {
          if (blist.size() == 1 && chlabs.size() == 2)
            bands = {blist[0], blist[0]};
          else if (blist.size() == chlabs.size())
            bands = blist;
          else
            Helper::halt("band= must give one value, or one per channel: " +
                         line);

          if (ftr == DPP_PLV && bands.size() == 2 && bands[0] != bands[1])
            Helper::halt(
                "PLV requires the same band for both channels -- cross-band "
                "phase-phase coupling is not supported (use PAC for "
                "cross-frequency "
                "phase-amplitude coupling): " +
                line);
        } else if (band_expandable) {
          bands = blist; // one spec per band value, expanded below
          expand_bands = true;
        }

        for (int k = 0; k < (int)bands.size(); k++)
          if (bands[k] != "" && !has_filter(bands[k]))
            Helper::halt(bands[k] + " not specified via FILTER: " + line);
      } else if (ftr != DPP_COH && ftr != DPP_PSI && ftr != DPP_PAC)
        bands.assign(chlabs.size(), "");

      if (band_env_plv)
        for (int k = 0; k < (int)bands.size(); k++)
          if (bands[k] == "")
            Helper::halt(ftr_lab + " requires band=: " + line);

      if (ftr == DPP_PSI && bands.size() == 0)
        Helper::halt("PSI requires band=: " + line);

      if (band_pac && bands.size() != 2)
        Helper::halt("PAC requires band=<phase-band>,<amplitude-band> (exactly "
                     "two values): " +
                     line);

      for (int w = 0; w < (int)windows.size(); w++) {
        if (expand_bands) {
          for (int b = 0; b < (int)bands.size(); b++) {
            dpp_spec_t spec;
            spec.block = block;
            spec.ftr = ftr;
            spec.ch = chlabs;
            spec.band = {bands[b]};
            spec.window_sec = windows[w];
            spec.arg = targs;
            specs.push_back(spec);
          }
        } else {
          dpp_spec_t spec;
          spec.block = block;
          spec.ftr = ftr;
          spec.ch = chlabs;
          spec.band = bands;
          spec.window_sec = windows[w];
          spec.arg = targs;
          specs.push_back(spec);
        }
      }
    }
  }

  IN1.close();
}

void dpp_specs_t::init_default(const std::vector<std::string> &channels) {
  init();

  for (int i = 0; i < (int)channels.size(); i++) {
    dpp_channel_t dc;
    dc.ch = channels[i];
    chs[channels[i]] = dc;
  }

  for (int i = 0; i < (int)channels.size(); i++) {
    const std::string &c = channels[i];

    dpp_spec_t s1;
    s1.block = "BASE";
    s1.ftr = DPP_PSD;
    s1.ch = {c};
    s1.window_sec = default_window_sec;
    dpp_spec_t s2;
    s2.block = "BASE";
    s2.ftr = DPP_SLOPE;
    s2.ch = {c};
    s2.window_sec = default_window_sec;
    dpp_spec_t s3;
    s3.block = "BASE";
    s3.ftr = DPP_HJORTH;
    s3.ch = {c};
    s3.window_sec = default_window_sec;

    specs.push_back(s1);
    specs.push_back(s2);
    specs.push_back(s3);
  }
}

void dpp_specs_t::apply_inline_overrides(const std::string &windows_arg,
                                         const std::string &filters_arg,
                                         const std::string &features_arg,
                                         const std::string &prefilter_arg) {
  // filters=name:lo-hi,name2:lo-hi2,...
  if (filters_arg != "") {
    std::vector<std::string> toks = Helper::parse(filters_arg, ',');
    for (int i = 0; i < (int)toks.size(); i++) {
      std::vector<std::string> kv = Helper::parse(toks[i], ':');
      if (kv.size() != 2)
        Helper::halt("bad filters= entry (expecting name:lo-hi): " + toks[i]);
      double lwr, upr;
      if (!parse_range(kv[1], &lwr, &upr))
        Helper::halt("bad filters= range: " + toks[i]);
      dpp_filter_t filt;
      filt.name = kv[0];
      filt.lwr = lwr;
      filt.upr = upr;
      filt.ripple = 0.02;
      filt.tw = 1;
      filters[filt.name] = filt;
    }
  }

  // prefilter=lo-hi  (applied to every already-declared channel)
  if (prefilter_arg != "") {
    double lwr, upr;
    if (!parse_range(prefilter_arg, &lwr, &upr))
      Helper::halt("bad prefilter= range: " + prefilter_arg);
    std::map<std::string, dpp_channel_t>::iterator ii = chs.begin();
    while (ii != chs.end()) {
      ii->second.has_prefilter = true;
      ii->second.prefilter_lwr = lwr;
      ii->second.prefilter_upr = upr;
      ++ii;
    }
  }

  // features=PSD,HJORTH,ENVELOPE,...  (add for every already-declared channel;
  // connectivity features are not addable this way -- they require explicit
  // channel-pairing via spec=)
  if (features_arg != "") {
    std::vector<std::string> toks = Helper::parse(features_arg, ',');
    const std::string default_band =
        filters.empty() ? "" : filters.begin()->first;

    for (int i = 0; i < (int)toks.size(); i++) {
      const std::string lab = Helper::toupper(toks[i]);
      if (lab2ftr.find(lab) == lab2ftr.end())
        Helper::halt("unrecognized feature in features=: " + toks[i]);
      const dpp_feature_t ftr = lab2ftr[lab];
      if (ftr == DPP_PLV || ftr == DPP_COH || ftr == DPP_PSI)
        continue;

      std::map<std::string, dpp_channel_t>::const_iterator ii = chs.begin();
      while (ii != chs.end()) {
        dpp_spec_t s;
        s.block = "EXTRA";
        s.ftr = ftr;
        s.ch = {ii->first};
        s.window_sec = default_window_sec;
        if (ftr == DPP_ENVELOPE) {
          if (default_band == "")
            Helper::halt("ENVELOPE requested via features= but no filters= "
                         "band declared");
          s.band = {default_band};
        }
        specs.push_back(s);
        ++ii;
      }
    }
  }

  // windows=w1,w2,...  (cross every existing spec across this window list)
  if (windows_arg != "") {
    std::vector<double> wl = parse_window_list(windows_arg);
    if (wl.size() > 0) {
      std::vector<dpp_spec_t> expanded;
      for (int i = 0; i < (int)specs.size(); i++)
        for (int w = 0; w < (int)wl.size(); w++) {
          dpp_spec_t s2 = specs[i];
          s2.window_sec = wl[w];
          expanded.push_back(s2);
        }
      specs = expanded;
      default_window_sec = wl[0];
    }
  }
}
