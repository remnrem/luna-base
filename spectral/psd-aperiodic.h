#ifndef LUNA_PSD_APERIODIC_H
#define LUNA_PSD_APERIODIC_H

#include <map>
#include <string>
#include <vector>
#include "spectral/aperiodic-fit.h"

struct param_t;
struct edf_t;

// Optional PSD sidecar. Does not modify the PSD, its bands, or existing tracks.
class psd_aperiodic_t {
public:
  psd_aperiodic_t( param_t &, double min, double max );
  void epoch( const std::vector<double> & f, const std::vector<double> & p, bool bands,
              const std::vector<double> & data, const std::string & channel,
              int epoch, double start, double stop );
  void summary( bool bands );
  void add( edf_t &, const std::string & label, int observed, int expected, int pad1, int pad2 );
private:
  bool enabled, db, silent, epoch_output, spectrum_output, epoch_spectrum,
    adding, adding_spectrum, any_variable = false;
  double lwr, upr, threshold, bin_fac;
  double pivot = 1;
  int iterations, count = 0, failed_count = 0, fallback_count = 0,
    flat_count = 0, gap = 0, max_gap = 0;
  std::string prefix, channel_label;
  std::vector<double> frequencies, sum, valid_windows, fallback_windows;
  std::map<std::string,std::vector<double> > tracks;
  std::vector<std::string> skip;
  void report( const std::vector<double> &, const std::vector<double> &,
               const psd_aperiodic::fit_t &, bool bands, bool spectrum );
  void values( double original, double background );
  bool band_in_range( int band ) const;
};
#endif
