
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


#ifdef HAS_LGBM

#include "pops/posteriors.h"

#include "edf/edf.h"
#include "helper/helper.h"
#include "helper/logger.h"

extern logger_t logger;

namespace {

std::vector<std::string> posterior_labels( const bool three_state )
{
  if ( three_state )
    return std::vector<std::string>{ "W" , "R" , "NR" };
  return std::vector<std::string>{ "W" , "R" , "N1" , "N2" , "N3" };
}

bool can_add_epoch_posteriors( edf_t * pedf , const int ne_total , const std::string & source_tag )
{
  if ( pedf == NULL )
    {
      logger << "  ** " << source_tag << ": no attached EDF, cannot add posterior channels\n";
      return false;
    }

  pedf->timeline.ensure_epoched();

  if ( pedf->is_actually_discontinuous() )
    {
      logger << "  ** " << source_tag << ": EDF is discontinuous, cannot add epoch posterior channels\n";
      return false;
    }

  if ( ! pedf->timeline.fixed_epoch_length() ||
       ! Helper::similar( pedf->timeline.epoch_length() , 30.0 , 1e-6 ) ||
       ! Helper::similar( pedf->timeline.epoch_inc() , 30.0 , 1e-6 ) )
    {
      logger << "  ** " << source_tag << ": EDF epochs are not fixed 30s non-overlapping epochs, skipping posterior channels\n";
      return false;
    }

  if ( ! Helper::similar( pedf->header.record_duration , 30.0 , 1e-6 ) )
    {
      logger << "  ** " << source_tag << ": EDF record duration is "
             << Helper::dbl2str( pedf->header.record_duration )
             << "s, but posterior channel export requires 30s records (one record = one epoch = one sample); skipping\n";
      return false;
    }

  if ( pedf->header.nr != ne_total )
    {
      logger << "  ** " << source_tag << ": EDF has " << pedf->header.nr
             << " records but POPS epoch map has " << ne_total
             << " total epochs; skipping posterior channels\n";
      return false;
    }

  return true;
}

}

std::string pops_posteriors::invalid_channel_label( const std::string & prefix )
{
  return prefix + "_INVALID";
}

bool pops_posteriors::add_edf_channels( edf_t * pedf ,
                                        const Eigen::MatrixXd & P ,
                                        const std::vector<int> & E ,
                                        const std::vector<bool> * flagged ,
                                        int ne_total ,
                                        bool three_state ,
                                        const std::string & prefix ,
                                        const std::string & source_tag )
{
  if ( prefix == "" )
    {
      logger << "  ** " << source_tag << ": empty posterior-prefix, skipping EDF channel export\n";
      return false;
    }

  if ( ! can_add_epoch_posteriors( pedf , ne_total , source_tag ) )
    return false;

  const std::vector<std::string> labels = posterior_labels( three_state );
  const int nc = labels.size();

  if ( P.cols() != nc )
    Helper::halt( source_tag + ": posterior column count does not match requested state space" );

  std::vector< std::vector<double> > data( nc , std::vector<double>( ne_total , -9.0 ) );
  std::vector<double> invalid( ne_total , 1.0 );
  int n_bad = 0;

  for (int e = 0; e < (int)E.size(); e++)
    {
      const int epoch = E[e];
      if ( epoch < 0 || epoch >= ne_total ) continue;

      const bool bad = flagged != NULL && e < (int)flagged->size() && (*flagged)[e];
      if ( bad )
        {
          ++n_bad;
          continue;
        }

      invalid[epoch] = 0.0;
      for (int j = 0; j < nc; j++)
        data[j][epoch] = P(e,j);
    }

  int added = 0;
  for (int j = 0; j < nc; j++)
    {
      const std::string ch = prefix + "_" + labels[j];
      if ( pedf->header.has_signal( ch ) )
        {
          logger << "  ** " << source_tag << ": channel " << ch
                 << " already exists, skipping\n";
          continue;
        }

      pedf->add_signal( ch , -1 , data[j] , -9.0 , 1.0 );
      ++added;
    }

  const std::string invalid_ch = invalid_channel_label( prefix );
  if ( pedf->header.has_signal( invalid_ch ) )
    logger << "  ** " << source_tag << ": channel " << invalid_ch
           << " already exists, skipping\n";
  else
    {
      pedf->add_signal( invalid_ch , -1 , invalid , 0.0 , 1.0 );
      ++added;
    }

  if ( added > 0 )
    logger << "  " << source_tag << ": added " << added
           << " channel(s) to EDF with prefix '" << prefix << "'"
           << ( n_bad ? " (" + Helper::int2str( n_bad ) + " flagged/skipped epochs written as -9)" : "" )
           << "\n";

  return added > 0;
}

bool pops_posteriors::add_edf_signals( edf_t * pedf ,
                                       const Eigen::MatrixXd & P ,
                                       const std::vector<int> & E ,
                                       const std::vector<bool> * flagged ,
                                       int ne_total ,
                                       bool three_state ,
                                       double row_duration_sec ,
                                       const std::string & prefix ,
                                       const std::string & source_tag )
{
  if ( prefix == "" )
    {
      logger << "  ** " << source_tag << ": empty posterior-prefix, skipping EDF signal export\n";
      return false;
    }

  if ( pedf == NULL )
    {
      logger << "  ** " << source_tag << ": no attached EDF, cannot add posterior signals\n";
      return false;
    }

  if ( pedf->is_actually_discontinuous() )
    {
      logger << "  ** " << source_tag << ": EDF is discontinuous, cannot add posterior signals\n";
      return false;
    }

  if ( row_duration_sec <= 0 )
    {
      logger << "  ** " << source_tag << ": invalid row duration, skipping EDF signal export\n";
      return false;
    }

  const double n_spr_d = pedf->header.record_duration / row_duration_sec;
  const int n_spr = (int)std::lround( n_spr_d );
  if ( std::fabs( n_spr_d - n_spr ) > 1e-6 || n_spr < 1 )
    {
      logger << "  ** " << source_tag << ": EDF record duration is "
             << Helper::dbl2str( pedf->header.record_duration )
             << "s, but posterior export at " << Helper::dbl2str( row_duration_sec )
             << "s requires an integer number of samples per record; skipping\n";
      return false;
    }

  const std::vector<std::string> labels = posterior_labels( three_state );
  const int nc = labels.size();
  if ( P.cols() != nc )
    Helper::halt( source_tag + ": posterior column count does not match requested state space" );

  std::vector< std::vector<double> > data( nc , std::vector<double>( ne_total , -9.0 ) );
  int n_bad = 0;
  for (int e = 0; e < (int)E.size(); e++)
    {
      const int epoch = E[e];
      if ( epoch < 0 || epoch >= ne_total ) continue;
      const bool bad = flagged != NULL && e < (int)flagged->size() && (*flagged)[e];
      if ( bad )
        {
          ++n_bad;
          continue;
        }
      for (int j = 0; j < nc; j++)
        data[j][epoch] = P(e,j);
    }

  int added = 0;
  const double Fs_code = -(double)n_spr;
  for (int j = 0; j < nc; j++)
    {
      const std::string ch = prefix + "_" + labels[j];
      if ( pedf->header.has_signal( ch ) )
        {
          logger << "  ** " << source_tag << ": channel " << ch
                 << " already exists, skipping\n";
          continue;
        }

      pedf->add_signal( ch , Fs_code , data[j] , -9.0 , 1.0 );
      ++added;
    }

  if ( added > 0 )
    logger << "  " << source_tag << ": added " << added
           << " posterior signal(s) to EDF with prefix '" << prefix << "'"
           << ( n_bad ? " (" + Helper::int2str( n_bad ) + " flagged/skipped rows written as -9)" : "" )
           << "\n";

  return added > 0;
}

#endif
