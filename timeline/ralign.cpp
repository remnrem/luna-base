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


//
// ALIGN-SCAN: scan for epoch/record/staging alignment edge cases
//
// Reports on:
//  - EDF record structure vs epoch size (divisibility, boundary alignment)
//  - Staging annotation starts vs epoch grid (drift, gaps)
//  - Per-epoch stage ambiguity (annotation boundary mid-epoch)
//  - Record-level ambiguity for MASK/RE (records shared between differently-staged epochs)
//  - EDF+D gap annotation issues
//  - Per-label verbose table (strata ANNOT): offset stats and epoch coverage per stage label
//
// Non-standard epoch modes (generic, contig, gap-spanning, overlapping) are detected
// and reported; inapplicable checks are skipped or flagged to avoid false positives.
//
// Parameters:
//   annots=N1,N2,...  custom annotation labels  (default: N1,N2,N3,R,W,?,L,U,M)
//   minimal           use N1,N2,N3,R,W only
//   annot[=T/F]       add derived issue annotations (default: on)
//   annot-prefix=P    prefix for derived annotations (default: align_)
//   annot-epoch[=T/F] also add epoch annotations as PEPOCHEPOCH
//   annot-rec[=T/F]   also add retained EDF record annotations as PRECREC
//   verbose           also emit per-epoch details (default: per-label ANNOT table only)
//

#include "timeline/timeline.h"
#include "edf/edf.h"
#include "db/db.h"
#include "helper/helper.h"
#include "helper/logger.h"
#include "annot/annot.h"
#include "defs/defs.h"

#include <set>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <limits>

extern writer_t writer;
extern logger_t logger;


// ---------------------------------------------------------------------------
// small helpers
// ---------------------------------------------------------------------------

static std::string frac( int n , int d )
{
  return Helper::int2str(n) + "/" + Helper::int2str(d);
}

static std::string pct( int n , int d )
{
  if ( d == 0 ) return "n/a";
  return Helper::dbl2str( 100.0 * n / d , 1 ) + "%";
}


// ---------------------------------------------------------------------------
// main function
// ---------------------------------------------------------------------------

void proc_review_alignment( edf_t & edf , param_t & param )
{

  // ------------------------------------------------------------------
  // 1. Annotation labels to check
  // ------------------------------------------------------------------

  std::vector<std::string> ann_labels;

  if ( param.has( "annots" ) )
    ann_labels = param.strvector( "annots" );
  else if ( param.has( "minimal" ) )
    ann_labels = { "N1" , "N2" , "N3" , "R" , "W" };
  else
    ann_labels = { "N1" , "N2" , "N3" , "R" , "W" , "?" , "L" , "U" , "M" };

  const bool emit_annots = param.yesno( "annot" , true , true );
  const bool emit_epoch_annots = param.has( "annot-epoch" ) ? param.yesno( "annot-epoch" , false , true ) : false;
  const bool emit_rec_annots = param.has( "annot-rec" ) ? param.yesno( "annot-rec" , false , true ) : false;
  const std::string annot_prefix = param.has( "annot-prefix" ) ? param.value( "annot-prefix" ) : "align_";


  // ------------------------------------------------------------------
  // 2. Ensure epochs exist; apply default 30s/30s if not
  // ------------------------------------------------------------------

  bool auto_epoched = false;

  if ( ! edf.timeline.epoched() )
    {
      edf.timeline.set_epoch( 30.0 , 30.0 );
      auto_epoched = true;
    }


  // ------------------------------------------------------------------
  // 3. Epoch mode flags
  //    Several modes make standard alignment checks inapplicable or
  //    misleading; detect them early so checks can be gated/flagged.
  // ------------------------------------------------------------------

  const bool is_generic      = edf.timeline.generic_epochs();
  const bool is_contig       = edf.timeline.contig_epochs();
  const bool is_fixed        = edf.timeline.fixed_epoch_length();
  const bool is_overlapping  = ! edf.timeline.exactly_contiguous_epochs() && ! is_generic && ! is_contig;
  const bool is_gap_spanning = edf.timeline.epochs_span_gaps();

  // Non-standard epoch mode: any mode other than fixed-length, non-overlapping, non-gap-spanning
  const bool nonstandard_epoch = is_generic || is_contig || is_gap_spanning;

  // Overlapping epochs: annotation-to-epoch mapping is many-to-many by design
  // Gap-spanning:      epoch boundaries deliberately cross segment gaps —
  //                    record-boundary alignment is less meaningful
  // Generic/contig:    variable-length or annotation-derived epochs —
  //                    per-epoch checks are not applicable


  // ------------------------------------------------------------------
  // 4. EDF structural info
  // ------------------------------------------------------------------

  const bool is_continuous = edf.header.continuous;
  const bool is_edfplus    = edf.header.edfplus;
  const int  nrecs         = edf.header.nr;
  const double rec_dur_sec = edf.header.record_duration;
  const uint64_t rec_dur_tp = edf.header.record_duration_tp;

  std::string edf_type;
  if      ( is_continuous ) edf_type = "EDF (continuous)";
  else if ( is_edfplus )    edf_type = "EDF+D (discontinuous)";
  else                      edf_type = "EDF+C";

  std::set<interval_t> segs = edf.timeline.segments();
  std::set<interval_t> gps  = edf.timeline.gaps( segs );

  int nsegs = (int)segs.size();
  int ngaps = is_continuous ? 0 : std::max( 0 , (int)gps.size() );

  double total_dur_sec = edf.timeline.total_duration_tp / (double)globals::tp_1sec;


  // ------------------------------------------------------------------
  // 5. Staging annotation presence + collect intervals
  // ------------------------------------------------------------------

  std::vector<std::string> found_labels;
  std::vector<std::string> missing_labels;

  // All annotation instances: (interval, label)
  std::vector<std::pair<interval_t,std::string>> all_anns;

  // Count and interval list per label
  std::map<std::string,int> ann_counts;
  std::map<std::string,std::vector<interval_t>> ann_intervals; // label → intervals

  for ( const auto & label : ann_labels )
    {
      annot_t * annot = edf.annotations->find( label );
      if ( annot == NULL )
        {
          missing_labels.push_back( label );
          continue;
        }
      found_labels.push_back( label );
      ann_counts[ label ] = 0;
      for ( const auto & kv : annot->interval_events )
        {
          all_anns.push_back( { kv.first.interval , label } );
          ann_intervals[ label ].push_back( kv.first.interval );
          ++ann_counts[ label ];
        }
    }

  // Unique annotation start time-points (across all labels)
  std::set<uint64_t> ann_starts_set;
  for ( const auto & ai : all_anns )
    ann_starts_set.insert( ai.first.start );

  const int total_ann_instances = (int)all_anns.size();
  const int total_ann_starts    = (int)ann_starts_set.size();


  // ------------------------------------------------------------------
  // 6. Epoch info
  // ------------------------------------------------------------------

  const int    nepochs         = edf.timeline.num_total_epochs();
  const double epoch_dur_sec   = edf.timeline.epoch_length();
  const double epoch_inc_sec   = edf.timeline.epoch_inc();
  const double epoch_off_sec   = edf.timeline.epoch_offset();  // effective
  const std::string align_str  = edf.timeline.align_string();
  const bool align_active      = ! align_str.empty();

  const uint64_t epoch_dur_tp = (uint64_t)std::round( epoch_dur_sec * globals::tp_1sec );
  const uint64_t epoch_inc_tp = (uint64_t)std::round( epoch_inc_sec * globals::tp_1sec );

  // Build set of all epoch start time-points
  std::set<uint64_t> epoch_starts_set;
  for ( int e = 0 ; e < nepochs ; ++e )
    epoch_starts_set.insert( edf.timeline.epoch(e).start );

  // Build set of all retained record start time-points
  std::set<uint64_t> rec_starts_set;
  for ( const auto & kv : edf.timeline.rec2tp )
    rec_starts_set.insert( kv.second );


  // ------------------------------------------------------------------
  // 7. Epoch → record boundary alignment
  //    Skip / note if nonstandard mode makes this less meaningful
  // ------------------------------------------------------------------

  bool rec_divides_epoch = false;
  int epochs_on_rec  = 0;
  int epochs_off_rec = 0;

  // In gap-spanning mode, epoch boundaries cross gaps by design —
  // the divisibility / boundary checks still matter for signal within each segment,
  // but the "epoch start on record boundary" check is less informative.
  // In generic/contig mode the epoch durations are variable — not applicable.

  if ( ! is_generic && ! is_contig && rec_dur_tp > 0 )
    {
      rec_divides_epoch = ( epoch_dur_tp % rec_dur_tp == 0 );

      for ( int e = 0 ; e < nepochs ; ++e )
        {
          uint64_t es = edf.timeline.epoch(e).start;
          if ( rec_starts_set.count(es) ) ++epochs_on_rec;
          else                            ++epochs_off_rec;
        }
    }


  // ------------------------------------------------------------------
  // 8. Annotation start → epoch grid
  //    For overlapping epochs: an annotation start may intentionally fall
  //    between two epoch starts by design.  Still informative but note it.
  // ------------------------------------------------------------------

  int ann_at_epoch  = 0;
  int ann_off_epoch = 0;

  // Per-label: offset from nearest epoch boundary (for verbose table)
  // offset = min |ann_start - nearest epoch_start|  in time-points
  struct LabelStats {
    int n_at_epoch   = 0;
    int n_off_epoch  = 0;
    double sum_off_tp = 0.0;   // sum of offsets (tp) for misaligned starts
    double max_off_tp = 0.0;
    int ep_single    = 0;      // epochs where this is the sole label
    int ep_multi     = 0;      // epochs where this appears with other labels
    int ep_dominant  = 0;      // epochs where this label dominates (most overlap)
  };
  std::map<std::string,LabelStats> lstats;
  for ( const auto & l : found_labels ) lstats[l] = LabelStats();

  struct AnnIssueInfo {
    interval_t interval;
    std::string label;
    bool on_epoch = false;
    double off_sec = 0.0;
  };
  std::vector<AnnIssueInfo> ann_issue_info;

  // Build sorted vector of epoch starts for binary-search nearest-boundary
  std::vector<uint64_t> epoch_starts_vec( epoch_starts_set.begin() , epoch_starts_set.end() );

  for ( const auto & ai : all_anns )
    {
      const uint64_t s = ai.first.start;
      const std::string & label = ai.second;
      AnnIssueInfo info;
      info.interval = ai.first;
      info.label = label;

      bool on_epoch = epoch_starts_set.count(s) > 0;
      info.on_epoch = on_epoch;
      if ( on_epoch ) { ++ann_at_epoch;  ++lstats[label].n_at_epoch;  }
      else            { ++ann_off_epoch; ++lstats[label].n_off_epoch; }

      if ( ! on_epoch && ! epoch_starts_vec.empty() )
        {
          // nearest epoch boundary
          auto it = std::lower_bound( epoch_starts_vec.begin() , epoch_starts_vec.end() , s );
          uint64_t nearest = std::numeric_limits<uint64_t>::max();
          if ( it != epoch_starts_vec.end() )
            nearest = std::min( nearest , *it - s );
          if ( it != epoch_starts_vec.begin() )
            {
              --it;
              nearest = std::min( nearest , s - *it );
            }
          if ( nearest != std::numeric_limits<uint64_t>::max() )
            {
              double off_sec = nearest / (double)globals::tp_1sec;
              info.off_sec = off_sec;
              lstats[label].sum_off_tp += off_sec;
              if ( off_sec > lstats[label].max_off_tp )
                lstats[label].max_off_tp = off_sec;
            }
        }

      ann_issue_info.push_back( info );
    }


  // ------------------------------------------------------------------
  // 9. Per-epoch stage coverage
  //    Not applicable for generic/contig epochs.
  // ------------------------------------------------------------------

  int ep_no_stage    = 0;
  int ep_one_stage   = 0;
  int ep_multi_stage = 0;

  std::map<int,std::set<std::string>> epoch_stage_labels;
  std::map<int,std::string> epoch_stage_str;

  if ( ! is_generic && ! is_contig )
    {
      for ( int e = 0 ; e < nepochs ; ++e )
        {
          const interval_t ep = edf.timeline.epoch(e);
          std::set<std::string> labs;

          for ( const auto & ai : all_anns )
            if ( ai.first.start < ep.stop && ai.first.stop > ep.start )
              labs.insert( ai.second );

          epoch_stage_labels[e] = labs;
          std::string stage_str;
          for ( const auto & l : labs ) stage_str += l + " ";
          epoch_stage_str[e] = stage_str.empty() ? "." : stage_str;

          if      ( labs.empty()     ) ++ep_no_stage;
          else if ( labs.size() == 1 ) ++ep_one_stage;
          else                         ++ep_multi_stage;
        }

      // Populate per-label epoch coverage stats
      for ( int e = 0 ; e < nepochs ; ++e )
        {
          const std::set<std::string> & labs = epoch_stage_labels[e];
          if ( labs.empty() ) continue;
          for ( const auto & l : labs )
            {
              if ( lstats.find(l) == lstats.end() ) continue;
              if ( labs.size() == 1 ) ++lstats[l].ep_single;
              else                    ++lstats[l].ep_multi;
            }
        }
    }


  // ------------------------------------------------------------------
  // 10. Record-level MASK/RE ambiguity
  //     Skip if gap-spanning (epochs cross records by design there too)
  //     or generic/contig.
  // ------------------------------------------------------------------

  int records_ambig = 0;
  std::map<int,std::string> rec_stage_str;

  if ( ! is_generic && ! is_contig && ! is_gap_spanning )
    {
      for ( const auto & kv : edf.timeline.rec2epoch )
        {
          if ( kv.second.size() <= 1 ) continue;
          std::set<std::string> stages_for_rec;
          for ( int e : kv.second )
            {
              auto it = epoch_stage_labels.find(e);
              if ( it != epoch_stage_labels.end() )
                stages_for_rec.insert( it->second.begin() , it->second.end() );
            }
          if ( stages_for_rec.size() > 1 )
            {
              ++records_ambig;
              std::string stage_str;
              for ( const auto & s : stages_for_rec ) stage_str += s + " ";
              rec_stage_str[ kv.first ] = stage_str.empty() ? "." : stage_str;
            }
        }
    }


  // ------------------------------------------------------------------
  // 11. EDF+D gap annotation issues
  // ------------------------------------------------------------------

  int anns_in_gap   = 0;
  int anns_span_gap = 0;

  if ( ! is_continuous )
    {
      for ( const auto & ai : all_anns )
        {
          const interval_t & iv = ai.first;
          for ( const auto & gap : gps )
            {
              const bool a_in_gap =
                ( iv.start >= gap.start && iv.stop <= gap.stop );
              const bool a_span_left =
                ( iv.start < gap.start && iv.stop > gap.start && iv.stop <= gap.stop );
              const bool a_span_right =
                ( iv.start >= gap.start && iv.start < gap.stop && iv.stop > gap.stop );

              if ( a_in_gap )                    { ++anns_in_gap;   break; }
              if ( a_span_left || a_span_right ) { ++anns_span_gap; break; }
            }
        }
    }


  // ------------------------------------------------------------------
  // 12. EPOCH align suggestion
  //     What offset would EPOCH align pick?
  //     How many annotation starts would be on-epoch after alignment?
  //     Compute for both full and partial fix scenarios.
  // ------------------------------------------------------------------

  const uint64_t align_offset_tp  = edf.annotations->first( ann_labels );
  const double   align_offset_sec = align_offset_tp / (double)globals::tp_1sec;

  // Simulate: how many annotation starts land on epoch boundaries if we use align_offset?
  // Works for standard, non-overlapping epochs only.
  int align_would_fix_n   = 0;  // number of ann_starts that would be on-epoch
  bool align_would_fix_all = false;

  if ( ! align_active && ann_off_epoch > 0 && epoch_inc_tp > 0
       && align_offset_tp > 0 && ! is_generic && ! is_contig )
    {
      for ( uint64_t s : ann_starts_set )
        {
          if ( s < align_offset_tp ) continue;
          uint64_t diff = s - align_offset_tp;
          if ( diff % epoch_inc_tp == 0 ) ++align_would_fix_n;
        }
      align_would_fix_all = ( align_would_fix_n == total_ann_starts );
    }
  else if ( align_active )
    {
      // already aligned — report current on-epoch count as the "fixed" count
      align_would_fix_n = ann_at_epoch;
      align_would_fix_all = ( ann_off_epoch == 0 );
    }

  // Suggest EPOCH align when: not already active, there is drift, and it helps
  const bool suggest_epoch_align =
    ! align_active && ann_off_epoch > 0 && align_would_fix_n > ann_at_epoch;


  // ------------------------------------------------------------------
  // 13. Overall warning flags
  // ------------------------------------------------------------------

  const bool warn_rec_not_div   = ! nonstandard_epoch && ! rec_divides_epoch;
  const bool warn_epoch_off_rec = ! nonstandard_epoch && epochs_off_rec > 0;
  const bool warn_stage_drift   = ann_off_epoch > 0;
  const bool warn_multi_stage   = ! is_generic && ! is_contig && ep_multi_stage > 0;
  const bool warn_gap_anns      = ( anns_in_gap + anns_span_gap ) > 0;
  const bool warn_nonstandard   = nonstandard_epoch;   // informational, not a hard error

  const bool any_warning =
    warn_rec_not_div || warn_epoch_off_rec || warn_stage_drift ||
    warn_multi_stage || warn_gap_anns;


  // ------------------------------------------------------------------
  // 14. Console output
  // ------------------------------------------------------------------

  const std::string SEP2 = "  =======================================================================";

  logger << "\n";

  // --- EDF structure ---
  logger << "\n  EDF structure\n";
  logger << "   * type        : " << edf_type << "\n";
  if ( ! is_continuous )
    logger << "                 " << nsegs << " segments, " << ngaps << " gaps\n";
  logger << "   * record size : " << Helper::dbl2str(rec_dur_sec) << "s"
         << "  (" << nrecs << " records, "
         << Helper::dbl2str(total_dur_sec / 60.0 , 1) << " min)\n";

  // --- Epoch definition ---
  logger << "\n  Epoch definition";
  if ( auto_epoched ) logger << "  [none set - applied default 30s]";
  logger << "\n";
  logger << "   * " << Helper::dbl2str(epoch_dur_sec) << "s epochs"
         << ", step " << Helper::dbl2str(epoch_inc_sec) << "s"
         << ", offset " << Helper::dbl2str(epoch_off_sec) << "s"
         << "  [" << nepochs << " epochs]\n";
  if ( is_generic )
    logger << "   * mode: GENERIC (variable-length, annotation-derived)\n"
           << "     -> record-boundary and per-epoch stage checks are not applicable\n";
  else if ( is_contig )
    logger << "   * mode: CONTIG (annotation-contiguous blocks)\n"
           << "     -> record-boundary and per-epoch stage checks are not applicable\n";
  else if ( is_gap_spanning )
    logger << "   * mode: gap-spanning (epochs deliberately cross segment boundaries)\n"
           << "     -> record-boundary alignment and MASK/RE ambiguity checks are skipped\n";
  else if ( is_overlapping )
    logger << "   * mode: overlapping epochs (step " << Helper::dbl2str(epoch_inc_sec)
           << "s < dur " << Helper::dbl2str(epoch_dur_sec) << "s)\n"
           << "     -> annotation-to-epoch mapping is intentionally many-to-many\n";

  if ( align_active )
    logger << "   * annotation alignment ACTIVE (align=" << align_str << ")\n";
  else
    logger << "   * annotation alignment: none\n";

  // --- Staging annotations ---
  logger << "\n  Staging annotations (" << (int)ann_labels.size() << " labels requested)\n";
  if ( found_labels.empty() )
    {
      logger << "   * NONE of the requested labels are present in this file\n";
    }
  else
    {
      std::ostringstream oss;
      for ( const auto & l : found_labels )
        oss << l << "(" << ann_counts[l] << ") ";
      logger << "   * found  : " << oss.str()
             << "  [" << total_ann_instances << " instances, "
             << total_ann_starts << " unique start times]\n";
    }
  if ( ! missing_labels.empty() )
    {
      std::ostringstream oss;
      for ( const auto & l : missing_labels ) oss << l << " ";
      logger << "   * absent : " << oss.str() << "\n";
    }

  // --- Epoch → record alignment ---
  if ( ! is_generic && ! is_contig )
    {
      logger << "\n  Epoch -> record alignment\n";
      bool ok_div = rec_divides_epoch;
      bool ok_bnd = ( epochs_off_rec == 0 );

      logger << "   * epoch (" << Helper::dbl2str(epoch_dur_sec)
             << "s) exact multiple of record (" << Helper::dbl2str(rec_dur_sec) << "s) : "
             << ( ok_div ? "YES" : "NO  *** WARN ***" ) << "\n";
      if ( ! ok_div )
        logger << "     -> epoch boundaries will always cross record boundaries;\n"
               << "       RE masking at record granularity will be imprecise\n";

      if ( is_gap_spanning )
        logger << "   * epoch-start-on-record-boundary: skipped (gap-spanning mode)\n";
      else
        {
          logger << "   * epoch starts on record boundaries : "
                 << frac(epochs_on_rec, nepochs) << " (" << pct(epochs_on_rec,nepochs) << ")"
                 << ( ok_bnd ? "  OK\n" : "  *** WARN ***\n" );
          if ( ! ok_bnd )
            logger << "     -> epoch grid is offset from the record grid\n";
        }
    }

  // --- Annotation → epoch alignment ---
  logger << "\n  Annotation start -> epoch grid\n";
  if ( total_ann_starts == 0 )
    {
      logger << "   * no staging annotations found; skipping\n";
    }
  else
    {
      bool ok = ( ann_off_epoch == 0 );
      logger << "   * annotation starts on an epoch boundary : "
             << frac(ann_at_epoch, total_ann_starts)
             << " (" << pct(ann_at_epoch, total_ann_starts) << ")"
             << ( ok ? "  OK\n" : "  *** WARN ***\n" );

      if ( is_overlapping && ann_off_epoch > 0 )
        logger << "     -> overlapping epoch mode: some drift is expected by design\n";

      if ( ! ok )
        {
          logger << "   * " << ann_off_epoch
                 << " annotation start(s) fall between epoch boundaries\n"
                 << "     -> those stages will not map 1:1 to a single epoch\n";

          if ( align_active )
            {
              logger << "     -> annotation alignment is already active"
                     << " (EPOCH align was run)\n";
              if ( ann_off_epoch > 0 )
                logger << "     -> remaining drift is likely due to stages outside the\n"
                       << "       retained data (masked-out epochs now in gaps)\n";
            }
          else if ( suggest_epoch_align )
            {
              logger << "\n   * SUGGESTION: run  EPOCH align  before this command\n"
                     << "     -> would align " << align_would_fix_n << "/" << total_ann_starts
                     << " (" << pct(align_would_fix_n,total_ann_starts) << ") annotation starts"
                     << " to epoch boundaries\n"
                     << "     -> suggested offset: " << Helper::dbl2str(align_offset_sec,1) << "s"
                     << " (from first annotation start)\n"
                     << "     -> example: luna ... -s 'EPOCH align & MASK ifnot=N2 & RE'\n";
              if ( ! align_would_fix_all )
                logger << "     -> partial fix only: "
                       << (total_ann_starts - align_would_fix_n)
                       << " start(s) cannot be aligned (check for sub-epoch-length annotations)\n";
            }
          else if ( ! align_active && align_offset_tp == 0 )
            logger << "     -> no annotation start found to anchor EPOCH align\n";
        }
    }

  // --- Per-epoch stage coverage ---
  if ( ! is_generic && ! is_contig )
    {
      logger << "\n  Per-epoch stage coverage\n";
      logger << "   * single unambiguous stage : " << ep_one_stage   << " epoch(s)\n";
      logger << "   * no stage (unscored)      : " << ep_no_stage    << " epoch(s)\n";
      logger << "   * spans >1 stage label     : " << ep_multi_stage << " epoch(s)";
      if ( ep_multi_stage > 0 )
        logger << "  *** WARN: stage boundary falls mid-epoch ***";
      logger << "\n";
    }

  // --- MASK/RE record ambiguity ---
  if ( ! is_generic && ! is_contig && ! is_gap_spanning )
    {
      logger << "\n  MASK/RE record ambiguity\n";
      if ( records_ambig == 0 )
        logger << "   * no records shared between differently-staged epochs : OK\n";
      else
        logger << "   * " << records_ambig
               << " record(s) shared between epochs with different stages  *** WARN ***\n"
               << "     -> RE masking by stage may retain signal from the wrong stage\n";
    }

  // --- EDF+D gap issues ---
  if ( ! is_continuous )
    {
      logger << "\n  EDF+D gap issues\n";
      if ( anns_in_gap == 0 && anns_span_gap == 0 )
        logger << "   * no staging annotations fall in or span a gap : OK\n";
      else
        {
          if ( anns_in_gap > 0 )
            logger << "   * " << anns_in_gap
                   << " annotation(s) entirely within a gap (signal absent)  *** WARN ***\n";
          if ( anns_span_gap > 0 )
            logger << "   * " << anns_span_gap
                   << " annotation(s) span a segment/gap boundary  *** WARN ***\n"
                   << "     -> partial signal only; RE will retain an incomplete epoch\n";
        }
    }

  // --- Per-label summary ---
  if ( ! found_labels.empty() && ! is_generic && ! is_contig )
    {
      logger << "\n  Per-label detail\n";
      // header
      logger << "\n  "
             << "  label"
             << "   total"
             << "  on-epoch"
             << "  off-epoch"
             << "  mean-off(s)"
             << "  max-off(s)"
             << "  ep-single"
             << "  ep-multi\n"
             << "  " << std::string(75,'-') << "\n";

      for ( const auto & l : found_labels )
        {
          const LabelStats & ls = lstats.at(l);
          int n = ann_counts.at(l);
          double mean_off = ls.n_off_epoch > 0
                            ? ls.sum_off_tp / ls.n_off_epoch
                            : 0.0;

          logger << "  "
                 << "  " << l
                 << std::string( std::max(1, 8 - (int)l.size() ) , ' ' )
                 << Helper::int2str(n)
                 << std::string( std::max(1, 8 - (int)Helper::int2str(n).size()) , ' ' )
                 << Helper::int2str(ls.n_at_epoch)
                 << std::string( std::max(1, 10 - (int)Helper::int2str(ls.n_at_epoch).size()) , ' ' )
                 << Helper::int2str(ls.n_off_epoch)
                 << std::string( std::max(1, 11 - (int)Helper::int2str(ls.n_off_epoch).size()) , ' ' )
                 << Helper::dbl2str(mean_off, 2)
                 << std::string( std::max(1, 13 - (int)Helper::dbl2str(mean_off,2).size()) , ' ' )
                 << Helper::dbl2str(ls.max_off_tp, 2)
                 << std::string( std::max(1, 12 - (int)Helper::dbl2str(ls.max_off_tp,2).size()) , ' ' )
                 << Helper::int2str(ls.ep_single)
                 << std::string( std::max(1, 11 - (int)Helper::int2str(ls.ep_single).size()) , ' ' )
                 << Helper::int2str(ls.ep_multi)
                 << "\n";
        }
    }

  // --- Overall summary ---
  logger << "\n" << SEP2 << "\n";
  if ( ! any_warning && ! warn_nonstandard )
    logger << "  OVERALL: OK - no alignment issues detected\n";
  else
    {
      logger << "  OVERALL: " << ( any_warning ? "WARNINGS DETECTED" : "NOTE" ) << "\n";
      if ( warn_nonstandard   ) logger << "   [NOTE] non-standard epoch mode - some checks were skipped\n";
      if ( warn_rec_not_div   ) logger << "   [WARN] epoch size not an exact multiple of record size\n";
      if ( warn_epoch_off_rec ) logger << "   [WARN] some epoch starts are not on record boundaries\n";
      if ( warn_stage_drift   ) logger << "   [WARN] some annotation starts fall between epoch boundaries\n";
      if ( suggest_epoch_align) logger << "          -> suggest: EPOCH align\n";
      if ( warn_multi_stage   ) logger << "   [WARN] some epochs span >1 stage label\n";
      if ( warn_gap_anns      ) logger << "   [WARN] some annotations fall in or span EDF+D gaps\n";
    }
  logger << SEP2 << "\n\n";


  // ------------------------------------------------------------------
  // 15. Optional issue annotations for visualization
  // ------------------------------------------------------------------

  if ( emit_annots || emit_epoch_annots || emit_rec_annots )
    {
      if ( emit_annots )
        {
          annot_t * conf_annot      = edf.annotations->add( annot_prefix + "CONF" );
          annot_t * off_epoch_annot = edf.annotations->add( annot_prefix + "OFF_EPOCH" );
          annot_t * in_gap_annot    = edf.annotations->add( annot_prefix + "IN_GAP" );
          annot_t * span_gap_annot  = edf.annotations->add( annot_prefix + "SPAN_GAP" );
          annot_t * rec_ambig_annot = edf.annotations->add( annot_prefix + "REC_AMBIG" );

          for ( int e = 0 ; e < nepochs ; ++e )
            {
              const std::set<std::string> & labs = epoch_stage_labels[e];
              if ( labs.size() <= 1 ) continue;

              const interval_t ep = edf.timeline.epoch(e);
              instance_t * inst = conf_annot->add( "." , ep , "." );
              inst->set( "TYPE" , "CONF" );
              inst->set( "EPOCH" , e + 1 );
              inst->set( "N_STAGES" , (int)labs.size() );
              inst->set( "STAGE" , epoch_stage_str[e] );
            }

          for ( const auto & info : ann_issue_info )
            {
              if ( ! info.on_epoch )
                {
                  instance_t * inst = off_epoch_annot->add( "." , info.interval , "." );
                  inst->set( "TYPE" , "OFF_EPOCH" );
                  inst->set( "STAGE" , info.label );
                  inst->set( "OFF_SEC" , info.off_sec );
                }

              if ( ! is_continuous )
                {
                  for ( const auto & gap : gps )
                    {
                      const bool a_in_gap =
                        ( info.interval.start >= gap.start && info.interval.stop <= gap.stop );
                      const bool a_span_left =
                        ( info.interval.start < gap.start && info.interval.stop > gap.start && info.interval.stop <= gap.stop );
                      const bool a_span_right =
                        ( info.interval.start >= gap.start && info.interval.start < gap.stop && info.interval.stop > gap.stop );

                      if ( a_in_gap )
                        {
                          instance_t * inst = in_gap_annot->add( "." , info.interval , "." );
                          inst->set( "TYPE" , "IN_GAP" );
                          inst->set( "STAGE" , info.label );
                          break;
                        }
                      if ( a_span_left || a_span_right )
                        {
                          instance_t * inst = span_gap_annot->add( "." , info.interval , "." );
                          inst->set( "TYPE" , "SPAN_GAP" );
                          inst->set( "STAGE" , info.label );
                          break;
                        }
                    }
                }
            }

          if ( ! is_generic && ! is_contig && ! is_gap_spanning )
            {
              for ( const auto & kv : rec_stage_str )
                {
                  std::map<int,uint64_t>::const_iterator rt = edf.timeline.rec2tp.find( kv.first );
                  if ( rt == edf.timeline.rec2tp.end() ) continue;
                  interval_t rec_iv( rt->second , rt->second + rec_dur_tp );
                  instance_t * inst = rec_ambig_annot->add( "." , rec_iv , "." );
                  inst->set( "TYPE" , "REC_AMBIG" );
                  inst->set( "REC" , kv.first );
                  inst->set( "STAGE" , kv.second );
                }
            }
        }

      if ( emit_epoch_annots )
        {
          annot_t * epoch_annot = edf.annotations->add( annot_prefix + "EPOCH" );
          for ( int e = 0 ; e < nepochs ; ++e )
            {
              const interval_t ep = edf.timeline.epoch(e);
              instance_t * inst = epoch_annot->add( "." , ep , "." );
              inst->set( "TYPE" , "EPOCH" );
              inst->set( "EPOCH" , e + 1 );
              if ( ! is_generic && ! is_contig )
                {
                  const std::set<std::string> & labs = epoch_stage_labels[e];
                  inst->set( "N_STAGES" , (int)labs.size() );
                  inst->set( "STAGE" , epoch_stage_str[e] );
                }
              inst->set( "ON_REC" , (int)( rec_starts_set.count( ep.start ) > 0 ) );
            }
        }

      if ( emit_rec_annots )
        {
          annot_t * rec_annot = edf.annotations->add( annot_prefix + "REC" );
          for ( const auto & kv : edf.timeline.rec2tp )
            {
              interval_t rec_iv( kv.second , kv.second + rec_dur_tp );
              instance_t * inst = rec_annot->add( "." , rec_iv , "." );
              inst->set( "TYPE" , "REC" );
              inst->set( "REC" , kv.first );
            }
        }
    }


  // ------------------------------------------------------------------
  // 16. Structured writer() output — summary level
  // ------------------------------------------------------------------

  writer.value( "OK"               , (int)( ! any_warning ) );

  writer.value( "EDF_TYPE"         , edf_type );
  writer.value( "REC_DUR"          , rec_dur_sec );
  writer.value( "NRECS"            , nrecs );
  writer.value( "TOT_DUR"          , total_dur_sec );
  writer.value( "NSEGS"            , nsegs );
  writer.value( "NGAPS"            , ngaps );

  writer.value( "EPOCH_DUR"            , epoch_dur_sec );
  writer.value( "EPOCH_INC"            , epoch_inc_sec );
  writer.value( "EPOCH_OFFSET"         , epoch_off_sec );
  writer.value( "EPOCH_N"              , nepochs );
  writer.value( "EPOCH_AUTO"           , (int)auto_epoched );
  writer.value( "EPOCH_ALIGN_ACTIVE"   , (int)align_active );
  writer.value( "EPOCH_GENERIC"        , (int)is_generic );
  writer.value( "EPOCH_CONTIG"         , (int)is_contig );
  writer.value( "EPOCH_OVERLAPPING"    , (int)is_overlapping );
  writer.value( "EPOCH_GAP_SPANNING"   , (int)is_gap_spanning );
  writer.value( "EPOCH_NONSTANDARD"    , (int)nonstandard_epoch );

  {
    std::ostringstream oss1, oss2;
    for ( const auto & l : found_labels )   oss1 << l << " ";
    for ( const auto & l : missing_labels ) oss2 << l << " ";
    writer.value( "ANN_FOUND"   , oss1.str() );
    writer.value( "ANN_MISSING" , oss2.str() );
  }
  writer.value( "ANN_N"      , total_ann_instances );
  writer.value( "ANN_STARTS" , total_ann_starts );

  writer.value( "REC_DIVIDES_EPOCH" , (int)rec_divides_epoch );
  writer.value( "EPOCH_ON_REC"      , epochs_on_rec );
  writer.value( "EPOCH_OFF_REC"     , epochs_off_rec );
  writer.value( "EPOCH_ON_REC_PCT"  ,
                nepochs > 0 ? 100.0 * epochs_on_rec / nepochs : 0.0 );

  writer.value( "ANN_AT_EPOCH"          , ann_at_epoch );
  writer.value( "ANN_OFF_EPOCH"         , ann_off_epoch );
  writer.value( "ANN_AT_EPOCH_PCT"      ,
                total_ann_starts > 0 ? 100.0 * ann_at_epoch / total_ann_starts : 0.0 );
  writer.value( "ALIGN_OFFSET"          , align_offset_sec );
  writer.value( "ALIGN_WOULD_FIX_N"     , align_would_fix_n );
  writer.value( "ALIGN_WOULD_FIX_PCT"   ,
                total_ann_starts > 0 ? 100.0 * align_would_fix_n / total_ann_starts : 0.0 );
  writer.value( "ALIGN_WOULD_FIX_ALL"   , (int)align_would_fix_all );
  writer.value( "SUGGEST_EPOCH_ALIGN"   , (int)suggest_epoch_align );

  writer.value( "EPOCH_ONE_STAGE"   , ep_one_stage );
  writer.value( "EPOCH_NO_STAGE"    , ep_no_stage );
  writer.value( "EPOCH_MULTI_STAGE" , ep_multi_stage );

  writer.value( "REC_AMBIG" , records_ambig );

  writer.value( "ANN_IN_GAP"   , anns_in_gap );
  writer.value( "ANN_SPAN_GAP" , anns_span_gap );

  writer.value( "WARN_REC_NOT_DIV"   , (int)warn_rec_not_div );
  writer.value( "WARN_EPOCH_OFF_REC" , (int)warn_epoch_off_rec );
  writer.value( "WARN_STAGE_DRIFT"   , (int)warn_stage_drift );
  writer.value( "WARN_MULTI_STAGE"   , (int)warn_multi_stage );
  writer.value( "WARN_GAP_ANNS"      , (int)warn_gap_anns );


  // ------------------------------------------------------------------
  // 17. Per-label table  (strata: ANNOT)
  //     Skipped for generic/contig modes where per-epoch stats N/A
  // ------------------------------------------------------------------

  if ( ! found_labels.empty() )
    {
      for ( const auto & l : found_labels )
        {
          writer.level( l , "ANNOT" );

          const LabelStats & ls = lstats.at(l);
          int n = ann_counts.at(l);
          double mean_off = ls.n_off_epoch > 0
                            ? ls.sum_off_tp / ls.n_off_epoch
                            : 0.0;

          writer.value( "ANN_N"          , n );
          writer.value( "ANN_AT_EPOCH"   , ls.n_at_epoch );
          writer.value( "ANN_OFF_EPOCH"  , ls.n_off_epoch );
          writer.value( "ANN_AT_EPOCH_PCT",
                        n > 0 ? 100.0 * ls.n_at_epoch / n : 0.0 );
          writer.value( "MEAN_OFF_SEC"   , mean_off );
          writer.value( "MAX_OFF_SEC"    , ls.max_off_tp );

          // epoch coverage (only meaningful for standard epochs)
          if ( ! is_generic && ! is_contig )
            {
              writer.value( "EPOCH_SINGLE"  , ls.ep_single );
              writer.value( "EPOCH_MULTI"   , ls.ep_multi );
            }

          writer.unlevel( "ANNOT" );
        }
    }


  // ------------------------------------------------------------------
  // 18. Verbose per-epoch table  (strata: E)
  //     Only emitted when verbose param is set.
  //     Not applicable for generic/contig modes.
  // ------------------------------------------------------------------

  const bool do_verbose = param.has( "verbose" );

  if ( do_verbose && ! is_generic && ! is_contig )
    {
      clocktime_t starttime( edf.header.starttime );
      const bool valid_hms = starttime.valid;

      for ( int e = 0 ; e < nepochs ; ++e )
        {
          writer.level( e + 1 , "E" );

          const interval_t ep = edf.timeline.epoch(e);

          // Times in seconds from EDF start
          const double start_sec = ep.start / (double)globals::tp_1sec;
          const double stop_sec  = ep.stop  / (double)globals::tp_1sec;

          writer.value( "START" , start_sec );
          writer.value( "STOP"  , stop_sec  );

          if ( valid_hms )
            {
              clocktime_t t1 = starttime;  t1.advance_seconds( start_sec );
              clocktime_t t2 = starttime;  t2.advance_seconds( stop_sec  );
              writer.value( "START_HMS" , t1.as_string( ':' , true ) );
              writer.value( "STOP_HMS"  , t2.as_string( ':' , true ) );
            }

          // Stage label(s)
          const std::set<std::string> & labs = epoch_stage_labels[e];
          std::string stage_str;
          if ( labs.empty() ) stage_str = ".";
          else for ( const auto & l : labs ) stage_str += l + " ";
          writer.value( "STAGE"    , stage_str );
          writer.value( "N_STAGES" , (int)labs.size() );

          // On a record boundary?
          writer.value( "ON_REC" , (int)( rec_starts_set.count( ep.start ) > 0 ) );

          // First / last record spanned by this epoch
          int ra = -1 , rb = -1;
          edf.timeline.epoch_records( e , &ra , &rb );
          writer.value( "REC_START" , ra );
          writer.value( "REC_STOP"  , rb );
          if ( ra >= 0 && rb >= 0 )
            writer.value( "NRECS" , rb - ra + 1 );

          writer.unlevel( "E" );
        }
    }


  // ------------------------------------------------------------------
  // 19. Verbose per-annotation-event table  (strata: ANNOT x ANN_N)
  //     Only emitted when verbose param is set.
  // ------------------------------------------------------------------

  if ( do_verbose && ! found_labels.empty() )
    {
      clocktime_t starttime( edf.header.starttime );
      const bool valid_hms = starttime.valid;

      // Sort events by start time for consistent ordering
      std::vector<std::pair<interval_t,std::string>> sorted_anns = all_anns;
      std::sort( sorted_anns.begin() , sorted_anns.end() ,
                 []( const auto & a , const auto & b )
                 { return a.first.start < b.first.start; } );

      // Per-label instance counter (1-based)
      std::map<std::string,int> label_idx;
      for ( const auto & l : found_labels ) label_idx[l] = 0;

      for ( const auto & ai : sorted_anns )
        {
          const interval_t & iv    = ai.first;
          const std::string & label = ai.second;
          const int idx = ++label_idx[ label ];

          writer.level( label , "ANNOT" );
          writer.level( idx   , "ANN_N" );

          const double start_sec = iv.start / (double)globals::tp_1sec;
          const double stop_sec  = iv.stop  / (double)globals::tp_1sec;

          writer.value( "START" , start_sec );
          writer.value( "STOP"  , stop_sec  );

          if ( valid_hms )
            {
              clocktime_t t1 = starttime;  t1.advance_seconds( start_sec );
              clocktime_t t2 = starttime;  t2.advance_seconds( stop_sec  );
              writer.value( "START_HMS" , t1.as_string( ':' , true ) );
              writer.value( "STOP_HMS"  , t2.as_string( ':' , true ) );
            }

          // On an epoch boundary?
          const bool at_ep = epoch_starts_set.count( iv.start ) > 0;
          writer.value( "AT_EPOCH" , (int)at_ep );

          // Which epoch (1-based), -1 if not on a boundary
          int epoch_num = -1;
          if ( at_ep )
            {
              for ( int e = 0 ; e < nepochs ; ++e )
                {
                  if ( edf.timeline.epoch(e).start == iv.start )
                    { epoch_num = e + 1; break; }
                }
            }
          writer.value( "EPOCH" , epoch_num );

          // Offset (seconds) from nearest epoch boundary
          double off_sec = 0.0;
          if ( ! at_ep && ! epoch_starts_vec.empty() )
            {
              auto it = std::lower_bound( epoch_starts_vec.begin() ,
                                          epoch_starts_vec.end() , iv.start );
              uint64_t nearest = std::numeric_limits<uint64_t>::max();
              if ( it != epoch_starts_vec.end() )
                nearest = std::min( nearest , *it - iv.start );
              if ( it != epoch_starts_vec.begin() )
                { --it; nearest = std::min( nearest , iv.start - *it ); }
              if ( nearest != std::numeric_limits<uint64_t>::max() )
                off_sec = nearest / (double)globals::tp_1sec;
            }
          writer.value( "OFF_SEC" , off_sec );

          // First / last retained record overlapping this annotation
          // Use tp2rec (map<uint64_t,int>: record_start_tp -> record_id)
          // to find bounding records without scanning all records.
          int rec_first = -1 , rec_last = -1;
          if ( ! edf.timeline.tp2rec.empty() )
            {
              // Last record starting at or before iv.start
              auto lo = edf.timeline.tp2rec.upper_bound( iv.start );
              if ( lo != edf.timeline.tp2rec.begin() )
                { --lo; rec_first = lo->second; }

              // Last record starting strictly before iv.stop
              auto hi = edf.timeline.tp2rec.lower_bound( iv.stop );
              if ( hi != edf.timeline.tp2rec.begin() )
                { --hi; rec_last = hi->second; }

              // If annotation is shorter than one record, first == last
              if ( rec_last < rec_first ) rec_last = rec_first;
            }
          writer.value( "REC_START" , rec_first );
          writer.value( "REC_STOP"  , rec_last  );

          writer.unlevel( "ANN_N"  );
          writer.unlevel( "ANNOT" );
        }
    }

}
