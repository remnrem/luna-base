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

#ifndef __COMBINE_EMG_H__
#define __COMBINE_EMG_H__

struct edf_t;
struct param_t;

// COMBINE-EMG: build a single continuous EMG channel from 2+ candidate
// channels (raw or already bipolar-referenced) by scoring short windows
// for artifact/quality and switching to the best available candidate,
// with hysteresis/crossfade so the output doesn't itself look artifactual.
struct combine_emg_t {
  combine_emg_t( edf_t & edf , param_t & param );
};

#endif
