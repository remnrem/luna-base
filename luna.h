
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

#ifndef __LUNA_H__
#define __LUNA_H__

#include <cstddef>

#include "eval.h" 
#include "dummy.h"
#include "cmddefs.h"

#include "defs/defs.h"
#include "helper/helper.h"
#include "helper/zfile.h"
#include "helper/token-eval.h"
#include "helper/token.h"
#include "helper/logger.h"
#include "helper/xml-parser.h"
#include "helper/json.h"
#include "helper/mapper.h"
#include "helper/validate.h"

#include "tinyxml/tinyxml.h"
#include "tinyxml/xmlreader.h"

#include "stats/Eigen/Dense"
#include "stats/eigen_ops.h"
#include "stats/cluster.h"
#include "stats/kmeans.h"
#include "stats/cpt.h"
#include "stats/nmf.h"
#include "stats/gpa.h"

#include "miscmath/miscmath.h"
#include "miscmath/dynam.h"
#include "miscmath/crandom.h"

#include "stats/glm.h"
#include "stats/lda.h"
#include "stats/qda.h"

#include "dynamics/qdynam.h"
#include "dynamics/e-cycles.h"

#include "annot/annot.h"
#include "annot/annotate.h"
#include "annot/align-annots.h"

#include "edf/edf.h"
#include "edfz/edfz2.h"
#include "edf/canonical.h"
#include "edf/slice.h"
#include "edf/sedf.h"
#include "edf/freezer.h"
#include "edf/align-epochs.h"
#include "edf/insert.h"

#include "timeline/timeline.h"

#include "annot/nsrr-remap.h"

#include "dsp/dsp.h"
#include "dsp/libsamplerate/samplerate.h"

#include "spectral/welch.h"
#include "spectral/irasa.h"
#include "spectral/mtm/mtm.h"

#include "fftw/fftwrap.h"

#include "ica/ica.h"

#include "spindles/spindles.h"
#include "spindles/slow-waves.h"

#include "artifacts/artifacts.h"
#include "artifacts/correct.h"

#include "cwt/cwt.h"

#include "pdc/pdc.h"

#include "staging/staging.h"

#include "suds/suds.h"

#include "lgbm/lgbm.h"

#include "pops/pops.h"

#include "assoc/assoc.h"
#include "assoc/massoc.h"

#include "clocs/topo.h"

#include "clocs/clocs.h"

#include "db/db.h"

#include "db/retval.h"

#include "sstore/sstore.h"

#include "resp/hb.h"

#include "models/predict.h"

#include "lunapi/lunapi.h"
#include "lunapi/segsrv.h"

#include <iostream>
#include <unistd.h>
#include <dirent.h>

extern globals global;
extern writer_t writer;
extern logger_t logger;
extern freezer_t freezer;

#endif
