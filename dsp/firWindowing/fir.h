
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

#include <vector>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include "../../depends/include/fftw3.h"
#include <complex>


struct fir_t
{
  
  // class for FIR design, using implementations from "FIR filters by
  // Windowing" by A.Greensted - Feb 2010
  // http://www.labbookpages.co.uk
  
  enum filterType { LOW_PASS, HIGH_PASS, BAND_PASS, BAND_STOP };
  enum windowType { RECTANGULAR, BARTLETT, HANNING, HAMMING, BLACKMAN };
  
  // Prototypes
  std::vector<double> create1TransSinc( int windowLength, double transFreq, double sampFreq, enum filterType type);
  std::vector<double> create2TransSinc( int windowLength, double trans1Freq, double trans2Freq, double sampFreq, enum filterType type);
  
  std::vector<double> createWindow(double *in, double *out, int windowLength, enum windowType type);

  void calculateKaiserParams(double ripple, double transWidth, double sampFreq, int *windowLength, double *beta);
  double *createKaiserWindow(double *in, double *out, int windowLength, double beta);
  double modZeroBessel(double x);

  int outputFFT(char *filename, double *window, int windowLength, double sampFreq);

  void demo();

};

