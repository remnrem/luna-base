
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

// this file directly includes the GFT library by Robert Brown (see note below)

#ifndef __LUNA_GFT_H__
#define __LUNA_GFT_H__

namespace dsptools {
  
  struct gft_t {
    
  };
  
}


/*
 *  gft.h
 *  GFT Framework
 *
 *  Created by Robert Brown on 30/05/08.
 *  Copyright 2008 Robert Brown. All rights reserved.
 *
 */

/*
  
  Notice

  FST (Fast S-Transform) Software, Copyright © 2010 UTI Limited
  Partnership, Original Authors: Robert A. Brown, M. Louis Lauzon,
  Richard Frayne.

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License version 3, as
  published by the Free Software Foundation, subject to the Additional
  Terms set forth below.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License v3 and the Additional Terms set forth below
  for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see
  <http://www.gnu.org/licenses/>.

  For general publications, we request that you cite: Brown, R.A.,
  Lauzon, M.L., and Frayne, R., ìA General Description of Linear
  Time-Frequency Transforms and Formulation of a Fast, Invertible
  Transform That Samples the Continuous S-Transform Spectrum
  Nonredundantly,î IEEE Trans. Signal Process. vol. 58, pp. 281-290,
  Jan 2010.

  Non-free versions of FST are available under terms different from
  those of the GNU General Public License. For these alternative terms
  you must purchase a license from University Technologies
  International.  Users interested in such a license should contact
  UTI (fst@uti.ca) for more information.

  Additional Terms

  1. Disclaimer of Warranty

  This disclaimer supplements the one included in the General Public
  License.  The Program is provided ìAS ISî without warranties,
  conditions or representations of any kind, and the authors, UTI and
  the University of Calgary expressly disclaim, to the fullest extent
  permitted by applicable law, any warranty or condition, express or
  implied, statutory or otherwise, whether arising from trade or
  course of dealing, including, without limitation, any warranty that
  the Program shall correspond with a particular description, is
  error-free, or that the Program will be compatible with third-party
  software or that any errors in the Program will be corrected, or
  that the Program does not and will not infringe any patent,
  trade-mark, trade-secret or other intellectual property or other
  proprietary rights of any third party.  No oral or written advice
  provided by the authors, UTI or the University of Calgary or any
  authorized representative shall create a warranty.  You shall make
  no statements, representations or warranties whatsoever to any third
  party which are inconsistent with the foregoing disclaimer.

  2. Limitation of Liability

  This limitation supplements the one included in the General Public
  License.  The total liability of the authors, UTI and the University
  of Calgary, whether under the express or implied terms of the
  General Public License, in tort (including negligence or other duty
  of care) or at common law, for any loss or damage suffered by you or
  your end-users or any third parties to whom you convey the Program
  or any modifications, whether direct, indirect or special, or any
  other similar damage that may arise or does arise, is limited to the
  amount (if any) paid by you to download the Program.

  3. Legal Notices and Author Attribution

  You must acknowledge the authors in any conveyance or propagation of
  the Program.  All copyright notices must be retained and must
  include author attributions: ìCopyright © 2010 UTI Limited
  Partnership, Original Authors: Robert A. Brown, M. Louis Lauzon,
  Richard Frayneî

  4. Trade-mark or Naming Rights

  No trade-mark or publicity rights are granted.  This license does
  NOT give you any right, title or interest in any marks of UTI of the
  University of Calgary.  You may not convey any modification of this
  Program using any marks of UTI or the University of Calgary or claim
  any affiliation or association with the authors, UTI, or the
  University of Calgary.

  You may not misrepresent the origins of this Program.  Modified
  versions of the Program must be clearly marked as such and not
  identified as the original program.

  5. Indemnification

  You agree to indemnify and hold harmless the authors, UTI, the
  author(s), the University of Calgary, their partners, Board of
  Governors, officers, employees, faculty, students, invitees and
  agents and their successors and assigns (collectively, the
  ìIndemnified Partiesî), for any claims, damages, liability or
  losses, consequential or otherwise, (including all associated legal
  fees and disbursements) incurred by the Indemnified Parties arising
  in any manner at all out of the use of the Program (or any
  modifications of it) by you, your end users or any third-party; or
  if you convey this Program (or any modifications of it) and assume
  contractual liability for the Program to recipients of it, you agree
  to indemnify, and hold harmless the Indemnified Parties, for any
  claims, damages, liability or losses, consequential or otherwise,
  (including all associated legal fees and disbursements) that those
  contractual assumptions impose on the Indemnified Parties.

  END OF ADDITIONAL TERMS.
*/


typedef void (gft_windowFunction)(double*,int,int);

// Windows
void gft_gaussian(double *win, int N, int freq);
void gft_box(double *win, int N, int freq);

// GFT partition and window screen generators
int gft_1dSizeOfPartitions(unsigned int N);
int *gft_1dPartitions(unsigned int N);
int *gft_1dMusicPartitions(unsigned int N, float samplerate, int cents);
double *gft_windows(int N, gft_windowFunction *window);
double *gft_windowsFromPars(int N, gft_windowFunction *window, int *pars);

// 1D GFT Functions
void gft_1dComplex64(double *signal, unsigned int N, double *win, int *pars, int stride);

// 2D GFT Functions
void gft_2dComplex64(double *signal, unsigned int N, unsigned int M, gft_windowFunction *window);

// Interpolation functions
double *gft_1d_interpolateNN(double *signal, unsigned int N, unsigned int M);

#endif
