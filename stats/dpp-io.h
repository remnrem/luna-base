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

#ifndef __LUNA_DPP_IO_H__
#define __LUNA_DPP_IO_H__

#include <string>
#include <vector>
#include <fstream>

// Per-individual DPP feature-matrix binary I/O -- same headerless,
// self-delimiting, cat-concatenable technique as pops/io.cpp (bwrite:
// 1-byte length + string; ints/doubles as raw fixed-size fields; no
// whole-file header; loop-until-EOF reader). No hypnodensity or
// cycle-context columns (out of scope for DPP itself, see plan Stage 2 SS4/9).
//
// Record layout, written once per individual (and simply appended for
// each subsequent individual -- this is why 'cat a.dat b.dat > all.dat'
// already works, exactly as for POPS's own format):
//   ID (string) | N_ROWS (int) | N_FEATURES (int) |
//   N_ROWS x [ TIME (double) | N_FEATURES x VALUE (double, NaN = missing) ]

struct dpp_matrix_t
{
  std::string id;
  std::vector<double> time_sec;              // one per row
  std::vector<std::vector<double> > X;        // [row][col]
};

namespace dpp_io {

  // append (or create) one individual's matrix in 'f'; 'n_features' is
  // passed explicitly (not inferred from X) so an individual with zero
  // valid rows can still be written with the correct schema
  void save( const std::string & f , const dpp_matrix_t & m , int n_features , bool append );

  // read all individuals from 'f'; halts if any individual's column count
  // doesn't match 'expected_n_features' (pass -1 to skip the check)
  std::vector<dpp_matrix_t> load( const std::string & f , int expected_n_features = -1 );
  std::vector<dpp_matrix_t> load_files( const std::vector<std::string> & files , int expected_n_features = -1 );

  // low-level primitives (mirrors pops_indiv_t::bwrite/bread_*)
  void bwrite( std::ofstream & O , const std::string & s );
  void bwrite( std::ofstream & O , int i );
  void bwrite( std::ofstream & O , double d );
  std::string bread_str( std::ifstream & I );
  int bread_int( std::ifstream & I );
  double bread_dbl( std::ifstream & I );

}

#endif
