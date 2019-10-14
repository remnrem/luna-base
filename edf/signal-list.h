
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

#ifndef __LUNA_SIGLIST_H__
#define __LUNA_SIGLIST_H__

struct signal_list_t
{
    
  std::vector<int> signals;  

  std::vector<std::string> signal_labels;  

  int size() const { return signals.size(); } 

  int operator()(int i) const { return signals[i]; } 

  std::string label(const int i) const { return signal_labels[i]; } 
  
  int find( const std::string & label ) const
  { 
    for (int j=0;j<signal_labels.size();j++) 
      if ( signal_labels[j] == label ) return j;
    return -1;
  }

  void add(int i , const std::string & label ) 
  { 
    //    std::cout << "adding " << i << " to " << label << "\n";
    for (int j=0;j<signals.size();j++) if ( signals[j] == i ) return;
    signals.push_back(i); 
    signal_labels.push_back(label); 
  }

  void clear() { signals.clear(); signal_labels.clear(); }

   
  
  static bool match( const std::set<std::string> * inp_signals , std::string * l , const std::set<std::string> & slabels );
  
};

#endif
