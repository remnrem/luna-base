
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


// this tool implements an NSRR-specific merge routine
// not currently designed for external use

#include "merge.h"
#include "merge-helpers.h"

#include <ftw.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>

int fn_process_data_dictionary( const char * fpath, const struct stat *ptr, int type );

int fn_process_study_data( const char * fpath, const struct stat *ptr, int type );

dataset_t data;

struct options_t {
  options_t() {
    skip_folders.insert( "extra" );
    missing_data_symbol = "NA";
    domain_includes.clear();
    hms_delim = ":.";
    date_delim = "/-";
  }

  std::set<std::string> skip_folders;
  std::string missing_data_symbol; // for output
  std::map<std::string,std::set<std::string> > domain_includes;

  // hh:mm:ss delimiter characters (allowed up to two)
  std::string hms_delim;
  // date delimiter chars
  std::string date_delim;

  bool read_domain( const std::string & domain , const std::string & group ) const {

    if ( domain_includes.size() == 0 ) return true;
    std::map<std::string,std::set<std::string> >::const_iterator ii = domain_includes.find( domain );

    // domain not in include list
    if ( ii == domain_includes.end() ) return false;

    // domain in include list, and no groups specified
    if ( ii->second.size() == 0 ) return true;

    // groups specified, is the domain::group pairing in the include
    return ii->second.find( group ) != ii->second.end() ;

  }
    
  
};

options_t options;

int main( int argc , char ** argv )
{

  // merge -d derived/domains -s /studies/Study2 -o merged/study2.txt                                                                    
  if ( argc < 7 )
    {
      std::cerr << "usage: merge -d path/to/domains -s path/to/study -o path/to/outputfile {domains domain-groups}\n";
      std::exit(1);
    }

  std::string domain_dir = "" , study_dir = "" , outfile = "";

  if ( std::string( argv[1] ) != "-d" ) halt( "usage: merge -d path/to/domains -s path/to/study -o path/to/outputfile {domains domain-groups}");
  if ( std::string( argv[3] ) != "-s" ) halt( "usage: merge -d path/to/domains -s path/to/study -o path/to/outputfile {domains domain-groups}");
  if ( std::string( argv[5] ) != "-o" ) halt( "usage: merge -d path/to/domains -s path/to/study -o path/to/outputfile {domains domain-groups}");
  domain_dir = argv[2];
  study_dir = argv[4];
  outfile = argv[6];
  
  // any other args are to only include these domains
  for (int i=7;i<argc;i++)
    {
      std::string t = argv[i];
      std::vector<std::string> tok = parse( t , "-" );
      if ( tok.size() > 2 ) halt( "invalid domain-group specification" );
      if ( tok.size() == 1 ) {
	std::set<std::string> empty;
	options.domain_includes[ tok[0] ] = empty;
      }
      else
	options.domain_includes[ tok[0] ].insert( tok[1] );
    }
  
  
  //
  // Step 1. Read domain-based data dictionaries, check across domains
  //
  
  if ( ftw( expand( domain_dir ).c_str() , fn_process_data_dictionary, 10 ) != 0 )
    halt( "problem traversing folder " + domain_dir );

  // ensure uniqueness of variable names across data domains
  data.check_variables_across_domains();
  

  //
  // Step 2. Read data-files
  //
  
  if ( ftw( expand( study_dir ).c_str() , fn_process_study_data, 10 ) != 0 )
    halt( "problem traversing folder " + study_dir );

  
  //
  // Step 3. Output
  //

  // nothing to do?
  if ( data.xvars.size() == 0 || data.indivs.size() == 0 )
    halt( "no data (variables and/or individuals) available for output" );
  
  std::ofstream OUT1( outfile.c_str() , std::ios::out );
  data.write(OUT1);
  OUT1.close();

  //
  // All done
  //

  std::exit(0);
}





int fn_process_data_dictionary( const char * fpath, const struct stat *ptr, int type )
{

  if ( type == FTW_F ) 
    {

      std::string filename( fpath );

      // ignore back up files ending ~
      if ( filename[ filename.size() - 1 ] != '~' )
	{      
	  domain_t domain( filename );
	  data.add( domain );
	}
    }  
  return 0;
}



int fn_process_study_data( const char * fpath, const struct stat *ptr, int type )
{

  if ( type == FTW_F ) 
    {

      std::string filename( fpath );

      // ignore back up files ending ~
      if ( filename[ filename.size() - 1 ] != '~' )
	{      	  
	  data.read( filename );
	}
    }  
  return 0;
}



//
// Read data dictionary
//

int domain_t::read( const std::string & filename )
{

  //
  // get domain-group naming from filename
  //

  // 1) remove any .txt extension and preceding '/' folder separators

  std::vector<std::string> tok = parse( filename , "/" );
  if ( tok.size() == 0 ) halt( "invalid " + filename );  
  
  // 2) name should be in 'domain-group' form, i.e. two hypen-delimited words
  std::string domain_group = remove_extension( tok[ tok.size() - 1 ] , "txt" );
  tok = parse( domain_group , "-" );
  if ( tok.size() != 2 ) halt( "expected 'domain-group' naming for " + domain_group );

  // all good, assign:
  name = tok[0];
  group = tok[1];

  //
  // shoud we read this?
  //

  bool read_this = options.read_domain( name , group );

  if ( ! read_this ) return 0 ;
    
  //
  // load variable definitions from file
  //
  
  if ( ! fileExists( filename ) ) halt( "could not open " + filename );

  std::ifstream IN1( filename.c_str() , std::ios::in );
  
  while ( ! IN1.eof() )
    {
      std::string line;
      safe_getline( IN1 , line ) ;
      std::vector<std::string> tok = parse(  line , "\t" );
      if ( tok.size() == 0 ) continue;

      // make all upper case for
      std::string varname = toupper( tok[0] );

      // special missing data code?
      if ( tok.size() == 3 && iequals( tok[1] , "missing" ) )
	{
	  missing = tok[2];
	  continue;
	}

      // otherwise, process as a variable

      if ( tok.size() != 3 )
	halt( "error in " + filename
	      + "\n  -- expecting 3 tab-delimited columns\n  -- "
	      + line + "\n" );
      
      
      // check this is not already populated      
      if ( variables.find( varname ) != variables.end() )
	halt( "duplicate of " + varname + " in " + filename );

      // variabkle names cannot contain periods
      if ( varname.find(".") != std::string::npos )
	halt( "variable names cannot contain periods: " );
	      
      // assign variable
      var_t var( *this , varname , tok[1] , tok[2] ) ;
    
      // and store
      variables[ tok[0] ] = var;
      
    }
  
  IN1.close();
  
  return variables.size();
}


//
// Step 2. Read study data for a given individual
//

void dataset_t::read( const std::string & filename )
{
  // expecting 
  
  // figure out which individual this is, as folder name
  // you can have subfolders, but not within an indiv folder
  // i.e.   study1/domain1/indiv1/file.txt   OKAY
  // i.e.   study1/indiv1/domain1/file.txt   NOT OKAY
  
  if ( ! fileExists( filename ) ) halt( "could not open " + filename );

  //
  // get ID from subfolder 
  //
  
  std::vector<std::string> tok = parse( filename , "/" );  
  if (tok.size() < 2 ) halt( "problem, expecting study/indiv/file.txt structure" );
  std::string folder_indiv_id = tok[ tok.size() - 2 ];

  //
  // Skip subfolders that start with code 
  //

  if ( options.skip_folders.find( folder_indiv_id ) != options.skip_folders.end() )
    {
      std::cerr << " -- skipping " << filename << "\n";
      return;
    }
  
  //
  // Get domain and factors from filename
  //
  
  // filename format:  {domain}-{group}-{tag}{.fac1}{.fac2}{.f3_l3}{.txt}
  // the 'tag' and factors/lvls can have '-' characters in them; i.e.
  // we just delimit based on the first two '-' characters;  everything after is 
  // taken 'as is'

  // domain-group-LUNA-COMMAND.F1.F2
  // [ domain ]  [ group ]  [ LUNA-COMMAND.F1.F2 ] 

  std::string fname = remove_extension( tok[ tok.size() - 1 ] , "txt" );
    
  std::vector<std::string> tok3 = parse( fname , "-" );

  if ( tok3.size() < 3 )
    {
      std::cerr << "found " << tok3.size() << " '-'-delimited items, expecting at least 3: " << fname << "\n";
      halt( "err1: expecting {domain}-{group}-{tag-name}{.fac1}{.fac2}{.f3_l3}{.txt}\n" );
    }

  //
  // Get data dictionary 
  //

  std::string domain_name = tok3[0];
  std::string group_name = tok3[1];

  const domain_t * domain = data.domain( domain_name , group_name );
  if ( domain == NULL ) halt( "could not find a dictionary for " + domain_name + "-" + group_name );

  // i.e. strip 'domain-group-' off start of filename
  std::string remainder = fname.substr( domain_name.size() + group_name.size() + 2 );
  
  
  //
  // shoud we read this?
  //

  bool read_this = options.read_domain( domain_name , group_name );

  if ( ! read_this ) return;

  //
  // domain-specific missing data code?
  //
  
  bool missing_code = domain->missing != "";

  //
  // split off any factors
  //
  
  std::vector<std::string> tokb = parse( remainder , "." );

  std::string tag_name = tokb[0];
  
  // any factors [ and optionally, levels ] 
  std::vector<std::string> facs;
  std::map<std::string,int> fac2col;
  std::map<std::string,std::string> setfac;// fixed factors, set by filename, e.g. SS_N2
  for (int i=1;i<tokb.size();i++)
    {
      std::vector<std::string> tokfl = parse( tokb[i] , "_" );
      if ( tokfl.size() == 1 ) { facs.push_back( toupper( tokfl[0] ) ) ; setfac[ toupper( tokfl[0] ) ] = "" ; }
      else if ( tokfl.size() == 2 ) { facs.push_back( toupper( tokfl[0] ) ); setfac[ toupper( tokfl[0] ) ] = tokfl[1] ; }
      else halt( "expecting {domain}-{group}-{tag}{.fac1}{.fac2}{.f3_l3}{.txt}\n" + filename );
      
      // check that tokfl[0] was specified /as a factor/ in the data dictionary
      if ( ! domain->has( toupper( tokfl[0] ) , FACTOR ) )
	halt( tokfl[0] + " not specified as a factor" );
    }
  
  //
  // Read actual data
  //
  
  std::set<var_t> vars;

  indiv_t indiv( folder_indiv_id );
  
  std::ifstream IN1( filename.c_str() , std::ios::in );
  bool had_header = false;
  int cols = 0;
  int id_col = 0;
  std::vector<std::string> colvar; // track column names
  std::set<std::string> colcheck; // for dupes
  while ( ! IN1.eof() )
    {
      
      if ( ! had_header )
	{
	  std::string header;
	  safe_getline( IN1 , header );

	  if ( IN1.eof() ) break;
	  std::vector<std::string> tok = parse( header , "\t" );
	  std::string id = "";
	  for (int i=0;i<tok.size();i++)
	    {
	      if ( iequals( tok[i] , "ID" ) )
		{
		  if ( id != "" ) halt( "multiple ID columns in " + filename );
		  id = tok[i];
		  id_col = i;
		}
	      else
		{
		  // enfore upper-case
		  std::string varname = toupper( tok[i] );

		  // check it exists in the dictionary
		  if ( ! domain->has( varname ) )
		    halt( varname + " not specified in data-dictionary for " + filename );
		  		  
		  // for checking factors exist
		  colcheck.insert( varname );

		  if ( domain->has( varname , FACTOR ) )
		    fac2col[ varname ] = i;
		}

	    } // next header col

	  // done reading header


	  //
	  // Check that ID col was present
	  //

	  if ( id == "" ) halt( "no ID column specified" );

	  //
	  // Check all factors that should be present are present
	  //

	  for (int f=0;f<facs.size();f++)
	    {
	      // if not a pre-specified factor (i.e. from filename)
	      if ( setfac[ facs[f] ] == "" )
		{
		  if ( colcheck.find( facs[f] ) == colcheck.end() )
		    halt( "could not find necessary factor " + facs[f] + "\n        in file: " + filename );
		}
	    }

	  //
	  // Done processing header
	  //
	  
	  cols = tok.size();
	  colvar = tok;
	  had_header = true;

	}
      else
	{

	  //
	  // Read a data row
	  //

	  std::string line;
	  safe_getline( IN1 , line );
	  if ( IN1.eof() ) break;
	  std::vector<std::string> tok = parse( line , "\t" );
	  if ( tok.size() == 0 ) continue;
	  if ( tok.size() != cols ) halt( "inconsistent number of columns versus header in " + filename );

	  //
	  // check ID matches subfolder ID
	  //

	  if ( tok[ id_col ] != folder_indiv_id )
	    halt( "folder for [" + folder_indiv_id + "] contains different ID ["
		  + tok[id_col] + "]\n        in file: " + filename );

	  //
	  // get any factor levels
	  //
	  
	  std::vector<std::string> lvls( facs.size() , "" );
	  
	  for (int f=0;f<facs.size();f++)
	    {
	      if ( setfac[ facs[f] ] != "" )
		lvls[f] = setfac[ facs[f] ] ; // from file name
	      else
		lvls[f] = tok[ fac2col[ facs[f] ] ]; // from data
	    }
	  

	  //
	  // add variables
	  //
	  
	  for (int i=0;i<tok.size();i++)
	    {

	      // ignore ID column
	      if ( i == id_col ) continue;
	      
	      // ignore FACTOR columns (they will go to VARNAMES)
	      if ( setfac.find( colvar[i] ) != setfac.end() ) continue; 

	      // get (base) variable; as header was checked, if here, var will always exist
	      const var_t * var = domain->variable( colvar[ i ] );

	      // get expanded variable

	      var_t xvar = data.xvar( *var , facs , lvls );
	      
	      // missing value?

	      if ( missing_code && tok[i] == domain->missing ) continue;

	      // get value

	      value_t value( tok[i] );

	      // check value

	      if ( ! type_check( tok[i] , var->type ) )
		halt( "invalid value [" + tok[i] + "] for " + var->name
		      + " (type " + var->print_type() + ")\n        in: " + filename );
	      
	      // insert
	      
	      indiv.add( xvar , value );

	    }
	  
	}

      
      
    }
  IN1.close();

  //
  // merge new data w/ existing
  //

  if ( ! had_header ) halt( "no header read for " + filename );
  
  data.add( indiv );  
  
  std::cerr << " ++ added to " << domain_name << "::" << group_name << " : " << filename << "\n";


}


//
// Step 3. Write data back
//

void dataset_t::write( std::ofstream & OUT1 )
{

  //
  // Data dictionary goes to stdout
  //

  // header, listing all factors present in the data

  std::cout << "COL\tVAR\tBASE\tDOMAIN\tGROUP\tTYPE\tDESC";

  std::set<std::string> factors;
  std::set<domain_t>::const_iterator dd = domains.begin();
  while ( dd != domains.end() )
    {
      std::map<std::string,var_t>::const_iterator vv = dd->variables.begin();
      while ( vv != dd->variables.end() )
	{
	  if ( vv->second.type == FACTOR ) factors.insert( vv->first );
	  ++vv;
	}
      ++dd;
    }

  std::set<std::string>::const_iterator ff = factors.begin();
  while ( ff != factors.end() )
    {
      std::cout << "\t" << *ff;
      ++ff;
    }

  std::cout << "\n";
  
  //
  // Factor descriptions (not in dataset)
  //

  ff = factors.begin();
  while ( ff != factors.end() )
    {

      if ( faclabels.find( *ff ) == faclabels.end() ) 
	halt( "internal error, could not look up factor label" );
      
      std::cout << "0\t"
		<< *ff << "\t"
		<< ".\t.\t.\t"
		<< "Factor\t"
		<< faclabels[ *ff ];
      
      std::set<std::string>::const_iterator ff2 = factors.begin();
      while ( ff2 != factors.end() )
	{
	  std::cout << "\t.";
	  ++ff2;
	}

      std::cout << "\n";
      ++ff;
    }

  
  
  //
  // First ID row
  //
  
  std::cout << "1\tID\t.\t.\t.\tID\tIndividual ID";

  ff = factors.begin();
  while ( ff != factors.end() )
    {
      std::cout << "\t.";
      ++ff;
    }
  std::cout << "\n";


  //
  // Variable/data rows
  //

  int cnt = 1;
  
  std::set<var_t>::const_iterator vv = xvars.begin();
  while ( vv != xvars.end() )
    {
      std::cout << ++cnt << "\t"
		<< vv->name << "\t"
		<< ( vv->base != vv->name ? vv->base : "." ) << "\t"
		<< vv->domain_name << "\t"
		<< vv->domain_group << "\t"
		<< vv->print_type() << "\t"
		<< vv->label ;
      
      ff = factors.begin();
      while ( ff != factors.end() )
	{
	  std::map<std::string,std::string>::const_iterator ll =  vv->fac2lvl.find( *ff );
	  if ( ll != vv->fac2lvl.end() )
	    std::cout << "\t" << ll->second;
	  else
	    std::cout << "\t.";
	  ++ff;
	}
      
      std::cout << "\n";
      
      ++vv;
    }

  
  //
  // End of data dictionary output
  //
  
  
  
  //
  // Compiled individual-level goes to O1 stream
  //

  //
  // Header
  //

  OUT1 << "ID";

  vv = xvars.begin();
  while ( vv != xvars.end() )
    {
      OUT1 << "\t" << vv->name;
      ++vv;
    }
  OUT1 << "\n";

  //
  // Individual rows
  //
  
  std::set<indiv_t>::const_iterator ii = indivs.begin();
  while ( ii != indivs.end() )
    {
      OUT1 << ii->id;      
      vv = xvars.begin();
      while ( vv != xvars.end() )
	{

	  std::map<var_t,value_t>::const_iterator kk = ii->values.find( var_t( vv->name ) );
	  if ( kk == ii->values.end() )
	    OUT1 << "\t" << options.missing_data_symbol;
	  else
	    OUT1 << "\t" << kk->second.data;
	  ++vv;
	}

      OUT1 << "\n";
      ++ii;
    }
  
}


void dataset_t::check_variables_across_domains()
{

  std::set<std::string> varnames;
  std::set<domain_t>::const_iterator ii = domains.begin();
  while ( ii != domains.end() )
    {      
      std::map<std::string,var_t>::const_iterator jj = ii->variables.begin();
      while ( jj != ii->variables.end() )
	{

	  // skip FACTORS (these can be duplicated across dictionaries)
	  if ( ii->has( jj->first , FACTOR ) )
	    {
	      // save FACTOR description, checking this is consistent across data dictionaries
	      if ( faclabels.find( jj->first ) != faclabels.end() )
		{
		  if ( faclabels.find( jj->first )->second !=  jj->second.label )
		    halt( "inconsistent label for factor " + jj->first + " across data dictionaries" );
		}
	      else
		faclabels[ jj->first ] = jj->second.label;
	      
	      ++jj;
	      continue;
	    }
	    
	  if ( varnames.find( jj->first ) != varnames.end() )
	    halt( jj->first + " is duplicated across data dictionaries" );
	    
	  varnames.insert( jj->first );

	  ++jj;
	}      
      ++ii;
    }
}



var_t::var_t( const domain_t & domain ,
	      const std::string & _name ,
	      const std::string & t , 
	      const std::string & _label ,
	      const std::string & _base )
{

  name = _name;

  label = _label;

  base = _base;

  if ( base == "" ) base = name;
  
  domain_name = domain.name;

  domain_group = domain.group;

  if      ( imatch( t , "factor" ) ) type = FACTOR;
  else if ( imatch( t , "text" ) ) type = TEXT;
  else if ( imatch( t , "int" ) ) type = INT;
  else if ( imatch( t , "numeric" ) ) type = FLOAT;
  else if ( imatch( t , "yesno" ) ) type = YESNO;
  else if ( imatch( t , "yn" ) ) type = YESNO;
  else if ( imatch( t , "date" ) ) type = DATE;
  else if ( imatch( t , "time" ) ) type = TIME;
  
   
}



bool type_check( const std::string & value , type_t type )
{
  if ( type == TEXT || type == FACTOR ) return true;

  if ( value == options.missing_data_symbol ) return true;
  
  if ( type == FLOAT )
    {
      double d;
      return str2dbl( value , &d );
    }

  if ( type == INT )
    {
      int i;
      bool okay = str2int( value , &i );
      if ( ! okay ) return false;
      // allow negative numebrs, but no 'e' or '.' notation
      // i.e. scan for those: if found, must be a float
      for (int c=0;c<value.size();c++)
	{
	  if ( value[c] == 'e' ) return false;
	  if ( value[c] == 'E' ) return false;
	  if ( value[c] == '.' ) return false;
	}
      return true;
    }

  if ( type == YESNO )
    {
      bool yes = imatch( value , "y" ) || imatch( value , "t" ) || value == "1" ;
      bool no = imatch( value , "n" ) || imatch( value , "f" ) || value == "0" ;
      return yes || no;
    }

  if ( type == DATE )
    {
      // requires M/Y (or M-Y)
      //  or D/M/Y (or D-M-Y)
      // **does not check for exact ordering** e..g valid date 
      std::vector<std::string> tok = parse( value , options.date_delim );
      if ( tok.size() < 2 || tok.size() > 3 ) return false;
    }

  if ( type == TIME )
    {
      std::vector<std::string> tok = parse( value , options.hms_delim );
     if ( tok.size() < 2 || tok.size() > 3 ) return false;
    }

  return true;
}
