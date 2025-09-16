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

#include "luna.h"
#include "main.h"

// DUMMY : a generic placeholder/scratchpad for templating new things
//  -- tests here generally hard-coded, not intended to be used/reproducible

void proc_dummy( const std::string & p , const std::string & p2 )
{

  // circular SD calcs

  if ( p == "circ" )
    {

      MiscMath::running_stats_calc_t cs(5);

      std::vector<double> sourceValues = { 1000000,22.2,33.3,44.4,55.5,66.6,77.7,88.8,0.0,100.1 };
      const int n = sourceValues.size();
      
      for (int i=0; i<n; i++) { 
	cs.update( sourceValues[i] );

	std::cout << i << "\t" << cs.mean() << "\t" << cs.sampleStdev() << "\n";
      }
      
      std::exit(0);
    }

  if ( p == "interval-tree" )
    {

      // parent(p) , interval(i) , id(s) , ch_str( ch_str )
      std::vector<instance_idx_t> v = {
	{NULL, {10,20}, "1", "A"},
	{NULL, {10,20}, "2", "B"},   // same [10,20), different id/name
	{NULL, {15,25}, "3", "C"}
      };
      
      interval_tree_t T(v.begin(), v.end());
      
      auto hits = T.query_ptrs(12, 18); // overlaps [12,18)

      for (auto* p : hits) {
	// You'll see both A and B separately
	std::cout << p->id << " " << p->ch_str << " ["
		  << p->interval.start << "," << p->interval.stop << ")\n";
      }
            
      std::exit(0);
      
    }
  
  if ( p == "circ2" )
    {
      Eigen::VectorXd X = Eigen::VectorXd::Zero( 8 );
      X << 0, 1 , 2, 3, 4, 5, 6, 7;
      const int sr = 1;
      Eigen::VectorXd Z = eigen_ops::rolling_norm( X , 5 );
      std::cout << Z << "\n";
      std::exit(0);
    }

  
  // quantify peaks in a power spectrum
  
  if ( p == "peaks" )
    {
      std::vector<double> x(100);
      for (int i=0;i<100;i++) x[i] = i*i;
      x[9] += 2000;
      x[20] += 2000;
      double m1, m2;
      std::vector<double> s1, s2, s3;
      psd_shape_metrics( x , x , 5 , &m1, &m2, &s1, &s2, &s3 );
      for (int i=0;i<s1.size();i++)
	std::cout << x[i] << "\t"
		  << s1[i] << "\t"
		  << s2[i] << "\t"
		  << s3[i] << "\n";

      std::cout << "m1\t" << m1 << "\n"
		<< "m2\t" << m2 << "\n";
     
      std::exit(0);
    }


  // K means clustering
  
  if ( p == "kmeans" )
    {
      const int nc = 4;
      const int nr = 150;
      const int nk = 3;
      Data::Matrix<double> X( nr , nc );
      for (int r=0;r<nr; r++)
	for (int c=0;c<nc; c++)
	  std::cin >> X(r,c);
      
      std::cout << "X. " << X.print() << "\n";

      kmeans_t kmeans;
      
      std::vector<int> sol;
      kmeans.kmeans( X , 3 , &sol );
      std::cout << "SOL\n";
      for (int i=0;i<150;i++)
	std::cout << sol[i] << "\n";

      std::exit(1);
    }


  // test JSON library
  
  if ( p == "json" )
    {

      // store a string in a JSON value
      nlohmann::json j_string = "this is a string";
      
      // retrieve the string value
      auto cpp_string = j_string.get<std::string>();
      // retrieve the string value (alternative when an variable already exists)
      std::string cpp_string2;
      j_string.get_to(cpp_string2);
      
      // retrieve the serialized value (explicit JSON serialization)
      std::string serialized_string = j_string.dump();
      
      // output of original string
      std::cout << cpp_string << " == " << cpp_string2 << " == " << j_string.get<std::string>() << '\n';
      // output of serialized value
      std::cout << j_string << " == " << serialized_string << std::endl;
      
      
      std::exit(1);
      
    }

  // runs test

  if ( p == "runs" )
    {
      std::vector<std::string> d = { "S", "S", "S", "F", "S", "F", "F", "F", "F" };
      std::cout << "runs p = " << Statistics::runs_test( d ) << "\n";
      std::exit(1);
    }
     
  // canonical correlation

  if ( p == "cancor" )
    {
      const int nrows = 100;
      const int nvars = 10;

      Eigen::MatrixXd X( nrows , nvars );
      Eigen::MatrixXd Y( nrows , nvars );

      int i = 0 , j = 0;
      std::ifstream INX( Helper::expand( "~/x.txt" ).c_str() , std::ios::in );
      while ( ! INX.eof() )
        {
          double d;
          INX >> d;
          if ( INX.eof() ) break;
          X(i,j) = d;
	  ++j;
	  if ( j == nvars ) { ++i; j=0; }
	}
      INX.close();
      
      i = j = 0;
      std::ifstream INY( Helper::expand( "~/y.txt" ).c_str() , std::ios::in );
      while ( ! INY.eof() )
        {
          double d;
          INY >> d;
          if ( INY.eof() ) break;
          Y(i,j) = d;
	  ++j;
	  if ( j == nvars ) { ++i; j=0; }
	}
      INY.close();
      
      Eigen::VectorXd CCA = eigen_ops::canonical_correlation( X , Y );
      
      std::cout << " CCA \n"
       		<< CCA << "\n";
      
      std::exit(1);
    }


  // quadratic discriminant analysis
  
  if ( p == "qda" )
    {
      std::vector<std::string> y;
      const int nrows = 1257;
      const int nvars = 18;
      
      Eigen::MatrixXd X( nrows , nvars );
      
      std::ifstream INY( Helper::expand("~/y.txt" ).c_str() , std::ios::in );
      int k = 0;
      while ( ! INY.eof() )
	{
	  std::string s;
	  INY >> s;
	  if ( INY.eof() ) break;
	  y.push_back(s);
	}

      INY.close();
      int i = 0 , j = 0;
      std::ifstream INX( Helper::expand( "~/x.txt" ).c_str() , std::ios::in );
      while ( ! INX.eof() )
        {
          double d;
          INX >> d;
          if ( INX.eof() ) break;
          X(i,j) = d;
	  ++j;
	  if ( j == nvars ) { ++i; j=0; }
	}
      INX.close();

      qda_t qda( y , X );
      
      qda_model_t fit = qda.fit();

      qda_posteriors_t pp = qda.predict( fit , X );

      for (int i=0;i<pp.pp.rows() ;i++)
	{
	  for (int j=0;j<pp.pp.cols();j++)
	    std::cout << " " << pp.pp(i,j);
	  std::cout << "\t" << pp.cl[i] << "\n";

	}      
      
      std::exit(1);
      
    }

  // linear discriminant analysis
  
  if ( p == "lda" )
    {
      std::vector<std::string> y;

      Eigen::MatrixXd X(500,10);
      
      std::ifstream INY( Helper::expand("~/y.txt" ).c_str() , std::ios::in );
      int k = 0;
      while ( ! INY.eof() )
	{
	  std::string s;
	  INY >> s;
	  if ( INY.eof() ) break;
	  y.push_back(s);
	}

      INY.close();
      int i = 0 , j = 0;
      std::ifstream INX( Helper::expand( "~/x.txt" ).c_str() , std::ios::in );
      while ( ! INX.eof() )
        {
          double d;
          INX >> d;
          if ( INX.eof() ) break;
          X(i,j) = d;
	  ++j;
	  if ( j == 10 ) { ++i; j=0; }
	}
      INX.close();
      
      lda_t lda( y , X );

      lda_model_t fit = lda.fit();

      lda_posteriors_t pp = lda.predict( fit , X );

      for (int i=0;i<pp.pp.rows() ;i++)
	{
	  for (int j=0;j<pp.pp.cols();j++)
	    std::cout << " " << pp.pp(i,j);
	  std::cout << "\t" << pp.cl[i] << "\n";

	}      
      
      std::exit(1);
      
    }


  // test cache mechanism

  if ( p == "cache" )
    {
      ctest();
      std::exit(0);
    }

  // microstate kmer analysis
  
  if ( p == "kmer" )
    {
      std::vector<int> x;
      while ( ! std::cin.eof() )
	{
	  int i;
	  std::cin >> i;
	  if ( std::cin.eof() ) break;
	  x.push_back(i);
	}

      ms_kmer_t kmers( x , 2 , 6 , 1000 , 0 );
      
      std::map<std::string,double>::const_iterator pp = kmers.basic.pval.begin();
      while ( pp != kmers.basic.pval.end() )
	{
	  std::cout << pp->first << "\t"
		    << pp->first.size() << "\t"
		    << kmers.basic.obs[ pp->first ] << "\t"
		    << pp->second << "\n";
	  ++pp;
	}
      
      std::exit(1);
    }


  // test command-defintions syntax/logic
  
  if ( p == "cmddefs" ) 
    {
      
      globals::cmddefs().add_domain( "misc" , "misc" ,  "Misc" );
      
      globals::cmddefs().add_cmd( "misc" , "comm1" , "this is a dummy command" );
      globals::cmddefs().add_table( "comm1" , "XX" , "A simple table" , false );
      globals::cmddefs().add_var( "comm1" , "XX" , "X", "Var X" );
      globals::cmddefs().add_var( "comm1" , "XX" , "Y" , "Var Y" );
      
      globals::cmddefs().add_cmd( "misc" , "comm2" , "this is a dummy command" );
      globals::cmddefs().add_table( "comm2" , "CH,B" , "A nice table" , true );
      globals::cmddefs().add_var( "comm2" , "CH,B" , "V1" , "Variable 1" );
      globals::cmddefs().add_var( "comm2" , "CH,B" , "V2" , "Variable 2" );
      
      //   std::cout << globals::cmddefs().help( "comm1" , true )  << "\n\n\n";
      std::cout << globals::cmddefs().help( "comm2" , true )  << "\n";
      
      // add a dummy tag
      globals::cmddefs().add_tag( "Z" );
      
      zfiles_t files( "folder1" , "indiv1" ); 
      
      zfile_t * z1 = files.file( "comm1" , NULL , "XX" ) ; 
      
      param_t param2;
      
      zfile_t * z2 = files.file( "comm2" , &param2 , "CH,B,Z" ) ; 
      
      z1->write_header();
      
      z2->write_header();
      
      //z1->display() ;
      //z2->display() ;
      
      z1->set_stratum( "XX" , "L1" );
      z1->set_value( "X" , 22 );
      z1->set_value( "Y" , 23 );
      z1->write_buffer();
      z1->set_stratum( "XX" , "L2" );
      z1->set_value( "X" , 24 );
      z1->set_value( "Y" , 25 );
      z1->write_buffer();
      
      z2->set_stratum( "CH" , "C3" );
      z2->set_stratum( "B" , "ALPHA" );
      z2->set_stratum( "Z" , "R1" );
      z2->set_value( "V1" , 22 );
      z2->set_value( "V2" , 23 );
      z2->write_buffer();
      
      files.close();
      
      std::exit(1);
      
      globals::cmddefs().add_cmd( "misc"   , "NEWONE" , "A test command" );
      globals::cmddefs().add_table( "NEWONE" , "" , "Table 0, baseline" );
      globals::cmddefs().add_table( "NEWONE" , "CH" , "Table 1, by channel" );
      globals::cmddefs().add_table( "NEWONE" , "CH,X" , "Table 2, by channel and X" );
      globals::cmddefs().add_table( "NEWONE" , "CH,X,Y" , "Table 2a, by channel and X/Y" , true );
      globals::cmddefs().add_table( "NEWONE" , "CH,X,Z" , "Table 2b, by channel and X/Z"  );

      globals::cmddefs().add_var( "NEWONE" , "" , "V1" , "some var1" );
      globals::cmddefs().add_var( "NEWONE" , "" , "V2" , "some var2" );

      globals::cmddefs().add_var( "NEWONE" , "CH" , "V1" , "some var1" );
      globals::cmddefs().add_var( "NEWONE" , "CH" , "V2" , "some var2" );
      globals::cmddefs().add_var( "NEWONE" , "CH" , "V3" , "some var3" );

      globals::cmddefs().add_var( "NEWONE" , "CH,X" , "V2a" , "some var2" );
      globals::cmddefs().add_var( "NEWONE" , "CH,X" , "V3a" , "some var3" );

      globals::cmddefs().add_var( "NEWONE" , "CH,X,Y" , "V2a" , "some var2" );
      globals::cmddefs().add_var( "NEWONE" , "CH,X,Y" , "V3a" , "some var3" );

      globals::cmddefs().add_var( "NEWONE" , "CH,X,Z" , "V2a" , "some var2" );
      globals::cmddefs().add_var( "NEWONE" , "CH,X,Z" , "V3a" , "some var3" );

      //  std::cout << globals::cmddefs().help_domains() << "\n";
      //std::cout << globals::cmddefs().help_commands() << "\n";
      
      globals::cmddefs().set_compressed( "NEWONE" , tfac_t( "CH,X,Z" ) );
      globals::cmddefs().set_compressed( "NEWONE" , tfac_t( "CH,X,Y" ) , false );
      
      std::cout << globals::cmddefs().help( "NEWONE" , true ) << "\n";
            
      std::exit(0);

      param_t param;
      param.add( "epoch" );
      param.add( "ep" );
      std::set<std::string> unk;
      std::cout << "st = " << globals::cmddefs().check( "ANNOTS" , param.keys() , &unk ) << "\n";
      
      std::set<std::string>::const_iterator uu = unk.begin();
      while ( uu != unk.end() ) { std::cout << " bad param " << *uu << "\n"; ++uu; } 

    }


  // test LightGBM

  if ( p == "lgbm" )
    {
#ifdef HAS_LGBM
      
      lgbm_t lgbm( "train.conf" ) ;
      
      lgbm.load_training_data( "binary.train" );

      const int n1 = lgbm_t::rows( lgbm.training );
      const int n2 = lgbm_t::cols( lgbm.training );

      lgbm_label_t labels( "luna.wgt" );
      std::cout << " from luna.wgt " << labels.n << "\n";

      lgbm.add_label_weights( lgbm.training , &lgbm.training_weights, labels );

      lgbm.apply_weights( lgbm.training , &lgbm.training_weights );
      
      //lgbm.load_weights( lgbm.training , "RENAMED_binary.train.weight" );

      std::vector<int> l = lgbm_t::labels( lgbm.training );
      std::vector<double> w = lgbm_t::weights( lgbm.training );
      
      std::cout << " l = " << l.size() << " ... \n";
      for (int i=0; i<30; i++) std::cout << l[i] << "\t" << w[i] << "\n";
                 
      lgbm.load_validation_data( "binary.test" );
      
      lgbm.create_booster();
      
      lgbm.save_model( "my-model.1" );      
      
#endif
      std::exit(0);
    }



  if ( p == "lgbm2" )
    {
#ifdef HAS_LGBM

      Eigen::MatrixXd X = eigen_ops::load_mat( "binary.test" );

      // remove first col
      X = X.rightCols( X.cols() - 1 );

      lgbm_t lgbm;

      lgbm.load_model( "my-model.1" );
      //lgbm.load_model( "LightGBM_model.txt" );
      //lgbm.load_model( "R.bst" );

      lgbm.predict( X );
      
#endif
      std::exit(0);
    }


  // test date/time functions

  if ( p == "datetime" )
    {
      
      if ( 0 )
	{
	  for (int c=0; c<10000; c++)
	    {
	      std::string ds = date_t::datestring( c );
	      date_t dt( ds );
	      int c2 = date_t::count( dt );
	      std::cout << ( c != c2 ? "*****" : "" ) 
			<< c << "\t"
			<< c2 << "\t"
			<< ds << "\n";
	    }
	}

      if ( 1 )
	{
	  std::string inp1, inp2;
	  std::cin >> inp1 >> inp2;
	  clocktime_t t1( inp1 );
	  clocktime_t t2( inp2 );

	  std::cout << "t1: "
		    << t1.valid << "\t"
		    << t1.as_string( ) << "\t"
		    << t1.as_datetime_string( ) << "\n";
	  //<< t1.d << " - " << t1.h << ":" << t1.m << ":" << t1.s << "\n";
	  
	  std::cout << "t2: "
		    << t2.valid << "\t"
		    << t2.as_string( ) << "\t"
		    << t2.as_datetime_string( ) << "\n";
	  //<< t2.d << " - " << t2.h << ":" << t2.m << ":" << t2.s << "\n";

	  const int earlier = clocktime_t::earlier( t1 , t2 ) ; 
	  
	  const double difh = earlier == 0 ? 0
	    : ( earlier == 1 ? clocktime_t::difference_hours( t1 , t2 )
	      : clocktime_t::difference_hours( t2 , t1 )   );
	  const double difs = earlier == 0 ? 0
	    : ( earlier == 1 ? clocktime_t::difference_seconds( t1 , t2 )
		: clocktime_t::difference_seconds( t2 , t1 ) );

	  clocktime_t midpoint;
	  midpoint.midpoint( t1 , t2 );
	  
 	  std::cout << " earlier = " << earlier << "\n";
	  std::cout << " t1 - t2 (hours) = " << difh << "\n";
	  std::cout << " t1 - t2 (secs) = " << difs << "\t" << difs / 3600.0 << "\n";
	  std::cout << " midpoint = " << midpoint.as_datetime_string() << "\n";
	  std::cout << "\n";
	  clocktime_t nt = t1;
	  for (int i=0;i<48;i++)
	    {
	      //nt.advance_hrs( 1.222 );
	      nt.advance( clocktime_t( "+1:30" ) );
	      std::cout << "  --> " << nt.as_string() << "\t" << nt.as_datetime_string() << "\n";
	    }
	}
      
      std::exit(0);
    }
  

  // test 'eval-expression' TRANS command

  if ( p == "trans" )
    {      
      std::string line;
      Helper::safe_getline( std::cin , line );
      std::vector<std::string> hdr = Helper::parse( line );
      const int k = hdr.size();
      
      std::cerr << " expr [" << p2 << "]\n";
      
      std::map<std::string,std::vector<double> > inputs;

      int rows = 0;
      
      while ( 1 )
	{
	  std::string line;
	  Helper::safe_getline( std::cin , line );
	  if ( line == "" ) break;
	  if ( std::cin.eof() || std::cin.bad() ) break;
	  std::vector<std::string> tok = Helper::parse( line );
	  if ( tok.size() != k ) Helper::halt( "wrong numbr of columns" );
	  for (int i=0; i<k; i++)
	    {
	      double d;
	      if ( ! Helper::str2dbl( tok[i] , &d ) ) Helper::halt( "bad numeric value" );
	      inputs[hdr[i]].push_back( d );
	    }
	  ++rows;
	}
      
      std::cerr << "read " << rows << "\n";
      
      // output
      instance_t out;
      
      // expression
      Eval expr( p2 );
      
      // bind input/output data to token evaluator
      expr.bind( inputs , &out );

      // evaluate
      bool is_valid = expr.evaluate( );

      // returned a valid bool? (single value)
      bool retval;      
      bool is_valid_retval = true;
      if ( ! expr.value( retval ) ) is_valid_retval = false;

      std::cerr << "parsed as a valid expression : " << ( is_valid ? "yes" : "no" ) << "\n";
      if ( is_valid_retval ) 
	std::cerr << "boolean return value         : " << ( retval ? "true" : "false" ) << "\n";
      std::cerr << "assigned meta-data           : " << out.print() << "\n";  

      //
      // actual output
      //
      
      std::vector<double> rr = expr.value().as_float_vector();
      for (int i=0; i<rr.size();i++)
	std::cout << rr[i] << "\n";

      std::exit(1);
    }


  // Test lunapi_t

  if ( p == "lunapi" )
    {
      std::cout << "lunapi test\n";

      std::cout << " firing up ... \n";
      lunapi_t * lp = lunapi_t::inaugurate();
      std::cout << " done\n";

      lunapi_inst_ptr p = lp->inst( "id1" );

      std::cout << " attaching EDFs (local)\n";

      std::cout << "re-init...\n";
      lp->var( "alias" , "XXX|Light" );
      
      p->attach_edf( "~/tutorial/edfs/learn-nsrr01.edf" );
      p->attach_annot( "~/tutorial/edfs/learn-nsrr01-profusion.xml" );

      std::cout << p->get_id() << " .. " << p->get_edf_file() << "\n";
      
      std::vector<std::string> d = p->desc();
      for (int i=0; i<d.size(); i++) std::cout << " d["<<i<<"] = " << d[i] << "\n";


      lp->re_init();

      lunapi_inst_ptr p2 = lp->inst( "id2" );
      p2->attach_edf( "~/tutorial/edfs/learn-nsrr02.edf" );
      p2->attach_annot( "~/tutorial/edfs/learn-nsrr02-profusion.xml" );

      std::cout << p2->get_id() << " .. " << p2->get_edf_file() << "\n";
      d = p2->desc();
      for (int i=0; i<d.size(); i++) std::cout << " d["<<i<<"] = " << d[i] << "\n";

      return;
      
      // introduce some gaps
      //p->eval( "EPOCH & MASK epoch=5-50,500-600,800-850 & MASK flip & RE" );
      //p->eval( "MASK ifnot=N2 & RE" );
      
      segsrv_t segsrv( p );

      std::vector<std::string> chs = { "SaO2" , "PR" , "EEG", "EEG_sec" , "ECG" , "EMG" ,  "EOG_L" , "EOG_R" , "EEG",
				       "AIRFLOW" , "THOR_RES" , "ABDO_RES" , "POSITION" , "LIGHT" , "OX_STAT"   };

      std::vector<std::string> anns = {  "Arousal", "Hypopnea", "N1", "N2", "N3", "Obstructive_Apnea", "R", "SpO2_artifact","SpO2_desaturation","W"  };
      
      std::vector<std::string> chs_b = { "EEG" };
      std::vector<std::string> chs_h = { "EEG" , "AIRFLOW" };

      std::cout << " about to pull\n";
      
      float a = 0 ;
      float b = 61440 ;

      segsrv.populate( chs , anns ) ;

      a = 0 ;
      b = 40918;
      b = 1;

      segsrv.set_window( a, b );
      segsrv.compile_evts( anns );

      std::cout << " done populating...\n";

      for (int s=2;s<3;s++)
        {
          std::cout << " --> " << chs[s] << "\n";
          Eigen::VectorXf XX = segsrv.get_scaled_signal( chs[s] , s );	  
	  //Eigen::VectorXf TT = segsrv.get_timetrack( chs[s]);
          std::cout << XX << "\n";
        }

      std::exit(0);


      // bands, hjorth
      //      segsrv.calc_bands( chs_b );
      //      segsrv.calc_hjorths( chs_h );

      segsrv.calc_bands( chs_b );
      segsrv.calc_hjorths( chs_h );

      //      segsrv.populate( chs , anns ) ;

      //std::exit(1);
      
      //std::cout << "EEG B " << segsrv.get_bands( "EEG" ) << "\n\n";
      segsrv.set_window( 0 , 30 );
      segsrv.compile_evts( anns );
      Eigen::VectorXf XX1 = segsrv.get_signal( "EEG" );

      segsrv.set_window( 30 , 60  );
      segsrv.compile_evts( anns );
      XX1 = segsrv.get_signal( "EEG" );
      
      segsrv.set_window( 0 , 30 );
      segsrv.compile_evts( anns );
      XX1 = segsrv.get_signal( "EEG" );
      
      //      std::cout << "EEG B " << segsrv.get_bands( "EEG" ) << "\n\n";
      // segsrv.set_window( 3000, 5000 );
      // segsrv.compile_evts( anns );
      //     std::exit(0);

      a = 0;
      b = 61440;

      segsrv.set_window( a, b );
      segsrv.compile_evts( anns );

      for (int s=0;s<chs.size();s++)
	{
	  std::cout << " --> " << chs[s] << "\n";	  
	  Eigen::VectorXf XX = segsrv.get_scaled_signal( chs[s] , s );
	  //Eigen::VectorXf TT = segsrv.get_timetrack( chs[s]);
	  std::cout << XX << "\n";
	}

      std::exit(0);
      
      a += 30.0;
      b += 30.0;

      segsrv.set_window( a, b );
      segsrv.compile_evts( anns );
      for (int s=0;s<chs.size();s++)
	{
	  Eigen::VectorXf XX = segsrv.get_scaled_signal( chs[s] , s);
	  Eigen::VectorXf TT = segsrv.get_timetrack( chs[s]);
	}

      a -= 30.001;
      b -= 30.001;

      segsrv.set_window( a, b );
      segsrv.compile_evts( anns );
      for (int s=0;s<chs.size();s++)
	{
	  Eigen::VectorXf XX = segsrv.get_scaled_signal( chs[s], s);
	  Eigen::VectorXf TT = segsrv.get_timetrack( chs[s]);
	}
      std::exit(0);


      int a1 = 0, b1 = 30;
      bool okay;

      std::cout << "\n-------\n test " << a1 << " " << b1 << "\n";
      std::cout << segsrv.set_window( a1, b1 ) << "\n" << segsrv.is_window_valid() << "\n";

      a1=0; b1 = 3000;
      std::cout << "\n-------\n test " << a1 << " " << b1 << "\n";
      std::cout << segsrv.set_window( a1, b1 ) << "\n" << segsrv.is_window_valid() << "\n";

      a1=3000; b1 = 3300;
      std::cout << "\n-------\n test " << a1 << " " << b1 << "\n";
      std::cout << segsrv.set_window( a1, b1 ) << "\n" << segsrv.is_window_valid() << "\n";

      a1=3200; b1 = 3200;
      std::cout << "\n-------\n test " << a1 << " " << b1 << "\n";
      std::cout << segsrv.set_window( a1, b1 ) << "\n" << segsrv.is_window_valid() << "\n";

      std::exit(1);
      
      std::set<std::pair<double,double> > gaps = segsrv.get_gaps();
      std::cout << " found gaps in seg sz = " << gaps.size() << "\n";
      std::set<std::pair<double,double> >::const_iterator ww = gaps.begin();
      while ( ww != gaps.end() )
	{
	  std::cout << " gapped in seg " << ww->first << " " << ww->second << "\n";
	  ++ww;
	}
      
      Eigen::VectorXf X1 = segsrv.get_signal( "EEG" );
      std::cout << "X1 = \n" << X1 << "\n";
      
      std::exit(1);
      std::cout << "EEG B " << segsrv.get_bands( "EEG" ) << "\n\n";
      
      
      // void set_scaling( const int nchs , const int nanns ,
      // 			const double yscale , const double ygroup ,
      // 			const double yheader , const double yfooter ,
      // 			const double scaling_fixed_annot ,
      //                        const bool clip );

      if ( 0 )
	{
	  segsrv.set_scaling( 2 , 1 ,
			      1 , 0 ,
			      0.1 , 0.15 , 0.1 , true );
	  
	  // check yscaling 
	  for (int i = 0 ; i < 3; i++)
	    {
	      double lwr, upr;
	      bool okay = segsrv.get_yscale_signal( i , &lwr, &upr );
	      std::cout << "yparam " << i << " " << okay << " " ;
	      if ( okay ) std::cout << lwr <<  " -- " << upr ;
	      std::cout << "\n";
	    }
	  std::exit(1);
	}
      
      std::vector<std::pair<double,double> > r = segsrv.get_time_scale();
      std::cout << " viz/clok = " << r.size() << "\n";
      for (int i=0; i<r.size(); i++)
	std::cout << r[i].first << "\t" << r[i].second << "\n";

      //std::exit(1);

      segsrv.fix_physical_scale( "EEG" , -50 , 50 );
	      
      for (int t=0; t<10; t+= 10 )
	{
	  std::cout <<"\n\n";
	  bool okay = segsrv.set_window( t, t+30 );
	  	  
	  std::set<std::pair<double,double> > g = segsrv.get_gaps();

	  std::cout << " n gaps = " << g.size() << "\n";
	  std::set<std::pair<double,double> >::const_iterator gg = g.begin();
	  while ( gg != g.end() )
	    {
	      std::cout << " gap = " << gg->first << " .. " << gg->second << "\n";
	      ++gg;
	    }
	    
	  Eigen::VectorXf X1 = segsrv.get_signal( "EEG" );
	  Eigen::VectorXf X2 = segsrv.get_signal( "SaO2" );
	  Eigen::VectorXf X3 = segsrv.get_signal( "AIRFLOW" );
	  
	  Eigen::VectorXf Z1 = segsrv.get_scaled_signal( "EEG" , 0 );
	  Eigen::VectorXf Z2 = segsrv.get_scaled_signal( "SaO2" , 1 );
	  Eigen::VectorXf Z3 = segsrv.get_scaled_signal( "AIRFLOW" , 2 );

	  Eigen::VectorXf T1 = segsrv.get_timetrack( "EEG" );
	  Eigen::VectorXf T2 = segsrv.get_timetrack( "SaO2" );
	  Eigen::VectorXf T3 = segsrv.get_timetrack( "AIRFLOW" );
	  
	  std::cout << " cols = " << X1.size() << " " <<  X2.size() << " " << X3.size() << " ";
	  std::cout << " times  = " << T1.size() << " " <<  T2.size() << " " << T3.size() << "\n";
	  
	  for (int i=0; i<X1.size(); i++)
	    std::cout << i << "\t" << T1[i] << "\t" << X1[i] << "\t" << Z1[i] << "\n";
	  std::cout << "\n----------------------------------------\n\n";

	  for (int i=0; i<X2.size(); i++)
	    std::cout << i << "\t" << T2[i] << "\t" << X2[i] << "\t" << Z2[i] << "\n";
	  std::cout << "\n----------------------------------------\n\n";
	  
	  for (int i=0; i<X3.size(); i++)
	    std::cout << i << "\t" << T3[i] << "\t" << X3[i] << "\t" << Z3[i] << "\n";
	  std::cout << "\n----------------------------------------\n\n";

	}
      
      std::exit(0);
    }


  
  //
  // Multiple test functions that take input from stdin
  //

  std::vector<double> x;
  

  if ( p == "fir" || p == "decimate" || p == "fft"
       || p == "dfa" || p == "fft-test" || p == "mtm"
       || p == "tv" || p == "psi" || p == "cwt" 
       || p == "dynam" || p == "ica" || p == "robust"
       || p == "fip" || p == "sl" || p == "acf"
       || p == "otsu" || p == "desats" || p == "zpks"
       || p == "gc" || p == "detrend" || p == "emd"
       || p == "tri" || p == "ngaus" || p == "ssa"
       || p == "xcorr" ) 
    {


      // populate single vector 'x'
      int cnt= 0;
      while ( ! std::cin.eof() )
	{
	  double xx;
	  std::cin >> xx;
	  if ( std::cin.bad() ) { std::cerr << "bad input\n"; std::exit(1);  } 
	  if ( std::cin.eof() ) break;	  
	  x.push_back( xx );	  
	  if ( ++cnt % 100000  == 0 ) std::cerr << cnt << "\n";
	}
      std::cerr << x.size() << " values read\n";
      
    }


  // Gaussian bandpass filter
  
  if ( p == "ngaus" )
    {
      double f = x[0];
      double fwhm = x[1];
      int sr = x[2];

      std::vector<double> y( x.size() - 3 );
      for (int i=0; i<y.size(); i++) y[i] = x[i+3];
      
      std::vector<double> z = narrow_gaussian_t::filter( y , sr, f , fwhm ) ;

      Eigen::VectorXd yy = Eigen::VectorXd::Zero( y.size() );
      for (int i=0; i<y.size(); i++) yy[i] = y[i];
      Eigen::VectorXd zz = narrow_gaussian_t::filter( yy , sr, f , fwhm ) ;
      
      for (int i=0; i<y.size(); i++) std::cout << y[i] << "\t" << z[i] << "\t" << zz[i] << "\n";
      std::exit(0);
      
    }


  // find desats
  
  if ( p == "desats" )
    {
      hb_find_desats_t r = hb_t::find_desats( eigen_ops::copy_array( x ) , 32 , 1.5 );

      std::cout << r.MagDown << "\n\n";

      std::cout << r.dsatStEnd << "\n\n";

      std::exit(1);
    }

  // EMD
  
  if ( p == "emd" )
    {

      emd_t emd;
      
      const int nk = emd.proc( &x );
      const int nr = emd.residual.size();
      
      for (int i=0; i<nr; i++)
	{
	  std::cout << x[i] ;
	  for (int j=0; j<nk; j++)
	    std::cout << "\t" << emd.imf[j][i] ;
	  std::cout << "\t" << emd.residual[i] << "\n";
	}
      
      std::exit(1);
    }


  // CWT

  if ( p == "cwt" )
    {
      std::cout << "input lengh = " << x.size() << "\n";
  
      int Fs = 100;
      if ( p2 != "" )
        if ( ! Helper::str2int( p2 , &Fs ) )
	  Helper::halt( "expecting sample rate as second parameter" );

      std::vector<double> fc, fwhm;
      fc.push_back( 1 ); fwhm.push_back( CWT::pick_fwhm( 1 ) );
      double timelength = 20 ;
      const bool wrapped = true;
      std::vector<double> mag, phase;
      dsptools::alt_run_cwt( x , Fs , fc[0] , fwhm[0] , timelength , wrapped , &mag , &phase );

      for ( int i=0; i< mag.size(); i++)
	std::cout << mag[i] << "\n";
      std::exit(1);

      CWT cwt;
      cwt.set_sampling_rate( 400 );
      cwt.add_wavelets( 0.5 , 5 , 30 , 0.25 , 35 , 20 ) ;

      std::vector<dcomp> w1 = cwt.alt_wavelet( 0 );
      for (int i=0;i<w1.size();i++)
	std::cout << i << "\t" << w1[i] << "\n";
      
      std::exit(1);
    }


  // Detrend
  
  if ( p == "detrend" )
    {
      const int n = x.size();
      
      std::vector<double> x2 = x;
      double beta , intercept;
      MiscMath::detrend( &x2 , &intercept , &beta );
      
      std::cout << "m, b = " << intercept << " " << beta << "\n";
      for (int i=0; i<n; i++)
	std::cout << x[i] << "\t" << x2[i] << "\n";
      
      Eigen::VectorXd T( n );
      for (int i=0; i<n; i++) T[i] = x[i];
      
      std::cout << "orig T\n" << T << "\n";

      std::cout << "detrended T\n";
      
      eigen_ops::detrend( T ) ;

      std::cout << T << "\n";
      
    }


  // Granger causality
  
  if ( p == "gc" )
    {
      
      int order = 3;

      if ( p2 != "" )
        if ( ! Helper::str2int( p2 , &order ) ) Helper::halt( "expecting integer model order as second parameter" );

      Eigen::MatrixXd X( x.size() / 2 , 2 );
      int cnt = 0;
      for (int r=0; r< x.size() / 2 ; r++)
	{
	  X(r,0) = x[cnt++] ;
	  X(r,1) = x[cnt++] ; 
	}

      std::cerr << "read " << x.size() / 2 << " observations (pairs)\n";

      // fix
      const int Nr = 99;
      const int Nl = 51;
      
      order = 15; 

      // armorf_t ar1( X.col(0) , Nr , Nl , order );
      
      // std::cout << "ar1.coeff\n" << ar1.coeff << "\n"
      // 		<< "ar1.E\n" << ar1.E << "\n";

      // armorf_t ar2( X.col(1) , Nr , Nl , order );
      
      // std::cout << "ar2.coeff\n" << ar2.coeff << "\n"
      // 		<< "ar2.E\n" << ar2.E << "\n";

      // armorf_t ar12( X , Nr , Nl , order );

      // std::cout << "ar12.coeff\n" << ar12.coeff << "\n"
      // 		<< "ar12.E\n" << ar12.E << "\n";

      signal_list_t signals;
      signals.add( 0, "S1" );
      signals.add( 1, "S2" );
      const int sr = 256;
      // std::vector<double> frqs(3);
      // frqs[0] = 1; frqs[30] ; frqs[2] = 10;

      std::vector<double> frqs = MiscMath::logspace( 10 , 40 , 15 );

      gc_t gc( X , signals , sr , 200 , 60 , &frqs );

      gc.report( signals );
      std::exit(0);

    }


  if ( p == "zpks" ) 
    {
      std::vector<interval_t> ints;
      //std::vector<int> s = MiscMath::smoothedZ( x , 30*256 , 3 , 0 , 128 , 20 , 2 , 256 , true , NULL, &ints) ;
      std::vector<int> s = MiscMath::smoothedZ( x , 400 , 3 , 0 , 96 , 0 , 0 , 0 , true , &ints , NULL, true ) ;
      
      for (int i=0; i<ints.size(); i++)
	std::cout << i << "\t" << ints[i].start << " -- " << ints[i].stop << "  " << ints[i].stop - ints[i].start << " " << ( ints[i].stop - ints[i].start ) / 256.0 << "\n";


      //std::vector<int> s = MiscMath::smoothedZ( x , 30 , 5 ) ;
      // for (int i=0; i<s.size(); i++)
      // 	std::cout << x[i] << "\t" << s[i] << "\n";
       std::exit(1);
    }
 
  if ( p == "psi" )
    {
      const int n = x.size() / 2 ;
      Data::Matrix<double> data( n , 2 );
      int r = 0;
      for (int i=0;i<n;i++)
	{
	  data(i,0) = x[ r++ ];
	  data(i,1) = x[ r++ ];
	}

      psi_t psi( &data , 100 , 200 , 200 );

      psi.calc();

      signal_list_t signals;
      signals.add( 0, "S1" );
      signals.add( 1, "S2" );
      
      psi.report( signals );
	    
      std::exit(0);
    }
  
  if ( p == "robust" )
    {
      
      const int n = x.size();
      Eigen::MatrixXd m( n , 1 );
      for (int i=0;i<n;i++) m(i,0) = x[i];
  
      eigen_ops::robust_scale( m , true , true , 0.05 );
      std::cout << "\n" << m  << "\n";
      std::exit(0);
    }

  if ( p == "otsu" )
    {
      std::map<double,double> tvals, fvals;      
      double f;
      double th = MiscMath::threshold2( x , &f, 0 , &fvals , &tvals );
      
      std::cout << "best th = " << th << "\n";

      std::map<double,double>::const_iterator ii = tvals.begin();
      while ( ii != tvals.end() )
       	{
       	  std::cout << "th = " << ii->first << "\t varB = " << ii->second << "\t F = " << fvals[ ii->first ] << "\n";
       	  ++ii;
       	}
      
      std::exit(0);
    }


  if ( p == "acf" )
    {
      acf_t acf( x );
      std::vector<double> rr = acf.acf();
      for (int i=0;i<rr.size();i++)
	std::cout << "lag = " << i << "\t" << rr[i] << "\n";
      std::exit(0);
    }


  if ( p == "anova" )
    {
      std::vector<std::string> group;
      Data::Vector<double> x;
      while ( ! std::cin.eof() ) 
	{
	  std::string g;
	  double t;
	  std::cin >> g >> t;
	  if ( std::cin.eof() ) break;
	  group.push_back( g );
	  x.push_back( t );
	}

      std::cout << Statistics::anova( group , x );
      std::exit(0);
      
    }


  // Freq-interval plots
  
  if ( p == "fip" )
    {

      int sr = 256;

      if ( p2 != "" )
        if ( ! Helper::str2int( p2 , &sr ) )
	  Helper::halt( "expecting integer sample rate as second parameter" );

      const uint64_t fs = globals::tp_1sec / sr;
      std::vector<uint64_t> tp( x.size() );
      for (int i=0;i<tp.size();i++) tp[i] = i * fs;
      double th = 0;
      bool norm = false;
      bool logit = false;
      double t_lwr = 0.1; double t_upr = 5;  double t_inc = 0.1;
      double f_lwr = 1;   double f_upr = 20; double f_inc = 0.5;
      bool logspace = false;
      bool cycles = false;

      int num_cyc = 7;

      fiplot_t fp( x , &tp , sr , 
		   th , norm , logit , 
		   t_lwr, t_upr, t_inc , cycles , 
		   f_lwr, f_upr, f_inc , num_cyc , logspace );

      std::exit(0);
    }

  // Generalized eigendecomposition

  if ( p == "ged" )
    {
      Eigen::MatrixXd X = Eigen::MatrixXd::Random(5,5);
      Eigen::MatrixXd A = X + X.transpose();
      std::cout << "Here is a random symmetric matrix, A:" << "\n" << A << "\n";
      X = Eigen::MatrixXd::Random(5,5);
      Eigen::MatrixXd B = X * X.transpose();
      std::cout << "and a random positive-definite matrix, B:" << "\n" << B << "\n\n";

      ged_t ged;
      ged.covar( A, B );
      ged.calc();
      
      std::exit(0);
    }
  

  // FIR test
  
  if ( p == "fir" ) 
    {
     
      double ripple = 0.1;
      double tw = 3;
      double f1 = 2;
      double f2 = 15;
      double fs = 1000;

      std::vector<double> fc = dsptools::design_bandpass_fir( ripple , tw , fs , f1, f2 );

      std::cerr << "bandpass FIR order " << fc.size() << "\n";
      fir_impl_t fir_impl ( fc );
      x = fir_impl.filter( &x );
      for (int i=0;i<x.size();i++) std::cout << x[i] << "\n";
      
      std::exit(1);
    }

  // FFT tests: 1) real_FFT vs mt_get_spec()

  if ( p == "fft-test" )
    {

      // test 1 : equivalence w/ mt_get_spec() and real_FFT
      
      //mtm::mt_get_spec ( b, npoints, klen, amp);  
      int fs = 256;
      double dt = 1/(double)fs;
      double * series = &(x)[0];
      int inum = x.size();
      int npoints = inum;
      int klen = MiscMath::nextpow2( inum );
      int num_freqs = 1 + klen/2;
      std::vector<double> amp( klen , 0 );
      
      //void  mtm::mt_get_spec(double *series, int inum, int klength, double *amp)
      // series = input time series
      // inum   = length of time series
      // klength = number of elements in power spectrum (a power of 2)
      // amp = returned power spectrum

      int             i, j, isign = 1;  
      unsigned long   nn;
      double          tsv;    
      nn = klen;
      
      /* copy amp onto series and apply zero padding to  klength */
      amp = x;
      amp.resize( klen );
      // for (i = 0; i < inum; i++) { amp[i] = series[i];  }
      // zero_pad(amp, inum, klength);
        
      /*  Fast Fourier Transform Routine:  here we are using the Numerical Recipes
	  routine jrealft which returns the fft in the 1-D input array
	  packed as pairs of real numbers.
	  The jrealft routine requires the input array to start at index=1
	  so we must decrement the index of amp
      */

      void jrealft(double data[], unsigned long n, int isign);
      double * pamp = &(amp)[0];
      jrealft(pamp-1, nn, isign);
      double anrm = sqrt(npoints/dt);  
      double sum = 0.0;

      // get spectrum from real fourier transform      
      double norm = 1.0/(anrm*anrm);
      
      for(int i=1; i<num_freqs-1; i++)
	{
	  if( 2*i+1 > klen) Helper::halt( "mtm_t error in index");  
	  double sqramp = pow( amp[2*i+1] , 2 ) + pow( amp[2*i] , 2 );
	  std::cout << 2 * norm * (sqramp) << "\n";
	  
	  // sqr_spec[i+kf] = norm * (sqramp);	   == real_FFT()
	  // sum += sqr_spec[i+kf];
	}

      std::cout << "DC " << norm * pow(fabs(amp[0]),2) << "\n"
		<< "NQ " << norm * pow(fabs(amp[1]),2) << "\n";
      
      // sqr_spec[0+kf] = norm * pow(fabs(amp[0]),2);
      // sqr_spec[num_freqs-1+kf] = norm * pow(fabs(amp[1]),2);
  
      std::exit(0);

      // test 2 : real_FFT()
    
      int index_length = x.size();	    
      int my_Fs = 256; // arbitrary

      std::cout << index_length << " is size\n";

      int index_start = 0;
      
      //FFT fftseg( index_length , index_length , my_Fs , FFT_FORWARD , WINDOW_NONE );
      real_FFT fftseg( index_length , index_length , my_Fs , WINDOW_NONE );
      
      const int reps = 5000;

      for ( int i=0; i<reps; i++)
	{
	  std::cout << "i\t" << i << "\n";
	  fftseg.apply( &(x[index_start]) , index_length );
	  
	  int my_N = fftseg.cutoff;
	  
	  // for (int f=0;f<my_N;f++)
	  //   {	  
	  //     // t : raw transform
	  //     // t2 : scaled by 1/N
	  //     // 
	  //     std::cout << f << "\t"
	  // 		<< fftseg.frq[f] << "\t"		    
	  // 		<< fftseg.X[f] << "\n";

	  //   }
	}
      std::exit(1);
    }


  // Detrended fluctuation analysis
  
  if ( p == "dfa" )
    {
      std::vector<double> w(0);
      dfa_t dfa;
      int nw = 100;
      dfa.set_windows( 200 );
      
      dfa.proc( &x );
      
      for (int i=0; i<nw; i++)
	std::cout << dfa.w[i] << "\t"
		  << dfa.fluctuations[i] << "\t"
		  << dfa.slopes[i] << "\n";
      
      std::exit(1);
    }

  // moving average w/ triangular window

  if ( p == "tri" )
    {
      int n = x.size();
      int h = 7 ;
      double w = 0.05;

      Eigen::VectorXd Y = eigen_ops::copy_array( x );

      Eigen::VectorXd Y2 = eigen_ops::tri_moving_average( Y , h , w );
      Eigen::VectorXd Y3 = eigen_ops::moving_average( Y , h );

      for (int i=0; i<n; i++)
	std::cout << Y[i] << "\t" << Y2[i] << "\t" << Y3[i] << "\n";
      
    }


  // decimate signal
  
  if ( p == "decimate" )
    {
      int q = 8;
      int sr = 200;
      const int n = x.size();
      Eigen::VectorXf X = Eigen::VectorXf::Zero( n );
      for (int i=0; i<n; i++) X[i] = x[i];
      
      Eigen::VectorXf Y = segsrv_t::decimate( X , sr , q );

      std::cout << Y << "\n";
      std::exit(0);
    }


  // generic application of FFT 
  
  if ( p == "fft" )
    {
      
      int index_length = x.size();

      int my_Fs = 256; // arbitrary

      if ( p2 != "" )
	if ( ! Helper::str2int( p2 , &my_Fs ) ) Helper::halt( "expecting integer sample rate as second parameter" );

      int index_start = 0;

      FFT fftseg( index_length , index_length , my_Fs , FFT_FORWARD , WINDOW_NONE );
      
      fftseg.apply( &(x[index_start]) , index_length );

      // Extract the raw transform
      std::vector<std::complex<double> > t = fftseg.transform();

      // Extract the raw transform scaled by 1/n
      std::vector<std::complex<double> > t2 = fftseg.scaled_transform();
      
      int my_N = fftseg.cutoff;      
      
      std::cout << "N" << "\t"
		<< "F" << "\t"
		<< "RE" << "\t"
		<< "IM" << "\t"
		<< "UNNORM_AMP" << "\t"
		<< "NORM_AMP" << "\t"
		<< "PSD" << "\t"
		<< "log10(PSD)" << "\n";

      for (int f=0;f<my_N;f++)
	{	  
	  // t : raw transform
	  // t2 : scaled by 1/N
	  // 
	  std::cout << f << "\t"
		    << fftseg.frq[f] << "\t"		    
		    << std::real( t[f] ) << "\t"
		    << std::imag( t[f] ) << "\t"
		    << fftseg.mag[f] << "\t"
		    << ( f == 0 ? 1 : 2 ) * fftseg.mag[f] / (double)index_length << "\t"
		    << fftseg.X[f] << "\t"
		    << log10( fftseg.X[f] ) << "\n";
	} 
      
      std::exit(1);
    }

  // surface laplacian
  
  if ( p == "sl" ) 
    {

      // CLOCLS
      
      clocs_t clocs;
      
      clocs.load_cart( "/Users/shaun/dropbox/projects/ltest/clocs.eegbook" );

      int i = 0 ; 
      signal_list_t signals;
      
      std::ifstream II( "/Users/shaun/dropbox/projects/ltest/clocs.eegbook" );
      while ( ! II.eof() ) 
	{
	  std::string l;
	  double x, y, z;
	  II >> l >> x >> y >> z;
	  if ( II.eof() ) break;
	  signals.add( i++ , l );
	  
	}
      II.close();

      sl_t sl( clocs , signals );
      
      
      // assume 64 channels; rows = channels; cols = time-points
      const int ns = 64;
      const int np = x.size() / ns;
      Eigen::MatrixXd X = Eigen::MatrixXd::Zero( np , ns );
      
      i = 0;
      for (int c=0;c<ns;c++)
	for (int t=0;t<np;t++) 
	  X(t,c) = x[i++];
      
      Eigen::MatrixXd O;

      sl.apply( X , O );
            
    }


  // epoch-level dynamics
  
  if ( p == "dynam" )
    {
      
      dynam_t dynam( x );

      double beta, rsq = 0;

      dynam.linear_trend( &beta , &rsq );
      
      std::cout << "beta = " << beta << "\n";
      std::cout << "rsq = " << rsq  << "\n";
      
      std::exit(0);
    }


  // multiscale entropy
  
  if ( p == "mse" )
    {

      while ( ! std::cin.eof() )
	{
	  double xx;
	  std::cin >> xx;
	  if ( std::cin.eof() ) break;	  
	  x.push_back( xx );	  
	}
      std::cerr << x.size() << " values read\n";
      
      mse_t mse( 1,20,1,2,0.15 );
      
      std::map<int,double> mses = mse.calc( x );

      std::map<int,double>::const_iterator ii = mses.begin();
      while ( ii != mses.end() )
	{
	  std::cout << ii->first << "\t" << ii->second << "\n";
	  ++ii;
	}

      std::exit(1);
    }


  // test of db -> retval mechanism

  if ( p == "db" )
    {
      const std::string db = p2;      
      retval_t retval = writer_t::dump_to_retval( db );
      retval.dump(); 
      std::cout << "\n";
      std::exit(1);
    }


  // XCORR / TSYNC 

  if ( p == "xcorr" )
    {
      // sensor data:
      //  2621    2621   28864 s1.txt
      // 2271    2271   25530 s2.txt

      if ( 1 )
	{
	  if ( x.size() != 2621 + 2271 ) Helper::halt( "expecting matlab sensor data example: 2621 + 2271 elements" );
	  std::vector<double> a(2621) , b(2271);
	  for (int i=0;i<2621;i++) a[i] = x[i];
	  for (int i=0;i<2271;i++) b[i] = x[i+2621];
	  std::cout << " a.size()  " << a.size() <<   " " << b.size() << "\n";
	  
	  xcorr_t xcorr( b , a );

	  std::cout << "main " << xcorr.mx << " " << xcorr.lags[xcorr.mx] << "\t" << xcorr.C[xcorr.mx] << "\n";

	  int n = xcorr.lags.size();
          for (int i=0; i<n; i++)
            std::cout << xcorr.lags[i] << "\t" << xcorr.C[i] << "\n";

	}
      else
	{
	  if ( x.size() != 32 ) Helper::halt( "expecting 2*16-element example" );

	  std::vector<double> a(16) , b(16);
          for (int i=0;i<16;i++) a[i] = x[i];
          for (int i=0;i<16;i++) b[i] = x[i+16];

	  xcorr_t xcorr( a , b );
	  
	  int n = xcorr.lags.size();
	  for (int i=0; i<n; i++)
	    std::cout << xcorr.lags[i] << "\t" << xcorr.C[i] << "\n";
	  
	}
            
    }


  // SSA

  if ( p == "ssa" )
    {
      const int n = x.size();
      ssa_t ssa( &x , 8 ) ;
      std::exit(1);
      
    }
    
  
  // ICA
  
  if ( p == "ica" ) 
    {
      // asssume two signals for now
      const int ns = 2;
      
      int rows = x.size() / ns;
      int cols = ns;

      Eigen::MatrixXd X( rows , cols );

      int p = 0;
      for (int i=0;i<rows;i++) 	
	for (int j=0;j<ns;j++)
	  X(i,j) = x[p++];
      
      int compc = 2;
      
      std::cerr << "performing ICA on " << rows << " x " << cols << " matrix\n";
      
      eigen_ica_t ica( X , compc );

      std::cerr << "K\n" << ica.K << "\n";
      std::cerr << "W\n" << ica.W << "\n";
      std::cerr << "A\n" << ica.A << "\n";
      std::cout <<  ica.S << "\n";
      
      std::exit(1);
    }
  
  
  
  // retval test

  if ( p == "retval" )
    {

      annotation_set_t anns;
      edf_t edf(&anns);
      edf.attach( "/Users/shaun/my-dropbox/my-sleep/Purcell06072016.edf" , "smp" );

      // mimic R leval() behavior
      retval_t retval ;
      
      writer.use_retval( &retval );

      // set command string
      cmd_t cmd( "PSD epoch sig=EEG1 & SPINDLES fc=11,15 sig=EEG1" );
      
      cmd.eval( edf );

      writer.use_retval( NULL );

      retval.dump();
	
      std::exit(1);
    }
  

  // Windows

  if ( p == "windows" )
    {
      const int N = 100;
      std::vector<double> W1(N), W2(N), W3(N);
      W1 = MiscMath::tukey_window(N,0.5);      
      W2 = MiscMath::hann_window(N);    
      W3 = MiscMath::hamming_window(N);
      
      for (int i=0;i<N;i++)
	std::cout << W1[i] << "\t"
		  << W2[i] << "\t"
		  << W3[i] << "\n";
      
      std::exit(1);
    }

  
  // MTM
  
  if ( p == "mtm" )
    {
      const int npi = 5;
      const int nwin = 9;	    
      const double segment_sec = 5;
      const double segment_step = 1;
      
      mtm_t mtm( npi , nwin );
      
      mtm.apply( &x , 256 , 256 * segment_sec , 256 * segment_step , true );
      
      std::cout << mtm.f.size() << "\t" << mtm.spec.size() << "\n";
      
      for (int f=0;f<mtm.f.size();f++)
	std::cout << "MTM\t" << f << "\t" << mtm.f[f] << "\t" << mtm.spec[f] << "\n";
      
      std::exit(0);
    }
  
  
  // TV denoiser

  if ( p == "tv" )
    {
      
      double lambda = 10;

      std::vector<double> y = dsptools::TV1D_denoise_copy( x , lambda );
      
      for (int i=0;i<x.size();i++)
	std::cout << x[i] << "\t" << y[i] << "\n";
      
      std::exit(1);
    }


  // topo functions

  if ( p == "topo" )
    {

      // read map from 'example.topo'   CH  THETA  RAD 
      // read data from std::cin,       CH  VALUE
      
      topo_t topo;      
      int ch = topo.load( "example.topo" );      
      topo.max_radius( 0.55 );      
      topo.grid( 67 , 67 );
      
      std::map<std::string, double> data;  
      
      while ( ! std::cin.eof() ) 
	{
	  std::string l;
	  double z;
	  std::cin >> l;
	  if ( std::cin.eof() ) continue;
	  if ( l == "" ) continue;
	  std::cin >> z ;
	  data[l] = z;	  
	}

      std::cerr << "read topo for " << ch << " channels\n";
      std::cerr << "read data for " <<  data.size() << " channels\n";

      Data::Matrix<double> I = topo.interpolate( data );      
      std::cout << I.dump() << "\n";
      
      std::exit(0);
    }



  if ( p == "clocs" )
    {
      
      clocs_t clocs;
      clocs.load_cart( "ex.clocs" );

      // read data : 64 by 
      int ns = 64;
      int np = 63360;
      
      Eigen::MatrixXd X = Eigen::MatrixXd::Zero( ns , np );
      
      for (int c=0;c<ns;c++) // channel by data points
	for (int r=0;r<np;r++)
	  {
	    double x1;
	    std::cin >> x1;
	    if ( std::cin.eof() ) Helper::halt("prob");
	    X(r,c) = x1;
	  }
      
      signal_list_t good_signals;
      signal_list_t bad_signals;
      
      int si = 0;
      std::ifstream IN1("good.sigs",std::ios::in);
      while ( !IN1.eof() ) 
	{
	  std::string l;
	  IN1 >> l;
	  if ( IN1.eof() ) break;
	  good_signals.add( si++ , l );
	}
      IN1.close();
      
      // bad
      si = 0;
      std::ifstream IN2("bad.sigs",std::ios::in);
      while ( !IN2.eof() ) 
	{
	  std::string l;
	  IN2 >> l;
	  if ( IN2.eof() ) break;
	  bad_signals.add( si++ , l );
	}
      IN2.close();
      
      
      Eigen::MatrixXd invG;
      Eigen::MatrixXd Gi;
      //      clocs.make_interpolation_matrices( good_signals , bad_signals , &invG , &Gi );
      std::vector<int> gi;
      for (int i=11;i<=64;i++) gi.push_back(i-1);

      Eigen::MatrixXd interp = clocs.interpolate( X , gi , invG , Gi );
      
      std::exit(1);
      
    }


  //
  // other misc tests 
  //

  if ( p == "json" )
    {      
      std::ifstream ifs( p2 );
      nlohmann::json j = nlohmann::json::parse(ifs);
      std::cout << j.dump(4) << "\n";
      std::exit(0);
    }

  if ( p == "tps" )
    {
      annotation_set_t anns;
      edf_t edf(&anns);
      const int rs = 1;
      bool okay = edf.init_empty( "id1" , 10000 , 1 , "01.01.85" , "00:00:00" );
      int fs;
      double s1, s2;
      std::cin >> fs >> s1 >> s2;
      // nb. need to fix up nasty floating point issues. e.g. 2.01 --> 2.00999999999  etc
      
      //uint64_t start = s1 * globals::tp_1sec;
      //uint64_t stop = s2 * globals::tp_1sec;

      // nb. use sec2tp() and not direct multiplicaiton...
      // to handle floating point noise
      uint64_t start = Helper::sec2tp( s1 );
      uint64_t stop = Helper::sec2tp( s2 );

      // Helper::sec2tp( "2.00" , &start ) ;
      // Helper::sec2tp( "2.01" , &stop ) ;

      std::cout << " start/stop = " << start <<" " << stop << "\n";


      uint64_t tp1;
      std::string st1 = "2.00";
      if ( Helper::sec2tp( st1 , &tp1 ) ) std::cout << "[" << st1 << "] -> [" << tp1 << "]\n"; 	

      st1 = "2.01";
      if ( Helper::sec2tp( st1,  &tp1 ) ) std::cout << "[" << st1 << "] -> [" << tp1 << "]\n"; 	
      else std::cout <<" prob w/ [" << st1 << "]\n";
      
      st1 = "202012.0192818721";
      if ( Helper::sec2tp( st1,  &tp1 ) ) std::cout << "[" << st1 << "] -> [" << tp1 << "]\n"; 	
      else std::cout <<" prob w/ [" << st1 << "]\n";
      
      st1 = "-1";
      if ( Helper::sec2tp( st1 , &tp1 ) ) std::cout << "[" << st1 << "] -> [" << tp1 << "]\n"; 	
      else std::cout <<" prob w/ [" << st1 << "]\n";
      
      st1 = "";
      if ( Helper::sec2tp( st1 , &tp1 ) ) std::cout << "[" << st1 << "] -> [" << tp1 << "]\n"; 	
      else std::cout <<" prob w/ [" << st1 << "]\n";

      st1 = "A";
      if ( Helper::sec2tp( st1 , &tp1 ) ) std::cout << "[" << st1 << "] -> [" << tp1 << "]\n"; 	
      else std::cout <<" prob w/ [" << st1 << "]\n";
      
      interval_t interval( start , stop ) ;

      int start_record, start_sample;
      int stop_record, stop_sample;
      int n_samples_in_record = fs * rs;
      
      okay = edf.timeline.interval2records( interval_t( start , stop ) ,
					    n_samples_in_record , 
					    &start_record, &start_sample ,
					    &stop_record, &stop_sample );
      
      
      int d = 0;
      int sr = start_record;
      int ss = start_sample;
      while ( 1 )
	{
	  ++d;
	  if ( sr == stop_record && ss == stop_sample )
	    break;

	  ++ss;
	  if ( ss == n_samples_in_record )
	    {
	      ss = 0;
	      ++sr;
	    }
	  
	}

      std::cout << " okay=" << okay << "  out = " << start_record << " " << start_sample << " ... " << stop_record << " " << stop_sample << "\t" << d << "\n";

      std::exit(1);
	      
    }
  
  if ( p == "randomize-kmer" )
    {

      // e.g.
      // awk ' { print $2 } ' seq.1 | sed -e 's/\(.\)/\1\'$'\n/g' | awk ' NF>0 ' | luna -d randomize-kmer > seq.2
      // from an existing sequence.

      // echo "file=seq.1 k=4 nreps=100 w=5" | luna --kmer -o out-local500-v1.db
      // echo "file=seq.2 k=4 nreps=100 w=5" | luna --kmer -o out-local500-v2.db 
      
      std::vector<char> s;
      std::map<char,int> u;
      
      while ( 1 )
	{
	  std::string c;
	  std::cin >> c;
	  if ( c == "" || std::cin.eof() ) break;
	  if ( c.size() != 1 ) break;
	  s.push_back( c[0] );
	  ++u[c[0]];
	}

      const int n = s.size();
      
      std::cerr << " read " << n << " elements\n";
      std::string s1( n , '.' );
      for (int i=0; i<n; i++)
	s1[i] = s[i];
	  
      std::map<char,int>::const_iterator uu = u.begin();
      while ( uu != u.end() )
	{
	  std::cerr << " " << uu->first << " = " << uu->second << "\n";
	  ++uu;
	}

      ms_kmer_t ms1;

      int w = 0;
      if ( p2 != "" )
        if ( ! Helper::str2int( p2 , &w ) )
	  Helper::halt( "expecting integer w as second parameter" );

      std::cerr << " w = " << w << "\n";

      std::string s2 = ms1.modified_random_draw( s1 , w );
      
      std::cout << "ID1\t" << s2 << "\n";
      
      std::exit(0);
    }


  //
  // end of dummy()
  //
      
}
