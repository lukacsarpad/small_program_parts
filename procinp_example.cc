/*********************************************************************************************
 * Procinp example file: demonstrates use of procinp                                         *
 * (c) 2016, D. Berényi, Á. Lukács                                                           *
 *                                                                                           *
 *********************************************************************************************/


#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <memory>

#include <boost/program_options/option.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

#include "procinp.h"

namespace po=boost::program_options;

int main(int argc, char **argv){
  po::options_description disp_opts("Allowed options"), hidden_opts("Hidden options"), all_opts("All options");
    
  std::string input_file, output_file;
  int num;

  disp_opts.add_options()("help,h", "Produce help message")
    ("chelp,c", "Display input file format help mesage");

  hidden_opts.add_options()
    ("input-file,I" , po::value<std::string>(&input_file)->required(), "Input file name");

  all_opts.add(disp_opts);
  all_opts.add(hidden_opts);

  po::positional_options_description pos_opts;
  pos_opts.add("input-file",1);

  po::variables_map cmdline_vars;
  try{
    po::store( po::command_line_parser(argc, argv).options(all_opts).positional(pos_opts).run() , cmdline_vars);
  } catch(po::error &e){
    std::cerr << "Error: " << e.what() << std::endl;
    return -1;
  }

  if( cmdline_vars.count("help") ){
    std::cout << "Usage: proccl [options] <input.in>" << std::endl
	      << "Where" << std::endl
	      << "\t<input.in> is the name of the input file." << std::endl;
      //	      << "\t<output.dat> is the name of the  output file" << std::endl;
    std::cout << disp_opts << std::endl;
    return 1;
  }


  // Define input variables
  std::vector<std::unique_ptr<input_variable> > myvars;
  double beta;
  std::vector<double> myvec, fixvec;
  std::string mytext;
  myvars.push_back( std::unique_ptr< input_variable >( new input_variable_scalar< double >( beta ,  "beta" , "a double called beta" ) ) );
  myvars.push_back( std::unique_ptr< input_variable >( new input_variable_vector< double >( myvec , "v" ,    "a dynamic vector" ) ) );
  myvars.push_back( std::unique_ptr< input_variable >( new input_variable_vector_fixedlength< double >( fixvec , "fv" , "a vector of lenght 5" , 5 ) ) );
  myvars.push_back( std::unique_ptr< input_variable >( new input_variable_scalar< std::string >( mytext, "text", "a string" ) ) );

  if( cmdline_vars.count("chelp") ){
    procinp_help( std::cout , myvars);
    return 1;
  }

  try{
    po::notify( cmdline_vars );
  } catch(const po::required_option & e){
    std::cerr << "Error: " << e.what() << std::endl;
    return -1;
  }
    
  std::cout << "Input file: " << input_file << std::endl;

  std::ifstream input_file_stream( input_file.c_str() );

  if( ! input_file_stream ){
    std::cerr << "Error opening file " << input_file << std::endl;
    return -2;
  }

  // process the input file
  try{
    process_input_file(input_file_stream, myvars);
  } catch( procinp_exception_invalid_name & e ){
    std::cerr << "Error: " << e.what() << ": " << e.which() << " in file " << input_file << ", line no. " << e.where() << std::endl;
    return -2;
  } catch( procinp_exception_variable_not_set & e ){
    std::cerr << "Error: " << e.what() << ": " << e.which() << " in file " << input_file << std::endl;
    return -2;
  } catch( procinp_exception_readfail & e ){
    std::cerr << "Error: " << e.what() << ", " << e.which() << " from file " << input_file << std::endl;
    return -2;
  } catch( procinp_exception_set_again & e){
    std::cerr << "Error: " << e.what() << ", " << e.which() << " in file " << input_file << " in line " << e.where() << std::endl;
    return -2;
  }

  input_file_stream.close();

  // print the variables read
  std::cout << "The value of beta: " << beta << std::endl;
  std::cout << "The vector v: "; for( std::vector<double>::const_iterator i=myvec.cbegin(); i != myvec.cend(); i++ ) std::cout << *i << ( ( i == myvec.cend()-1 ) ? "" : ", " );
  std::cout << std::endl;
  std::cout << "The vector fv: "; for( std::vector<double>::const_iterator i=fixvec.cbegin(); i != fixvec.cend(); i++ ) std::cout << *i << ( ( i == fixvec.cend()-1 ) ? "" : ", ");
  std::cout << std::endl;
  std::cout << "The text variable: " << mytext << std::endl;

  return 0;
}
