/**********************************************************************************************
 * procinp.cc                                                                                 *
 * (c) 2016,  D. Berényi, Á. Lukács                                                           *
 * purpose: process an input file of structure key value                                      *
 *                                                                                            *
 **********************************************************************************************/

#define __PROCINP_SOURCE
#include "procinp.h"


// A variant of getline: ignore lines starting with a #, concatenate the next line to all lines 
// ending with a \ and at each getline, increment the lineno variable, used to display error messages
void getline_nocomment_continue(std::istream &file, std::string &line, int & lineno){
  std::string temp;
  line = "\\";
  while(  line[line.size()-1] == '\\' && getline(file, temp) ){
    lineno++;
    if( temp[0] != '#' ){
      line.erase(line.size()-1, 1);
      line = line + temp;
    }
  }
  if ( line=="\\" ) line="";
}


// Internal function to process a single line: read the first string from it as a variable name, and give the rest to the
// << operator of the input variable class
void process_line( std::string line, std::vector<std::unique_ptr<input_variable> > &input_variables, const int lineno){
  std::string varname;

  std::stringstream line_stream(line);
  line_stream >> varname;

  if( varname != ""){
    bool found = false;

    for( std::vector<std::unique_ptr<input_variable> >::const_iterator i = input_variables.cbegin() ; i != input_variables.cend() ; i++){
      if( (*i)->name==varname ) {
	if( (*i)->check_read() ) throw( procinp_exception_set_again( varname, lineno ) );
	(*i)->read_value(line_stream);
	found = true;
	break;
      }
    }
    if( ! found ){
      throw( procinp_exception_invalid_name( varname, lineno ) );
    }
  }
}


// Main subroutine to process a complete input file
// The input file is expected in the first variable, as an input stream, the second variable is the vector of unique_ptrs to the input variables
// Throws an exception of not all variables could be read
void process_input_file(std::istream &file, std::vector<std::unique_ptr<input_variable> > & input_variables ){
  int lineno=0;
  
  while( file.good() ) {
    std::string line;
    getline_nocomment_continue(file, line, lineno);

    process_line(line, input_variables, lineno);
  }

  for( std::vector<std::unique_ptr<input_variable> >::const_iterator i = input_variables.cbegin() ; i != input_variables.cend() ; i++){
    if( ! (*i)->check_read() ){
      throw( procinp_exception_variable_not_set( (*i)->name ) );
    }
  }
}



// This function is used to display a help message on the format of the output file
void procinp_help( std::ostream & out, std::vector<std::unique_ptr<input_variable> > & input_variables ){
  out << "Input file format:" << std::endl
      << '\t' << "Lines beginning with a # are ignored, used for comments" << std::endl
      << '\t' << "To all lines ending whith a \\, the next line is appended" << std::endl
      << '\t' << "<variable_name> <value> pairs denote variables to be loaded" << std::endl
      << "Allowed variables are:" << std::endl;
  for( std::vector<std::unique_ptr<input_variable> >::const_iterator i = input_variables.cbegin() ; i != input_variables.cend() ; i++){
    out << '\t' << (*i)->name << ": " << (*i)->print_description() << std::endl;
  }
}

