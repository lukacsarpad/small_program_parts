/**********************************************************************************************
 * procinp.h                                                                                  *
 * (c) 2016,  D. Berényi, Á. Lukács                                                           *
 * purpose: process an input file of structure key value                                      *
 *                                                                                            *
 **********************************************************************************************/

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <memory>

#ifndef __PROCINP_H

#define __PROCINP_H


#ifdef __PROCINP_SOURCE
void getline_nocomment_continue(std::istream &, std::string &, int &);

#endif


// generic procinp exception class, only used for derived classes
class procinp_exception{
 protected:
  const std::string name;
 public:
  std::string what(){
    return name;
  }
 procinp_exception( const std::string &namex ) : name( namex ) {}
};

// procinp exception that is thrown when a line in the file starts with
// an invalid variable name
class procinp_exception_invalid_name : public procinp_exception{
 protected: 
  const int lineno;
  const std::string invalid_name;
 public:
  int where(){
    return lineno;
  }
  std::string which(){
    return invalid_name;
  }
 procinp_exception_invalid_name( const std::string & invalid_namex, const int linenox ) : procinp_exception( "Invalid variable name" ) , invalid_name( invalid_namex) , lineno( linenox ) { }
};


// procimp exception thrown when a variable was not set in the input file
class procinp_exception_variable_not_set : public procinp_exception {
 protected: 
  const std::string varname;
 public:
  std::string which(){
    return varname;
  }
 procinp_exception_variable_not_set(const std::string varnamex ) : procinp_exception ( "Variable not set" ) , varname( varnamex ) { }
};


// procinp exception thrown when a fixed lenght vector could not be loaded with the desired number of elements (too few/too many)
class procinp_exception_vector_length : public procinp_exception {
 protected: 
  const std::string varname;
 public:
  std::string which(){
    return varname;
  }
 procinp_exception_vector_length( const std::string varnamex ) : procinp_exception ( "Vector length invalid" ) , varname( varnamex ) { }
};

// procinp exception thrown if there is a read error
class procinp_exception_readfail : public procinp_exception {
 protected: 
  const std::string varname;
 public:
  std::string which(){
    return varname;
  }
 procinp_exception_readfail( const std::string varnamex ) : procinp_exception ( "Variable read failed" ) , varname( varnamex ) { }  
};

// procinp exception thrown if a variable is already set, and appears again
class procinp_exception_set_again: public procinp_exception {
 protected: 
  const int lineno;
  const std::string varname;
 public:
  int where(){
    return lineno;
  }
  std::string which(){
    return varname;
  }
 procinp_exception_set_again( const std::string & invalid_namex, const int linenox ) : procinp_exception( "Variable set again" ) , varname( invalid_namex) , lineno( linenox ) { }
};

// generic input variable class, only used from deriving specific input variable classes from it
class input_variable{
 public:
  const std::string name;
  const std::string description;
  virtual void read_value(std::stringstream &) = 0 ;
 input_variable( const std::string & namex ) : name( namex ), description( "" ) {}
 input_variable( const std::string & namex , const std::string & descriptionx ) : name( namex ), description( descriptionx ) {}  

  bool check_read(){
    return loaded;
  }
  std::string print_description(){
    return description;
  }
 protected:
  bool loaded=false;
};


// an input variable class template for loading a scalar input variable
// the type T should have an operator << for loading
template <typename T> class input_variable_scalar : public input_variable {
 public:
  //  using input_variable::name;
  T & value;
  void read_value(std::stringstream & src){
    src >> value;
    if( src.fail() ) throw( procinp_exception_readfail( name ) );
    loaded = true;
  }
 input_variable_scalar( T & target , const std::string & namex ) : input_variable( namex ), value( target ) { }
 input_variable_scalar( T & target , const std::string & namex , const std::string & descriptionx ) : input_variable( namex , descriptionx ), value(target) { }
};


// an input variable class template for loading a vector of dynamic length
// of the type T; T shall have an operator << for loading
template <typename T> class input_variable_vector : public input_variable {
 public:
  std::vector<T> & value;
  void read_value( std::stringstream & src ){
    while( src.good() ){
      T temp;
      src >> temp;
      if( src.fail() ) throw( procinp_exception_readfail( name ) );
      value.push_back( temp );
    }
    loaded = true;
  }
 input_variable_vector( std::vector<T> & target , const std::string & namex ) : input_variable( namex ), value( target ) { }
 input_variable_vector( std::vector<T> & target , const std::string & namex , const std::string & descriptionx ) : input_variable( namex , descriptionx ), value(target) { }  
};

// an input variable class template for loading a vector of given length
// of the type T; T shall have an operator << for loading
template <typename T> class input_variable_vector_fixedlength : public input_variable_vector<T> {
 public:
  using input_variable_vector<T>::value;
  using input_variable_vector<T>::name;
  void read_value( std::stringstream &src ){
    input_variable_vector<T>::read_value( src );
    if( value.size() != length ) throw( procinp_exception_vector_length( name ) );
  }
 protected:
  int length;
 public:
 input_variable_vector_fixedlength( std::vector<T> & target , const std::string & namex , const int lengthx ) : input_variable_vector<T>( target , namex ), length( lengthx ) { }
 input_variable_vector_fixedlength( std::vector<T> & target , const std::string & namex , const std::string & descriptionx , const int lengthx ) : input_variable_vector<T>( target , namex , descriptionx ), length( lengthx ) { }  
};

// Main subroutine to process a complete input file
// The input file is expected in the first variable, as an input stream, the second variable is the vector of unique_ptrs to the input variables
// Throws an exception of not all variables could be read
void process_input_file(std::istream &, std::vector<std::unique_ptr<input_variable> > &);

// This function is used to display a help message on the format of the output file
void procinp_help( std::ostream &, std::vector<std::unique_ptr<input_variable> > & );

#endif
