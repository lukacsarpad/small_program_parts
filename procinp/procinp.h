/**********************************************************************************************
 * procinp.h                                                                                  *
 * (c) 2016,  D. Berényi, Á. Lukács                                                           *
 * purpose: process an input file of structure key value                                      *
 *                                                                                            *
 **********************************************************************************************/
#pragma once
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <memory>
#include <exception>
#include <iterator>
#include <algorithm>

#ifndef __PROCINP_H
#define __PROCINP_H

// generic procinp exception class, only used for derived classes
class procinp_exception : public std::exception
{
protected:
	const std::string name;
public:
	char const* what() const noexcept { return name.c_str(); }
	procinp_exception( const std::string& name_in ) : name( name_in ) {}
};

// procinp exception that is thrown when a line in the file starts with
// an invalid variable name
class procinp_exception_invalid_name : public procinp_exception
{
protected:
	const std::string invalid_name;
	const int lineno;
public:
	int         where() const { return lineno; }
	std::string which() const { return invalid_name; }
	procinp_exception_invalid_name( const std::string& invalid_name_in, int lineno_in ) : procinp_exception( "Invalid variable name" ), invalid_name( invalid_name_in ), lineno( lineno_in ) { }
};

// procimp exception thrown when a variable was not set in the input file
class procinp_exception_variable_not_set : public procinp_exception
{
protected: 
	const std::string varname;
public:
	std::string which() const { return varname; }
	procinp_exception_variable_not_set(const std::string& varname_in ) : procinp_exception ( "Variable not set" ), varname( varname_in ) { }
};

// procinp exception thrown when a fixed lenght vector could not be loaded with the desired number of elements (too few/too many)
class procinp_exception_vector_length : public procinp_exception
{
protected: 
	const std::string varname;
public:
	std::string which() const { return varname; }
	procinp_exception_vector_length( const std::string& varname_in ) : procinp_exception ( "Vector length invalid" ), varname( varname_in ) { }
};

// procinp exception thrown if there is a read error
class procinp_exception_readfail : public procinp_exception
{
protected: 
	const std::string varname;
public:
	std::string which() const { return varname; }
	procinp_exception_readfail( const std::string& varname_in ) : procinp_exception ( "Variable read failed" ) , varname( varname_in ) { }
};

// procinp exception thrown if a variable is already set, and appears again
class procinp_exception_set_again: public procinp_exception
{
protected:
	const std::string varname;
	const int lineno;
public:
	int         where() const { return lineno; }
	std::string which() const { return varname; }
	procinp_exception_set_again( const std::string& invalid_name_in, const int lineno_in ) : procinp_exception( "Variable set again" ) , varname( invalid_name_in ) , lineno( lineno_in ) { }
};

// generic input variable class, only used from deriving specific input variable classes from it
class input_variable_base
{
public:
	const std::string name;
	const std::string desc;
protected:
	bool loaded;
public:
	virtual void read_value( std::stringstream& ) = 0;
	input_variable_base( const std::string& name_in ) : name( name_in ), desc(), loaded(false) {}
	input_variable_base( const std::string& name_in , const std::string& description_in ) : name( name_in ), desc( description_in ), loaded(false) { }
	
	bool               is_loaded()   const { return loaded; }
	const std::string& description() const { return desc; }
};

// Holder class for the references to read the data file values into
class variable_binder
{
protected:
	std::vector< std::unique_ptr<input_variable_base> > bindings;
public:
	void add( input_variable_base* data ){ bindings.push_back( std::unique_ptr<input_variable_base>(data) ); }
	input_variable_base* find( const std::string& name ) 
	{
		auto it = std::find_if( bindings.cbegin(), bindings.cend(), [&](std::unique_ptr<input_variable_base> const& ptr){ return ptr->name == name; } );
		if( it == bindings.cend() ){ return nullptr; }
		return it->get();
	}

	void check_all_loaded() const
	{
		std::for_each( bindings.cbegin(), bindings.cend(), [](std::unique_ptr<input_variable_base> const& ptr)
		{
			if( !ptr->is_loaded() ){ throw( procinp_exception_variable_not_set( ptr->name ) ); }
		} );
	}

	void print_variables( std::ostream& sout ) const
	{
		std::for_each( bindings.cbegin(), bindings.cend(), [&](std::unique_ptr<input_variable_base> const& ptr)
		{
			sout << '\t' << ptr->name << ": " << ptr->description() << std::endl;
		} );
	}
};

// an input variable class template for loading a scalar input variable
// the type T should have an operator << for loading
template <typename T> class input_variable_scalar : public input_variable_base
{
public:
  //  using input_variable::name;
	T& value;
	void read_value(std::stringstream& src)
	{
		src >> value;
		if( src.fail() ){ throw( procinp_exception_readfail( name ) ); }
		loaded = true;
	}
	input_variable_scalar( T& target, const std::string& name_in )                                    : input_variable_base( name_in ),                 value( target ) { }
	input_variable_scalar( T& target, const std::string& name_in, const std::string& description_in ) : input_variable_base( name_in, description_in ), value( target ) { }
};

// an input variable class template for loading a vector of dynamic length
// of the type T; T shall have an operator << for loading
template <typename T> class input_variable_vector : public input_variable_base
{
public:
	std::vector<T>& value;
	void read_value( std::stringstream& src )
	{
		while( src.good() )
		{
			T temp;
			src >> temp;
			if( src.fail() ){ throw( procinp_exception_readfail( name ) ); }
			value.push_back( temp );
		}
		loaded = true;
	}
	input_variable_vector( std::vector<T>& target, const std::string& name_in )                                    : input_variable_base( name_in ),                  value( target ) { }
	input_variable_vector( std::vector<T>& target, const std::string& name_in, const std::string& description_in ) : input_variable_base( name_in , description_in ), value( target ) { }
};

// an input variable class template for loading a vector of given length
// of the type T; T shall have an operator << for loading
template <typename T> class input_variable_vector_fixedlength : public input_variable_vector<T>
{
protected:
	unsigned int length;
public:
	using input_variable_vector<T>::value;
	using input_variable_vector<T>::name;
	void read_value( std::stringstream &src )
	{
		input_variable_vector<T>::read_value( src );
		if( value.size() != length ){ throw( procinp_exception_vector_length( name ) ); }
	}
	
	input_variable_vector_fixedlength( std::vector<T>& target , const std::string& name_in,                                    const unsigned int length_in ) : input_variable_vector<T>( target, name_in ),                 length( length_in ) { }
	input_variable_vector_fixedlength( std::vector<T>& target , const std::string& name_in, const std::string& description_in, const unsigned int length_in ) : input_variable_vector<T>( target, name_in, description_in ), length( length_in ) { }
};

namespace procinp_impl
{
	template<typename T> struct input_variable
	{
		T&                data;
		const std::string name;
		const std::string desc;
		int               size;

		input_variable( T& data_in, const std::string& name_in, const std::string& desc_in = std::string(), int size_in = 0 ) : data(data_in), name(name_in), desc(desc_in), size(size_in) { }
	};
}

template<typename T> procinp_impl::input_variable<T> variable( T& data, const std::string& name )                                   { return procinp_impl::input_variable<T>(data, name); }
template<typename T> procinp_impl::input_variable<T> variable( T& data, const std::string& name, const std::string& desc )          { return procinp_impl::input_variable<T>(data, name, desc); }
template<typename T> procinp_impl::input_variable<T> variable( T& data, const std::string& name, const std::string& desc, int size ){ return procinp_impl::input_variable<T>(data, name, desc, size); }

template<typename T>
variable_binder& operator<<( variable_binder& binder, procinp_impl::input_variable<T>&& var ){ binder.add( new input_variable_scalar<T>(var.data, var.name, var.desc) ); return binder; }

template<typename T>
variable_binder& operator<<( variable_binder& binder, procinp_impl::input_variable<std::vector<T>>&& var )
{
	if( var.size != 0 ){ binder.add( new input_variable_vector_fixedlength<T>(var.data, var.name, var.desc, var.size) ); }
	else{                binder.add( new input_variable_vector<T>            (var.data, var.name, var.desc)           ); }
	return binder; 
}

namespace procinp_impl
{
	// A variant of getline: ignore lines starting with a #, concatenate the next line to all lines 
	// ending with a \ and at each getline, increment the lineno variable, used to display error messages
	void getline_nocomment_continue(std::istream& file, std::string& line, int& lineno)
	{
		std::string temp;
		line = "\\";
		while( line[line.size()-1] == '\\' && getline(file, temp) )
		{
			lineno++;
			auto ihash = std::find(temp.begin(), temp.end(), '#');
			if( ihash != temp.end() )
			{
				if( ihash == temp.begin()){}//skip full line
				else
				{
					line.erase(line.size()-1, 1);
					std::copy(temp.begin(), --ihash, std::back_inserter(line));
				}
			}
			else if( temp.size() == 0 ){}//empty line was read
			else
			{
				line.erase(line.size()-1, 1);
				line = line + temp;
			}
		}
		if( line=="\\" ){ line=""; }
	}
	
	// Internal function to process a single line: read the first string from it as a variable name, and give the rest to the
	// << operator of the input variable class
	void process_line( const std::string& line, variable_binder& bindings, const int lineno )
	{
		std::string varname;
		std::stringstream line_stream(line);
		line_stream >> varname;
		
		if( varname != "")
		{
			auto var = bindings.find(varname);
			if( var )
			{
				if( var->is_loaded() ){ throw( procinp_exception_set_again( varname, lineno ) ); }
				var->read_value(line_stream);
			}
			else{ throw( procinp_exception_invalid_name( varname, lineno ) ); }
		}
	}
}

// Main subroutine to process a complete input file
// The input file is expected in the first variable, as an input stream, the second variable is the vector of unique_ptrs to the input variables
// Throws an exception of not all variables could be read
void operator>>(std::istream& file, variable_binder& bindings )
{
	int lineno=0;
	while( file.good() )
	{
		std::string line;
		procinp_impl::getline_nocomment_continue(file, line, lineno);
		procinp_impl::process_line(line, bindings, lineno);
	}

	bindings.check_all_loaded();
}

// This function is used to display a help message on the format of the output file
void procinp_help( std::ostream& sout, const variable_binder& bindings )
{
	sout << "Input file format:" << std::endl
		 << '\t' << "Lines beginning with a # are ignored, used for comments" << std::endl
		 << '\t' << "To all lines ending whith a \\, the next line is appended" << std::endl
		 << '\t' << "<variable_name> <value> pairs denote variables to be loaded" << std::endl
		 << "Allowed variables are:" << std::endl;
	bindings.print_variables(sout);
}

#endif
