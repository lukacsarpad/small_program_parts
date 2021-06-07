/****************************************************************************************
 * Code to save scalar, vector, and matrix data to HDF5 files                           *
 * Uses HDF5 C++ API                                                                    *
 *                                                                                      *
 * (c) 2020, Á. Lukács                                                                  *
 *                                                                                      *
 ****************************************************************************************/

#include <iostream>
#include <string>
#include <H5Cpp.h>
#include "boost/numeric/ublas/matrix.hpp"
using namespace H5;

// create hdf5 file -- we always add a group "properties", and a string "description"
void HDF5_create(const std::string &filename, const std::string &description){
  Exception::dontPrint();

  H5File outfile(filename, H5F_ACC_TRUNC);

  Group properties = outfile.createGroup("/properties");

  hsize_t description_dims[] = {description.length()};
  DataSpace description_dataspace(1, description_dims);

  DataSet description_dataset = outfile.createDataSet("/properties/description", PredType::C_S1, description_dataspace);

  description_dataset.write(description.c_str(), PredType::C_S1);
}

void HDF5_create_group(const std::string &filename, const std::string &groupname){
  H5File outfile(filename, H5F_ACC_RDWR);

  Group newgroup = outfile.createGroup(groupname);
}

// save: scalar types

void HDF5_save(const std::string &filename, const std::string &varname, const std::string &data){
  H5File outfile(filename, H5F_ACC_RDWR);

  hsize_t data_dims[] = {data.length()};
  DataSpace data_dataspace(1, data_dims);

  DataSet data_dataset = outfile.createDataSet(varname, PredType::C_S1, data_dataspace);

  data_dataset.write(data, PredType::C_S1);
}

void HDF5_save(const std::string &filename, const std::string &varname, const double data){
  H5File outfile(filename, H5F_ACC_RDWR);

  hsize_t data_dims[] = {1};
  DataSpace data_dataspace(1, data_dims);

  DataSet data_dataset = outfile.createDataSet(varname, PredType::NATIVE_DOUBLE, data_dataspace);

  data_dataset.write(&data, PredType::NATIVE_DOUBLE);
}

void HDF5_save(const std::string &filename, const std::string &varname, const std::complex<double> data){
  H5File outfile(filename, H5F_ACC_RDWR);

  hsize_t data_dims[] = {1};
  DataSpace data_dataspace(1, data_dims);

  CompType NATIVE_COMPLEX( sizeof(std::complex<double>) );
  NATIVE_COMPLEX.insertMember( "r", 0, PredType::NATIVE_DOUBLE);
  NATIVE_COMPLEX.insertMember( "i", sizeof(double), PredType::NATIVE_DOUBLE);
  NATIVE_COMPLEX.pack();

  DataSet data_dataset = outfile.createDataSet(varname, NATIVE_COMPLEX, data_dataspace);

  data_dataset.write(&data, NATIVE_COMPLEX);
}

void HDF5_save(const std::string &filename, const std::string &varname, const int data){
  H5File outfile(filename, H5F_ACC_RDWR);

  hsize_t data_dims[] = {1};
  DataSpace data_dataspace(1, data_dims);

  DataSet data_dataset = outfile.createDataSet(varname, PredType::NATIVE_INT, data_dataspace);

  data_dataset.write(&data, PredType::NATIVE_INT);
}

// save: double vector, matrix

void HDF5_save(const std::string &filename, const std::string &varname, const std::vector<double> & data){
  H5File outfile(filename, H5F_ACC_RDWR);

  hsize_t data_dims[] = {data.size()};
  DataSpace data_dataspace(1, data_dims);

  DataSet data_dataset = outfile.createDataSet(varname, PredType::NATIVE_DOUBLE, data_dataspace);

  data_dataset.write(&data[0], PredType::NATIVE_DOUBLE);
}

void HDF5_save(const std::string &filename, const std::string &varname, const boost::numeric::ublas::vector<double> & data){
  H5File outfile(filename, H5F_ACC_RDWR);

  hsize_t data_dims[] = {data.size()};
  DataSpace data_dataspace(1, data_dims);

  DataSet data_dataset = outfile.createDataSet(varname, PredType::NATIVE_DOUBLE, data_dataspace);

  data_dataset.write(&data[0], PredType::NATIVE_DOUBLE);
}

void HDF5_save(const std::string &filename, const std::string &varname, const boost::numeric::ublas::matrix<double> & data){
  H5File outfile(filename, H5F_ACC_RDWR);

  hsize_t data_dims[] = {data.size1(), data.size2()};
  DataSpace data_dataspace(2, data_dims);

  DataSet data_dataset = outfile.createDataSet(varname, PredType::NATIVE_DOUBLE, data_dataspace);

  data_dataset.write(&data(0, 0), PredType::NATIVE_DOUBLE);
}


// save: int vector, matrix

void HDF5_save(const std::string &filename, const std::string &varname, const std::vector<int> & data){
  H5File outfile(filename, H5F_ACC_RDWR);

  hsize_t data_dims[] = {data.size()};
  DataSpace data_dataspace(1, data_dims);

  DataSet data_dataset = outfile.createDataSet(varname, PredType::NATIVE_INT, data_dataspace);

  data_dataset.write(&data[0], PredType::NATIVE_INT);
}

void HDF5_save(const std::string &filename, const std::string &varname, const boost::numeric::ublas::vector<int> & data){
  H5File outfile(filename, H5F_ACC_RDWR);

  hsize_t data_dims[] = {data.size()};
  DataSpace data_dataspace(1, data_dims);

  DataSet data_dataset = outfile.createDataSet(varname, PredType::NATIVE_INT, data_dataspace);

  data_dataset.write(&data[0], PredType::NATIVE_INT);
}

void HDF5_save(const std::string &filename, const std::string &varname, const boost::numeric::ublas::matrix<int> & data){
  H5File outfile(filename, H5F_ACC_RDWR);

  hsize_t data_dims[] = {data.size1(), data.size2()};
  DataSpace data_dataspace(2, data_dims);

  DataSet data_dataset = outfile.createDataSet(varname, PredType::NATIVE_INT, data_dataspace);

  data_dataset.write(&data(0, 0), PredType::NATIVE_INT);
}

// save: complex vector, matrix

void HDF5_save(const std::string &filename, const std::string &varname, const std::vector<std::complex<double> > & data){
  H5File outfile(filename, H5F_ACC_RDWR);

  hsize_t data_dims[] = {data.size()};
  DataSpace data_dataspace(1, data_dims);

  CompType NATIVE_COMPLEX( sizeof( std::complex<double> ) );
  NATIVE_COMPLEX.insertMember( "r", 0, PredType::NATIVE_DOUBLE);
  NATIVE_COMPLEX.insertMember( "i", sizeof(double), PredType::NATIVE_DOUBLE);
  NATIVE_COMPLEX.pack();

  DataSet data_dataset = outfile.createDataSet(varname, NATIVE_COMPLEX, data_dataspace);

  data_dataset.write(&data[0], NATIVE_COMPLEX);
}

void HDF5_save(const std::string &filename, const std::string &varname, const boost::numeric::ublas::vector<std::complex<double> > & data){
  H5File outfile(filename, H5F_ACC_RDWR);

  hsize_t data_dims[] = {data.size()};
  DataSpace data_dataspace(1, data_dims);

  CompType NATIVE_COMPLEX( sizeof( std::complex<double> ) );
  NATIVE_COMPLEX.insertMember( "r", 0, PredType::NATIVE_DOUBLE);
  NATIVE_COMPLEX.insertMember( "i", sizeof(double), PredType::NATIVE_DOUBLE);
  NATIVE_COMPLEX.pack();

  DataSet data_dataset = outfile.createDataSet(varname, NATIVE_COMPLEX, data_dataspace);

  data_dataset.write(&data[0], NATIVE_COMPLEX);
}

void HDF5_save(const std::string &filename, const std::string &varname, const boost::numeric::ublas::matrix<std::complex<double> > & data){
  H5File outfile(filename, H5F_ACC_RDWR);

  hsize_t data_dims[] = {data.size1(), data.size2()};
  DataSpace data_dataspace(2, data_dims);

  CompType NATIVE_COMPLEX( sizeof( std::complex<double> ) );
  NATIVE_COMPLEX.insertMember( "r", 0, PredType::NATIVE_DOUBLE);
  NATIVE_COMPLEX.insertMember( "i", sizeof(double), PredType::NATIVE_DOUBLE);
  NATIVE_COMPLEX.pack();

  DataSet data_dataset = outfile.createDataSet(varname, NATIVE_COMPLEX, data_dataspace);

  data_dataset.write(&data(0, 0), NATIVE_COMPLEX);
}

// read: scalar types

double HDF5_read(const std::string &filename, const std::string &varname){
  H5File infile(filename, H5F_ACC_RDONLY);

  DataSet data_dataset = infile.openDataSet(varname);

  if( data_dataset.getTypeClass() != H5T_FLOAT ) throw("Not a double.");
  DataSpace data_dataspace = data_dataset.getSpace();
  if( data_dataspace.getSimpleExtentNdims() != 1 ) throw("Multidimensional");
  hsize_t dims[1];
  data_dataspace.getSimpleExtentDims( dims );
  if( dims[0] != 1 ) throw("Vector");

  double retval;
  data_dataset.read(&retval, PredType::NATIVE_DOUBLE);

  return retval;
}

int HDF5_read_int(const std::string &filename, const std::string &varname){
  H5File infile(filename, H5F_ACC_RDONLY);

  DataSet data_dataset = infile.openDataSet(varname);

  if( data_dataset.getTypeClass() != H5T_INTEGER ) throw("Not an integer.");
  DataSpace data_dataspace = data_dataset.getSpace();
  if( data_dataspace.getSimpleExtentNdims() != 1 ) throw("Multidimensional");
  hsize_t dims[1];
  data_dataspace.getSimpleExtentDims( dims );
  if( dims[0] != 1 ) throw("Vector");

  int retval;
  data_dataset.read(&retval, PredType::NATIVE_INT);

  return retval;
}

std::complex<double> HDF5_read_complex(const std::string &filename, const std::string &varname){
  H5File infile(filename, H5F_ACC_RDONLY);

  DataSet data_dataset = infile.openDataSet(varname);
  if( data_dataset.getTypeClass() != H5T_COMPOUND ) throw("Type error");


  CompType NATIVE_COMPLEX( sizeof( std::complex<double> ) );
  NATIVE_COMPLEX.insertMember( "r", 0, PredType::NATIVE_DOUBLE);
  NATIVE_COMPLEX.insertMember( "i", sizeof(double), PredType::NATIVE_DOUBLE);
  NATIVE_COMPLEX.pack();

  DataSpace data_dataspace = data_dataset.getSpace();
  if( data_dataspace.getSimpleExtentNdims() != 1 ) throw("Multidimensional");
  hsize_t dims[1];
  data_dataspace.getSimpleExtentDims( dims );
  if( dims[0] != 1 ) throw("Vector");

  std::complex<double> retval;
  data_dataset.read(&retval, NATIVE_COMPLEX);

  return retval;
}


void HDF5_read(const std::string &filename, const std::string &varname, std::string & str){
  H5File infile(filename, H5F_ACC_RDONLY);

  DataSet data_dataset = infile.openDataSet(varname);

  if( data_dataset.getTypeClass() != H5T_STRING ) throw("Not a string.");
  DataSpace data_dataspace = data_dataset.getSpace();
  if( data_dataspace.getSimpleExtentNdims() != 1 ) throw("Not a string");
  hsize_t dims[1];
  data_dataspace.getSimpleExtentDims( dims );

  str = std::string(dims[0], ' ' );
  data_dataset.read(&str[0], PredType::C_S1);
}

// read operation: double vector, matrix

void HDF5_read(const std::string &filename, const std::string &varname, std::vector<double>  & vec){
  H5File infile(filename, H5F_ACC_RDONLY);

  DataSet data_dataset = infile.openDataSet(varname);

  if( data_dataset.getTypeClass() != H5T_FLOAT ) throw("Not a double.");
  DataSpace data_dataspace = data_dataset.getSpace();
  if( data_dataspace.getSimpleExtentNdims() != 1 ) throw("Multidimensional");
  hsize_t dims[1];
  data_dataspace.getSimpleExtentDims( dims );

  vec = std::vector<double>(dims[0]);
  data_dataset.read(&vec[0], PredType::NATIVE_DOUBLE);
}

void HDF5_read(const std::string &filename, const std::string &varname, boost::numeric::ublas::vector<double>  & vec){
  H5File infile(filename, H5F_ACC_RDONLY);

  DataSet data_dataset = infile.openDataSet(varname);

  if( data_dataset.getTypeClass() != H5T_FLOAT ) throw("Not a double.");
  DataSpace data_dataspace = data_dataset.getSpace();
  if( data_dataspace.getSimpleExtentNdims() != 1 ) throw("Multidimensional");
  hsize_t dims[1];
  data_dataspace.getSimpleExtentDims( dims );

  vec = boost::numeric::ublas::vector<double>(dims[0]);
  data_dataset.read(&vec[0], PredType::NATIVE_DOUBLE);
}

void HDF5_read(const std::string &filename, const std::string &varname, boost::numeric::ublas::matrix<double> & mat){
  H5File infile(filename, H5F_ACC_RDONLY);

  DataSet data_dataset = infile.openDataSet(varname);

  if( data_dataset.getTypeClass() != H5T_FLOAT ) throw("Not a double.");
  DataSpace data_dataspace = data_dataset.getSpace();
  if( data_dataspace.getSimpleExtentNdims() != 2 ) throw("Not a matrix");
  hsize_t dims[2];
  data_dataspace.getSimpleExtentDims( dims );

  mat = boost::numeric::ublas::matrix<double>(dims[0], dims[1]);
  data_dataset.read( &mat(0, 0), PredType::NATIVE_DOUBLE);
}


// read operation: integer vector, matrix

void HDF5_read(const std::string &filename, const std::string &varname, std::vector<int> & vec){
  H5File infile(filename, H5F_ACC_RDONLY);

  DataSet data_dataset = infile.openDataSet(varname);

  if( data_dataset.getTypeClass() != H5T_INTEGER ) throw("Not an integer.");
  DataSpace data_dataspace = data_dataset.getSpace();
  if( data_dataspace.getSimpleExtentNdims() != 1 ) throw("Multidimensional");
  hsize_t dims[1];
  data_dataspace.getSimpleExtentDims( dims );

  vec = std::vector<int>(dims[0]);
  data_dataset.read(&vec[0], PredType::NATIVE_INT);
}

void HDF5_read(const std::string &filename, const std::string &varname, boost::numeric::ublas::vector<int> & vec){
  H5File infile(filename, H5F_ACC_RDONLY);

  DataSet data_dataset = infile.openDataSet(varname);

  if( data_dataset.getTypeClass() != H5T_INTEGER ) throw("Not an integer.");
  DataSpace data_dataspace = data_dataset.getSpace();
  if( data_dataspace.getSimpleExtentNdims() != 1 ) throw("Multidimensional");
  hsize_t dims[1];
  data_dataspace.getSimpleExtentDims( dims );

  vec = boost::numeric::ublas::vector<int>(dims[0]);
  data_dataset.read(&vec[0], PredType::NATIVE_INT);
}

void HDF5_read(const std::string &filename, const std::string &varname, boost::numeric::ublas::matrix<int> & mat){
  H5File infile(filename, H5F_ACC_RDONLY);

  DataSet data_dataset = infile.openDataSet(varname);

  if( data_dataset.getTypeClass() != H5T_INTEGER ) throw("Not an integer.");
  DataSpace data_dataspace = data_dataset.getSpace();
  if( data_dataspace.getSimpleExtentNdims() != 2 ) throw("Not a matrix");
  hsize_t dims[2];
  data_dataspace.getSimpleExtentDims( dims );

  mat = boost::numeric::ublas::matrix<int>(dims[0], dims[1]);
  data_dataset.read( &mat(0, 0), PredType::NATIVE_INT);
}


// read operations: complex vector, matrix

void HDF5_read(const std::string &filename, const std::string &varname, std::vector<std::complex<double> > & vec){
  H5File infile(filename, H5F_ACC_RDONLY);

  DataSet data_dataset = infile.openDataSet(varname);

  CompType NATIVE_COMPLEX( sizeof( std::complex<double> ) );
  NATIVE_COMPLEX.insertMember( "r", 0, PredType::NATIVE_DOUBLE);
  NATIVE_COMPLEX.insertMember( "i", sizeof(double), PredType::NATIVE_DOUBLE);
  NATIVE_COMPLEX.pack();

  if( data_dataset.getTypeClass() != H5T_COMPOUND ) throw("Type error.");
  DataSpace data_dataspace = data_dataset.getSpace();
  if( data_dataspace.getSimpleExtentNdims() != 1 ) throw("Multidimensional");
  hsize_t dims[1];
  data_dataspace.getSimpleExtentDims( dims );

  vec = std::vector<std::complex<double> >(dims[0]);
  data_dataset.read(&vec[0], NATIVE_COMPLEX);
}

void HDF5_read(const std::string &filename, const std::string &varname, boost::numeric::ublas::vector<std::complex<double> > & vec){
  H5File infile(filename, H5F_ACC_RDONLY);

  DataSet data_dataset = infile.openDataSet(varname);

  CompType NATIVE_COMPLEX( sizeof( std::complex<double> ) );
  NATIVE_COMPLEX.insertMember( "r", 0, PredType::NATIVE_DOUBLE);
  NATIVE_COMPLEX.insertMember( "i", sizeof(double), PredType::NATIVE_DOUBLE);
  NATIVE_COMPLEX.pack();

  if( data_dataset.getTypeClass() != H5T_COMPOUND ) throw("Type error.");
  DataSpace data_dataspace = data_dataset.getSpace();
  if( data_dataspace.getSimpleExtentNdims() != 1 ) throw("Multidimensional");
  hsize_t dims[1];
  data_dataspace.getSimpleExtentDims( dims );

  vec = boost::numeric::ublas::vector<std::complex<double> >(dims[0]);
  data_dataset.read(&vec[0], NATIVE_COMPLEX);
}

void HDF5_read(const std::string &filename, const std::string &varname, boost::numeric::ublas::matrix<std::complex<double> > & mat){
  H5File infile(filename, H5F_ACC_RDONLY);

  DataSet data_dataset = infile.openDataSet(varname);

  CompType NATIVE_COMPLEX( sizeof( std::complex<double> ) );
  NATIVE_COMPLEX.insertMember( "r", 0, PredType::NATIVE_DOUBLE);
  NATIVE_COMPLEX.insertMember( "i", sizeof(double), PredType::NATIVE_DOUBLE);
  NATIVE_COMPLEX.pack();

  if( data_dataset.getTypeClass() != H5T_COMPOUND ) throw("Type error.");
  DataSpace data_dataspace = data_dataset.getSpace();
  if( data_dataspace.getSimpleExtentNdims() != 2 ) throw("Not a matrix");
  hsize_t dims[2];
  data_dataspace.getSimpleExtentDims( dims );

  mat = boost::numeric::ublas::matrix<std::complex<double> >(dims[0], dims[1]);
  data_dataset.read( &mat(0, 0), NATIVE_COMPLEX);
}


// replace operations

void HDF5_replace(const std::string & filename, const std::string &varname, const double d){
  H5File outfile(filename, H5F_ACC_RDWR);

  DataSet data_dataset = outfile.openDataSet(varname);

  if( data_dataset.getTypeClass() != H5T_FLOAT ) throw("Not a double.");
  DataSpace data_dataspace = data_dataset.getSpace();
  if( data_dataspace.getSimpleExtentNdims() != 1 ) throw("Multidimensional");
  hsize_t dims[1];
  data_dataspace.getSimpleExtentDims( dims );
  if( dims[0] != 1 ) throw("Vector");

  data_dataset.write(&d, PredType::NATIVE_DOUBLE);
}

void HDF5_replace(const std::string & filename, const std::string &varname, const int i){
  H5File outfile(filename, H5F_ACC_RDWR);

  DataSet data_dataset = outfile.openDataSet(varname);

  if( data_dataset.getTypeClass() != H5T_INTEGER ) throw("Not an integer.");
  DataSpace data_dataspace = data_dataset.getSpace();
  if( data_dataspace.getSimpleExtentNdims() != 1 ) throw("Multidimensional");
  hsize_t dims[1];
  data_dataspace.getSimpleExtentDims( dims );
  if( dims[0] != 1 ) throw("Vector");

  data_dataset.write(&i, PredType::NATIVE_INT);
}

void HDF5_replace(const std::string & filename, const std::string &varname, const std::complex<double> c){
  H5File outfile(filename, H5F_ACC_RDWR);

  DataSet data_dataset = outfile.openDataSet(varname);

  CompType NATIVE_COMPLEX( sizeof( std::complex<double> ) );
  NATIVE_COMPLEX.insertMember( "r", 0, PredType::NATIVE_DOUBLE);
  NATIVE_COMPLEX.insertMember( "i", sizeof(double), PredType::NATIVE_DOUBLE);
  NATIVE_COMPLEX.pack();

  if( data_dataset.getTypeClass() != H5T_COMPOUND ) throw("Type error.");
  DataSpace data_dataspace = data_dataset.getSpace();
  if( data_dataspace.getSimpleExtentNdims() != 1 ) throw("Multidimensional");
  hsize_t dims[1];
  data_dataspace.getSimpleExtentDims( dims );
  if( dims[0] != 1 ) throw("Vector");

  data_dataset.write(&c, NATIVE_COMPLEX);
}

void HDF5_replace(const std::string & filename, const std::string &varname, const std::vector<double> &v){
  H5File outfile(filename, H5F_ACC_RDWR);

  DataSet data_dataset = outfile.openDataSet(varname);

  if( data_dataset.getTypeClass() != H5T_FLOAT ) throw("Not a double.");
  DataSpace data_dataspace = data_dataset.getSpace();
  if( data_dataspace.getSimpleExtentNdims() != 1 ) throw("Multidimensional");
  hsize_t dims[1];
  data_dataspace.getSimpleExtentDims( dims );
  if( dims[0] != v.size() ) throw("Size mismatch.");

  data_dataset.write(&v[0], PredType::NATIVE_DOUBLE);
}

void HDF5_replace(const std::string & filename, const std::string &varname, const boost::numeric::ublas::vector<double> &v){
  H5File outfile(filename, H5F_ACC_RDWR);

  DataSet data_dataset = outfile.openDataSet(varname);

  if( data_dataset.getTypeClass() != H5T_FLOAT ) throw("Not a double.");
  DataSpace data_dataspace = data_dataset.getSpace();
  if( data_dataspace.getSimpleExtentNdims() != 1 ) throw("Multidimensional");
  hsize_t dims[1];
  data_dataspace.getSimpleExtentDims( dims );
  if( dims[0] != v.size() ) throw("Size mismatch.");

  data_dataset.write(&v[0], PredType::NATIVE_DOUBLE);
}

void HDF5_replace(const std::string & filename, const std::string &varname, const boost::numeric::ublas::matrix<double> &m){
  H5File outfile(filename, H5F_ACC_RDWR);

  DataSet data_dataset = outfile.openDataSet(varname);

  if( data_dataset.getTypeClass() != H5T_FLOAT ) throw("Not a double.");
  DataSpace data_dataspace = data_dataset.getSpace();
  if( data_dataspace.getSimpleExtentNdims() != 2 ) throw("Not a matrix.");
  hsize_t dims[2];
  data_dataspace.getSimpleExtentDims( dims );
  if( dims[0] != m.size1() || dims[1] != m.size2() ) throw("Size mismatch.");

  data_dataset.write(&m(0, 0), PredType::NATIVE_DOUBLE);
}

void HDF5_replace(const std::string & filename, const std::string &varname, const std::vector<int> &v){
  H5File outfile(filename, H5F_ACC_RDWR);

  DataSet data_dataset = outfile.openDataSet(varname);

  if( data_dataset.getTypeClass() != H5T_INTEGER ) throw("Not an integer.");
  DataSpace data_dataspace = data_dataset.getSpace();
  if( data_dataspace.getSimpleExtentNdims() != 1 ) throw("Multidimensional");
  hsize_t dims[1];
  data_dataspace.getSimpleExtentDims( dims );
  if( dims[0] != v.size() ) throw("Size mismatch.");

  data_dataset.write(&v[0], PredType::NATIVE_INT);
}

void HDF5_replace(const std::string & filename, const std::string &varname, const boost::numeric::ublas::vector<int> &v){
  H5File outfile(filename, H5F_ACC_RDWR);

  DataSet data_dataset = outfile.openDataSet(varname);

  if( data_dataset.getTypeClass() != H5T_INTEGER ) throw("Not an integer.");
  DataSpace data_dataspace = data_dataset.getSpace();
  if( data_dataspace.getSimpleExtentNdims() != 1 ) throw("Multidimensional");
  hsize_t dims[1];
  data_dataspace.getSimpleExtentDims( dims );
  if( dims[0] != v.size() ) throw("Size mismatch.");

  data_dataset.write(&v[0], PredType::NATIVE_INT);
}

void HDF5_replace(const std::string & filename, const std::string &varname, const boost::numeric::ublas::matrix<int> &m){
  H5File outfile(filename, H5F_ACC_RDWR);

  DataSet data_dataset = outfile.openDataSet(varname);

  if( data_dataset.getTypeClass() != H5T_INTEGER ) throw("Not an integer.");
  DataSpace data_dataspace = data_dataset.getSpace();
  if( data_dataspace.getSimpleExtentNdims() != 2 ) throw("Not a matrix.");
  hsize_t dims[2];
  data_dataspace.getSimpleExtentDims( dims );
  if( dims[0] != m.size1() || dims[1] != m.size2() ) throw("Size mismatch.");

  data_dataset.write(&m(0, 0), PredType::NATIVE_INT);
}


void HDF5_replace(const std::string & filename, const std::string &varname, const std::vector<std::complex<double> > &v){
  H5File outfile(filename, H5F_ACC_RDWR);

  DataSet data_dataset = outfile.openDataSet(varname);

  CompType NATIVE_COMPLEX( sizeof( std::complex<double> ) );
  NATIVE_COMPLEX.insertMember( "r", 0, PredType::NATIVE_DOUBLE);
  NATIVE_COMPLEX.insertMember( "i", sizeof(double), PredType::NATIVE_DOUBLE);
  NATIVE_COMPLEX.pack();

  if( data_dataset.getTypeClass() != H5T_COMPOUND ) throw("Type error.");
  DataSpace data_dataspace = data_dataset.getSpace();
  if( data_dataspace.getSimpleExtentNdims() != 1 ) throw("Multidimensional");
  hsize_t dims[1];
  data_dataspace.getSimpleExtentDims( dims );
  if( dims[0] != v.size() ) throw("Size mismatch.");

  data_dataset.write(&v[0], NATIVE_COMPLEX);
}

void HDF5_replace(const std::string & filename, const std::string &varname, const boost::numeric::ublas::vector<std::complex<double> > &v){
  H5File outfile(filename, H5F_ACC_RDWR);

  DataSet data_dataset = outfile.openDataSet(varname);

  CompType NATIVE_COMPLEX( sizeof( std::complex<double> ) );
  NATIVE_COMPLEX.insertMember( "r", 0, PredType::NATIVE_DOUBLE);
  NATIVE_COMPLEX.insertMember( "i", sizeof(double), PredType::NATIVE_DOUBLE);
  NATIVE_COMPLEX.pack();

  if( data_dataset.getTypeClass() != H5T_COMPOUND ) throw("Type error.");
  DataSpace data_dataspace = data_dataset.getSpace();
  if( data_dataspace.getSimpleExtentNdims() != 1 ) throw("Multidimensional");
  hsize_t dims[1];
  data_dataspace.getSimpleExtentDims( dims );
  if( dims[0] != v.size() ) throw("Size mismatch.");

  data_dataset.write(&v[0], NATIVE_COMPLEX);
}

void HDF5_replace(const std::string & filename, const std::string &varname, const boost::numeric::ublas::matrix<std::complex<double> > &m){
  H5File outfile(filename, H5F_ACC_RDWR);

  DataSet data_dataset = outfile.openDataSet(varname);

  CompType NATIVE_COMPLEX( sizeof( std::complex<double> ) );
  NATIVE_COMPLEX.insertMember( "r", 0, PredType::NATIVE_DOUBLE);
  NATIVE_COMPLEX.insertMember( "i", sizeof(double), PredType::NATIVE_DOUBLE);
  NATIVE_COMPLEX.pack();

  if( data_dataset.getTypeClass() != H5T_COMPOUND ) throw("Type error.");
  DataSpace data_dataspace = data_dataset.getSpace();
  if( data_dataspace.getSimpleExtentNdims() != 2 ) throw("Not a matrix.");
  hsize_t dims[2];
  data_dataspace.getSimpleExtentDims( dims );
  if( dims[0] != m.size1() || dims[1] != m.size2() ) throw("Size mismatch.");

  data_dataset.write(&m(0, 0), NATIVE_COMPLEX);
}

