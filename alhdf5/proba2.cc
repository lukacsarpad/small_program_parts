#include <iostream>
#include <iomanip>
#include "alhdf5.h"
#include <cmath>
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/matrix.hpp"
#include "boost/numeric/ublas/io.hpp"


// Opens an hdf5 file, reads data from it

int main(){
  std::string s;
  HDF5_read("testfile.h5", "properties/description", s);
  std::cout << "properties/description: " << s << std::endl;

  double d = HDF5_read("testfile.h5", "my_data/inch");
  std::cout << std::scientific << "variable my_data/inch: " << d << std::endl;

  std::cout << "variables/a_complex: " << HDF5_read_complex("testfile.h5", "my_data/a_complex") << std::endl;

  std::vector<double> v;
  HDF5_read("testfile.h5", "my_data/a_vector", v);
  std::cout << "variable my_data/a_vector: "; for(unsigned int i = 0; i < v.size(); i++) std::cout << v[i] << ", "; std::cout <<std::endl;

  boost::numeric::ublas::vector<double> v2;
  HDF5_read("testfile.h5", "my_data/another_vector", v2);
  std::cout << "variable my_data/another_vector: " << v2 << std::endl;

  std::vector<std::complex<double> > w;
  HDF5_read("testfile.h5", "my_data/a_complex_vector", w);
  std::cout << "variable my_data/a_complex_vector: "; for(unsigned int i = 0; i < w.size(); i++) std::cout << w[i] << ", "; std::cout <<std::endl;

  boost::numeric::ublas::vector<std::complex<double> > w2;
  HDF5_read("testfile.h5", "my_data/a_complex_vector", w2);
  std::cout << "variable my_data/a_complex_vector: " << w2 << std::endl;


  boost::numeric::ublas::matrix<double> m;
  HDF5_read("testfile.h5", "my_data/a_matrix", m);
  std::cout << "variable my_data/a_matrix: " << m << std::endl;

  std::vector<int> iv;
  HDF5_read("testfile.h5", "my_data/intvec", iv);
  std::cout << "variable my_data/intvec: "; for(unsigned int i = 0; i < iv.size(); i++) std::cout << iv[i] << ", "; std::cout <<std::endl;

  return 0;
}
