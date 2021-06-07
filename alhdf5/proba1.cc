#include "alhdf5.h"
#include <cmath>
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/matrix.hpp"

// Creates a hdf5 file, deletes it if it exists
// Also writes data to it

int main(){
  HDF5_create("testfile.h5", "hdf5 file written by proba1");

  double d = M_PI, e = M_E;
  HDF5_save("testfile.h5", "properties/pi", d);

  std::complex<double> ii(0, 1);

  HDF5_create_group("testfile.h5", "/my_data");
  HDF5_save("testfile.h5", "my_data/inch", e);
  HDF5_save("testfile.h5", "my_data/a_complex", ii);
  HDF5_save("testfile.h5", "my_data/name", "Juli");

  std::vector<double> vec1(5);
  vec1[0] = 0.0;
  vec1[1] = M_PI;
  vec1[2] = 2.54;
  vec1[3] = M_SQRT2;
  vec1[4] = M_E;

  HDF5_save("testfile.h5", "my_data/a_vector", vec1);

  std::vector<std::complex<double> > vec2(2);
  vec2[0] = std::complex<double>(1, -1);
  vec2[1] = std::complex<double>(M_PI, M_E);

  HDF5_save("testfile.h5", "my_data/a_complex_vector", vec2);

  boost::numeric::ublas::vector<double> vec3(2);
  vec3[0] = 0.0;
  vec3[1] = 12.7;

  HDF5_save("testfile.h5", "my_data/another_vector", vec3);

  boost::numeric::ublas::matrix<double>  m(2, 2);
  m(0, 0) = 0.0;
  m(0, 1) = M_PI;
  m(1, 0) = -1.25;
  m(1, 1) = -M_E;

  HDF5_save("testfile.h5", "my_data/a_matrix", m);

  std::vector<int> iv(3);
  iv[0] = 1;
  iv[1] = 2;
  iv[2] = 3;

  HDF5_save("testfile.h5", "my_data/intvec", iv);
  return 0;
}
