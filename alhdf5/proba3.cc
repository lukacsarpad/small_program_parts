#include "alhdf5.h"
#include <cmath>
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/matrix.hpp"

// Opens a hdf5 file, replaces data in it

int main(){
  double d = -M_PI, e = -2.54;
  HDF5_replace("testfile.h5", "properties/pi", d);

  HDF5_replace("testfile.h5", "my_data/inch", e);


  std::vector<double> vec1(5);
  vec1[0] = 0.0;
  vec1[1] = M_PI;
  vec1[2] = -2.54;
  vec1[3] = -M_SQRT2;
  vec1[4] = M_E;

  HDF5_replace("testfile.h5", "my_data/a_vector", vec1);

  std::vector<std::complex<double> > vec2(2);
  vec2[0] = std::complex<double>(-1, -1);
  vec2[1] = std::complex<double>(-M_PI, M_E);

  HDF5_replace("testfile.h5", "my_data/a_complex_vector", vec2);

  boost::numeric::ublas::matrix<double> m(2, 2);
  m(0, 0) = 0.0;
  m(0, 1) = -M_PI;
  m(1, 0) = 1.25;
  m(1, 1) = -M_E;

  HDF5_replace("testfile.h5", "my_data/a_matrix", m);

  std::vector<int> iv(3);
  iv[0] = 3;
  iv[1] = 2;
  iv[2] = 1;

  HDF5_replace("testfile.h5", "my_data/intvec", iv);

  return 0;
}
