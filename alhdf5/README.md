## ALhdf5

Routines to save matrices and vectors (from std:: and from boost::numeric::ublas) into HDF5 files. At present, they can be used to save and load

 * `int`
 * `double`
 * `std::complex<double>`

 * `std::vector<int>`
 * `std::vector<double>`
 * `std::vector<std::complex<double> >`

 * `boost::numeric::ublas::matrix<int>`
 * `boost::numeric::ublas::matrix<double>`
 * `boost::numeric::ublas::matrix<std::complex<double> >`

 * `std::string`

# DOCUMENTATION

To see functions implemented see alhdf5.h and alhdf5.txt. The latter one is automatically generated form the comments of the source file alhdf5.cpp.

# Requierements

Compilation requieres `g++` (another C++ compiler may also work), `libhdf5-dev` and `libboost-dev`.
