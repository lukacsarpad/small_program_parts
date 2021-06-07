#include<string>
#include<vector>
#include "boost/numeric/ublas/matrix.hpp"


#ifndef ALHDF5_ALHDF5_H_
#define ALHDF5_ALHDF5_H_

void HDF5_create(const std::string &, const std::string &);

void HDF5_create_group(const std::string &, const std::string &);

void HDF5_save(const std::string &, const std::string &, const double);
void HDF5_save(const std::string &, const std::string &, const int);
void HDF5_save(const std::string &, const std::string &, const std::complex<double>);

void HDF5_save(const std::string &, const std::string &, const std::vector<double> &);
void HDF5_save(const std::string &, const std::string &, const boost::numeric::ublas::vector<double> &);
void HDF5_save(const std::string &, const std::string &, const boost::numeric::ublas::matrix<double> &);

void HDF5_save(const std::string &, const std::string &, const std::vector<int> &);
void HDF5_save(const std::string &, const std::string &, const boost::numeric::ublas::vector<int> &);
void HDF5_save(const std::string &, const std::string &, const boost::numeric::ublas::matrix<int> &);

void HDF5_save(const std::string &, const std::string &, const std::vector<std::complex<double> > &);
void HDF5_save(const std::string &, const std::string &, const boost::numeric::ublas::vector<std::complex<double> > &);
void HDF5_save(const std::string &, const std::string &, const boost::numeric::ublas::matrix<std::complex<double> > &);

void HDF5_save(const std::string &, const std::string &, const std::string&);


int HDF5_read_int(const std::string &, const std::string &);
double HDF5_read(const std::string &, const std::string &);
std::complex<double> HDF5_read_complex(const std::string &, const std::string &);


void HDF5_read(const std::string &, const std::string &, std::string&);
void HDF5_read(const std::string &, const std::string &, std::vector<double> &);
void HDF5_read(const std::string &, const std::string &, boost::numeric::ublas::vector<double> &);
void HDF5_read(const std::string &, const std::string &, boost::numeric::ublas::matrix<double> &);
void HDF5_read(const std::string &, const std::string &, std::vector<int> &);
void HDF5_read(const std::string &, const std::string &, boost::numeric::ublas::vector<int> &);
void HDF5_read(const std::string &, const std::string &, boost::numeric::ublas::matrix<int> &);
void HDF5_read(const std::string &, const std::string &, std::vector<std::complex<double> >&);
void HDF5_read(const std::string &, const std::string &, boost::numeric::ublas::vector<std::complex<double> > &);
void HDF5_read(const std::string &, const std::string &, boost::numeric::ublas::matrix<std::complex<double> > &);


void HDF5_replace(const std::string &, const std::string &, const double);
void HDF5_replace(const std::string &, const std::string &, const int);
void HDF5_replace(const std::string &, const std::string &, const std::complex<double>);

void HDF5_replace(const std::string &, const std::string &, const std::vector<double> &);
void HDF5_replace(const std::string &, const std::string &, const boost::numeric::ublas::vector<double> &);
void HDF5_replace(const std::string &, const std::string &, const boost::numeric::ublas::matrix<double> &);

void HDF5_replace(const std::string &, const std::string &, const std::vector<int> &);
void HDF5_replace(const std::string &, const std::string &, const boost::numeric::ublas::vector<int> &);
void HDF5_replace(const std::string &, const std::string &, const boost::numeric::ublas::matrix<int> &);

void HDF5_replace(const std::string &, const std::string &, const std::vector<std::complex<double> > &);
void HDF5_replace(const std::string &, const std::string &, const boost::numeric::ublas::matrix<std::complex<double> > &);

#endif  // ALHDF5_ALHDF5_H_
