#include<string>
#include<vector>
#include "boost/numeric/ublas/matrix.hpp"

#pragma once
#ifndef ALHDF5_ALHDF5_H_
#define ALHDF5_ALHDF5_H_

void HDF5_create(const std::string &, const std::string &);

void HDF5_create_group(const std::string &, const std::string &);

void HDF5_save(const std::string &, const std::string &, const double);
void HDF5_save(const std::string &, const std::string &, const int);
void HDF5_save(const std::string &, const std::string &, const unsigned int);
void HDF5_save(const std::string &, const std::string &, const uint64_t);
void HDF5_save(const std::string &, const std::string &, const std::complex<double>);

void HDF5_save(const std::string &, const std::string &, const std::vector<double> &);
void HDF5_save(const std::string &, const std::string &, const boost::numeric::ublas::vector<double> &);
void HDF5_save(const std::string &, const std::string &, const boost::numeric::ublas::matrix<double> &);

void HDF5_save(const std::string &, const std::string &, const std::vector<int> &);
void HDF5_save(const std::string &, const std::string &, const boost::numeric::ublas::vector<int> &);
void HDF5_save(const std::string &, const std::string &, const boost::numeric::ublas::matrix<int> &);

void HDF5_save(const std::string &, const std::string &, const std::vector<unsigned int> &);
void HDF5_save(const std::string &, const std::string &, const boost::numeric::ublas::vector<unsigned int> &);
void HDF5_save(const std::string &, const std::string &, const boost::numeric::ublas::matrix<unsigned int> &);

void HDF5_save(const std::string &, const std::string &, const std::vector<uint64_t> &);
void HDF5_save(const std::string &, const std::string &, const boost::numeric::ublas::vector<uint64_t> &);
void HDF5_save(const std::string &, const std::string &, const boost::numeric::ublas::matrix<uint64_t> &);

void HDF5_save(const std::string &, const std::string &, const std::vector<std::complex<double> > &);
void HDF5_save(const std::string &, const std::string &, const boost::numeric::ublas::vector<std::complex<double> > &);
void HDF5_save(const std::string &, const std::string &, const boost::numeric::ublas::matrix<std::complex<double> > &);

void HDF5_save(const std::string &, const std::string &, const std::string&);


int HDF5_read_int(const std::string &, const std::string &);
unsigned int HDF5_read_unsigned_int(const std::string &, const std::string &);
uint64_t HDF5_read_uint64(const std::string &, const std::string &);
double HDF5_read(const std::string &, const std::string &);
std::complex<double> HDF5_read_complex(const std::string &, const std::string &);


void HDF5_read(const std::string &, const std::string &, std::string&);
void HDF5_read(const std::string &, const std::string &, std::vector<double> &);
void HDF5_read(const std::string &, const std::string &, boost::numeric::ublas::vector<double> &);
void HDF5_read(const std::string &, const std::string &, boost::numeric::ublas::matrix<double> &);
void HDF5_read(const std::string &, const std::string &, std::vector<int> &);
void HDF5_read(const std::string &, const std::string &, boost::numeric::ublas::vector<int> &);
void HDF5_read(const std::string &, const std::string &, boost::numeric::ublas::matrix<int> &);
void HDF5_read(const std::string &, const std::string &, std::vector<unsigned int> &);
void HDF5_read(const std::string &, const std::string &, boost::numeric::ublas::vector<unsigned int> &);
void HDF5_read(const std::string &, const std::string &, boost::numeric::ublas::matrix<unsigned int> &);
void HDF5_read(const std::string &, const std::string &, std::vector<uint64_t> &);
void HDF5_read(const std::string &, const std::string &, boost::numeric::ublas::vector<uint64_t> &);
void HDF5_read(const std::string &, const std::string &, boost::numeric::ublas::matrix<uint64_t> &);
void HDF5_read(const std::string &, const std::string &, std::vector<std::complex<double> >&);
void HDF5_read(const std::string &, const std::string &, boost::numeric::ublas::vector<std::complex<double> > &);
void HDF5_read(const std::string &, const std::string &, boost::numeric::ublas::matrix<std::complex<double> > &);


void HDF5_replace(const std::string &, const std::string &, const double);
void HDF5_replace(const std::string &, const std::string &, const int);
void HDF5_replace(const std::string &, const std::string &, const unsigned int);
void HDF5_replace(const std::string &, const std::string &, const uint64_t);
void HDF5_replace(const std::string &, const std::string &, const std::complex<double>);

void HDF5_replace(const std::string &, const std::string &, const std::vector<double> &);
void HDF5_replace(const std::string &, const std::string &, const boost::numeric::ublas::vector<double> &);
void HDF5_replace(const std::string &, const std::string &, const boost::numeric::ublas::matrix<double> &);

void HDF5_replace(const std::string &, const std::string &, const std::vector<int> &);
void HDF5_replace(const std::string &, const std::string &, const boost::numeric::ublas::vector<int> &);
void HDF5_replace(const std::string &, const std::string &, const boost::numeric::ublas::matrix<int> &);

void HDF5_replace(const std::string &, const std::string &, const std::vector<unsigned int> &);
void HDF5_replace(const std::string &, const std::string &, const boost::numeric::ublas::vector<unsigned int> &);
void HDF5_replace(const std::string &, const std::string &, const boost::numeric::ublas::matrix<unsigned int> &);

void HDF5_replace(const std::string &, const std::string &, const std::vector<uint64_t> &);
void HDF5_replace(const std::string &, const std::string &, const boost::numeric::ublas::vector<uint64_t> &);
void HDF5_replace(const std::string &, const std::string &, const boost::numeric::ublas::matrix<uint64_t> &);

void HDF5_replace(const std::string &, const std::string &, const std::vector<std::complex<double> > &);
void HDF5_replace(const std::string &, const std::string &, const boost::numeric::ublas::matrix<std::complex<double> > &);

bool HDF5_exists(const std::string &, const std::string &);

#endif  // ALHDF5_ALHDF5_H_
