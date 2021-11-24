#include <numeric>

#pragma once
#ifndef __KAHANSUM_H_
#define __KAHANSUM_H_

template<typename T> struct KahanAccumulation{
  T sum;
  T correction;
};

template<typename T> KahanAccumulation<T> KahanSum(KahanAccumulation<T> accumulation, const T& value){
  KahanAccumulation<T> result;
  T y = value - accumulation.correction;
  T t = accumulation.sum + y;
  result.correction = (t - accumulation.sum) - y;
  result.sum = t;
  return result;
}

#endif  // __KAHANSUM_H_
