#include <numeric>

#pragma once
#ifndef __NEUMAIERSUM_H_
#define __NEUMAIERSUM_H_

template<typename T> struct NeumaierAccumulation{
  T tempsum;
  T sum;
  T correction;
};

template<typename T> NeumaierAccumulation<T> NeumaierSum(NeumaierAccumulation<T> accumulation, const T& value){
  NeumaierAccumulation<T> result;
  result.tempsum = accumulation.sum + value;
  result.correction = (abs(accumulation.tempsum) > abs(value)) ? accumulation.correction + (accumulation.tempsum - result.tempsum) + value
                                                               : accumulation.correction + (value - result.tempsum) + accumulation.tempsum;
  result.sum = result.tempsum + result.correction;
  return result;
}

#endif  // __NEUMAIERSUM_H_
