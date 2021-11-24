#include <iostream>
#include <iomanip>
#include <vector>
#include <numeric>
#include <cmath>

#include <boost/math/constants/constants.hpp>

#include "neumaiersum.h"

int main(){
  const unsigned int n = 1000;
  std::vector<float> x(n);

  for(unsigned int i = 0; i < x.size(); i++) x[i] = 1.0 / (i + 1);

  NeumaierAccumulation<float> init = {0.0, 0.0};
  NeumaierAccumulation<float> res = std::accumulate(x.begin(), x.end(), init, NeumaierSum<float>);
  float bfres = 0.0;
  for(unsigned int i = 0; i < x.size(); i++) bfres += x[i];

  std::cout << "Neumaier summation, res = " << std::setprecision(15) << std::scientific << res.sum << ", error = " << res.sum - log(static_cast<float>(n)) - boost::math::constants::euler<float>() << std::endl;
  std::cout << "normal addition, res =    " << std::setprecision(15) << std::scientific << bfres << ", error = " << bfres - log(static_cast<float>(n)) - boost::math::constants::euler<float>() << std::endl;

  return 0;
}
