#include <iostream>
#include <array>

#include "pw85/pw85.hpp"

using Vec = std::array<double, pw85::dim>;
using SymMat = std::array<double, pw85::sym>;

int main() {
  Vec x1{-0.5, 0.4, -0.7};
  Vec n1{0., 0., 1.};
  double a1 = 10.;
  double c1 = 0.1;

  Vec x2{0.2, -0.3, 0.4};
  Vec n2{1., 0., 0.};
  double a2 = 0.5;
  double c2 = 5.;

  Vec r12;
  for (int i = 0; i < pw85::dim; i++) r12[i] = x2[i] - x1[i];

  SymMat q1, q2;
  pw85::spheroid(a1, c1, n1.data(), q1.data());
  pw85::spheroid(a2, c2, n2.data(), q2.data());

  std::array<double, 2> out;
  pw85::contact_function(r12.data(), q1.data(), q2.data(), out.data());
  std::cout << "mu^2 = " <<  out[0] << std::endl;
  std::cout << "lambda = " << out[1] << std::endl;
}
