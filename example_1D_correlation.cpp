/*
 * =====================================================================================
 *
 *       Filename:  example_1D_correlation.cpp
 *
 *    Description:  Example program that generates many lattices in one dimension
 *                  with long-range correlated continuous variables according to
 *                  different correlation functions and measures it.
 *
 *        Created:  04/18/2017
 *
 *         Author:  Niklas Fricke, Johannes Zierenberg, 
 *
 * =====================================================================================
 */

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <random>
#include <cassert>
#include "correlated_disorder_continuous.hpp"
#include "correlated_disorder_discrete.hpp"


int main()
{
  int l         = 21; 
  double a      = 0.5;
  std::stringstream filebase;
  
  Disorder disorder(1,l);
  disorder.init_correlation_power_law_euclidean(a);

  std::vector<int> lag;
  lag.push_back(0);
  lag.push_back(1);
  double n=2.;
  while(std::round(n)<=disorder.L/2.0){
    lag.push_back(std::round(n));
    n*=sqrt(2.);
  }
  std::vector<double> C(lag.size(), 0.0);
  std::vector<double> C2(lag.size(), 0.0);
  std::vector<double> err(lag.size(), 0.0);
 
  int sample_size = 100;
  for(int seed=1000; seed<1000+sample_size; seed++){
    // generate continuous correlated variables on lattice
    std::vector<double> lattice1(disorder.N,0.0);
    std::vector<double> lattice2(disorder.N,0.0);
    disorder.generate_correlated_continuous(lattice1, lattice2, seed);

    std::vector<double> c_1 = disorder.calc_correlation_x(lattice1, lag, 0);
    std::vector<double> c_2 = disorder.calc_correlation_x(lattice2, lag, 0);
    for(unsigned i=0; i<lag.size(); i++){
      C[i]  += c_1[i]        + c_2[i];
      C2[i] += c_1[i]*c_1[i] + c_2[i]*c_2[i];
    }
  }

  for(unsigned i=0; i<lag.size(); i++){
    C[i]  /= 2*sample_size;
    C2[i] /= 2*sample_size;
    err[i] = std::sqrt((C2[i]-C[i]*C[i])/2/sample_size);
  }
  
  std::ofstream out;
  filebase << "l" << l << "_a" << a;
  out.open("./correlation_1D_"+filebase.str()+".dat");
  out << std::scientific;
  out << "l\t C_x(l)\t\terr(C_x)\n";
  for(unsigned i=0; i<lag.size(); i++){
    out << lag[i] << "\t" << C[i] << "\t" << err[i] << "\n";
  }
  out.close();


  return EXIT_SUCCESS;
}
