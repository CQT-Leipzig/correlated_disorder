/*
 * =====================================================================================
 *
 *       Filename:  example_generate_configuration.cpp
 *
 *    Description:  Example program that generates a lattice in two dimension
 *                  with long-range correlated continuous variables according to
 *                  different correlation functions. 
 *                  Afterwards these are mapped to discrete variables to produce
 *                  a desired density.
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

/**
 * Writes @a lattice into file @a filename
 *
 * Assumes thath the values of @a lattice originates from
 * a two dimensional lattice of linear @a length.
 */
static void write2D(const std::vector<double> &lattice, unsigned length,
                    const std::string &filename) {
  assert(length * length == lattice.size());

  std::ofstream out(filename);
  int index=0;
  for (unsigned j = 0; j < length; j++) {
    for (unsigned i = 0; i < length; i++) {
      index = length * j + i;
      out << lattice[index] << " ";
    }
    out << "\n";
  }
}


int main()
{
  int l         = 8; 
  int seed      = 1000;
  double a      = 0.5;
  double sigma2 = 1.0;
  double p_free = 0.5210; // size-dependent percolation threshold
  std::stringstream filebase;
  filebase << "l" << l << "_a" << a << "_p" << p_free << "_seed" << seed;
 
  {
  Disorder disorder(2,l);
  disorder.init_correlation_power_law_euclidean(a);
 
  // generate continuous correlated variables on lattice
  std::vector<double> lattice1(disorder.N,0.0);
  std::vector<double> lattice2(disorder.N,0.0);
  disorder.generate_correlated_continuous(lattice1, lattice2, seed);
  write2D(lattice1, disorder.L, "continuous_"+filebase.str()+".dat");
 
  // map to discrete variables on lattice
  discretize( lattice1, calc_theta_global(p_free,sigma2) );
  write2D(lattice1, disorder.L, "discrete_"+filebase.str()+".dat");
  }

  return EXIT_SUCCESS;
}
