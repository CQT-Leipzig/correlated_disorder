/*
 * =====================================================================================
 *
 *       Filename:  correlated_disorder_discrete.hpp
 *
 *    Description:  class for mapping correlated disorder (generated with help
 *                  of correlated_disorder_continuous.hpp) to discrete variables
 *                  on a hypercubic lattice.
 *                  Includes functions to measure percolation, correlation etc.
 *
 *        Version:  1.0
 *        Created:  11/18/2014
 *
 *         Author:  Niklas Fricke, Martin Marenz, Johannes Zierenberg 
 *
 * =====================================================================================
 */
#pragma once

// macro for free site and defect site on lattice
#define FREE   0
#define DEFECT 1

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>


/** Calculate occupation probability p from the mapping threshold theta  */
inline double calc_p_from_theta(double theta, double sigma2){
  return 0.5 * std::erfc((-theta) / std::sqrt(2.0 * sigma2));
}

/**
 * Calculate the global theta value (as opposed to local theta), defined as
 * fixed disorder threshold theta such that $p(\theta) = p_{\testrm{free}}$, which
 * can be determined by bisection.
 */
double calc_theta_global(double p_free, double sigma2){
  double p_i    = 0;
  double theta  = 0;
  double sigma  = sqrt(sigma2);
  double theta_min=-sigma*5;
  double theta_max= sigma*5;
  // if the free case is to be simulated: set it to free
  if (p_free == 1) {
    theta = sigma * 10;
  } else {
    while (std::fabs(p_free - p_i) > 1e-12) {
      p_i = calc_p_from_theta(theta, sigma2);
      if (p_i < p_free) {
        theta_min = theta;
      }
      if (p_i > p_free) {
        theta_max = theta;
      }
      theta = theta_min + (theta_max - theta_min) / 2.0;
    }
  }
  return theta;
}

/**
 * Measure a local theta (as opposed to global theta), defined as a threshold
 * per disorder realization such that p_disorder(theta)=p_free.
 *  
 * This requires the correlated disorder realization, sorts all values and returns the
 * threshold such that p_free*N sites have values lower than theta.
 */
template <typename T>
double calc_theta_local(const T& lattice, double p_free) {
  std::vector<double> sorted_list(lattice);
  std::sort(sorted_list.begin(), sorted_list.end());
  const int index = static_cast<int>(p_free * lattice.size());
  const double theta = sorted_list[index];
  return theta;
}

/**
 * Map a lattice with continuous (Gaussian) variables to discrete variables 
 * by given threshold theta.
 * 
 * @param [out] lattice Filled with discrete values.
 * @return Density in terms of Free sites (number of free sites devided by lattice size).
 */
template <typename T>
double discretize(T& lattice, double theta) {
  double p_i = 0; 
  for(unsigned i = 0; i< lattice.size(); i++){
    if( lattice[i] > theta){
      lattice[i] = DEFECT;
    }
    else {
      lattice[i] = FREE;
      p_i += 1;
    }
  }
  return p_i/lattice.size();
}


////////////////////////////////////////////////////////////////////
/// cluster identification routines
////////////////////////////////////////////////////////////////////

/**
 * Nearest neighbor in d dimension (2*dim neighbors)
 * 
 * @param [in] index Lattice index of which neighbor should be calculated
 * @param [in] n     Number of neighbor $n \elem [0, 2L)$
 */
inline int get_neighbor(const int index, const int n, const int dim,
                        const int L) {
  assert(n<2*dim);
  // get d-dim indices [0,L)
  int d=n/2;
  int L_d = 1;

  for(int i=0;i<d;i++){
    L_d *= L;
  }

  int coord = (index/L_d) % L;
  int step  = (-1+2*(n%2));
  if(coord+step > L-1) return index + (step - L)*L_d;
  if(coord+step <   0) return index + (step + L)*L_d;
  return index + step*L_d;
}

/** returns the size of the marked cluster  -> good for mass of percolating cluster */
template <typename T>
int mark_cluster(T& lattice, int index, int value, int dim, int L) {
  int reference = lattice[index];
  lattice[index] = value;
  int size = 1;
  for (int n = 0; n < 2 * dim; n++) {
    int index_n = get_neighbor(index, n, dim, L);
    if (lattice[index_n] == reference) {
      size += mark_cluster(lattice, index_n, value, dim, L);
    }
  }
  return size;
}

/** Identify largest cluster and return its size as well as one index
 * 
 * @param [out] lattice
 */
template<typename T>
std::pair<int, int> identify_clusters(T& lattice, int dim, int L){
  int size  = -1;
  int index = -1;
  int id    = -1;
  for(unsigned i=0; i<lattice.size(); i++){
    if(lattice[i] > -1){
      int size_local = mark_cluster(lattice, i, id, dim, L);
      if(size_local > size){
        size  = size_local;
        index = i;
      }
      id--;
    }
  }
  return std::make_pair(size, index);
}

////////////////////////////////////////////////////////////////////
/// percolation routines
////////////////////////////////////////////////////////////////////

/** check crossing only in x-dimension! */
bool check_percolation_here(const std::vector<double>& lattice, int index,
                            int crossed, std::vector<int>& crossed_boundaries,
                            const double threshold, const int dim,
                            const int L) {
  crossed_boundaries[index] = crossed;
  for (int n = 0; n < 2*dim; n++) {
    int index_n = get_neighbor(index, n, dim, L);
    int crossed_n = crossed;
    if (index % L - index_n % L == L - 1)
      ++crossed_n;
    else if (index % L - index_n % L == - L + 1)
      --crossed_n;
    if (lattice[index_n] < threshold) {
      if (crossed_boundaries[index_n] == 10000) {
        if (check_percolation_here(lattice, index_n, crossed_n, crossed_boundaries, threshold, dim, L)) {
          return true;
        }
      } else if (crossed_n != crossed_boundaries[index_n]) {
        return true;
      }
    }  
  }    
  return false;
}

/**
 * check if percolates; 
 * 
 * @return -1 if lattice does not perculate, else a lattice index which
 *         belongs to perculation cluster.
 */
int check_percolation(const std::vector<double>& lattice, double threshold,
                      const int dim, const int L) {
  int perc = 0;
  std::vector<int> crossed_boundaries; 
  crossed_boundaries.resize(lattice.size(), 10000);
  for (size_t i = 0; i < lattice.size(); ++i){
    if ((lattice[i] < threshold) and (crossed_boundaries[i] == 10000)) {
      perc = check_percolation_here(lattice, i, 0, crossed_boundaries,
                                    threshold, dim, L);
      if (perc){
        return i;
      }
    }
  }
  return -1;
}
  
/**
 * Bisection to find percolation threshold of 2D sequence
 * first value, theta_c, second p_c
 */
double find_percolation_treshold(std::vector<double>& lattice, const int dim,
                                 const int L) {
  std::vector<double> sorted_Lattice(lattice);
  std::sort(sorted_Lattice.begin(), sorted_Lattice.end());

  const int N = lattice.size();
  assert(check_percolation(lattice, sorted_Lattice[0], dim, L) == -1);
  assert(check_percolation(lattice, sorted_Lattice[N - 1], dim, L) != -1);

  size_t start_index, end_index;
  size_t check_index = N/2;
  if (check_percolation(lattice, sorted_Lattice[check_index], dim, L) != -1) {
    start_index = 0;
    end_index = check_index;
  }else{
    start_index= check_index;
    end_index = N-1;
  }
  check_index = (end_index + start_index)/2;
  while (not (start_index + 1 == end_index)) {
    assert(start_index < end_index);
    if (check_percolation(lattice, sorted_Lattice[check_index], dim, L) != -1) {
       end_index = check_index;
     } else {
       start_index = check_index;
     }
    check_index = (end_index + start_index) / 2;
  }
  assert(check_percolation(lattice, sorted_Lattice[start_index], dim, L) == -1);
  assert(check_percolation(lattice, sorted_Lattice[end_index], dim, L) != -1);
  double theta_p = sorted_Lattice[end_index];
  return theta_p;
}
