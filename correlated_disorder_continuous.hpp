/*
 * =====================================================================================
 *
 *       Filename:  correlated_disorder_continuous.hpp
 *
 *    Description:  class for correlated disorder with continuous variables on a 
 *                  hypercubic lattice
 *
 *        Version:  1.0
 *        Created:  11/18/2014
 *
 *         Author:  Johannes Zierenberg, Niklas Fricke 
 *
 * =====================================================================================
 */
#pragma once

#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <random>
#include <functional>
#include <algorithm>
#include "fft_nr.h"

/**
 * Class to generated long-range correlated disorder
 *
 * The usage of a class allows to reuse precomputed spectral densities and
 * hence speed up the sequential generation of disorder averages. Before the
 * correlated disorder can be generated via generate_correlated_continous the
 * correlation function must be initialized via init_correlation_power_law_euclidean,
 * init_correlation_power_law_manhattan, or init_correlation_custom.
 */
class Disorder{

 public:
  unsigned dim;                //!< Dimension of the used lattice
  unsigned l;                  //!< Linear lattice size as power of two $(L=2^l)$
  unsigned L;                  //!< Linear lattice size (the lattice contains L*L sites)
  unsigned N;                  //!< Number of lattice sites $(N = L \times L)$
  double inv_N;                //!< Inverse lattice size $( = 1/N )$
  std::vector<int> fft_nn;     //!< fast-fourier transform vector (see Numerical Recipies)
  std::vector<double> Cx;      /**< Correlation function in real-space; will be
                                    "cleared" after sqrtSk computed */
  std::vector<double> sqrtSk;  /**< square-root of the spectral density, the
                                    fourier transform of Cx in frequency space */

 protected:
  bool flag_Sk;  //!< flag that specifies if Sk has been precomputed
  bool flag_Cx;  //!< flag that specifies if Cx has been defined

 public:

  /** Initialize correlation function @a Cx $C(x)=(1+x^2)^(-a/2)$ in euclidean metric */
  void init_correlation_power_law_euclidean(double a) {
    const double minus_a_half = -0.5 * a;

#ifndef SILENT
    std::cout << "... init C(x)=(1+x^2)^(-" << a
              << "/2) in real space with Euclidean distance\n";
#endif

    Cx.resize(N, 0.0);
    for (unsigned index = 0; index < N; index++) {
      // get coordinates in d-dim
      double r2 = 0;
      unsigned index_ = index;

      for (unsigned d = 0; d < dim; d++) {
        // get linear r_i in [0,L)
        double r_i = index_ % L;
        // make symmetric around x=0, i.e. (-L/2,L/2]; x=L/2 is stored in the
        // middle
        r_i = convert_linear_symmetric(r_i);
        // add to squared distance
        r2 += r_i * r_i;
        // modify index_ for next coordinate
        index_ /= L;
      }
      Cx[index] = std::pow(1 + r2, minus_a_half);
    }
    flag_Cx = true;
    flag_Sk = false;
  }

  /** Initialize correlation function @a Cx $C(x)=(1+x^2)^(-a/2)$ in Manhattan metric */
  void init_correlation_power_law_manhattan(double a) {
    const double minus_a_half = -0.5 * a;

#ifndef SILENT
    std::cout << "... init C(x)=(1+x^2)^(-" << a
              << "/2) in real space with Manhattan distance\n";
#endif

    Cx.resize(N, 0.0);
    for (unsigned index = 0; index < N; index++) {
      // get coordinates in d-dim
      double r = 0;
      unsigned index_ = index;
      for (unsigned d = 0; d < dim; d++) {
        // get linear r_i in [0,L)
        double r_i = index_ % L;
        // make symmetric around x=0, i.e. (-L/2,L/2]; x=L/2 is stored in the
        // middle
        r_i = convert_linear_symmetric(r_i);
        // add to squared distance
        r += std::fabs(r_i);
        // modify index_ for next coordinate
        index_ /= L;
      }
      Cx[index] = pow(1 + r * r, minus_a_half);
    }
    flag_Cx = true;
    flag_Sk = false;
  }

  /** Initialize correlation function @a Cx $C(x)=(1+x^2)^(-a/2)$ in euclidean metric*/
  void init_correlation_power_law_euclidean(double a, double alpha) {
    const double alpha_half = 0.5 * alpha;
    const double minus_a_div_alpha = -1.0 * a / alpha;

#ifndef SILENT
    std::cout << "... init C(x)=(1+|x|^" << alpha << ")^(-" << a << "/" << alpha
              << ") in real space with Euclidean distance\n";
#endif

    Cx.resize(N, 0.0);
    for (unsigned index = 0; index < N; index++) {
      // get coordinates in d-dim
      double r2 = 0;
      unsigned index_ = index;

      for (unsigned d = 0; d < dim; d++) {
        // get linear r_i in [0,L)
        double r_i = index_ % L;
        // make symmetric around x=0, i.e. (-L/2,L/2]; x=L/2 is stored in the
        // middle
        r_i = convert_linear_symmetric(r_i);
        // add to squared distance
        r2 += r_i * r_i;
        // modify index_ for next coordinate
        index_ /= L;
      }
      Cx[index] = std::pow(1 + std::pow(r2, alpha_half), minus_a_div_alpha);
    }
    flag_Cx = true;
    flag_Sk = false;
  }

  /** Sets correlation function to @a Cx_ */
  void init_correlation_custom(const std::vector<double>& Cx_) {
    assert(Cx_.size() == Cx.size());

    Cx = Cx_;

    flag_Cx = true;
    flag_Sk = false;
  }

  ////////////////////////////////////////////////////////////////
  ///////////// generate different types of correlated disorder
  ////////////////////////////////////////////////////////////////

  /** Generate uncorrelated disorder */
  void generate_uncorrelated_continuous(std::vector<double>& lattice1,
                                        std::vector<double>& lattice2,
                                        int seed) {
    assert(lattice1.size() == N);
    assert(lattice2.size() == N);

    std::mt19937_64 random(seed);
    std::uniform_real_distribution<double> distribution(0, 1);
    for (unsigned i = 0; i < N; i++) {
      lattice1[i] = distribution(random);
      lattice2[i] = distribution(random);
    }
  }

  /** Generate correlated disorder */
  void generate_correlated_continuous(std::vector<double>& lattice1,
                                      std::vector<double>& lattice2,
                                      int seed) {
    assert(lattice1.size() == N);
    assert(lattice2.size() == N);

    ///////////////////////////////////////////////////////////
    /// precompute spectral density if not already computed
    if (not flag_Sk){
      precompute_sqrt_Sk();
    }

    ////////////////////////////////////////////////////////////
    /// Random number generator
    // requirement on random: <r^2> = N
    // i.e. Gaussian sigma = sqrt(N)
    // uniform (-a,a) a = sqrt(3N)
    double mean = 0.0;
    double sigma = 1.0;
    std::mt19937_64 mt(seed);
    std::normal_distribution<double> gauss(mean, sigma * sqrt(N));
    auto random = std::bind(gauss, mt);

////////////////////////////////////////////////////////////
/// correlation in Fourier space
#ifndef SILENT
    std::cout << "... correlation in Fourier space\n";
#endif

    std::vector<double> complex_lattice(2 * N, 0.0);
    for (unsigned i = 0; i < N; i++) {
      complex_lattice[2 * i] = sqrtSk[i] * random();
      complex_lattice[2 * i + 1] = sqrtSk[i] * random();
    }

///////////////////////////////////////////////////////////
/// Inverse FFT
#ifndef SILENT
    std::cout << "... iFFT and normalization into real number arrays\n";
#endif

    fourn<double>(complex_lattice, fft_nn, -1);

    /// normalize the IFFT and assign to lattices
    for (unsigned i = 0; i < N; i++) {
      lattice1[i] = complex_lattice[2 * i] * inv_N;
      lattice2[i] = complex_lattice[2 * i + 1] * inv_N;
    }
  }

  ////////////////////////////////////////////////////////////////
  ///////////// measurements
  ////////////////////////////////////////////////////////////////

  /** Calculate x-correaltion */
  template <typename T>
  std::vector<double> calc_correlation_x(const T& lattice,
                                         std::vector<int>& lag, double mean) {
    // correlation for lag
    unsigned x_p;
    std::vector<double> C_x(lag.size(), 0.0);
    for (unsigned index = 0; index < N; index++) {
      // get d-dim indices [0,L)
      int x = index % L;
      // calculate correlation for different distances along x-axes
      for (unsigned l = 0; l < lag.size(); l++) {
        x_p = x + lag[l];
        if (x_p >= L) {
          x_p -= L;
        }
        unsigned index_n = index + (x_p - x);

        C_x[l] += (lattice[index] - mean) * (lattice[index_n] - mean);
      }
    }
    for (unsigned l = 0; l < lag.size(); l++) {
      C_x[l] /= static_cast<double>(N);
    }
    return C_x;
  }

  /** Calculate diagonal correlation on only one diagonal direction */
  template <typename T>
  std::vector<double> calc_correlation_diagonal(const T& lattice,
                                                std::vector<int>& lag,
                                                double mean) {
    std::vector<double> C(lag.size(), 0.0);
    std::vector<int> Ld(dim, 1);
    for (unsigned d = 1; d < dim; d++) Ld[d] = Ld[d - 1] * L;
    for (unsigned index = 0; index < lattice.size(); index++) {
      // get d-dim indices [0,L)
      std::vector<int> r(dim, 0);
      int index_ = index;
      for (unsigned d = 0; d < dim; d++) {
        r[d] = index_ % L;
        index_ /= L;
      }
      std::vector<int> r_n(dim, 0);
      // calculate correlation for different diagonal distances (+lag in each
      // dimension)
      for (unsigned l = 0; l < lag.size(); l++) {
        for (unsigned d = 0; d < dim; d++) {
          r_n[d] = r[d] + lag[l];
          if (r_n[d] > (int)L - 1) r_n[d] -= L;
        }
        int index_n = 0;
        for (unsigned d = 0; d < dim; d++) {
          index_n += r_n[d] * Ld[d];
        }
        C[l] += (lattice[index] - mean) * (lattice[index_n] - mean);
      }
    }
    double norm = N;
    for (unsigned l = 0; l < lag.size(); l++) {
      C[l] /= norm;
    }
    return C;
  }

  ///////////////////////////////////////////////////////////////////

  /**
   * Constructs a disorder class
   *
   * @param [in] dim_ Dimension of the lattice
   * @param [in] l_   Lattice size on power of 2 (linear lattice size $L = 2^l$
   */
  Disorder(int dim_, int l_)
      : dim(dim_),
        l(l_),
        L(std::pow(2, l)),
        N(std::pow(L, dim)),
        inv_N(1.0 / static_cast<double>(N)),
        fft_nn(dim, L),
        Cx(N, 0.0),
        sqrtSk(N, 0.0),
        flag_Sk(false),
        flag_Cx(false) {}

  /** Size of the lattice in number of total lattice sites */
  unsigned size() const { return N; }

  /* Clear the correlation function and spectral density */
  void clear() {
    Cx.clear();
    sqrtSk.clear();
  }

  /**
   * Convert @a coord from linear convention $x_i \elem [0,L)$
   * into symmetric convention $r_i \elem [-L/2, L/2]$
   */
  int convert_linear_symmetric(int coord) const {
    return (coord < L / 2. + 1) ? coord : coord - L;
  }

  /**
   * Convert @a coord from symmetric convention $r_i \elem [-L/2, L/2]$
   * into linear convention  $x_i \elem [0,L)$
   */
  int convert_symmetric_linear(int coord) const {
    return (coord >= 0) ? coord : L + coord;
  }


  /** Calculate the square root of the spectral density */
  double precompute_sqrt_Sk(bool error_message = true) {

    assert(flag_Cx);

    // fft complex Cx to complex Sk (only generated here!)
    std::vector<double> Sk(2 * N, 0.0);
    assert(Sk.size() == 2 * Cx.size());

    for (unsigned i = 0; i < N; i++) {
      Sk[2 * i] = Cx[i];
    }
    fourn<double>(Sk, fft_nn, 1);

    // check for negative Sk values and imaginary contributions
    if (error_message) {
      double max_Sk = 0;
      double min_negative = 0;
      double max_imaginary = 1e-7;
      for (unsigned i = 0; i < N; i++) {
        if (Sk[2 * i] > max_Sk) max_Sk = Sk[2 * i];
        if (Sk[2 * i] < min_negative) min_negative = Sk[2 * i];
        if (fabs(Sk[2 * i + 1]) > max_imaginary)
          max_imaginary = fabs(Sk[2 * i + 1]);
      }
      if (min_negative < 0) {
        std::cerr << "WARNING: negative values in S(k) not valid for spectral "
                     "density; min(negative)/max_value = "
                  << min_negative / max_Sk << "\n";
        std::cerr << "         method sets S(k)=0 for S(k)<0 instead;\n";
        std::cerr << "WARNING: Correlation function not numerically exact\n";
      }
      if (max_imaginary > 1e-7) {
        std::cerr << "ERROR: imaginary contribution to S(k) not suitable for "
                     "algorithm but max(imaginary) = "
                  << max_imaginary << ".\n";
        std::cerr << "       try symmmetric C(x) around x=0!!!\n";
        std::cerr << "Abort\n";
        exit(1);
      }
    }
    // take square root and compute true c0=sum(Sk) value
    double c0 = 0;
    for (unsigned i = 0; i < N; i++) {
      if (Sk[2 * i] < 0) {
        sqrtSk[i] = 0;
      } else {
        sqrtSk[i] = sqrt(Sk[2 * i]);
        c0 += Sk[2 * i];
      }
    }

    // clean Cx memory
    Cx.clear();
    std::vector<double>().swap(Cx);
    // set pre-compute flag to true
    flag_Sk = true;
    // return true c0 value
    return c0 / N;
  }
};
