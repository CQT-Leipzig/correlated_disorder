# Generating correlated disorder on hypercubic lattices
 
This is a general implementation of the Fourier Filter Method for generating
correlated disorder on a hypercubic lattice. We provide all relevant functions
in header-only c++ files that can be easily included in any existing code. An
example file shows how to generate long-range power-law correlated disorder. 

This implementation was used to produce the data in [1] for C(r)=(1-r^2)^{-a/2}.
Consequently, long-range correlated disorder generated according to the
respective function (see below) should show percolation thresholds and fractal
dimensions according to [1]. However, we do not take responsibility for any
deviations.

Authors: Johannes Zierenberg, Niklas Fricke, Martin Marenz, F. Paul Spitzner,
Viktoria Blavatska and Wolfhard Janke

If you use parts of this code please cite [1].

## Requirements
- C++11 compatible compiler (tested with GNU gcc)
- FFTW3 installed (e.g. on debian 'sudo apt-get install libfftw3 libfftw3-dev')

## Usage
The header-only functions need to be included in a c++ program. The following
brief manual gives a rough idea of the steps 

1.  Steps to generate a d-dimensional lattice of linear size L=2^l with
    long-range power-law correlated site variables:

    1. construct a Disorder object, e.g.
        ```
        Disorder myDisorder(d,l)  
        ```

    1. choose a correlation function:
        * init_correlation_power_law_euclidean(double a): 
        ```
        C(r) = (1+|r|^2)^(-a/2); |r|=sqrt(\sum r_i^2)
        ```
        * init_correlation_power_law_manhattan(double a): 
        ```
        C(r) = (1+|r|^2)^(-a/2); |r|=\sum r_i          
        ```
        * init_correlation_power_law_euclidean(double a, double alpha):
        ```
        C(r) = (1+|r|^{alpha})^(-a/alpha); euclidian 
        ```
        * init_correlation_custom(std::vector<double>& Cx_)
        ```
        C(r) = custom
        ```
    1. if several disorder are required sequentially, then precompute the
         spectral density by calling 
         ```code
         myDisorder.precompute_sqrt_Sk();
         ```
         
    1. each call to generate a correlated disorder yields two independent
       disordered lattices.  Predefine lattices
         ```
         std::vector<double> lattice1(disorder.N,0.0);
         std::vector<double> lattice2(disorder.N,0.0);
         ```
         And generated them
         ```
         myDisorder.generate_correlated_continuous(lattice1, lattice2, seed);
         ```

1.  Steps to generate discrete lattices with p=density of non-defect sites:
    1. Repeat 1)
    1. Be aware that sigma2=C(0) and call
         ```
         discretize( lattice1, calc_theta_global(p,sigma2) );
         ```
         The lattice1 now holds the discretized lattice

## Example

For an example first build via
```bash
make
```
and then executed
```bash
./example
```
The output can be compared to the data provided in 
```bash
./output/
```
where we also provide a gnuplot file to generate the corresponding disorder
illustrations. Again, we provide a Makefile to generate the figures via
```bash
make
```

## File description

### FourierTransfromWrapper.hpp
Wrapper arround fftw.

### correlated_disorder_continuous.hpp
All functions needed to generated continuous variables on a hypercubic lattice
correlated according to predefined and customized correlation functions.

### correlated_disorder_discrete.hpp
Functions needed to map the continuous site variables to discrete ones. This also
includes functions we used to measure the percolation threshold.

### example_generate_configuration.cpp
Example test program which produces the same files as provided in 
```bash
./output/
```

## References
[1] J. Zierenberg, N. Fricke, M. Marenz, F.P. Spitzner, V. Blavatska, and W. Janke,
    "Percolation thresholds and fractal dimensions for square and cubic lattices
    with long-range correlated defects",
    submitted (2017).
    
[2] W. H. Press, S. A. Teukolsky, W. T. Vetterling, and B. P.  Flannery,
    Numerical Recipes 3rd edition: The Art of Scientific Computing 
    (Cambridge University Press, Cambridge, 2007).
