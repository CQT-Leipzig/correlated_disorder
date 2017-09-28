/*
 * =====================================================================================
 *
 *       Filename:  FourierTransfromWrapper.hpp
 *
 *    Description:  Compatibility layer between correlation code and FFTW
 *
 *        Version:  1.0
 *        Created:  09/09/2017 09:22:30 AM
 *       Compiler:  gcc
 *
 *         Author:  Martin Marenz 
 *
 * =====================================================================================
 */

#pragma once

#include <fftw3.h>
#include <vector>
#include <algorithm>
#include <functional>
#include <cmath>

// calls fftw, heavily unoptimized
void fourn(std::vector<double> &data, std::vector<int> &nn, const int isign)
{
  assert(data.size() > 0);
  assert(data.size()%2 == 0);

  const int rank = nn.size();
  const size_t num = data.size()/2;

  fftw_complex * fftw_data;
  fftw_complex * fftw_res;
  
  fftw_data = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num);
  fftw_res = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num);

  // plan should be resused if parameter are identical
  // ATENTION: plan must be created before transforming the data
  fftw_plan plan = fftw_plan_dft(rank, nn.data(), fftw_data, fftw_res, isign, FFTW_MEASURE);

  // transform data
  for (size_t i = 0; i < num; ++i){
    const size_t index_r = 2*i;
    const size_t index_i = index_r+1;
    fftw_data[i][0] = data[index_r];
    fftw_data[i][1] = data[index_i];
  }

  fftw_execute(plan);

  // back transformation
  for (size_t i = 0; i < num; ++i){
    const size_t index_r = 2*i;
    const size_t index_i = index_r+1;
    data[index_r] = fftw_res[i][0];
    data[index_i] = fftw_res[i][1];
  }

  fftw_destroy_plan(plan);
  fftw_free(fftw_data);
  fftw_free(fftw_res);
}
