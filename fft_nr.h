/************************************************
* FFT code from the book Numerical Recipes in C *
* Visit www.nr.com for the licence.             *
************************************************/

#pragma once

#include <vector>

// The following line must be defined before including math.h to correctly define M_PI
#define _USE_MATH_DEFINES
#include <math.h>
#define PI	M_PI	/* pi to machine precision, defined in math.h */
#define TWOPI	(2.0*PI)

template <typename Type>
void fourn(std::vector<Type> &data, std::vector<int> &nn, const int isign) {
  int idim,i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
  int ibit,k1,k2,n,nprev,nrem,ntot=1,ndim=nn.size();
  double tempi,tempr,theta,wi,wpi,wpr,wr,wtemp;
  for (idim=0;idim<ndim;idim++) ntot *= nn[idim];
  if (ntot<2 || ntot&(ntot-1)) throw("must have powers of 2 in fourn");
  nprev=1;
  for (idim=ndim-1;idim>=0;idim--) {
    n=nn[idim];
    nrem=ntot/(n*nprev);
    ip1=nprev << 1;
    ip2=ip1*n;
    ip3=ip2*nrem;
    i2rev=0;
    for (i2=0;i2<ip2;i2+=ip1) {
      if (i2 < i2rev) {
        for (i1=i2;i1<i2+ip1-1;i1+=2) {
          for (i3=i1;i3<ip3;i3+=ip2) {
            i3rev=i2rev+i3-i2;
            tempr = data[i3];   data[i3]   = data[i3rev];   data[i3rev]   = tempr;
            tempr = data[i3+1]; data[i3+1] = data[i3rev+1]; data[i3rev+1] = tempr;
          }
        }
      }
      ibit=ip2 >> 1;
      while (ibit >= ip1 && i2rev+1 > ibit) {
        i2rev -= ibit;
        ibit >>= 1;
      }
      i2rev += ibit;
    }
    ifp1=ip1;
    while (ifp1 < ip2) {
      ifp2=ifp1 << 1; 
			theta=isign* 2.0 * PI / (ifp2/ip1);
			wpr = cos(theta);
			wpi = sin(theta);
			wr=1.0;
			wi=0.0;
      for (i3=0;i3<ifp1;i3+=ip1) {
        for (i1=i3;i1<i3+ip1-1;i1+=2) {
          for (i2=i1;i2<ip3;i2+=ifp2) {
            k1=i2;
            k2=k1+ifp1;
            tempr=wr*data[k2]-wi*data[k2+1];
            tempi=wr*data[k2+1]+wi*data[k2];
            data[k2]=data[k1]-tempr;
            data[k2+1]=data[k1+1]-tempi;
            data[k1] += tempr;
            data[k1+1] += tempi;
          }
        }
        wtemp=wr;
				wr=wtemp*wpr - wi*wpi;
				wi=wi*wpr    + wtemp*wpi;
      }
      ifp1=ifp2;
    }
    nprev *= n;
  }
}
