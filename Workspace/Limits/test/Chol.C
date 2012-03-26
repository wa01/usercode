#include "TFile.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompChol.h"
#include "TRandom2.h"
#include <vector>
#include <iostream>
//
// return Cholesky decomposition of a vector of std.dev. / correlations
//

TMatrixD Chol (TVectorD& covV)
{
  int nCov = covV.GetNrows();
  int n = int((sqrt(8*nCov+1.)-1.)/2.+0.5);
  if ( nCov != n*(n+1)/2. ) {
    std::cout << "Chol: length of vector " << nCov << " is not n*(n+1)/2" << std::endl;
    return TMatrixD();
  }
  

  // get diagonal elements
  int ind(0);
  TVectorD sigmas(n);
  for ( int i=0; i<n; ++i ) {
    for ( int j=0; j<=i; ++j ) {
      if ( j == i )  sigmas[i] = covV(ind);
      ++ind;
    }
  }
  // fill cov matrix (could be more elegant ...)
  ind = 0;
  TMatrixDSym covMatrix(n);
  for ( int i=0; i<n; ++i ) {
    for ( int j=0; j<=i; ++j ) {
      if ( j == i )
	covMatrix(i,i) = sigmas(i)*sigmas(i);
      else
	covMatrix(i,j) = covMatrix(j,i) = covV(ind)*sigmas(i)*sigmas(j);
      ++ind;
    }
  }
  covMatrix.Print();
  
  TDecompChol tdc(covMatrix);
  bool worked = tdc.Decompose();
  if ( !worked ) {
    std::cout << "Decomposition failed" << std::endl;
    return TMatrixD();
  }
  
  TMatrixD matU = tdc.GetU();
  return matU;

//   //
//   // cross check with random generation
//   //  
//   double sum0(0.);
//   TVectorD sum1(n);
//   TMatrixDSym sum2(n);


//   TRandom2 rgen;
//   TVectorD xrnd(n);
//   TVectorD xrndRot(n);
//   for ( unsigned int i=0; i<1000000; ++i ) {
//     for ( unsigned int j=0; j<n; ++j )  xrnd(j) = 0.;
//     for ( unsigned int j=0; j<n; ++j ) {
//       TVectorD aux(n);
//       for ( int k=0; k<n; ++k )  aux(k) = matU(j,k);
//       xrnd += rgen.Gaus(0.,1.)*aux;
//     }
// //       xrnd *= matUT;
//     sum0 += 1.;
//     for ( unsigned int j0=0; j0<n; ++j0 ) {
//       sum1(j0) += xrnd(j0);
//       for ( unsigned int j1=0; j1<n; ++j1 ) {
// 	sum2(j0,j1) += xrnd(j0)*xrnd(j1);
//       }
//     }
//   }
//   for ( unsigned int j0=0; j0<n; ++j0 ) {
//     printf("%10.3g",sum1(j0)/sum0);
//   }
//   printf("  sum1 \n");
//   printf("\n");
//   for ( unsigned int j0=0; j0<n; ++j0 ) {
//     for ( unsigned int j1=0; j1<n; ++j1 ) {
//       printf("%10.3g",sum2(j0,j1)/sum0);
//     }
//     printf(" sum2 \n");
//   }
//   return matU;

}

TMatrixD Chol (TVectorD& covV, TVectorD& newSig)
{
  int nCov = covV.GetNrows();
  int n = newSig.GetNrows();
  std::cout << nCov << " " << n << std::endl;
  if ( nCov != n*(n+1)/2. ) {
    std::cout << "vecTest: mismatch in inputs" << std::endl;
    return TMatrixD();
  }
  //
  // create modified vector (replacing std.dev.s)
  //
  TVectorD newCovV(covV);
  int ind(0);
  for ( int i=0; i<n; ++i ) {
    for ( int j=0; j<=i; ++j ) {
      if ( j==i )  newCovV[ind] = newSig(i);
      ++ind;
    }
  }
  return Chol(newCovV,newSig);
}
