//
// Quasiperiodic Automated Transit Search (QATS) Algorithm
//
// C++ implementation by J. Carter (2012).
//
// If you use this code, please cite Carter & Agol (2012)
// *AND* Kel'manov and Jeon 2004:
//
// \bibitem[Kel'Manov \& Jeon (2004)]{Kel'Manov522004}Kel'Manov, 
//	A.~V.~\& B.~Jeon,  A Posteriori Joint Detection and Discrimination of Pulses in a Quasiperiodic Pulse Train.
//	{\it IEEE Transactions on Signal Processing} {\bf 52}, 645--656 (2004).
//

// 
// Refer to `qats.h' for more detailed function descriptions 
//

#include <math.h>
#include "qats.h"
using namespace std;


#define ROWMINUS1(r)	((r-1)%3 < 0 ? 3+(r-1)%3 : (r-1)%3)


#define PARITY -1 // Parity of box (down = -1, up = 1)

void shConvol(double *  y,int N,int q, double *  & d) {
  // Box convolution (computes D in Eqn. 15 of Carter & Agol 2012).  Assumed d is sized as d[N-q+1]
  // and y sized as y[N].
  for (int n = 0; n <= N-q; n++) {
    d[n] = 0;
    for (int j = 0; j < q; j++) {
      d[n] += PARITY*y[j+n]; 
    }
  }
}



double maximum(double * v, int lb, int rb, int & index) {
  // Maximum of vector v between indices lb and rb.  Index of maximum is index.
  double m = v[lb];
  int i;
  index=lb;
  for (i=lb; i<=rb;i++) {
    if (v[i] > m) {
      m = v[i];
      index = i;
    }
  }
  return m;
}

void omegaBounds(int m, int M, int DeltaMin, int DeltaMax, int N, int q, int & lb, int & rb) {
  //Returns bounds for omega set as listed in paper (Eqns. 16--18 of Carter & Agol 2012).
  lb = ((m-1)*DeltaMin > N-DeltaMax-(M-m)*DeltaMax ? (m-1)*DeltaMin : N-DeltaMax-(M-m)*DeltaMax);
  rb = (DeltaMax-q+(m-1)*DeltaMax < N-q-(M-m)*DeltaMin ? DeltaMax-q+(m-1)*DeltaMax : N-q-(M-m)*DeltaMin);
}

void gammaBounds(int m, int n, int DeltaMin, int DeltaMax, int q, int & lb, int & rb) {
  //Returns bounds for gamma set as listed in paper (Eqns. 19--21 of Carter & Agol 2012).
  lb = ((m-1)*DeltaMin > n-DeltaMax ? (m-1)*DeltaMin : n-DeltaMax);
  rb = (DeltaMax-q+(m-1)*DeltaMax < n-DeltaMin ? DeltaMax-q+(m-1)*DeltaMax : n-DeltaMin);
}
  

void computeSmnRow(double *  d, int M, int m, int DeltaMin, int DeltaMax, int N, int q, double **  & Smn) {
  // Recursive generation of Smn (Eqns. 22,23 of Carter & Agol (2012)) called by computeSmn
  int omegaLb, omegaRb, gammaLb, gammaRb,index;
  if ( m != 1 ) computeSmnRow(d,M,m-1,DeltaMin,DeltaMax,N,q,Smn);
  omegaBounds(m,M,DeltaMin,DeltaMax,N,q,omegaLb,omegaRb);
  for (int n = omegaLb; n <= omegaRb; n++) {
    if ( m == 1 ) {
      Smn[m-1][n] = d[n]; 
    } else {
      gammaBounds(m-1,n,DeltaMin,DeltaMax,q,gammaLb,gammaRb);
      Smn[m-1][n] = d[n]+maximum(Smn[m-2],gammaLb,gammaRb,index);
    }
  }
}

void computeSmn(double *  d, int M, int DeltaMin, int DeltaMax, int N, int q, double **  & Smn, double & Sbest) {
  // Assumed Smn has been allocated as M x N-q+1, d is data, others are from paper (Eqns. 22,23 of Carter & Agol (2012)). 
  // Sbest holds merit for this set of parameters. Memory Intensive; refer to computeSmnInPlace.
  int omegaLb, omegaRb,index;
  computeSmnRow(d,M,M,DeltaMin,DeltaMax,N,q,Smn);
  omegaBounds(M,M,DeltaMin,DeltaMax,N,q,omegaLb,omegaRb);
  Sbest = maximum(Smn[M-1],omegaLb,omegaRb,index)/sqrt(((double) M)*q);
} 

void computeSmnInPlace(double *  d, int M, int DeltaMin, int DeltaMax, int N, int q, double **  & Smn, double & Sbest) {
  // Assumed Smn  has been allocated as 3 x N-q+1. Use this algorithm to compute `Sbest' from Smn 
  // (Eqns. 22,23 of Carter & Agol (2012)) when returned indices are not important.
  // This computation should be used to produce the "spectrum."  Much less memory intensive than computing Smn in full.
  int omegaLb, omegaRb,gammaLb, gammaRb, index,m,n,r = 1;
  for (m = 1; m <= M; m++) {
    omegaBounds(m,M,DeltaMin,DeltaMax,N,q,omegaLb,omegaRb);
    
    if (m != 1) {
      for (n = omegaLb; n <= omegaRb; n++) {
	gammaBounds(m-1,n,DeltaMin,DeltaMax,q,gammaLb,gammaRb);
	Smn[r][n] = d[n]+maximum(Smn[ROWMINUS1(r)],gammaLb,gammaRb,index);
      }
    } else {
      for (n = omegaLb; n <= omegaRb; n++) {
	Smn[r][n] = d[n];
      }
    }
    r = (r+1) % 3;
    
  }
  
  omegaBounds(M,M,DeltaMin,DeltaMax,N,q,omegaLb,omegaRb);
  Sbest = maximum(Smn[ROWMINUS1(r)],omegaLb,omegaRb,index)/sqrt(((double) M)*q);
}

void optIndices(int M, int DeltaMin, int DeltaMax, int N, int q, double ** & Smn, int * & indices) {
  // Optimal starting indices for M, DeltaMin, DeltaMax, q (according to Eqns. 25,26 of Carter & Agol 2012).
  // Given Smn (Eqns. 22,23 of Carter & Agol (2012)).
  // Assumed indices has been allocated as having M values. Called by qats_indices.
  int omegaLb, omegaRb, gammaLb, gammaRb, index;
  omegaBounds(M,M,DeltaMin,DeltaMax,N,q,omegaLb,omegaRb);
  maximum(Smn[M-1],omegaLb,omegaRb,index); indices[M-1] = index;
  for (int m = M-1; m >= 1; m--) {
    gammaBounds(m,indices[m],DeltaMin,DeltaMax,q,gammaLb,gammaRb);
    maximum(Smn[m-1],gammaLb,gammaRb,index); indices[m-1] = index;
  }
}

void qats(double *  d, int DeltaMin, int DeltaMax, int N, int q, double & Sbest, int & MBest) {
  // Determine highest likelihood for this set of input parameters 
  // (Algorithm 1, for a single value of q in Carter & Agol (2012)).  Start indices are not returned.
  // Useful for calculating the QATS "spectrum."
  // Assumed d is pre-convolved data (D in Eqn. 15 of Carter & Agol -- see shConvol above)
  int MMin = floor((N+q-1.0)/DeltaMax);
  int MMax = floor((N-q-0.0)/DeltaMin)+1;

  Sbest = 0; double c_Sbest;
  
  for (int M = MMin; M <= MMax; M++) {
    double **  Smn = new double*[3]; for (int i = 0; i < 3; i++) { Smn[i] = new double[N-q+1]; }
    computeSmnInPlace(d, M, DeltaMin, DeltaMax, N, q, Smn, c_Sbest);
    if (c_Sbest > Sbest) {
      Sbest = c_Sbest;
      MBest = M;
    }
    for (int i=0; i < 3; i++) { delete[] Smn[i];}
    delete [] Smn;
  }
}

void qats_indices(double * d, int M, int DeltaMin, int DeltaMax, int N, int q, double & Sbest, int * & indices) {
  // For a specific collection of M, DeltaMin,DeltaMax, and q find the optimal starting indices of transit.
  // Use after finding the most likely parameters with detectMqNoTimes or even wider search on DeltaMin, DeltaMax or q.
  double ** Smn = new double*[M]; for (int i = 0; i < M; i++) { Smn[i] = new double[N-q+1]; }
  computeSmn(d, M, DeltaMin, DeltaMax, N, q, Smn, Sbest);	    
  indices = new int[M];
  optIndices(M, DeltaMin, DeltaMax, N, q, Smn, indices);
  
  for (int i = 0; i < M; i++) { delete[] Smn[i]; }
  delete [] Smn;
  delete [] d;
}

