//
// Quasiperiodic Automated Transit Search (QATS) Algorithm/ Example QATS call from C++.
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

#include <iostream>
#include <fstream>
#include <cstdlib>
#include "qats.h"
using namespace std;

int main(int argc,char * argv[]) {
  if (argc != 6) {
    cerr << "Call as: call_qats_indices [data sequence filename] [M] [DeltaMin] [DeltaMax] [q]\n";
    exit(1);
  }
 
  ifstream in;
  int N = 0, M, DeltaMin, DeltaMax, q;
  int * indices;
  
  double SBest;

  M = atoi(argv[2]);
  DeltaMin = atoi(argv[3]);
  DeltaMax = atoi(argv[4]);
  q = atoi(argv[5]);  

  in.open(argv[1],ios_base::in);

  char input[100];
  if (in.is_open()) {
    while (!in.eof()) {
      in >> input;
      N++;
    }
  } else { cerr << "File: " << argv[1] << " not found." << endl; exit(1); }
  in.close();
  in.open(argv[1],ios_base::in);
  double * y = new double[N];

  for (int i = 0; i < N; i++) {
    in >> input;
    y[i] = atof(input);
  }
  in.close();

  double * d = new double[N-q+1];

  shConvol(y,N,q, d); // Calculate D_n
  qats_indices(d,M,DeltaMin,DeltaMax,N,q,SBest,indices);
  
  cout << "# QATS; Optimal Indices " << endl;
  cout << "# Input sequence from " << argv[1] << endl;
  cout << "# DeltaMin = " << DeltaMin << endl;
  cout << "# DeltaMax = " << DeltaMax << endl;
  cout << "# q = " << q << endl;
  cout << "# M = " << M << endl;
  cout << "# S_Best = "<< SBest << endl;
  cout << "# Optimal Indices (transit #, index): " << endl;
  for (int i = 0; i < M; i++) {
    cout << i << "\t" << indices[i] << endl;
  }

  return 0;
}
