#ifndef _QATS_H
#define _QATS_H

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
// shConvol: convolves the data, y, comprised of N datum, with a box of unit
// depth (or height -- see PARITY in qats.cpp) and returns the result of size
// N-q+1 in d.
//  
// inputs:	y - data
//		N - size of y
//		q - size of box
// outputs:	d - convolved data of size N-q+1 
// pre:		d - is allocated with appropriate size
//
void shConvol(double *  y,int N,int q, double *  & d);

//
// qats: determines the maximum metric Sbest over the set of transits instants
// for minimum interval DeltaMin and maximum interval DeltaMax for transits of
// duration (box-width) q by searching the box-convolved data d (of size N). 
// d may be computed with shConvol prior to this call (or by another 
// appropriate filter).
//
// inputs:	d - convolved data
//		DeltaMin - minimum interval (in cadences)
//		DeltaMax - maximum interval (in cadences)
//		N - size of convolved data
//		q - duration (width) of box signal
// outputs:	Sbest - maximum of obj. function for optimal transit instants
//		Mbest - number of transits in d corresponding to Sbest
//
void qats(double *  d, int DeltaMin, int DeltaMax, int N, int q, double & Sbest, int & MBest);

//
// qats_indices: determines the maximum metric Sbest over the set of transit
// instants for minimum interval DeltaMin and maximum interval DeltaMax 
// for M transits of duration (box-width) q by searching the box-convolved 
// data d (of size N). d may be computed with shConvol prior to this call 
// (or by another appropriate filter). The indices of the best transit instants
// are also returned.
//
// inputs:	d - convolved data
//		M - number of transits
//		DeltaMin - minimum interval (in cadences)
//		DeltaMax - maximum interval (in cadences)
//		N - size of convolved data
//		q - duration (width) of box signal
// outputs:	Sbest - maximum of obj. function for optimal transit instants
//		indices - array of M optimal transit instants (cadence #)
//
void qats_indices(double * d, int M, int DeltaMin, int DeltaMax, int N, int q, double & Sbest, int * & indices);

#endif

