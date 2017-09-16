QATS C++ Implementation
C++ implementation by J. Carter (2012).

If you use this code, please cite Carter & Agol (2012)
 *AND* Kel'manov and Jeon 2004:

 \bibitem[Kel'Manov \& Jeon (2004)]{Kel'Manov522004}Kel'Manov, 
	A.~V.~\& B.~Jeon,  A Posteriori Joint Detection and Discrimination of Pulses in a Quasiperiodic Pulse Train.
	{\it IEEE Transactions on Signal Processing} {\bf 52}, 645--656 (2004).

INSTALL.

Type make in unpacked directory.

TEST OUT EXAMPLE Kepler-36 data.

kep36.dat is a single column of detrended, median subtracted, normalized by sigma Kepler data for Kepler-36 (as shown in paper).  Missing cadences have been filled with zeros. 

- Get Sbest for DeltaMin = 794, DeltaMax = 794, q = 14

./call_qats kep36.dat 794 794 14
# QATS; Maximum S_{eta} 
# Input sequence from kep36.dat
# DeltaMin = 794
# DeltaMax = 794
# q = 14
# Output (S_Best, M_Best)
82.6976 55

- Get Sbest for DeltaMin = 794, DeltaMax = 797, q = 14

./call_qats kep36.dat 794 797 14
# QATS; Maximum S_{eta} 
# Input sequence from kep36.dat
# DeltaMin = 794
# DeltaMax = 797
# q = 14
# Output (S_Best, M_Best)
127.442 55

- Get optimal starting indices for the 55 transits for circumstance above

./call_qats_indices kep36.dat 55 794 797 14
# QATS; Optimal Indices 
# Input sequence from kep36.dat
# DeltaMin = 794
# DeltaMax = 797
# q = 14
# M = 55
# S_Best = 127.442
# Optimal Indices (transit #, index): 
0	110
1	904
2	1699
3	2494
4	3290
5	4085
6	4880
7	5676
8	6470
9	7265
10	8059
11	8853
12	9647
13	10441
14	11235
15	12029
16	12823
17	13617
18	14411
19	15205
20	15999
21	16793
22	17587
23	18381
24	19175
25	19969
26	20763
27	21557
28	22351
29	23146
30	23941
31	24736
32	25532
33	26327
34	27122
35	27917
36	28713
37	29507
38	30302
39	31096
40	31890
41	32684
42	33478
43	34272
44	35066
45	35860
46	36654
47	37448
48	38242
49	39036
50	39830
51	40624
52	41418
53	42212
54	43006
