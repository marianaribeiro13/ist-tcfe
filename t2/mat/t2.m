close all
clear all

%%SYMBOLIC COMPUTATIONS

pkg load symbolic

R1 = sym ('1.02604398505');
R2 = sym ('2.03674964113');
R3 = sym ('3.09851037606');
R4 = sym ('4.10284190715');
R5 = sym ('3.13763431318');
R6 = sym ('2.01979218864');
R7 = sym ('1.03837216927');

Vs = sym ('5.23334661904');
C = sym ('1.0462067838');
Kb = sym ('7.03552682478');
Kd = sym ('8.12635172903');

printf("\n\n \\subsection{Nodal Analysis}\n");
syms V1 V2 V3 V5 V6 V7 V8;
N1= V1==Vs;
N2= (V2-V1)/R1 + (V2-V3)/R2 + (V2-V5)/R3==0;
N3= V5/R4 + (V5-V6)/R5 + (V8-V7)/R7 + (V5-V2)/R3==0;
N5= V2 - V5 == (V3- V2)/(R2*Kb);
N6= (V3-V2)/R2==(V5-V6)/R5;
N7= V7==R6*((V8-V7)/R7);
N8= Kd*((V7-V8)/R7)==V5-V8;
ns= solve(N1,N2,N3,N5,N6,N7,N8, [V1 V2 V3 V5 V6 V7 V8]);

V1 = vpa(ns.V1)
V2 = vpa(ns.V2)
V3 = vpa(ns.V3)
V5 = vpa(ns.V5)
V6 = vpa(ns.V6)
V7 = vpa(ns.V7)
V8 = vpa(ns.V8)
  
syms I1 I2 I3;
M1=(R4+R3+R1)*I1-R3*I2-R4*I3==-Vs;
M2=-I1*Kb*R3+I2*(Kb*R3-1)==0;
M3=-I1*R4+I3*(R4+R6+R7-Kd)==0;
ms= solve(M1,M2,M3,[I1 I2 I3]);
I1 = double(ms.I1)
I2 = double(ms.I2)
I3 = double(ms.I3)

Vx=V6-V8
  
syms V1x V2x V3x V5x V6x V7x V8x;
N2x= (V2x)/R1 + (V2x-V3x)/R2 + (V2x-V5x)/R3==0;
N3x= V2x - V5x == (V3x- V2x)/(R2*Kb);
N5x= V5x/R4 + (V3x-V2x)/R2 + (V8x-V7x)/R7 + (V5x-V2x)/R3==0;
N6x= V6x-V8x==Vx;
N7x= V7x==R6*((V8x-V7x)/R7);
N8x= Kd*((V7x-V8x)/R7)==V5x-V8x;
ns= solve(N2x,N3x,N5x,N6x,N7x,N8x, [V2x V3x V5x V6x V7x V8x]);

V2x = vpa(ns.V2x)
V3x = vpa(ns.V3x)
V5x = vpa(ns.V5x)
V6x = vpa(ns.V6x)
V7x = vpa(ns.V7x)
V8x = vpa(ns.V8x)
