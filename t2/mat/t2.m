close all
clear all
format long

%%SYMBOLIC COMPUTATIONS
#include "oct-stream.h"
pkg load symbolic
pkg load control
digits(11)
output_precision(11)
fid=fopen('../data.txt', 'r');
m_p = textscan(fid,'%s %s %s %f','delimiter', ' ', 'HeaderLines', 8);
fclose(fid);
fid=fopen('../data.txt', 'r');
m_s = textscan(fid,' %s %s %f','delimiter', ' ', 'HeaderLines', 9)
A=cell2mat(m_p(1,4));
B=cell2mat(m_p(1,3));
C=cell2mat(m_s(1,3));
fclose(fid);
D=C(9);
E=C(10);
Vsd=C(7);
R1d = A(1);
R2d = C(1);
R3d = C(2);
R4d = C(3);
R5d = C(4);
R6d = C(5);
R7d = C(6);
Cd =  C(8);

R1 = sym (A(1));
R2 = sym (R2d);
R3 = sym (R3d);
R4 = sym (R4d);
R5 = sym (R5d);
R6 = sym (R6d);
R7 = sym (R7d);

Vs = sym (Vsd);
C = sym (Cd);
Kb = sym (D);
Kd = sym (E);

Kbd =  D;
Kdd =  E;

filename='octave.txt';
fp=fopen('octave.txt', 'w');
fprintf(fp, " V1 1 0 dc %.11f; \n R1 2 1 %.11fk; \n R2 3 2 %.11fk; \n R3 2 5 %.11fk; \n R4 5 0 %.11fk; \n R5 6 5 %.11fk; \n R6 7 0 %.11fk; \n R7 8 9 %.11fk; \n G1 6 3 2 5 %.11fm; \n v2 7 9 dc 0;\n H1 5 8 v2 %.11fk; \n", Vsd, R1d, R2d, R3d, R4d, R5d, R6d, R7d, Kbd, Kdd);
fclose(fp);


printf("\n\n \\subsection{Nodal Analysis}\n");
syms V1 V2 V3 V5 V6 V7 V8;
N1= V1==Vs;
N2= (V2-V1)/R1 + (V2-V3)/R2 + (V2-V5)/R3==0;
N3= V2 - V5 == (V3- V2)/(R2*Kb);
N5= V5/R4 + (V5-V6)/R5 + (V8-V7)/R7 + (V5-V2)/R3==0;
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

  Vm1=double(V1);
 Vm2=double(V2);
 Vm3=double(V3);
 Vm5=double(V5);
 Vm6=double(V6);
 Vm7=double(V7);
Vm8=double(V8);
Ib=(Vm3-Vm2)/double(R2);

Im1 = (Vm2 -Vm1)/double(R1);
Im3 = (Vm2 -Vm5)/double(R3);
Im4 = Vm5/double(R4);
Im5 = (Vm6 -Vm5)/double(R5);
Im6 = Vm7/double(R6);
Im7 = (Vm8 -Vm7)/double(R7);

filename='minor.tex';
minor=fopen('minor.tex', 'w');
 fprintf(minor, "$ \\left(\\begin{array}{c} V_1 \\\\ V_2 \\\\ V_3 \\\\ V_5 \\\\ V_6 \\\\ V_7 \\\\ V_8 \\end{array}\\right)= \\left(\\begin{array}{c} %f \\\\ %f \\\\ %f \\\\ %f \\\\ %f \\\\ %f \\\\ %f \\end{array}\\right) $", Vm1, Vm2, Vm3, Vm5, Vm6, Vm7, Vm8);

  
fprintf(minor, "\n \\begin{table}[H]\n \\footnotesize\n \\centering\n \\caption{Nodal Analysis results for t<0}\n \\label{tab:tables}\n \\begin{center}\n \\begin{tabular}{cccc}\n\\hline \n Node & Voltage (V) & R & Current(mA) \\\\ \n \\hline \n 1 & %f& $R_1$ & %f \\\\ \n \\hline \n 2& %f & $R_2$ & %f\\\\ \n \\hline \n 3 & %f& $R_3$ & %f \\\\ \n \\hline \n 4 & 0 &$R_4$ & %f \\\\ \n \\hline \n 5 &%f& $R_5$ & %f\\\\ \n \\hline \n 6 & %f & $R_6$ & %f\\\\ \n \\hline \n 7 & %f & $R_7$ & %f\\\\ \n \\hline \n 8 & %f & $I_b$ & %f \\\\ \n \\hline \\end{tabular} \n \\end{center} \n \\end{table}", Vm1, Im1, Vm2, Ib, Vm3, Im3, Im4, Vm5, Im5,Vm6, Im6, Vm7, Im7,Vm8, Ib);
fclose(minor);
  
Vx=V6-V8
vx = double (Vx);
filename='oc0.txt';
f0=fopen('oc0.txt', 'w');
fprintf(f0, " V3 6 8 dc %.11f; \n R1 1 2 %.11fk; \n R2 2 3 %.11fk; \n R3 2 5 %.11fk; \n R4 0 5 %.11fk; \n R5 6 5 %.11fk; \n R6 0 7 %.11fk; \n R7 9 8 %.11fk; \n G1 6 3 2 5 %.11fm; \n v2 7 9 dc 0; \n H1 5 8 v2 %.11fk; \n", vx, R1d, R2d, R3d, R4d, R5d, R6d, R7d, Kbd, Kdd);
fclose(f0);

syms V2x V3x V5x V6x V7x V8x;
N2x= (V2x)/R1 + (V2x - V3x)/R2 + (V2x - V5x)/R3==0;
N3x= V2x - V5x == (V3x - V2x)/(R2*Kb);
N5x= V5x/R4 + (V3x - V2x)/R2 + (V8x - V7x)/R7 + (V5x - V2x)/R3==0;
N6x= V6x-V8x==Vx;
N7x= V7x==R6*((V8x - V7x)/R7);
N8x= Kd*((V7x - V8x)/R7)==V5x - V8x;
ns= solve(N2x,N3x,N5x,N6x,N7x,N8x, [V2x V3x V5x V6x V7x V8x]);

V2x = vpa(ns.V2x)
V3x = vpa(ns.V3x)
V5x = vpa(ns.V5x)
V6x = vpa(ns.V6x)
V7x = vpa(ns.V7x)
V8x = vpa(ns.V8x)

Ix=(V2x-V3x)/R2 + (V5x-V6x)/R5
Req=-Vx/Ix
Tau=Req*C/1000


  Vx1=0;
 Vx2=double(V2x);
 Vx3=double(V3x);
 Vx5=double(V5x);
 Vx6=double(V6x);
 Vx7=double(V7x);
Vx8=double(V8x);
Ib=(Vx3-Vx2)/double(R2);
Ix1 = (Vx2 -Vx1)/double(R1);
Ix3 = (Vx2 -Vx5)/double(R3);
Ix4 = Vx5/double(R4);
Ix5 = (Vx6 -Vx5)/double(R5);
Ix6 = Vx7/double(R6);
Ix7 = (Vx8 -Vx7)/double(R7);

filename='equal.tex';
equal=fopen('equal.tex', 'w');
 fprintf(equal, "$ \\left(\\begin{array}{c} V_1 \\\\ V_2 \\\\ V_3 \\\\ V_5 \\\\ V_6 \\\\ V_7 \\\\ V_8 \\end{array}\\right)= \\left(\\begin{array}{c} %f \\\\ %f \\\\ %f \\\\ %f \\\\ %f \\\\ %f \\\\ %f \\end{array}\\right) $", Vx1, Vx2, Vx3, Vx5, Vx6, Vx7, Vx8);

   
fprintf(equal, "\n \\begin{table}[H]\n \\footnotesize\n \\centering\n \\caption{Nodal Analysis results}\n \\label{tab:tables}\n \\begin{center}\n \\begin{tabular}{cccc}\n\\hline \n Node & Voltage (V) & R & Current(mA) \\\\ \n \\hline \n 1 & %f& $R_1$ & %f \\\\ \n \\hline \n 2& %f & $R_2$ & %f\\\\ \n \\hline \n 3 & %f& $R_3$ & %f \\\\ \n \\hline \n 4 & 0 &$R_4$ & %f \\\\ \n \\hline \n5 &%f& $R_5$ & %f\\\\ \n \\hline \n 6 & %f & $R_6$ & %f\\\\ \n \\hline \n 7 & %f & $R_7$ & %f\\\\ \n \\hline \n 8 & %f & $I_b$ & %f \\\\ \n \\hline \n \\end{tabular} \n \\end{center} \n \\end{table}", Vx1, Ix1, Vx2, Ib, Vx3, Ix3, Ix4, Vx5, Ix5, Vx6, Ix6, Vx7, Ix7,Vx8, Ib);

fprintf(equal, "\n $Ix=%f mA$  $Req=%f \\Omega$ $\\tau=%f$", double(Ix), double(Req),double(Tau));
fclose(equal);


  
printf("\n\nNatural solution:\n");

A=double(Vx)
%time axis: 0 to 20ms with 2us steps
t=0:2e-6:20e-03; %s
TC=double(Tau)
v6=A*exp(-t/TC);

hf = figure ();
plot (t*1000, v6, "r");

xlabel ("t[ms]");
ylabel ("v6n(t) [V]");
print (hf, "natural.eps", "-depsc");
C1=double(C);
zc=1/(2*pi*1000*C*i)
syms V1p V2p V3p V5p V6p V7p V8p;
N1p= V1p==1;
N2p= (V2p-V1p)/R1 + (V2p-V3p)/R2 + (V2p-V5p)/R3==0;
N3p= V2p - V5p == (V3p- V2p)/(R2*Kb);
N5p= V5p/R4 + (V5p-V6p)/R5 + (V8p-V7p)/R7 + (V8p-V6p)/zc + (V5p-V2p)/R3==0;
N6p= (V8p-V6p)/zc + (V5p-V6p)/R5==(V3p-V2p)/R2;
N7p= V7p==R6*((V8p-V7p)/R7);
N8p= Kd*((V7p-V8p)/R7)==V5p-V8p;
ns= solve(N1p,N2p,N3p,N5p,N6p,N7p,N8p, [V1p V2p V3p V5p V6p V7p V8p]);
  
V1p = vpa(ns.V1p)
  V2p = vpa(ns.V2p)
  V3p = vpa(ns.V3p)
  V5p = vpa(ns.V5p)
  V6p = vpa(ns.V6p)
  V7p = vpa(ns.V7p)
  V8p = vpa(ns.V8p)


 Vp1=double(V1);
 Vp2=double(V2p);
 Vp3=double(V3p);
 Vp5=double(V5p);
 Vp6=double(V6p);
 Vp7=double(V7p);
Vp8=double(V8p);

filename='major.tex';
major=fopen('major.tex', 'w');
 fprintf(major, "$ \\left(\\begin{array}{c} V_1 \\\\ V_2 \\\\ V_3 \\\\ V_5 \\\\ V_6 \\\\ V_7 \\\\ V_8 \\end{array}\\right)= \\left(\\begin{array}{c} %f \\\\ %f \\\\ %f \\\\ %f \\\\ %f \\\\ %f \\\\ %f \\end{array}\\right) $", Vp1, Vp2, Vp3, Vp5, Vp6, Vp7, Vp8);


fprintf(major, "\n \\begin{table}[H]\n \\footnotesize\n \\centering\n \\caption{Nodal Analysis results for t>0}\n \\label{tab:tables}\n \\begin{center}\n \\begin{tabular}{ccc} \n & Voltage (V)\\\\ \n \\hline \n\n\n \\hline \n 1 & %f \\\\ \n \\hline \n 2 & %f \\\\ \n \\hline \n 3 & %f \\\\ \n \\hline \n 4 & 0 \\\\ \n \\hline \n 5 & %f \\\\ \n \\hline \n 6 & %f \\\\ \n \\hline \n 7 & %f \\\\ \n \\hline \n 8 &  %f \\\\ \n \\hline \n \\\\ \n \\hline \n \\end{tabular} \n \\end{center} \n \\end{table}", Vp1, Vp2,Vp3,Vp5,Vp6,Vp7,Vp8);
fclose(major);


V6xa = double (V6x);
V8xa = double (V8x);
filename='oc1.txt';
f1=fopen('oc1.txt', 'w');
fprintf(f1, "R1 2 1 %.11fk; \n R2 2 3 %.11fk \n R3 2 5 %.11fk; \n R4 5 0 %.11fk; \n R5 6 5 %.11fk; \n R6 7 0 %.11fk; \n R7 9 8 %.11fk; \n G1 6 3 2 5 %.11fm; \n v2 7 9 dc 0; \n H1 5 8 v2 %.11fk; \n C1 6 8 %.11fuF;\n .ic v(6) = %f v(8) = %f;\n", R1d, R2d, R3d, R4d, R5d, R6d, R7d, Kbd, Kdd, Cd, V6xa, V8xa);
fclose(f1);

V6ra = real(V6p);
V6ima = imag (V6p);

V6r = double(V6ra);
V6im = double (V6ima);

w=2*pi*1000;
v6f=V6r*sin(w*t)+V6im*cos(w*t);
hff = figure ();
plot (t*1000, v6f, "b");

xlabel ("t[ms]");
ylabel ("v6f(t) [V]");
print (hff, "forced.eps", "-depsc");

v6t=V6r*sin(w*t)+V6im*cos(w*t)+A*exp(-t/TC);
hfff = figure ();
hold on;
plot ([-5,0], [double(V6), double(V6)], "r");
plot ([0,0], [double(V6), double(Vx)], "r");
plot (t*1000, v6t, "r- ;v6(t);");
plot (t*1000, v6f, "b- ;v6f(t);");
plot (t*1000, v6, "g- ;v6n(t);");
plot ([-5,0], [double(Vs), double(Vs)], "m- ;vs(t);");
plot ([0,0], [double(Vs), 0], "m");
plot (t*1000, sin(2*pi*1000*t), "m", "linewidth", 0.5);
hold off;
xlabel ("t[ms]");
ylabel ("[V]");
axis ([-5 20 -1.5 10]) 
print (hfff, "total.eps", "-depsc");

R2f=double(R2);
R5f=double(R5);
V1f = double (V1p);
V3f=double(V3p);
V2f=double(V2p);
V5f=double(V5p);
V8f=double(V8p);
f=logspace(-1,6,200);

T = 1 ./(1+i*2*pi*f*double(Tau));
T6 = 1 ./(1+i*2*pi*f*double(Tau)) + V8f;

ht = figure ();
hold on;
plot (log10(f), arg(T)*180/pi, "g- ;Vc(f);");
plot (log10(f), arg(T6)*180/pi,"r- ;V6(f);");
plot ([-1,6], [arg(V1f)*180/pi,arg(V1f)*180/pi],"b- ;Vs(f);");
hold off;
xlabel ("Log_{10}(f)[Hz]");
ylabel ("Phase(degrees)");
print (ht, "Phase.eps", "-depsc");

hm = figure ();
hold on;
plot (log10(f), mag2db(abs(T)), "g- ;Vc(f);");
plot (log10(f), mag2db(abs(T6)),"r- ;V6(f);");
plot ([-1,6], [mag2db(abs(V1f)),mag2db(abs(V1f))],"b- ;Vs(f);");
hold off;
xlabel ("Log_{10}(f)[Hz]");
ylabel ("Magnitude(dB)");
legend ('location', 'west');
print (hm, "Magnitude.eps", "-depsc");

numer = [0, 1];
denom = [double(Tau), 1];
sys = tf (numer, denom);
figure

bode(sys, f*2*pi);
print ("RC_bode.png", "-dpng");
