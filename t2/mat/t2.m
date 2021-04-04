close all
clear all
format long
%%SYMBOLIC COMPUTATIONS
#include "oct-stream.h"
pkg load symbolic
pkg load control
fid=fopen('../data.txt', 'r');
m_p = textscan(fid,'%s %s %s %f','delimiter', ' ', 'HeaderLines', 8);
m_s = textscan(fid,' %s %s %f','delimiter', ' ')
A=cell2mat(m_p(1,4));
B=cell2mat(m_p(1,3));
C=cell2mat(m_s(1,3))
fclose(fid);
D=C(8)
E=C(9)
R1 = sym (A(1));
R2 = sym (B(2,1));
R3 = sym (C(1));
R4 = sym (C(2));
R5 = sym (C(3));
R6 = sym (C(4));
R7 = sym (C(5));

Vs = sym (C(6));
C = sym (C(7));
Kb = sym (D);
Kd = sym (E);
filename='octave.txt';
fp=fopen('octave.txt', 'w')
fprintf(fp, "
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
 
Vx=V6-V8
  
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

Ix=-V6x/R5
Req=-Vx/Ix
Tau=Req*C/1000

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
zc=1/(2*pi*100*C*i)
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

V6r=-0.57979607293241249067088011934093;
V6im=- 0.000081810508973260330098219669037616;

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
V3f=double(V3p);
V2f=double(V2p);
V5f=double(V5p);
V8f=double(V8p);
f=logspace(-1,6,35);
phase = 1:35;
magn=1:35;
for j=1:35,
	zcf=1/(2*pi*f(j)*C1*i);
%V6f=((R5f*(V3f-V2f)-V5f*R2f)/(2*pi*f(j)*C1*i) -R2f*R5f*V8f)/(-(R2f*(R5f+1/(2*pi*f(j)*C1*i))))
  V6f = - (R5f*zcf*(V3f-V2f) - V8f*R5f*R2f - zcf*V5f*R2f)/(R2f*(R5f+zcf));
aux = V6f-V8f;
aux2=arg(aux);
 phase(j)=aux2 * 180/pi;
aux3 = abs(aux);
magn(j)=mag2db(aux3);
printf("V6\n")
disp(V6f)
printf("\naux\n")
disp(aux)
end;

hp = figure ();
plot (log10(f), phase, "r");
xlabel ("f[Hz]");
ylabel ("phase(degrees)");
print (hp, "phasev6.eps", "-depsc");

hp = figure ();
plot (log10(f), magn, "r");
xlabel ("f[Hz]");
ylabel ("magnitude[V]");
print (hp, "magnitudev6.eps", "-depsc");

T = 1 ./(1+i*2*pi*f*double(Tau));
numer = [0, 1];
denom = [double(Tau), 1];
sys = tf (numer, denom);
figure

bode(sys, f*2*pi);
hold on;
plot (f, phase, "r");
print ("RC_bode.png", "-dpng");
