pkg load symbolic;
%gain stage


VIN=exp(j*pi/2);
R1=1000;
ZR1=1/(2/R1);
R2=1000;
C1=3*1E-6+2*0.22E-6;
C2=0.22E-6;
R3=1/(3/10000);
R4=100000*3;


f=1000;
w=2*pi*f;
C3=1E-9;
ZC3=1/(j*w*C3);
ZC1=1/(j*w*(3*C1));
ZR1=(1/(3/R1));
Z1=ZR1+ZC1;
Z2=1/(1/R2+1/ZC3);
N=[1,0,0;
    -1/Z1,1/Z2+1/Z1, -1/Z2;
    0,Z2/Z1,-1
   ];
b=[VIN;0;0];
V=N\b;
Gain=20*log10(abs(V(3)/V(1)))
  I1=(V(1)-V(2))/ZC1;
ZIN=(V(1)/I1)

  
%--------------------------------------------
G=zeros(1,70);
f=logspace(1,8,70);
for i=1:70
w=2*pi*f(i);
ZC1=1/(j*w*C1);
ZC2=1/(j*w*C2);
N=[1,0,0,0,0;
    -1/ZC1,1/ZR1+1/ZC1, 0, 0, 0;
    0,0, 1/R3+1/R4, -1/R4,0;
    0, -1-R4/R3, 1+R4/R3,-1,0;
    0,0,0,-1/R2, 1/R2+1/ZC2
   ];
b=[VIN;0;0;0;0];
V=N\b
  G(i)=20*log10((V(5)/V(1)));
endfor


fig=figure();
plot(log10(f),G);
xlabel ("Frequency[Hz]")
ylabel ("vo/vi[dB]")
print(fig, "gain.eps", "-depsc");
