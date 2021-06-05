pkg load symbolic;
%gain stage


VIN=exp(j*pi/2);
R1=1000;
R2=1000;
ZR2=1/(2/1000);
C1=0.22E-6;
C2=0.22E-6;
R3=1/(3/10000);
R4=100000*3;


G=@(f) (1+R4/R3)*1/(1/R1+(j*2*pi*f*C1))*1/((2*pi*j*f*C2)+1/ZR2)/(1/(C1*2*pi*j*f)*ZR2);
fc=1000;
Gain=20*log10(abs(G(fc)))
Phase=arg(G(fc))

Iin=-VIN/(R1+1/(2*pi*fc*j*C1));
ZIN=-VIN/Iin
  V2=VIN/(1/(2*pi*fc*C1))*(1/(1/R1+2*pi*j*fc*C1));
  V5=(1+R4/R3)*V2;
Iout=1/(2*pi*fc*j*C2)+(1-V5)/ZR2;
Zout=1/Iout

  
%--------------------------------------------
G=zeros(1,70);
Gain=zeros(1,70);
Phase=zeros(1,70);
f=logspace(1,8,70);
for i=1:70
w=2*pi*f(i);
ZC1=1/(j*w*C1);
ZC2=1/(j*w*C2);
G(i)=(1+R4/R3)*1/(1/R1+1/ZC1)*1/(1/ZC2+1/ZR2)/(ZC1*ZR2);
Gain(i)=20*log10(abs(G(i)));
Phase(i)=angle(G(i));
endfor

MaxGain=max((Gain))
Cutoff=MaxGain-3
 
flowcutoff=1/(2*pi*C1*R1)
fhighcutoff=1/(2*pi*C2*ZR2)
centralf=sqrt(flowcutoff*fhighcutoff)

  
H=@(f) 20*log10(abs((1+R4/R3)*1/(1/R1+(j*2*pi*f*C1))*1/((2*pi*j*f*C2)+1/ZR2)/(1/(C1*2*pi*j*f)*ZR2)))-32.659;
S=[10, 1000];
flow=fzero(H,S)
M=[1000, 10000];
fhigh=fzero(H,M)
centralf2=sqrt(flow*fhigh)
Max=H(centralf)+32.659
  
fig=figure();
plot(log10(f),Gain);
xlabel ("Frequency[Hz]")
ylabel ("vo/vi[dB]")
print(fig, "gain.eps", "-depsc");



figg=figure();
plot(log10(f),Phase);
xlabel ("Frequency[Hz]")
ylabel ("Phase[rad]")
print(figg, "phase.eps", "-depsc");
