pkg load symbolic;
%gain stage

VT=25e-3
BFN=178.7
VAFN=69.7
RE1=100
RC1=1000
RB1=80000
RB2=20000
VBEON=0.7
VCC=12
RS=100

RB=1/(1/RB1+1/RB2)
VEQ=RB2/(RB1+RB2)*VCC
IB1=(VEQ-VBEON)/(RB+(1+BFN)*RE1)
IC1=BFN*IB1
IE1=(1+BFN)*IB1
VE1=RE1*IE1
VO1=VCC-RC1*IC1
VCE=VO1-VE1
  
%voltage nodes gain stage
V1gs=VEQ
V2gs=VEQ-RB*IB1
V3gs=VO1
V4gs=VCC
V5gs=RE1*IE1


gm1=IC1/VT
rpi1=BFN/gm1
ro1=VAFN/IC1

RSB=RB*RS/(RB+RS)

AV1 = RSB/RS * RC1*(RE1-gm1*rpi1*ro1)/((ro1+RC1+RE1)*(RSB+rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)
AVI_DB = 20*log10(abs(AV1))
AV1simple = RB/(RB+RS) * gm1*RC1/(1+gm1*RE1)
AVIsimple_DB = 20*log10(abs(AV1simple))

RE1=0
AV1 = RSB/RS * RC1*(RE1-gm1*rpi1*ro1)/((ro1+RC1+RE1)*(RSB+rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)
AVI_DB = 20*log10(abs(AV1))
AV1simple =  - RSB/RS * gm1*RC1/(1+gm1*RE1)
AVIsimple_DB = 20*log10(abs(AV1simple))

RE1=100
ZI1 = 1/(1/RB+1/(((ro1+RC1+RE1)*(rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1)))
ZX = ro1*((RSB+rpi1)*RE1/(RSB+rpi1+RE1))/(1/(1/ro1+1/(rpi1+RSB)+1/RE1+gm1*rpi1/(rpi1+RSB)))
ZX = ro1*(   1/RE1+1/(rpi1+RSB)+1/ro1+gm1*rpi1/(rpi1+RSB)  )/(   1/RE1+1/(rpi1+RSB) ) 
ZO1 = 1/(1/ZX+1/RC1)

RE1=0
ZI1 = 1/(1/RB+1/(((ro1+RC1+RE1)*(rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1)))
ZO1 = 1/(1/ro1+1/RC1)

%ouput stage
BFP = 227.3
VAFP = 37.2
RE2 = 100
VEBON = 0.7
VI2 = VO1
IE2 = (VCC-VEBON-VI2)/RE2
IC2 = BFP/(BFP+1)*IE2
VO2 = VCC - RE2*IE2

%voltage nodes and currents output stage
V1os=VI2
V2os=VO2
V3os=VCC
IBos=IE2-IC2
ICos=IC2  
IEos=IE2

  
gm2 = IC2/VT
go2 = IC2/VAFP
gpi2 = gm2/BFP
ge2 = 1/RE2
RL=8;

AV2 = gm2/(gm2+gpi2+go2+ge2)
ZI2 = (gm2+gpi2+go2+ge2)/gpi2/(gpi2+go2+ge2)
ZO2 = 1/(gm2+gpi2+go2+ge2)

RE1=100
%frequency response
gain=zeros(1,70);
f=logspace(1,8,70);
for i=1:70
VI=0.010*exp(j*pi/2);
w=2*pi*f(i);
C1=1.7e-3;
CE=1.7e-3;
ZS=RS+1./(j*w*C1);
ZE=1/(1/RE1+j*w*CE);
C3=1.7e-3;
Z2=RL+1./(j*w*C3);
Z3=1./(go2+1/Z2);
ZT=1./(ge2+1/Z3);
N=[1,0,0,0,0;
    -1/ZS,1/ZS+1/RB+1/rpi1,-1/rpi1,0,0;
    0,-gm1-1/rpi1,gm1+1/ZE+1/rpi1+1/ro1,-1/ro1,0;
    0,gm1,-gm1-1/ro1,1/ro1+1/RC1+gpi2,-gpi2;
    0,0,0,-gpi2-gm2,1/ZT+gpi2+gm2
   ];
b=[VI;0;0;0;0];
V=N\b;
Vout=V(5)*j*w*C3/(j*w*C3+1/RL);
Voutr=abs(Vout);
VIr=abs(VI);
gain(i)=20*log10(abs(Vout/VI));
endfor


fig=figure();
plot(log10(f),gain);
xlabel ("Frequency[Hz]")
ylabel ("vo/vi[dB]")
print(fig, "gain.eps", "-depsc");
%total
gB = 1/(1/gpi2+ZO1)
AV = (gB+gm2/gpi2*gB)/(gB+ge2+go2+gm2/gpi2*gB)*AV1
AV_DB = 20*log10(abs(AV))
ZI=ZI1
ZO=1/(go2+gm2/gpi2*gB+ge2+gB)
