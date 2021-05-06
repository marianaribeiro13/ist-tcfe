close all 
clear all
vT=25e-3;
f=50;
w=2*pi*f;
R=3e3;
C=3e-6;
vON=0.65;

t=linspace(0, 0.2, 10000);
A = 20;
vS=A*cos(w*t);
vOhr= zeros(1,length(t));
for i=1:length(t)
  if(vS(i) >= 2*vON)
    vOhr(i) = vS(i) - 2*vON;
  elseif (vS(i)<= -2*vON)
    vOhr(i) = -vS(i) - 2*vON;
  else
    vOhr(i) = 0;
  endif 
endfor

figure
plot(t*1000, vOhr)
title("Output voltage v_o(t)")
xlabel ("t[ms]")
ylabel ("v_o[V]")
print ("vo.eps", "-depsc");

%envelope detector
tOFF=atan(1/(w*R*C))/w
	 
vOnexp = (A*cos(w*tOFF)- 2*vON)*exp(-(t-tOFF)/(R*C));

figure

plot(t*1000, vOhr)
hold
T=1/f
for i=1:length(t)
if t(i) < tOFF
	vO(i) = vOhr(i);
  elseif vOnexp(i) > vOhr(i)
   vO(i) = vOnexp(i);
  else 
    vO(i) = vOhr(i);
  endif
endfor

plot(t*1000, vO)
axis ([0 10])
title("Output voltage v_o(t)")
xlabel ("t[ms]")
legend("rectified","envelope")
print ("venvlope.eps", "-depsc");



%lpf response
Is=1e-12;
vT=25e-3;
vON=0.65;
Vd = 12/20;
vlim =20*vON;
rd=vT/(Is*exp(Vd/vT))
vo =20*rd/(R+20*rd) *A* cos(w*t);
figure;
title("Output voltage v_o(t)")
xlabel ("t[ms]")
ylabel ("v_{lpf}[V]")
%limit
for i=1:length(t)
  if vO(i) > vlim
    vO(i) = vlim;
  elseif vO(i) < -vlim
    vO(i) = -vlim;
  endif
endfor
  VO=3*vON;
vO=VO+vo;
plot(t*1000, vo);
print ("vrec.eps", "-depsc");
