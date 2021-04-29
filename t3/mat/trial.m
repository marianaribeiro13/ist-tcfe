close all 
clear all
pkg load symbolic
vT=25e-3;
f=1e3;
w=2*pi*f;
RON=80
R=5e3
C=5e-6
vON=0.65;

t=linspace(0, 5e-3, 1000);
A = 5;
vS=A*cos(w*t);
vO= zeros(1,length(t));
for i=1:length(t)
  if(vS(i) - vON >=0)
    vO(i) = vS(i) - vON;
  else
    vO(i) = 0;
  endif 
endfor

figure
plot(t*1000, vO)
title("Output voltage v_o(t)")
xlabel ("t[ms]")
ylabel ("v_o[V]")
print ("vo.eps", "-depsc");

%envelope detector
A=5
t=linspace(0, 5e-3, 1000);
f=1000;
w=2*pi*f;
vS = A * cos(w*t);
vOhr = zeros(1, length(t));
vO = zeros(1, length(t));
a= sym (A);
c= sym (C);
r= sym (R);
W= sym (w);
syms to;
n=C*A*R*w*sin(w*to)==A*cos(w*to);
tof=solve(n);
  for i=1:2
	  toff=tof(i);
  if toff>=0
    tOFF=double(toff)
endif
endfor
	 
vOnexp = A*cos(w*tOFF)*exp(-(t-tOFF)/R/C);

figure
for i=1:length(t)
if (vS(i)-vON >= 0)
    vOhr(i) = vS(i) - vON;
  else
    vOhr(i) = 0;
  endif
endfor

  for i=1:length(t)
  if(vS(i) - vON >=0)
    vO(i) = vS(i) - vON;
  else
    vO(i) = 0;
  endif 
endfor
  
plot(t*1000, vOhr)
hold
T=1/f
for i=1:length(t)
if 0<t(i) < tOFF
	vO(i) = vOhr(i);
elseif vOnexp(i)-vON> vOhr(i)
   vO(i) = vOnexp(i)-vON;
  else 
    vO(i) = vOhr(i);
  endif
endfor

plot(t*1000, vO)
title("Output voltage v_o(t)")
xlabel ("t[ms]")
legend("rectified","envelope")
print ("venvlope.eps", "-depsc");








%lpf response
t=linspace(0, 50e-3, 1000);
Is=1e-9;
VON=0.65
vlim =3*VON
  rd=2*vT/(Is*exp(3*vON/2/vT));
f=400;
w=2*pi*f;
vo =3*rd/(R+3*rd) *A* cos(w*t);
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
  VO=3*VON;
vO=VO+vo;
plot(t*1000, vo);
print ("vrec.eps", "-depsc");
