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
plot(t*1000, vS, "r-;vS;")
hold on;
plot(t*1000, vOhr, "c-;Rectified;")
xlabel ("t[ms]")
ylabel ("v[V]")
print ("vo.eps", "-depsc");

%envelope detector
tOFF=atan(1/(w*R*C))/w

vOnexp = (A*cos(w*tOFF)- 2*vON)*exp(-(t-tOFF)/(R*C));

figure
plot(t*1000, vOhr, "c-; Rectified;")
hold on;
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

plot(t*1000, vO, "b-;Envelope;")
hold off;
axis ([0 10])
xlabel ("t[ms]")
ylabel ("v[V]")
legend ('location', 'west');
print ("venvlope.eps", "-depsc");


%voltage regulator
Is=1e-12;
vT=25e-3;
vON=0.65;
Vd = 12/18;
vlim =18*vON;
rd=vT/(Is*exp(Vd/vT))
vo= zeros(1,length(t));


figure;
%limit
%for i=1:length(t)
 % if vO(i) > vlim
  %  vO(i) = vlim;
 % elseif vO(i) < -vlim
  %  vO(i) = -vlim;
 % endif
%endfor

for i=1:length(t)
vo(i) =18*rd/(R+18*rd) *vO(i);
endfor

  VO=12;
  vOr= zeros(1,length(t));
for i=1:length(t)
vOr(i)=VO+vo(i);
endfor

plot(t*1000, vOr, 'y');
axis ([0 10])
xlabel ("t[ms]")
ylabel ("V[V]")
print ("vregulator.eps", "-depsc");

plot(t*1000, vo, 'g');
axis ([0 10])
xlabel ("t[ms]")
ylabel ("v[V]")
print ("vfinal.eps", "-depsc");


%DC output level
t2 = linspace(0, 0.2, 50);
aux= zeros(1,length(t2));
for i=1:length(t2)
  aux(i) = vOr(i);
endfor

vav = mean(aux, "a")


%Ripple
vmax = max(aux);
vmin = min(aux);

rip = vmax - vmin
