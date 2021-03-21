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

Va = sym ('5.23334661904');
Id = sym ('1.0462067838');
Kb = sym ('7.03552682478');
Kc = sym ('8.12635172903');

filename='octaveresults.tex';
fp=fopen('octaveresults.tex', 'w');
fprintf(fp, "\n\n \\subsection{Mesh Analysis}\n");
syms I1 I2 I3;
M1=(R4+R3+R1)*I1-R3*I2-R4*I3==-Va;
M2=-I1*Kb*R3+I2*(Kb*R3-1)==0;
M3=-I1*R4+I3*(R4+R6+R7-Kc)==0;
ms= solve(M1,M2,M3,[I1 I2 I3]);
fprintf(fp, "\nSolution (mA):\n");
I1 = double(ms.I1)
I2 = double(ms.I2)
I3 = double(ms.I3)
fprintf(fp, "$ \\left(\\begin{array}{c} I_1 \\\\ I_2 \\\\ I_3  \\end{array}\\right) = \\left(\\begin{array}{c} %f \\\\ %f \\\\ %f \\end{array}\\right) $", I1, I2, I3);
fprintf("\nCurrents in resistors:\n");
IM1=I1
IM2=I2
IM3=I2-I1
IM4=I1-I3
IM5=I2-double(Id)
IM6=I3
IM7=I3
printf("\nVoltages in resistors:\n");
VM1=IM1*double(R1)
VM2=IM2*double(R2)
VM3=IM3*double(R3)
VM4=IM4*double(R4)
VM5=IM5*double(R5)
VM6=IM6*double(R6)
VM7=IM7*double(R7)
printf("\nDependent Sources Currents and Voltages:\n");
Ib=I2
Ic=I3
Vc = double(Kc)*Ic
Vb = Ib/double(Kb)
  fprintf(fp, "\n \\begin{table}[H]\n \\footnotesize\n \\centering\n \\caption{Mesh Analysis results}\n \\label{tab:tables}\n \\begin{center}\n \\begin{tabular}{ccc} \n & Voltage (V) & Current (mA) \\\\ \n \\hline \n\n\n \\hline \n $R_1$ & %f & %f \\\\ \n \\hline \n $R_2$ & %f & %f \\\\ \n \\hline \n $R_3$ & %f & %f \\\\ \n \\hline \n $R_4$ & %f & %f \\\\ \n \\hline \n $R_5$ & %f & %f \\\\ \n \\hline \n $R_6$ & %f & %f \\\\ \n \\hline \n $R_7$ & %f & %f \\\\ \n \\hline \n $V_c$ & %f & %f \\\\ \n \\hline \n $I_b$ & %f & %f \\\\ \n \\hline \n \\end{tabular} \n \\end{center} \n \\end{table}", VM1, IM1, VM2, IM2, VM3, IM3, VM4, IM4, VM5, IM5, VM6, IM6, VM7, IM7, Vc, Ic, Vb, Ib)
fprintf(fp, "\n\n \\subsection{Nodal Analysis}\n");
syms V1 V2 V3 V4 V5 V6 V7;
N1= V1==Va;
N2= (V1-V2)/R1 + (V3-V2)/R2 + (V4-V2)/R3==0;
N3= V2 - V4 + (V2-V3)/(R2*Kb)==0;
N4=Id+(V7-V6)/R7 + (V4-V5)/R5 + (V4-V2)/R3 + V4/R4==0;
N5=(V3-V2)/R2 + (V5-V4)/R5==Id;
N6=(V7-V6)/R7 -V6/R6==0;
N7=V4-V7+(Kc*V6)/R6==0;
ns= solve(N1,N2,N3,N4,N5,N6,N7, [V1 V2 V3 V4 V5 V6 V7]);
fprintf(fp, "\nSolution (V):\n");

V1 = double(ns.V1)
V2 = double(ns.V2)
V3 = double(ns.V3)
V4 = double(ns.V4)
V5 = double(ns.V5)
V6 = double(ns.V6)
V7 = double(ns.V7)
  
  fprintf(fp, "$ \\left(\\begin{array}{c} V_1 \\\\ V_2 \\\\ V_3 \\\\ V_4 \\\\ V_5 \\\\ V_6 \\\\ V_7 \\end{array}\\right)= \\left(\\begin{array}{c} %f \\\\ %f \\\\ %f \\\\ %f \\\\ %f \\\\ %f \\\\ %f \\end{array}\\right) $", V1, V2, V3, V4, V5, V6, V7);
printf("\nCurrents in resistors:\n");
IN1=(V2-V1)/double(R1)
IN2=(V3-V2)/double(R2)
IN3=(V2-V4)/double(R3)
IN4=-V4/double(R4)
IN5=(V4-V5)/double(R5)
IN6=-V6/double(R6)
IN7=(V6-V7)/double(R7)
printf("\nVoltages in resistors:\n");
VN1=IN1*double(R1)
VN2=IN2*double(R2)
VN3=IN3*double(R3)
VN4=IN4*double(R4)
VN5=IN5*double(R5)
VN6=IN6*double(R6)
VN7=IN7*double(R7)
printf("\nDependent Sources Currents and Voltages:\n");
ICN=IN6
IBN=IN2
VBN=IBN/double(Kb)
VCN=ICN*double(Kc)
fprintf(fp, "\n \\begin{table}[H]\n \\footnotesize\n \\centering\n \\caption{Nodal Analysis results}\n \\label{tab:tables}\n \\begin{center}\n \\begin{tabular}{ccc} \n & Voltage (V) & Current (mA)\\\\ \n \\hline \n\n\n \\hline \n $R_1$ & %f & %f \\\\ \n \\hline \n $R_2$ & %f & %f \\\\ \n \\hline \n $R_3$ & %f & %f \\\\ \n \\hline \n $R_4$ & %f & %f \\\\ \n \\hline \n $R_5$ & %f & %f \\\\ \n \\hline \n $R_6$ & %f & %f \\\\ \n \\hline \n $R_7$ & %f & %f \\\\ \n \\hline \n $V_c$ & %f & %f \\\\ \n \\hline \n $I_b$ & %f & %f \\\\ \n \\hline \n \\end{tabular} \n \\end{center} \n \\end{table}", VN1, IN1, VN2, IN2, VN3, IN3, VN4, IN4, VN5, IN5, VN6, IN6, VN7, IN7, VCN, ICN, VBN, IBN)
  fclose(fp);
