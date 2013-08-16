function new_baseline_average
tic;                        %start timing
clear all;                  %clear memory
close all;                  %close existing plots

global SimValues sol NumVars;
global para;                %declare the temporary variable para to read the parameters from the input file

fid = fopen('para_Fe_dynamics.txt','r'); %reading the parameters from para.in
C = textscan(fid,'%s%f');   %names of parameters, values of parameters
para = C{2};                %we only need the values of parameters
fclose(fid);                %close file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BLOCK ONE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%***** USER DEFINED *****%
x = linspace(0,15,511);         % definition of the spatial domain, x=[0 15] cm with resolution of 511
t = linspace(0,99,100);        % definition of the time domain, t=[0 400] years with resolution of 401
% defining the species
VarNames = {'O2(aq)',...        %u1
            'Fe(OH)3(s)',...    %u2
            'FeOOH(s)', ...     %u3
            'SO4(2-)(aq)', ...  %u4
            'Fe(2+)(aq)', ...   %u5
            'S(-II)(aq)',...    %u6
            'FeS(s)',...        %u7
            'S(0)(aq)',...      %u8
            'OM1(s)',...        %u9
            'OM2(s)',...        %u10
			'PO4(aq)',...		%u11
			'PO4adsa(s)',...	%u12
			'PO4adsb(s)',...	%u13
			'Ca2(aq)',...		%u14
			'CaPO4(s)',...		%u15
			'NO3(aq)',...		%u16
			'NH4(aq)'};%u17
NumVars = int16(length(VarNames)) ;

SimValues = cell(NumVars,1);                      %creates an NumVars-by-1 cell array of empty matrices
                                                    %to store the concentration
rt_time = clock;                                    %display starting time
fprintf('The starting time is %.0f/%.0f/%.0f %.0f:%.0f:%.0f\n', rt_time(1),...
    rt_time(2), rt_time(3), rt_time(4), rt_time(5), rt_time(6));

disp('spin up ');
sol = pdepe(0,@pdex14pde,@pdex1ic,@pdex1bc,x,t);    %calls the solver and stores the result in 'sol'
toc;                                                %stop timing

for j = 1:NumVars
    SimValues{j,1} = sol(:,:,j);                     %store the solution in the cell,
end                                                 %cleaner as compared to the sol matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BLOCK TWO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,f,s] = pdex14pde(x,~,u,DuDx)
global BC0_O2 BC0_FeOH3 BC0_FeOOH BC0_SO4 BC0_Fe BC0_H2S BC0_FeS BC0_S0...
       BC0_OM1 BC0_OM2 BC0_PO4 BC0_PO4adsa...
       BC0_PO4adsb BC0_Ca2 BC0_CaPO4 BC0_NO3 BC0_NH4 w F NumVars;
global para;
%***** USER DEFINED *****%
%boundary conditions at sediment-water interface---------------------------
BC0_O2 = (0.05+0.275)/2;
BC0_FeOH3 = 20;
BC0_FeOOH = 10;
BC0_SO4 = 0.625;
BC0_Fe = 1.9e-3; 
BC0_H2S = 0; 
BC0_FeS = 0;
BC0_S0 = 0;
BC0_OM1 = 10;
BC0_OM2 = BC0_OM1*0.5;
BC0_PO4 = (3.8e-3+4.2e-4)/2;
fPFe = 0.25; %0.08
BC0_PO4adsa = BC0_FeOH3*fPFe;
BC0_PO4adsb = BC0_FeOOH*fPFe;
BC0_Ca2 = 1.125;
BC0_CaPO4 =1;
BC0_NO3 = 0.003;
BC0_NH4 = 0.004;
%bioturbatiuon coef--------------------------------------------------------
D_bio = 0.0694;
%molecular diffusion coefs-------------------------------------------------
D_O2 = 375; 
D_SO4 = 175;
D_Fe = 118;
D_H2S = 284;
D_S0 = 100;
D_PO4 = 232;
D_Ca2 = 212;
D_NO3 = 508;
D_NH4 = 530;
%xyz ratio-----------------------------------------------------------------
Cx1 = para(1);              %each value has been replaced with the    
Ny1 = para(2);
Pz1 = para(3);              %corresponding entry in 'para'
Cx2 = para(4);
Ny2 = para(5);
Pz2 = para(6);
%half saturation coefs-----------------------------------------------------
KO2 = para(7);
KNO3 = para(8);
KFeOH3 = para(9);               
KFeOOH = para(10);
KSO4 = para(11);
%inhibition coefs----------------------------------------------------------
kinO2 = para(12);
kinNO3 = para(13);
kinFeOH3 = para(14);
kinFeOOH = para(15);
%saturation constants------------------------------------------------------
KFeS = para(16);
%Secondary reaction constants----------------------------------------------
kOM1 = 1;                           %degredation rate constant for OM1, Rc1
kOM2 = 0.1;                        %degredation rate constant for OM2, Rc2
ktsox = para(17);                   %S(II) oxidation by O2, R5
ktsfe = para(18);                   %Fe(OH)3 reduction by S2-, R6
kfeox = para(19);                   %Fe(II) oxidation by O2, R7
kfedis = para(20);                  %dissolution rate of FeS, R8
kfepre = para(21);                  %precipitation rate of FeS, R_8
ksorg = para(22);                   %Sulfidization of OM, R9
k_psorb_a = para(23);               %PO4 sorption onto Fe(OH)3, R22
k_psorb_b = para(24);               %PO4 sorption onto FeOOH, R23
kapa = para(25);                    %Apatite precipitation, R25
k_apa = para(26);
knh4ox = para(27);
Ksorp_Fe=para(28);
kFeOOH=para(29);
kMag=para(30);
KMag=para(31);
SA=para(32);
SD=para(33);
MW=para(34);
SAR=para(35);


%other parameters----------------------------------------------------------
rho=1.5;       % solid phase density in gr/cm^3
fi=0.95;     %porosity
F= rho*(1-fi)/fi;    % conversion factor
w = 0.17; %0.23 or SAR/F;     %time-dependant burial rate
alfa0 = 10;                           %bioirrigation constant at sediment-water interface
alfax = alfa0*exp(-0.25*x);             %depth-dependant bioirrigation
h_plus = 1e-5;% %umol cm-3                 %[H+] concentration, equals 10^(-pH)

%indices-------------------------------------------------------------------
SFeOH3=u(2)*SA*SD*MW;%total surface available for sorption of Fe(II) (Umol/gr)
Fe_sorb=Ksorp_Fe*SFeOH3*u(5)/(h_plus+Ksorp_Fe*u(5));
p=1;%p=(SFeOH3-Fe_sorb)/SFeOH3;%passivation of ferrihydrite due to sorption of Fe(II)
Rret=.06*Ksorp_Fe*SFeOH3*h_plus/(h_plus+Ksorp_Fe*u(5))^2+1;%Retardation Factor
fO2 = u(1)/(KO2+u(1));
fNO3 = u(16)/(KNO3+u(16))*kinO2/(kinO2+u(1));
fFeOH3 = u(2)/(KFeOH3+u(2))*kinO2/(kinO2+u(1))*kinNO3/(kinNO3+u(16))*p;
fFeOOH = u(3)/(KFeOOH+u(3))*kinO2/(kinO2+u(1))*kinNO3/(kinNO3+u(16))*kinFeOH3/(kinFeOH3+u(2));
fSO4 = u(4)/(KSO4+u(4))*kinO2/(kinO2+u(1))*kinNO3/(kinNO3+u(16))*kinFeOH3/(kinFeOH3+u(2))*kinFeOOH/(kinFeOOH+u(3));
Sat_FeS = u(5)*u(6)/(KFeS*h_plus^2);    %saturation index of FeS
P_index1 = Pz1/Cx1;
P_index2 = Pz2/Cx2;
N_index1 = Ny1/Cx1;
N_index2 = Ny2/Cx2;


%reaction rates------------------------------------------------------------
Rc1 = kOM1*u(9);            %first pool rate
R1a = Rc1*fO2;           %OM oxidation
R2a = Rc1*fNO3;             %OM oxidation by NO3
R3a = Rc1*fFeOH3;           %OM oxidation by Fe(OH)3
R4a = Rc1*fFeOOH;           %OM oxidation by FeOOH
R5a = Rc1*fSO4;             %OM oxidation by SO4

Rc2 = kOM2*u(10);           %second pool rate
R1b = Rc2*fO2;           %OM oxidation
R2b = Rc2*fNO3;             %OM oxidation by NO3
R3b = Rc2*fFeOH3;           %OM oxidation by Fe(OH)3
R4b = Rc2*fFeOOH;           %OM oxidation by FeOOH
R5b = Rc2*fSO4;             %OM oxidation by SO4

R1 = R1a+R1b;
R2 = R2a+R2b;
R3 = R3a+R3b;
R4 = R4a+R4b;
R5 = R5a+R5b;

Ra = R1a+R2a+R3a+R4a+R5a;                   %degradation of OM1
Rb = R1b+R2b+R3b+R4b+R5b;                   %degradation of OM2

R6 = ktsox*u(1)*u(6);                   %S(II) oxidation by O2
R7 = ktsfe*u(2)*u(6);                   %Fe(OH)3 reduction by S2-
R8 = kfeox*u(1)*u(5);                   %Fe(II) oxidation by O2

if (Sat_FeS>=1)
    R9 = 0;
    R_9 = kfepre*(Sat_FeS-1);           %precipitation rate of FeS
else
    R9 = kfedis*u(7)*(1-Sat_FeS);       %dissolution rate of FeS
    R_9 = 0;
end

R10 = ksorg*u(6)*(u(9)+u(10));           %Sulfidization of OM
R11 = k_psorb_a*u(2)*u(11);             %PO4 sorption onto Fe(OH)3
R12 = k_psorb_b*u(3)*u(11);             %PO4 sorption onto FeOOH
R13 = fPFe*(4*R3+2*R7);                 %PO4 release from pool a
R14 = fPFe*4*R4;                        %PO4 release from pool b

if u(11)>kapa
    R15 = k_apa*(u(11)-kapa);           %Apatite precipitation
else
    R15 = 0;
end

R16 = knh4ox*u(1)*u(17);
% R17=kFeOOH*u(5)*(u(2)>0);%Goethite/lepidocrocite formation
% R18a=Rc1*u(18)/(u(18)+KMag)*kinO2/(kinO2+u(1));%Magnetite respiration 1st pool of OM
% R18b=Rc2*u(18)/(u(18)+KMag)*kinO2/(kinO2+u(1));%Magnetite respiration 2nd pool of OM
% R18=Rc1+Rc2;%Magnetite respiration
% R19=kMag*u(5)*(u(2)>0);%magnetite formation
R17=0;R18=0;R19=0;


%constant porosity---------------------------------------------------------
c = ones(NumVars,1);%c

%Transport: constant porosity----------------------------------------------
f = [(D_bio+D_O2)*DuDx(1)-w*u(1);...
    D_bio*DuDx(2)-w*u(2);...
    D_bio*DuDx(3)-w*u(3);...
    (D_bio+D_SO4)*DuDx(4)-w*u(4);...
    (D_bio+D_Fe)*DuDx(5)-w*u(5);...
    (D_bio+D_H2S)*DuDx(6)-w*u(6);...
    D_bio*DuDx(7)-w*u(7);...
    (D_bio+D_S0)*DuDx(8)-w*u(8);...
    D_bio*DuDx(9)-w*u(9);...
    D_bio*DuDx(10)-w*u(10);...
	(D_bio+D_PO4)*DuDx(11)-w*u(11);...
	D_bio*DuDx(12)-w*u(12);...
	D_bio*DuDx(13)-w*u(13);...
	(D_bio+D_Ca2)*DuDx(14)-w*u(14);...
	D_bio*DuDx(15)-w*u(15);...
	(D_bio+D_NO3)*DuDx(16)-w*u(16);...
	(D_bio+D_NH4)*DuDx(17)-w*u(17)];%f
 
%Reaction: constant porosity----------------------------------------------- 
s = [(BC0_O2-u(1))*alfax + (-2*R6-0.25*R8-2*R16) - R1*F;...               %the first part are all in unit of umol/cm3/yr
    R8/F + (-4*R3-2*R7)-(R17+2/3*R18)/F;...                                         %while the second part are all in unit of umol/g/yr
    -4*R4+R17/F;...                                                       %this is the basis of using the conversion factor F
    (BC0_SO4-u(4))*alfax + R6 - 0.5*R5*F;...                    
    (BC0_Fe-u(5))*alfax - R8 + (4*R3+4*R4+2*R7+R9-R_9+6*R18)*F-1/3*R19;...  
    (BC0_H2S-u(6))*alfax - R6 + (0.5*R5-R7+R9-R_9-R10)*F;... 
    -R9+R_9;... 
    (BC0_S0-u(8))*alfax + R7*F;...
    -Ra;...
    -Rb;...
	(BC0_PO4-u(11))*alfax - 2*R15 + (-R11-R12+P_index1*Ra+P_index2*Rb+R13+R14)*F;...
	+R11-R13;...
	+R12-R14;...
	(BC0_Ca2-u(14))*alfax - 3*R15;...
	+R15/F;...
	(BC0_NO3-u(16))*alfax + R16 - 0.8*R2*F;...
	(BC0_NH4-u(17))*alfax - R16 + (N_index1*Ra+N_index2*Rb)*F];%s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BLOCK THREE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u0 = pdex1ic(~)                            %initial condition
global NumVars;
u0 = zeros(NumVars,1);%u0

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BLOCK FOUR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pl,ql,pr,qr] = pdex1bc(~,ul,~,ur,~)       %boundary condition
global BC0_O2 BC0_FeOH3 BC0_FeOOH BC0_SO4 BC0_Fe BC0_H2S BC0_FeS BC0_S0...
       BC0_OM1 BC0_OM2 BC0_PO4 BC0_PO4adsa...
       BC0_PO4adsb BC0_Ca2 BC0_CaPO4 BC0_NO3 BC0_NH4 w F NumVars;

%Upper bounday: constant porosity------------------------------------------
pl = [ul(1)-BC0_O2;...
      BC0_FeOH3/F;...
      BC0_FeOOH/F;...
      ul(4)-BC0_SO4;...
      ul(5)-BC0_Fe;...
      ul(6)-BC0_H2S;...
      BC0_FeS/F;...
      ul(8)-BC0_S0;...
      BC0_OM1/F;...
      BC0_OM2/F;...
	  ul(11)-BC0_PO4;...
	  BC0_PO4adsa/F;...
	  BC0_PO4adsb/F;...
	  ul(14)-BC0_Ca2;...
	  BC0_CaPO4/F;...
	  ul(16)-BC0_NO3;...
	  ul(17)-BC0_NH4];%pl

%1 for solid, 0 for solute for both constant and depth-dependant porosity
ql = [0;1;1;0;0;0;1;0;1;1;0;1;1;0;1;0;0];%ql

%Lower bounday: constant porosity------------------------------------------
pr = [w*ur(1);...
      w*ur(2);...
      w*ur(3);...
      w*ur(4);...
      w*ur(5);...
      w*ur(6);...
      w*ur(7);...
      w*ur(8);...
      w*ur(9);...
      w*ur(10);...
      w*ur(11);...
      w*ur(12);...
      w*ur(13);...
      w*ur(14);...
	  w*ur(15);...
	  w*ur(16);...
	  w*ur(17)];%pr
  
qr = ones(NumVars,1);%qr
%----------------------------------------------------------------