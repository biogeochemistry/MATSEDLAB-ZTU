function MATSEDLAB
tic;                        %start timing
clear all;                  %clear memory
close all;                  %close existing plots

global SimValues sol NumVars;
global para;                %declare the temporary variable para to read the parameters from the input file

fid = fopen('para.in','r'); %reading the parameters from para.in
C = textscan(fid,'%s%f');   %names of parameters, values of parameters
para = C{2};                %we only need the values of parameters
fclose(fid);                %close file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BLOCK ONE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%***** USER DEFINED *****%
x = linspace(0,15,511);         % definition of the spatial domain, x=[0 15] cm with resolution of 511
t = linspace(0,200,201);        % definition of the time domain, t=[0 400] years with resolution of 401
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
			'NO3(aq)',...		%u12
			'NH4(aq)',...       %u13
             'NO2(aq)',...      %u14 AML
             'N2O(aq)',...      %u15 AML
             'N2(aq)'};         %u16 AML
NumVars = int16(length(VarNames)) ;

SimValues = cell(NumVars, 201);                      %creates an NumVars-by-1 cell array of empty matrices
                                                    %to store the concentration
disp('solving PDE ');
sol = pdepe(0,@pdex14pde,@pdex1ic,@pdex1bc,x,t);    %calls the solver and stores the result in 'sol'
toc;                                                %stop timing

rt_time = clock;                                    %display ending time
fprintf('The ending time is %.0f/%.0f/%.0f %.0f:%.0f:%.0f\n', rt_time(1),...
    rt_time(2), rt_time(3), rt_time(4), rt_time(5), rt_time(6));

for j = 1:NumVars
    SimValues{j,1} = sol(:,:,j);                     %store the solution in the cell,
end                                                 %cleaner as compared to the sol matrix

save('Result.mat','SimValues','para');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BLOCK TWO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,f,s] = pdex14pde(x,~,u,DuDx)
global BC0_O2 BC0_FeOH3 BC0_FeOOH BC0_SO4 BC0_Fe BC0_H2S BC0_FeS BC0_S0...
       BC0_OM1 BC0_OM2 BC0_PO4...
       BC0_NO3 BC0_NH4 BC0_NO2 BC0_N2O BC0_N2 w F NumVars;
global para;
%***** USER DEFINED *****%
%boundary conditions at sediment-water interface---------------------------
BC0_O2 = 0.238;
BC0_FeOH3 = 14.7;%from Canavan et al AML
BC0_FeOOH = 0; %4.6;
BC0_SO4 = 0.638; %from Canavan et al AML
BC0_Fe = 0; 
BC0_H2S = 0; 
BC0_FeS = 0;
BC0_S0 = 0;
BC0_OM1 = 630; %from Canavan et al AML
BC0_OM2 = 315; %from Canavan et al AML
BC0_PO4 = 3e-4;
BC0_NO3 = 0.154; 
BC0_NH4 = 0;
BC0_NO2 = 0;
BC0_N2O = 0;
BC0_N2 = 0; %<VALUE TO BE MODIFIED>
%bioturbatiuon coef--------------------------------------------------------
D_bio = 5;%Canavan et al D_bio between 0-5, 5 in the top layers
%molecular diffusion coefs-------------------------------------------------
D_O2 = 375; 
D_SO4 = 175;
D_Fe = 118;
D_H2S = 284;
D_S0 = 100;
%D_AsO4 = 160;
D_PO4 = 232;
%D_Ca2 = 212;
D_NO3 = 539; %AML
D_NH4 = 555; %AML
D_NO2 = 520; %AML
D_N2O = 666; %AML
D_N2 = 550; %AML
%xyz ratio-----------------------------------------------------------------
Cx1 = para(1);              %each value has been replaced with the    
Ny1 = para(2);
Pz1 = para(3);              %corresponding entry in 'para'
Cx2 = para(4);
Ny2 = para(5);
Pz2 = para(6);
%half saturation coefs-----------------------------------------------------
KmO2 = 0.01;%Km O2 for aerobic degradation of OM
KmNO3 = 0.005;%Km for nitrate reduction of OM
KmFeOH3 = para(9);               
%KFeOOH = para(10);
KmSO4 = para(11);
KmNO2 = 0.01; %Km for nitrite reduction AML 
KmNO2DNRA = 0.01; %Km for nitrite reduction to ammonium AML 
KmN2O = 0.005; %Km for n2o reduction AML 
KmNH4ao = 0.005; %Km for NH4 ammonia oxidation AML 
KmNO2no = 0.015; %Km for NO2 nitrite oxidation AML 
KmNO2amx = 0.005; %Km for nitrite anammox AML 
KmNH4amx = 0.005; %Km for ammonium anammox AML 
KmO2ao = 0.005; %Km for O2 ammonia oxidation AML 
KmO2no = 0.001; %Km for O2 nitrite oxidation AML 

%inhibition coefs----------------------------------------------------------
kinO2 = para(12);
kinNO3 = 0.01;
kinFeOH3 = para(14);
kinNO2 = 0.01; % inhibition coefficient of NO2 for N2O reduction AML 

%Secondary reaction constants----------------------------------------------
kOM1 = 1;                           %degredation rate constant for OM1, Rc1 Canavan et al
kOM2 = 0.01;                        %degredation rate constant for OM2, Rc2 Canavan et al
kfeox = para(19);                   %Fe(II) oxidation by O2, R7
knh4ox = 20; %max ammonium oxidation rate AML 
kno2ox = 20; %max nitrite oxidation rate AML 
kamx = 20; %max anammox rate AML 

%other parameters----------------------------------------------------------
w = 5;     %time-dependant burial rate from Canavan et al AML
alfax = 10;             %fixed  bioirrigation
%h_plus = 3.4e-4;                        %[H+] concentration, equals 10^(-pH)
F = 0.26;                               %convertion factor=rhob*(1-fi)/fi; where fi=porosity and rhob=solid phase density
                                        %2.1*(1-0.89)/0.89


%indices-------------------------------------------------------------------
fO2 = u(1)/(KmO2+u(1));
fNO3 = (u(12)/(KmNO3+u(12)))*(kinO2/(kinO2+u(1)));%nitrate reduction to nitrite
fNO2 = (u(14)/(KmNO2+u(14)))*(kinO2/(kinO2+u(1))); %nitrite reduction to N2O AML
fNO2DNRA = (u(14)/(KmNO2DNRA+u(14)))*(kinO2/(kinO2+u(1))); %nitrite reduction to ammonium AML
fN2O = (u(15)/(KmN2O+u(15)))* (kinO2/(kinO2+u(1))); %nitrous oxide reduction to N2 AML
fFeOH3 = u(2)/(KmFeOH3+u(2))*kinO2/(kinO2+u(1))*kinNO3/(kinNO3+u(12));
fSO4 = u(4)/(KmSO4+u(4))*kinO2/(kinO2+u(1))*kinNO3/(kinNO3+u(12))*kinFeOH3/(kinFeOH3+u(2));
%Sat_FeS = 0;% u(5)*u(6)/(KFeS*h_plus^2);    %saturation index of FeS
P_index1 = Pz1/Cx1;
P_index2 = Pz2/Cx2;
N_index1 = Ny1/Cx1;
N_index2 = Ny2/Cx2;

%reaction rates------------------------------------------------------------
Rc1 = kOM1*u(9);            %first pool rate
R1a = Rc1*fO2*25;           %OM oxidation
R2a = Rc1*fNO3*25;          %OM oxidation by NO3
R3a = Rc1*fFeOH3;           %OM oxidation by Fe(OH)3
%R4a = 0; %Rc1*fFeOOH;      %OM oxidation by FeOOH
R5a = Rc1*fSO4;             %OM oxidation by SO4
R6a = Rc1*fNO2*25;          %OM oxidation by NO2 AML
R7a = Rc1*fN2O*25;          %OM oxidation by N2O AML
R9a = Rc1*fNO2DNRA*25;      %OM oxidation by NO2 to NH4 AML

Rc2 = kOM2*u(10);           %second pool rate
R1b = Rc2*fO2*25;           %OM oxidation
R2b = Rc2*fNO3*25;          %OM oxidation by NO3
R3b = Rc2*fFeOH3;           %OM oxidation by Fe(OH)3
%R4b = 0; % Rc2*fFeOOH;     %OM oxidation by FeOOH
R5b = Rc2*fSO4;             %OM oxidation by SO4
R6b = Rc2*fNO2*25;          %OM oxidation by NO2 AML
R7b = Rc2*fN2O*25;          %OM oxidation by N2O AML
R9b = Rc2*fNO2DNRA*25;      %OM oxidation by NO2 to NH4 AML


R1 = R1a+R1b;
R2 = R2a+R2b;
R3 = R3a+R3b;
%R4 = R4a+R4b;
R5 = R5a+R5b;
R6 = R6a+R6b;
R7 = R7a+R7b;
R9 = R9a+R9b;

Ra = R1a+R2a+R3a+R5a+R6a+R9a;     %degradation of OM1
Rb = R1b+R2b+R3b+R5b+R6b+R9b;     %degradation of OM2

R8 = kfeox*u(1)*u(5);             %Fe(II) oxidation by O2
R21 = knh4ox*(u(1)/(KmO2ao+u(1)))*(u(13)/(KmNH4ao+u(13))); %ammonia oxidation with MM kinetics AML
R22 = kno2ox*(u(1)/(KmO2no+u(1)))*(u(14)/(KmNO2no+u(14))); %nitrite oxidation with MM kinetics AML
R23 = kinO2/(kinO2+u(1))*kamx*u(13)/(KmNO2amx+u(13))*u(14)/(KmNH4amx+u(14)); %anammox with MM kinetics AML

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
	(D_bio+D_NO3)*DuDx(12)-w*u(12);...
	(D_bio+D_NH4)*DuDx(13)-w*u(13);...
    (D_bio+D_NO2)*DuDx(14)-w*u(14);...
    (D_bio+D_N2O)*DuDx(15)-w*u(15);...
    (D_bio+D_N2)*DuDx(16)-w*u(16)];%f
 
%Reaction: constant porosity----------------------------------------------- 
s = [(BC0_O2-u(1))*alfax + (-0.25*R8-2*R21) - R1*F;...               %the first part are all in unit of umol/cm3/yr
    R8/F + (-4*R3);...                                         %while the second part are all in unit of umol/g/yr
    0;...                                                       %this is the basis of using the conversion factor F
    (BC0_SO4-u(4))*alfax - 0.5*R5*F;...                    
    (BC0_Fe-u(5))*alfax - R8 + (4*R3)*F;...  
    (BC0_H2S-u(6))*alfax + (0.5*R5)*F;... 
    0;... 
    (BC0_S0-u(8))*alfax *F;...
    -Ra;...
    -Rb;...
	(BC0_PO4-u(11))*alfax + (P_index1*Ra+P_index2*Rb)*F;...
	(BC0_NO3-u(12))*alfax + R22 - 2*R2*F;...
	(BC0_NH4-u(13))*alfax - R21 - R23 + 0.67*R9*F + (N_index1*Ra+N_index2*Rb)*F;...
    (BC0_NO2-u(14))*alfax + 2*R2*F + R21 - R22 - R23 - 0.67*R9*F - 2*R6*F;... % NO2 produced by AO (R21), NO3 red (R2)
    %NO2 consumed by NO2 oxidation (R22), anammox (R23)dnra (R9) and NO2
    %red to N2O (R6)- AML
    (BC0_N2O-u(15))*alfax + R6*F - 2*R7*F;... % N2O produced by NO2 red to N2O (R6)reduced by R7
    (BC0_N2-u(16))*alfax + 2*R7*F];%s

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
       BC0_OM1 BC0_OM2 BC0_PO4 ...
       BC0_NO3 BC0_NH4 BC0_NO2 BC0_N2O BC0_N2 w F NumVars;

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
	  ul(12)-BC0_NO3;...
	  ul(13)-BC0_NH4;...
      ul(14)-BC0_NO2;...
      ul(15)-BC0_N2O;...
      ul(16)-BC0_N2];%pl

%1 for solid, 0 for solute for both constant and depth-dependant porosity
ql = [0;1;1;0;0;0;1;0;1;1;0;0;0;0;0;0];%ql

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
      w*ur(16)];%pr
  
qr = ones(NumVars,1);%qr
%----------------------------------------------------------------