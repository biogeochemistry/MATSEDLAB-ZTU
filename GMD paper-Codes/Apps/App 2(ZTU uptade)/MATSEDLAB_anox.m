function MATSEDLAB_anox
tic;                        %start timing

global sol sol1 SimValues u_ox NumVars;
global n_segment;

%***** USER DEFINED *****%
x = linspace(0,15,511);         % definition of the spatial domain, x=[0 15] cm with resolution of 511
t = linspace(0,0.5,3);        % definition of the time domain, t=[0 400] years with resolution of 401

u_ox=zeros(NumVars,511);
for i1=1:NumVars
    u_ox(i1,:) = sol1(end,:,i1);
end

disp('solving PDE ');                               %calling pdepe solver by passing the spatial-temporal domain to it
sol = pdepe(0,@pdex14pde,@pdex1ic,@pdex1bc,x,t);    %Extract each species concentration at each time and depth
toc;                                                %stop timing

%rt_time=clock;
%fprintf('The ending time is %.0f/%.0f/%.0f %.0f:%.0f:%.0f\n', rt_time(1),...
%    rt_time(2), rt_time(3), rt_time(4), rt_time(5), rt_time(6));

for j=1:NumVars
    SimValues{j,n_segment} = sol(:,:,j);
end

n_segment = n_segment+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BLOCK TWO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,f,s] = pdex14pde(x,~,u,DuDx)
global BC0_O2 BC0_FeOH3 BC0_FeOOH BC0_SO4 BC0_Fe BC0_H2S BC0_FeS BC0_S0...
       BC0_OM1 BC0_OM2 BC0_AsFea BC0_AsFeb BC0_AsO4 BC0_PO4...
       BC0_PO4adsa BC0_PO4adsb BC0_Ca2 BC0_CaPO4 w F NumVars;
global para;
%***** USER DEFINED *****%
%boundary conditions at sediment-water interface---------------------------
BC0_O2 = 0.152;
BC0_FeOH3 = 6.7;
BC0_FeOOH = 4.6;
BC0_SO4 = 0.02;
BC0_Fe = 0; 
BC0_H2S = 0; 
BC0_FeS = 0;
BC0_S0 = 0;
BC0_OM1 = 50;
BC0_OM2 = BC0_OM1*0.5;
fAsFe = 1e-6;                         %Amount of As associated with Fe
BC0_AsFea = BC0_FeOH3*fAsFe;
BC0_AsFeb = BC0_FeOOH*fAsFe;
BC0_AsO4 = 3e-4;
BC0_PO4 = 3e-4;
fPFe = 1e-6;
BC0_PO4adsa = BC0_FeOH3*fPFe;
BC0_PO4adsb = BC0_FeOOH*fPFe;
BC0_Ca2 = 4e-2;
BC0_CaPO4 = 0;
%bioturbatiuon coef--------------------------------------------------------
D_bio = 0.0694;
%molecular diffusion coefs-------------------------------------------------
D_O2 = 375; 
D_SO4 = 175;
D_Fe = 118;
D_H2S = 284;
D_S0 = 100;
D_AsO4 = 160;
D_PO4 = 232;
D_Ca2 = 212;
%xyz ratio-----------------------------------------------------------------
Cx1 = para(1);              %each value has been replaced with the         
Pz1 = para(2);              %corresponding entry in 'para'
Cx2 = para(3);
Pz2 = para(4);
%half saturation coefs-----------------------------------------------------
KO2 = para(5);                  
KFeOH3 = para(6);               
KFeOOH = para(7);
KSO4 = para(8);
%inhibition coefs----------------------------------------------------------
kinO2 = para(9);
kinFeOH3 = para(10);
kinFeOOH = para(11);
%saturation constants------------------------------------------------------
KFeS = para(12);
%Secondary reaction constants----------------------------------------------
kOM1 = 1;                           %degredation rate constant for OM1, Rc1
kOM2 = 0.01;                        %degredation rate constant for OM2, Rc2
ktsox = para(13);                   %S(II) oxidation by O2, R5
ktsfe = para(14);                   %Fe(OH)3 reduction by S2-, R6
kfeox = para(15);                   %Fe(II) oxidation by O2, R7
kfedis = para(16);                  %dissolution rate of FeS, R8
kfepre = para(17);                  %precipitation rate of FeS, R_8
ksorg = para(18);                   %Sulfidization of OM, R9
kAsO4_adsa = para(19);              %As sorption onto Fe(OH)3, R15
kAsO4_adsb = para(20);              %As sorption onto FeOOH, R16
kAs_FeS = para(21);                 %As sorption onto FeS, R17
k_psorb_a = para(22);               %PO4 sorption onto Fe(OH)3, R22
k_psorb_b = para(23);               %PO4 sorption onto FeOOH, R23
kapa = para(24);                    %Apatite precipitation, R25
k_apa = para(25);

%other parameters----------------------------------------------------------
w = 0.1;     %time-dependant burial rate
alfa0 = 14.4;                           %bioirrigation constant at sediment-water interface
alfax = alfa0*exp(-0.25*x);             %depth-dependant bioirrigation
h_plus = 3.4e-4;                        %[H+] concentration, equals 10^(-pH)
F = 0.06;                               %convertion factor=rhob*(1-fi)/fi; where fi=porosity and rhob=solid phase density

%indices-------------------------------------------------------------------
fO2 = u(1)/(KO2+u(1));
fFeOH3 = u(2)/(KFeOH3+u(2))*kinO2/(kinO2+u(1));
fFeOOH = u(3)/(KFeOOH+u(3))*kinO2/(kinO2+u(1))*kinFeOH3/(kinFeOH3+u(2));
fSO4 = u(4)/(KSO4+u(4))*kinO2/(kinO2+u(1))*kinFeOH3/(kinFeOH3+u(2))*kinFeOOH/(kinFeOOH+u(3));
Sat_FeS = u(5)*u(6)/(KFeS*h_plus^2);    %saturation index of FeS
P_index1 = Pz1/Cx1;
P_index2 = Pz2/Cx2;

%reaction rates------------------------------------------------------------
Rc1 = kOM1*u(9);         %first pool rate
R1a = Rc1*fO2*25;         %OM oxidation
R2a = Rc1*fFeOH3;         %OM oxidation by Fe(OH)3
R3a = Rc1*fFeOOH;         %OM oxidation by FeOOH
R4a = Rc1*fSO4;           %OM oxidation by SO4

Rc2 = kOM2*u(10);         %second pool rate
R1b = Rc2*fO2*25;         %OM oxidation
R2b = Rc2*fFeOH3;         %OM oxidation by Fe(OH)3
R3b = Rc2*fFeOOH;         %OM oxidation by FeOOH
R4b = Rc2*fSO4;           %OM oxidation by SO4

R1 = R1a+R1b;
R2 = R2a+R2b;
R3 = R3a+R3b;
R4 = R4a+R4b;

Ra = R1a+R2a+R3a+R4a;                   %degradation of OM1
Rb = R1b+R2b+R3b+R4b;                   %degradation of OM2

R5 = ktsox*u(1)*u(6);                   %S(II) oxidation by O2
R6 = ktsfe*u(2)*u(6);                   %Fe(OH)3 reduction by S2-
R7 = kfeox*u(1)*u(5);                   %Fe(II) oxidation by O2

if (Sat_FeS>=1)
    R8 = 0;
    R_8 = kfepre*(Sat_FeS-1);           %precipitation rate of FeS
else
    R8 = kfedis*u(7)*(1-Sat_FeS);       %dissolution rate of FeS
    R_8 = 0;
end

R9 = ksorg*u(6)*(u(9)+u(10));           %Sulfidization of OM
R10 = kAsO4_adsa*u(2)*u(13);            %As sorption onto Fe(OH)3
R11 = kAsO4_adsb*u(3)*u(13);            %As sorption onto FeOOH
R12 = kAs_FeS*u(7)*u(13);               %As sorption onto FeS
R13 = fAsFe*(4*R2+2*R6);                %As release from pool a
R14 = fAsFe*4*R3;                       %As release from pool b
R15 = k_psorb_a*u(2)*u(14);             %PO4 sorption onto Fe(OH)3
R16 = k_psorb_b*u(3)*u(14);             %PO4 sorption onto FeOOH
R17 = fPFe*(4*R2+2*R6);                 %PO4 release from pool a
R18 = fPFe*4*R3;                        %PO4 release from pool b

if u(14)>kapa
    R19 = k_apa*(u(14)-kapa);           %Apatite precipitation
else
    R19 = 0;
end

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
    D_bio*DuDx(11)-w*u(11);...
    D_bio*DuDx(12)-w*u(12);...
	(D_bio+D_AsO4)*DuDx(13)-w*u(13);...
	(D_bio+D_PO4)*DuDx(14)-w*u(14);...
	D_bio*DuDx(15)-w*u(15);...
	D_bio*DuDx(16)-w*u(16);...
	(D_bio+D_Ca2)*DuDx(17)-w*u(17);...
	D_bio*DuDx(18)-w*u(18)];%f
 
%Reaction: constant porosity----------------------------------------------- 
s = [(BC0_O2-u(1))*alfax + (-2*R5-0.25*R7) - R1*F;...               %the first part are all in unit of umol/cm3/yr
    R7/F + (-4*R2-2*R6);...                                         %while the second part are all in unit of umol/g/yr
    -4*R3;...                                                       %this is the basis of using the conversion factor F
    (BC0_SO4-u(4))*alfax + R5 - 0.5*R4*F;...                    
    (BC0_Fe-u(5))*alfax - R7 + (4*R2+4*R3+2*R6+R8-R_8)*F;...  
    (BC0_H2S-u(6))*alfax - R5 + (0.5*R4-R6+R8-R_8-R9)*F;... 
    -R8+R_8;... 
    (BC0_S0-u(8))*alfax + R6*F;...
    -Ra;...
    -Rb;...
	+R10+R12-R13;...
    +R11-R14;...
	(BC0_AsO4-u(13))*alfax + (-R10-R11-R12+R13+R14)*F;...
	(BC0_PO4-u(14))*alfax - 2*R19 + (-R15-R16+P_index1*Ra+P_index2*Rb+R17+R18)*F;...
	+R15-R17;...
	+R16-R18;...
	(BC0_Ca2-u(17))*alfax - 3*R19;...
	+R19/F];%s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BLOCK THREE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u0 = pdex1ic(x) 
global u_ox
u0=u_ox(:,round(x/15*510+1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BLOCK FOUR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pl,ql,pr,qr] = pdex1bc(~,ul,~,ur,~)       %boundary condition
global BC0_O2 BC0_FeOH3 BC0_FeOOH BC0_SO4 BC0_Fe BC0_H2S BC0_FeS BC0_S0...
       BC0_OM1 BC0_OM2 BC0_AsFea BC0_AsFeb BC0_AsO4 BC0_PO4...
       BC0_PO4adsa BC0_PO4adsb BC0_Ca2 BC0_CaPO4 w F NumVars;

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
	  BC0_AsFea/F;...
	  BC0_AsFeb/F;...
	  ul(13)-BC0_AsO4;...
	  ul(14)-BC0_PO4;...
	  BC0_PO4adsa/F;...
	  BC0_PO4adsb/F;...
	  ul(17)-BC0_Ca2;...
	  BC0_CaPO4/F];%pl

%1 for solid, 0 for solute for both constant and depth-dependant porosity
ql = [0;1;1;0;0;0;1;0;1;1;1;1;0;0;1;1;0;1];%ql

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
	  w*ur(17);...
	  w*ur(18)];%pr
  
qr = ones(NumVars,1);%qr
%----------------------------------------------------------------