
function MATSEDLAB_app_01
tic;
clear all;
close all;
clc;
% By defining the concentrations of the species, x and t as 'global' 
% variables there will be accessible through workspace
global sol x t fT

% 1-D diagenesis problem
m = 0;
% definition of the spatial domain, x=[0 15] cm with resolution of 300
%***** USER DEFINED *****%
x = linspace(0,15,500);
% definition of the spatial domain, t=[0 50] years with resolution of 155
%***** USER DEFINED *****%
t = linspace(0,200,255);   

disp('solving PDE ');        
%calling pdepe solver by passing the spatial-temporal domain to it
 sol = pdepe(m,@pdex14pde,@pdex1ic,@pdex1bc,x,t);
%Extract each species concentration at each time and depth
u1=sol(:,:,1);      % O2(aq)
u2=sol(:,:,2);      % Fe(OH)3(s)
u3=sol(:,:,3);      % SO4(2-)(aq)
u4=sol(:,:,4);      % Fe(2+)(aq)
u5=sol(:,:,5);      % S(-II)(aq)
u6=sol(:,:,6);      % FeS(s)
u7=sol(:,:,7);      % As(s)
u8=sol(:,:,8);      % As(III,V)(aq)

%Arsenic Flux calculation
noT=length(t);
noX=length(x);
fD=160*.97*(u8(noT,2)-u8(noT,1))/x(2);
    for j=1:noX
        Irr(j)=14.4*exp(-.25*x(j))*(u8(155,j)-u8(1,j))*.97;
    end
    fIrr=trapz(x,Irr);
    fT=fIrr+fD;
toc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BLOCK TWO%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,f,s] = pdex14pde(x,t,u,DuDx)
%The boundary values, sedimentation rate, bioturbation coef and conversion 
%factor (F) must be defined as global since they will be used in BLOCK FOUR
%If the problem is run for depth-dependant porosity then phi and Dphi have
%to set as global varibles too.
global BC0_FeOH3 BC0_O2 BC0_SO4 BC0_Fe BC0_H2S BC0_FeS BC0_AsFeOx BC0_AsO4 ...
       D_bio w F

% All the reaction and transport parameters and boundary values are defined 
% in this block
%***** USER DEFINED *****%
%boundary conditions at sediment-water interface
BC0_O2=.152;
BC0_Fe=0; 
BC0_H2S=0; 
BC0_FeOH3=6.7;
BC0_FeS=0;
BC0_AsO4=1e-6;
%for Steady-state simulation use the constant backgound concentrations
% BC0_SO4=.033; % Present day conditions 
% BC0_AsFeOx=BC0_FeOH3*0.32e-3;% Present day conditions 
%for non-steady state make the following uncomment
BC0_SO4 = 0.022 + 0.06*exp(-0.5*((t-182)/10)^2);
BC0_AsFeOx=2.14e-3 + 1.9e-3*exp(-0.5*((t-152)/6)^2);
D_bio=0.0694;% biturbation coef
%molecular diffusion coefs
D_O2=375;
D_SO4=175;
D_Fe=118;
D_H2S=284;
D_AsO4=160;
%half saturation coefs
KSO4=0.05;
KFeOH3=2000;
KO2=0.004;
%inhibition coefs
kinO2=3.2e-6;%3.2e-6
kinFeOH3=200;%200
%Secondary reaction constants
ktsox=1e3; %Range 1e2 - 1.6e6 validated 1e3
kfeox=4e4; %4e4
ktsfe=2.5; %Range 2.5-95 validated 2.5
kfedis=1e-3; 
kfepre=1500;
KFeS=1.78e3;%
kAsO4_ads=1.35; %Range 2.9-4.7 validated 1.35
kAs_FeS=1; %Range 1-10 validated 1
w=(.131*(t>62)+.095*(t<=62));%time-dependant burial rate
alfa0=14.4; % bioirrigation constant at sediment-water interface
alfax=alfa0*exp(-.25*x);% depth-dependant bioirrigation
h_plus=3.4e-4;%[H+] concentration equals 10^(-pH)
fAsFe=3.2e-4;% amount of As associated with Fe
%uncomment if porosity is depth-dependant
% phi=.9*exp(-0.2*x);
% Dphi=-.18*exp(-.2*x);
F=.06;%convertion factor=rhob*(1-fi)/fi; where fi=porosity and rhob=solid phase density
%contribution of each mineralization pathway
fO2=u(1)/(KO2+u(1));
fFeOH3=u(2)/(KFeOH3+u(2))*kinO2/(kinO2+u(1));
fSO4=u(3)/(KSO4+u(3))*kinO2/(kinO2+u(1))*kinFeOH3/(kinFeOH3+u(2));
%Saturation index for FeS precipitation
Sat_FeS=u(4)*u(5)/(KFeS*h_plus^2);
Rc=400*exp(-.1831*x);% depth dependant OM degredation
R1=Rc*fO2*26.5; % OM oxiation by O2, acceleration factor=26.5
R2=Rc*fFeOH3;%OM oxiation by Fe(OH)3
R3=Rc*fSO4;%OM oxiation by SO4
R4=ktsox*u(5)*u(1);%S(II) oxidation by O2
R5=kfeox*u(1)*u(4);%Fe(II) oxidation by O2
R6=ktsfe*u(2)*u(5);%Fe(OH)3 reduction by S(II)
if (Sat_FeS>=1)
    R7=0;
    R_7=kfepre*(Sat_FeS-1);%recipitation rate of FeS
else
    R7=kfedis*u(6)*(1-Sat_FeS);%dissolotion rate of FeS
    R_7=0;
end
R8=kAsO4_ads*u(2)*F*u(8);%As sorption on Fe(OH)3
R9=kAs_FeS*u(6)*u(8)*F;%As sorption on FeS
R10=117/F*u(5); %Sulfidization of OM
%adding nitrate to the reaction network
% KNO3=10;%NO3 half saturation
% kinNO3=10;
% fNO3=u(10)/(KNO3+u(10))*kinO2/(kinO2+u(1));
% R11=Rc*fNO3;  %OM oxiation by NO3 (denitrification)
% finNO3=kinNO3/(kinNO3+u(10)); % inhibition factor of NO3
% R12=knh4ox*u(1)*u(11);%Nitrification

%constant porosity
c=ones(1,8);
%depth dependant porosity
% c = [ phi;...
%       1-phi;...
%       phi;...
%       phi;...
%       phi;...
%       1-phi;...
%       phi;...
%       1-phi];

% Transport:constant porosity
f = [(D_bio+D_O2)*DuDx(1)-w*u(1);...
    D_bio*DuDx(2)-w*u(2);...
    (D_bio+D_SO4)*DuDx(3)-w*u(3);...
    (D_bio+D_Fe)*DuDx(4)-w*u(4);...
    (D_bio+D_H2S)*DuDx(5)-w*u(5);...
    D_bio*DuDx(6)-w*u(6);...
    D_bio*DuDx(7)-w*u(7);...
    (D_bio+D_AsO4)*DuDx(8)-w*u(8)];
% Transport:depth-dependant porosity
% f = [((D_bio+D_O2)*DuDx(1)-w*u(1))*phi+D_bio*Dphi*u(1);...
%     (D_bio*DuDx(2)-w*u(2))*(1-phi)-D_bio*Dphi*u(2);...
%     ((D_bio+D_SO4)*DuDx(3)-w*u(3))*phi+D_bio*Dphi*u(3);...
%     ((D_bio+D_Fe)*DuDx(4)-w*u(4))*phi+D_bio*Dphi*u(4);...
%     ((D_bio+D_H2S)*DuDx(5)-w*u(5))*phi+D_bio*Dphi*u(5);...
%     (D_bio*DuDx(6)-w*u(6))*(1-phi)-D_bio*Dphi*u(6);...
%     ((D_bio+D_AsO4)*DuDx(7)-w*u(7))*phi+D_bio*Dphi*u(7);...
%     (D_bio*DuDx(8)-w*u(8))*(1-phi)-D_bio*Dphi*u(8)];
 
% Reaction: constant porosity 
s = [(BC0_O2-u(1))*alfax-F*R1-.25*R5-2*R4;...
    -4*R2+R5/F-8*R6;...
    (BC0_SO4-u(3))*alfax+F*(R6-.5*R3)+R4;...
    (BC0_Fe-u(4))*alfax+F*(4*R2+8*R6+R7-R_7)-R5;...
    (BC0_H2S-u(5))*alfax+F*(.5*R3-R6+R7-R_7-R10);...
    -R7+R_7;...
    (-4*fAsFe*R2)+(-8*fAsFe*R6)+(R8+R9)/F;...
    (BC0_AsO4-u(8))*alfax+F*fAsFe*(4*R2+8*R6)-R8-R9];
% Reaction: depth-dependant porosity
% s = [(BC0_O2-u(1))*alfax-F*R1-.25*R5-2*R4;...
%     -4*R2+R5/F-8*R6;...
%     (BC0_SO4-u(3))*alfax+F*(R6-.5*R3)+R4;...
%     (BC0_Fe-u(4))*alfax+F*(4*R2+8*R6+R7-R_7)-R5;...
%     (BC0_H2S-u(5))*alfax+F*(.5*R3-R6+R7-R_7-R10);...
%     -R7+R_7;...
%     (-4*fAsFe*R2)+(-8*fAsFe*R6)+(R8+R9)/F;...
%     (BC0_AsO4-u(8))*alfax+F*fAsFe*(4*R2+8*R6)-R8-R9].*c;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BLOCK THREE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u0 = pdex1ic(x) 
%***** USER DEFINED *****%
% the default is zero initial concentrations for all the species. it can be
% modified to assign different background concentrations for any of them.
u0=zeros(1,8);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%BLOCK FOUR%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
  
global BC0_FeOH3 BC0_O2 BC0_SO4 BC0_Fe BC0_H2S BC0_FeS BC0_AsFeOx BC0_AsO4 w F

%Upper bounday: constant porosity
pl = [ul(1)-BC0_O2;...
      BC0_FeOH3/F;...
      ul(3)-BC0_SO4;...
      ul(4)-BC0_Fe;...
      ul(5)-BC0_H2S;...
      BC0_FeS/F;...
      BC0_AsFeOx/F;...
      ul(8)-BC0_AsO4];
  
 %Upper bounday: depth-dependant porosity
 %rhob=2;% dry density of the sediment
%   pl = [ul(1)-BC0_O2;...
%       BC0_FeOH3/rhob;...
%       ul(3)-BC0_SO4;...
%       ul(4)-BC0_Fe;...
%       ul(5)-BC0_H2S;...
%       BC0_FeS/rhob;...
%       BC0_AsFeOx/rhob;...
%       ul(8)-BC0_AsO4];

%1 for solid, 0 for solute for both constant and depth-dependant porosity
ql = [0;1;0;0;0;1;1;0];

%Lower bounday: constant porosity
pr = [w*ur(1);...
      w*ur(2);...
      w*ur(3);...
      w*ur(4);...
      w*ur(5);...
      w*ur(6);...
      w*ur(7);...
      w*ur(8)];
  
%Lower bounday: depth-dependant porosity  
% pr = [w*ur(1)*phi-D_bio*Dphi*ur(1);...
%       w*ur(2)*(1-phi)+D_bio*Dphi*ur(2);...
%       w*ur(3)*phi-D_bio*Dphi*ur(3);...
%       w*ur(4)*phi-D_bio*Dphi*ur(4);...
%       w*ur(5)*phi-D_bio*Dphi*ur(5);...
%       w*ur(6)*(1-phi)+D_bio*Dphi*ur(6);...
%       w*ur(7)*(1-phi)+D_bio*Dphi*ur(7);...
%       w*ur(8)*phi-D_bio*Dphi*ur(8)];
  
qr=ones(1,8);
%---------------------------------------------------------------