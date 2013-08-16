%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   MATSEDLAB
%
% In script MATSEDLAB_app_02.m, the description of organic matter degradation
% is based on the classic G-model (Westrich and Berner, 1984). Two pools
% of reactive organic carbon are assumed. OM1 (reactive) and OM2 (refractory) 
% are assigned first-order degradation constants kOM1 and kOM2 of 0.1 
% and 0.001 yr-1, respectively (Canavan et al., 2006). The deposition fluxes 
% of the two pools of organic carbon, JOM1 and JOM2, are specified as upper
% boundary conditions. (Note: here, the reactions involving As have been removed.) 
% The MATSEDLAB_app_02.m script presents two options for JOM1 and JOM2, which
% can be switched on or off by the user (see lines 104-105). In the first
% option, JOM1 and JOM2 are assigned constant values of 300 and 150
% µmol C cm-2 yr-1. In the second option, the deposition fluxes JOM1 and 
% JOM2 are assumed to peak in early summer with the maximum value of 600 and
% 300 µmol C cm-2 yr-1 respectively, according to the sinusoidal boundary
% conditions, JOM1=300[cos?(2?t)+1] and JOM2=150[cos?(2?t)+1], where t=0 
% corresponds to July 1, and the depositional fluxes drop to zero in January 
%
% Further documentation can be found at: http://tinyurl.com/matsedlab
%
function MATSEDLAB_app_02
tic;
% clear all;
close all;
clc;

 global u1 u2 u3 u4 u5 u6 u7 u8 x t sol

% 1-D diagenesis problem
m = 0;
% definition of the spatial domain, x=[0 15] cm
x = linspace(0,15,250);
% definition of the spatial domain, t=[0 50] years
t = linspace(0,20,200);
     
% %calling pdepe solver by passing the spatial-temporal domain to it
sol = pdepe(m,@pdex14pde,@pdex1ic,@pdex1bc,x,t);
% % Extract the first solution component as u.
u1=sol(:,:,1);      % OM1(s)
u2=sol(:,:,2);      % OM2(s)
u3=sol(:,:,3);      % O2(aq)
u4=sol(:,:,4);      % Fe(OH)3(s)
u5=sol(:,:,5);      % SO4(2-)(aq)
u6=sol(:,:,6);      % Fe(2+)(aq)
u7=sol(:,:,7);      % S(-II)(aq)
u8=sol(:,:,8);      % FeS(s)

VarNames = {'OM1',...
            'OM2',...
            'O2(aq)',...
            'Fe(OH)3(s)', ... 
            'SO4(2-)(aq)', ... 
            'Fe(2+)(aq)', ... 
            'S(-II)(aq)',...
            'FeS(s)'};
         
ql = [1;1;0;1;0;0;0;1]; 
%Plotting all the species vs. depth
for i=1:8
    MatValues = sol(:,:,i);
    [m, n]= size(MatValues);
    subplot(2,4,i);plot(MatValues(m,:),x);set(gca,'YDir','reverse');
    title(char(VarNames(i)), 'FontSize', 14, 'FontWeight','bold');
     if(ql(i)== 0)
            %aq phase
            xlabel('Concentration (\mumol/cm^3)','FontSize', 14);
        else
            %s phase
            xlabel('Concentration (\mumol/g)','FontSize', 14);
     end
     if ((i==1)||(i==5))
        ylabel('Depth (cm)','FontSize', 14);
     end   
end    
% Plotting FeS profiles for 4 yrs+Figure of FeS and SO4
% FeSplot(u8(m/4,:), x, u8(m/2,:), u8(m/4*3,:),u8(m,:),u5(m/4,:),u5(m,:));

%time evolution of FeS
subplot(4,1,1);plot(t,u8(:,1)); ylabel('FeS (\mumol/g)','fontsize',14);title('X=0','fontsize',12);
subplot(4,1,2);plot(t,u8(:,84)); ylabel('FeS (\mumol/g)','fontsize',14);title('X=5 cm','fontsize',12);
subplot(4,1,3);plot(t,u8(:,167));ylabel('FeS (\mumol/g)','fontsize',14);title('X=10 cm','fontsize',12);
subplot(4,1,4);plot(t,u8(:,250));xlabel('time (yr)','fontsize',14); ylabel('FeS (\mumol/g)','fontsize',14);title('X=15 cm','fontsize',12);

%Aquouse phases vs time at X=5cm from Jan to Jul 
subplot(4,1,1);plot(t(5:15),u3(5:15,84)); ylabel('FeS (\mumol/g)','fontsize',14);title('X=0','fontsize',12);
subplot(4,1,2);plot(t(5:15),u5(5:15,84)); ylabel('FeS (\mumol/g)','fontsize',14);title('X=5 cm','fontsize',12);
subplot(4,1,3);plot(t(5:15),u6(5:15,84));ylabel('FeS (\mumol/g)','fontsize',14);title('X=10 cm','fontsize',12);
subplot(4,1,4);plot(t(5:15),u7(5:15,84));xlabel('time (yr)','fontsize',14); ylabel('FeS (\mumol/g)','fontsize',14);title('X=15 cm','fontsize',12)

subplot(1,4,1);plot(u3(9,:),x,'g',u3(6,:),x,'r'); set(gca,'YDir','reverse');ylabel('depth(cm))','fontsize',14);xlabel('O_2 (\mumol/cm^3)','fontsize',14);
subplot(1,4,2);plot(u5(9,:),x,'g',u5(6,:),x,'r'); set(gca,'YDir','reverse');xlabel('SO_4^-^2 (\mumol/cm^3)','fontsize',14);
subplot(1,4,3);plot(u6(9,:),x,'g',u6(6,:),x,'r');set(gca,'YDir','reverse');xlabel('S(II) (\mumol/cm^3)','fontsize',14);
subplot(1,4,4);plot(u7(9,:),x,'g',u7(6,:),x,'r');set(gca,'YDir','reverse');xlabel('Fe(II) (\mumol/cm^3)','fontsize',14);
% --------------------------------------------------------------
function [c,f,s] = pdex14pde(x,t,u,DuDx)

global BC0_OM1 BC0_OM2 BC0_FeOH3 BC0_O2 BC0_SO4 BC0_Fe BC0_H2S BC0_FeS  ...
       D_bio D_O2 D_SO4 D_Fe D_H2S  w

%boundary conditions at sediment-water interface
BC0_OM1=(300*cos(t*(2*pi))+300);
BC0_OM2=BC0_OM1*0.5;
% BC0_OM1=300;
% BC0_OM2=BC0_OM1*0.5;
BC0_O2=.152;
BC0_Fe=0; % 2e-3
BC0_H2S=0; % 1e-5
BC0_FeOH3=6.7;
BC0_FeS=0;
BC0_SO4=.033; % Present day conditions 

D_bio=0.0694;% biturbation coef
%molecular diffusion coefs
D_O2=375;
D_SO4=175;
D_Fe=118;
D_H2S=284;

kOM1=10^-1;%degredation rate constant OM1
kOM2=10^-3;%degredation rate constant OM2

%half saturation coefs
KSO4=.05;%0.05
KFeOH3=2000;
KO2=.004;%.004
%inhibition coefs
kinO2=3.2e-6;%3.2e-6
kinFeOH3=200;%200
%eaction constants
ktsox=1e3; %Range 1e2 - 1.6e6 validated 1e3
kfeox=4e4; %4e4
ktsfe=2.5; %Range 2.5-95 validated 2.5
kfedis=1e-3; 
kfepre=1500;
KFeS=1.78e3;%

w=0.131*8;%burial rate
alfa0=14.4*.1; % bioirrigation constant at sediment-water interface
alfax=alfa0*exp(-.25*x);% depth-dependant bioirrigation
F=.06;%convertion factor=rhob*(1-fi)/fi; where fi=porosity and rhob=solid phase density
h_plus=3.4e-4;%[H+] concentration equals 10^(-pH)
%contribution of each mineralization pathway
fO2=u(3)/(KO2+u(3));
fFeOH3=u(4)/(KFeOH3+u(4))*kinO2/(kinO2+u(3));
fSO4=u(5)/(KSO4+u(5))*kinO2/(kinO2+u(3))*kinFeOH3/(kinFeOH3+u(4));
%Saturation index for FeS precipitation
Sat_FeS=u(6)*u(7)/(KFeS*h_plus^2);
% depth dependant OM degredation
Rc1=kOM1*u(1);% first pool rate
Rc2=kOM1*u(2);% 2nd pool rate

R1a=Rc1*fO2*26.5; % OM oxiation by O2, acceleration factor=26.5
R1b=Rc2*fO2*26.5;
R2a=Rc1*fFeOH3;%OM oxiation by Fe(OH)3
R2b=Rc2*fFeOH3;%OM oxiation by Fe(OH)3
R3a=Rc1*fSO4;%OM oxiation by SO4
R3b=Rc2*fSO4;%OM oxiation by SO4

R4=ktsox*u(7)*u(3);%S(II) oxidation by O2
R5=kfeox*u(3)*u(6);%Fe(II) oxidation by O2
R6=ktsfe*u(4)*u(7);%Fe(OH)3 reduction by S(II)

if (Sat_FeS>=1)
    R7=0;
    R_7=kfepre*(Sat_FeS-1);%precipitation rate of FeS
else
    R7=kfedis*u(8)*(1-Sat_FeS);%dissolotion rate of FeS
    R_7=0;
end

c=ones(1,8);

% Transport
f = [D_bio*DuDx(1)-w*u(1);...
    D_bio*DuDx(2)-w*u(2);...
    (D_bio+D_O2)*DuDx(3)-w*u(3);...
    D_bio*DuDx(4)-w*u(4);...
    (D_bio+D_SO4)*DuDx(5)-w*u(5);...
    (D_bio+D_Fe)*DuDx(6)-w*u(6);...
    (D_bio+D_H2S)*DuDx(7)-w*u(7);...
    D_bio*DuDx(8)-w*u(8)];
 
% Reaction 
s = [-(R1a+R2a+R3a);...
    -(R1b+R2b+R3b);...
    (BC0_O2-u(3))*alfax-0.06*(R1a+R1b)-.25*R5-2*R4;...
    -4*(R2a+R2b)+R5/0.06-8*R6;...
    (BC0_SO4-u(5))*alfax+0.06*(R6-.5*(R3a+R3b))+R4;...
    (BC0_Fe-u(6))*alfax+0.06*(4*(R2a+R2b)+8*R6+R7-R_7)-R5;...
    (BC0_H2S-u(7))*alfax+0.06*(.5*(R3a+R3b)-R6+R7-R_7);...
    -R7+R_7];
    
% --------------------------------------------------------------
function u0 = pdex1ic(x) 
u0=zeros(1,8);
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
  
global BC0_OM1 BC0_OM2 BC0_FeOH3 BC0_O2 BC0_SO4 BC0_Fe BC0_H2S BC0_FeS  w
F=0.06;
pl = [BC0_OM1/F
      BC0_OM2/F
      ul(3)-BC0_O2;...
      BC0_FeOH3/F;...
      ul(5)-BC0_SO4;...
      ul(6)-BC0_Fe;...
      ul(7)-BC0_H2S;...
      BC0_FeS/F];

%1 for solid, 0 for solute
ql = [1;1;0;1;0;0;0;1];

pr = [w*ur(1);...
      w*ur(2);...
      w*ur(3);...
      w*ur(4);...
      w*ur(5);...
      w*ur(6);...
      w*ur(7);...
      w*ur(8)];
  
qr=ones(1,8);
