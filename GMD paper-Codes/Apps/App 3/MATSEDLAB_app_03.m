%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   MATSEDLAB
%    
% In MATASEDLAB_app_03, we modify the baseline code to represent Fe(III)-bearing mineral 
% phases. The phases selected are Fe(OH)3, goethite/lepidocrocite (lumped in one phase)
% and magnetite. The reactions involving the Fe(III) minerals are summarized in Table 7. 
% The formulation are taken from Pallud et al. (2010), who studied the reductive 
% dissolution of ferrihydrite and the subsequent formation of secondary 
% Fe minerals in soil aggregates. The reaction stoichiometries and reaction 
% rate expressions are implemented in MATSEDLAB_app_03.m on lines 162-168.
% respectively.    
%
% Further documentation can be found at: http://tinyurl.com/matsedlab
%
function MATSEDLAB_app_03
tic;
clear all;
close all;
clc;
global u1 u2 u3 u4 u5 u6 u7 u8 x t
m = 0;
x = linspace(0,15,500);
t = linspace(0,100,450);%300 for fT calculations 100 for baseline

sol = pdepe(m,@pdex14pde,@pdex1ic,@pdex1bc,x,t);

u1=sol(:,:,1);      % O2(aq)
u2=sol(:,:,2);      % Ferrihydrite
u3=sol(:,:,3);      % SO4(2-)(aq)
u4=sol(:,:,4);      % Fe(2+)(aq)
u5=sol(:,:,5);      % S(-II)(aq)
u6=sol(:,:,6);      % FeS(s)
u7=sol(:,:,7);      %Goethte
u8=sol(:,:,8);      %Magnetite

subplot(2,2,1);plot(u2(450,:),x);set(gca,'YDir','reverse');xlabel('Concentration (\mumol/gr)','fontsize',12); ylabel('Depth (cm)','fontsize',12);title('Ferrihydrite','fontsize',12);
subplot(2,2,2);plot(u7(450,:),x);set(gca,'YDir','reverse');xlabel('Concentration (\mumol/gr)','fontsize',12); ylabel('Depth (cm)','fontsize',12);title('Goethite','fontsize',12);
subplot(2,2,3);plot(u8(450,:),x);set(gca,'YDir','reverse');xlabel('Concentration (\mumol/gr)','fontsize',12); ylabel('Depth (cm)','fontsize',12);title('Magnetite','fontsize',12);
subplot(2,2,4);plot(u4(450,:),x);set(gca,'YDir','reverse');xlabel('Concentration (\mumol/cm3)','fontsize',12); ylabel('Depth (cm)','fontsize',12);title('Fe(II)','fontsize',12);
% subplot(3,2,5);plot(u3(450,:),x);set(gca,'YDir','reverse');xlabel('Concentration (\mumol/cm3)'); ylabel('Depth (cm)');title('Sulfate');
% subplot(3,2,6);plot(u5(450,:),x);set(gca,'YDir','reverse');xlabel('Concentration (\mumol/cm3)'); ylabel('Depth (cm)');title('S(II)');
h_plus=3.4e-4;

Ksorp=10^-2.5;%effective sorption constant
SA=600;%surface area (m2/gFe(OH)3)
SD=3.5*10^-6;%site density(mol sites/m2)
MW=168.70;%gr/mol
SFeOH3=u2(:,:)*SA*SD*MW;
Fe_sorb=Ksorp*SFeOH3.*u4(:,:)./(h_plus+Ksorp*u4(:,:));
Fe_tot=u2(:,:)+u7(:,:)+u8(:,:);
fer_per=u2(:,:)./Fe_tot;
goe_per=u7(:,:)./Fe_tot;
magn_per=u8(:,:)./Fe_tot;
Fe_per=u4(:,:)./Fe_tot;
Fe_sorb_per=Fe_sorb./Fe_tot;
figure;
% subplot(3,2,1);plot(fer_per(450,:),x);set(gca,'YDir','reverse');xlabel('percentage of Iron'); ylabel('Depth (cm)');title('Ferrihydrite');
% subplot(3,2,2);plot(goe_per(450,:),x);set(gca,'YDir','reverse');xlabel('percentage of Iron'); ylabel('Depth (cm)');title('Goethite');
% subplot(3,2,3);plot(magn_per(450,:),x);set(gca,'YDir','reverse');xlabel('percentage of Iron'); ylabel('Depth (cm)');title('Magnetite');
% subplot(3,2,4);plot(Fe_per(450,:),x);set(gca,'YDir','reverse');xlabel('percentage of Iron'); ylabel('Depth (cm)');title('Fe(II)');
% subplot(3,2,5);plot(Fe_sorb_per(450,:),x);set(gca,'YDir','reverse');xlabel('percentage of Iron'); ylabel('Depth (cm)');title('Fe_s_o_r_b(II)');
subplot(2,2,[1 3]);plot(fer_per(450,:),x);set(gca,'YDir','reverse');xlabel('percentage of Iron','fontsize',12); ylabel('Depth (cm)','fontsize',12);title('Ferrihydrite','fontsize',12);box off;
subplot(2,2,2);plot(goe_per(450,:),x);set(gca,'YDir','reverse');ylabel('Depth (cm)','fontsize',12);title('Goethite','fontsize',12);xlim([0 max(max(goe_per(:,:)))]);ylim([0 15]);box off;
subplot(2,2,4);plot(magn_per(450,:),x);set(gca,'YDir','reverse'); ylabel('Depth (cm)','fontsize',12);title('Magnetite','fontsize',12);xlim([0 .01]);ylim([0 15]);box off;

figure;
% subplot(2,4,1:4);plot(t,fer_per(:,500));ylabel('percentage of Iron','fontsize',12); xlabel('Time (yr)','fontsize',12);title('Ferrihydrite','fontsize',12);box off;
% subplot(2,4,5);plot(t,goe_per(:,500));ylabel('percentage of Iron','fontsize',12); xlabel('Time (yr)','fontsize',12);title('Goethite','fontsize',12);box off;
% subplot(2,4,6);plot(t,magn_per(:,500));ylabel('percentage of Iron','fontsize',12); xlabel('Time (yr)','fontsize',12);title('Magnetite','fontsize',12);box off;
% subplot(2,4,7);plot(t,Fe_per(:,500));ylabel('percentage of Iron','fontsize',12); xlabel('Time (yr)','fontsize',12);title('Fe(II)','fontsize',12);box off;
% subplot(2,4,8);plot(t,Fe_sorb_per(:,500));ylabel('percentage of Iron','fontsize',12); xlabel('Time (yr)','fontsize',12);title('Fe_s_o_r_b(II)','fontsize',12);box off;
subplot(2,2,1:2);plot(t,fer_per(:,500));ylabel('percentage of Iron','fontsize',12); xlabel('Time (yr)','fontsize',12);title('Ferrihydrite','fontsize',12);box off;
subplot(2,2,3);plot(t,goe_per(:,500));ylabel('percentage of Iron','fontsize',12); xlabel('Time (yr)','fontsize',12);title('Goethite','fontsize',12);box off;
subplot(2,2,4);plot(t,magn_per(:,500));ylabel('percentage of Iron','fontsize',12); xlabel('Time (yr)','fontsize',12);title('Magnetite','fontsize',12);box off;

toc;
% --------------------------------------------------------------
function [c,f,s] = pdex14pde(x,t,u,DuDx)

global BC0_FeOH3 BC0_O2 BC0_SO4 BC0_Fe BC0_H2S BC0_FeS ...
       D_bio D_O2 D_SO4 D_Fe D_H2S  w
      
BC0_O2=.152;
BC0_Fe=0; % 2e-3
BC0_H2S=0; % 1e-5
BC0_FeOH3=6.7;
BC0_FeS=0;

% BC0_SO4=.033; % Steady-state sulfate 
 BC0_SO4 = 0.022 + 0.06*exp(-0.5*((t-182)/10)^2);

D_bio=0.0694;
D_O2=375;
D_SO4=175;
D_Fe=118;
D_H2S=284;

KSO4=0.05;%0.05
KFeOH3=2000;
KO2=0.004;

kinO2=3.2e-6;
kinFeOH3=200;

ktsox=1e3; 
kfeox=4e4; 
ktsfe=2.5; 
kfedis=1e-3; 
kfepre=1500;
KFeS=1.78e3;

w=.131*(t>62)+.095*(t<=62);
alfa0=14.4; % 
alfax=alfa0*exp(-.25*x);
F=.06;
h_plus=3.4e-4;

% Km_lac=10;%half saturation constant for lactate (UM)
% Km_FeOH3=10^-4;%half saturation constant for ferrihydrite (mol Fe/gr solid)
Ksorp=10^-2.5;%effective sorption constant
SA=600;%surface area (m2/gFe(OH)3)
SD=3.5*10^-6;%site density(mol sites/m2)
MW=168.70;%gr/mol
% rho_solid=2.83;%solid phase density (gr/cm3)
% kresp=4.8*10^-19;%mol/cell/s
% krespmag=4.8*10^-19;%mol/cell/s;
% fi=.90;%porosity
kFeOOH=7.2*10^-5*365*24*3600/10000;%1/yr
kMag=4.5*10^-6*365*24*3600/10000;%1/yr
Km_mag=2000;%(Umol/gr)

%total surface available for sorption of Fe(II) (Umol/gr)
SFeOH3=u(2)*SA*SD*MW;
Fe_sorb=Ksorp*SFeOH3*u(4)/(h_plus+Ksorp*u(4));
%passivation of ferrihydrite due to sorption of Fe(II)
p=(SFeOH3-Fe_sorb)/SFeOH3;
R_retFe=.06*Ksorp*SFeOH3*h_plus/(h_plus+Ksorp*u(4))^2+1;


fO2=u(1)/(KO2+u(1));
fFeOH3=u(2)/(KFeOH3+u(2))*kinO2/(kinO2+u(1))*p;
fSO4=u(3)/(KSO4+u(3))*kinO2/(kinO2+u(1))*kinFeOH3/(kinFeOH3+u(2));
Sat_FeS=u(4)*u(5)/(KFeS*h_plus^2);

Rc=400*exp(-.1831*x);

R1=Rc*fO2*26.5; 
R2=Rc*fFeOH3;
R3=Rc*fSO4;
R4=ktsox*u(5)*u(1);
R5=kfeox*u(1)*u(4);
R6=ktsfe*u(2)*u(5);

if (Sat_FeS>=1)
    R7=0;
    R_7=kfepre*(Sat_FeS-1);
else
    R7=kfedis*u(6)*(1-Sat_FeS);
    R_7=0;
end

R_OMS=40/.06*u(5); 
%Goethite/lepidocrocite formation
Rformgoe=kFeOOH*u(4)*(u(2)>0);
%Magnetite respiration
Rrespmag=Rc*u(8)/(u(8)+Km_mag)*kinO2/(kinO2+u(1));
%magnetite formation
Rformag=kMag*u(4)*(u(2)>0);

% c=ones(1,12);
c=[ones(1,3),R_retFe,ones(1,4)];

% Transport
f = [(D_bio+D_O2)*DuDx(1)-w*u(1);...
    D_bio*DuDx(2)-w*u(2);...
    (D_bio+D_SO4)*DuDx(3)-w*u(3);...
    (D_bio+D_Fe)*DuDx(4)-w*u(4);...
    (D_bio+D_H2S)*DuDx(5)-w*u(5);...
    D_bio*DuDx(6)-w*u(6);...
    D_bio*DuDx(7)-w*u(7);...
    D_bio*DuDx(8)-w*u(8)];
 
% Reaction 
s = [(BC0_O2-u(1))*alfax-0.06*R1-.25*R5-2*R4;...
    -4*R2+R5/0.06-8*R6-(Rformgoe+2/3*Rformag)/.06;...
    (BC0_SO4-u(3))*alfax+0.06*(R6-.5*R3)+R4;...
    (BC0_Fe-u(4))*alfax+0.06*(4*R2+8*R6+R7-R_7+6*Rrespmag)-R5-1/3*Rformag-w*(R_retFe-1)*DuDx(4);...
    (BC0_H2S-u(5))*alfax+0.06*(.5*R3-R6+R7-R_7-R_OMS);...
    -R7+R_7;...
     Rformgoe/.06;...
    1/3*Rformag/.06-2*Rrespmag];
    
% --------------------------------------------------------------
function u0 = pdex1ic(x) 

% u0=zeros(1,12);
u0=[0,10,zeros(1,6)];
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
  
global BC0_FeOH3 BC0_O2 BC0_SO4 BC0_Fe BC0_H2S BC0_FeS w

% pl = [ul(1)-BC0_O2;...
%       BC0_FeOH3/.05961;...
%       ul(3)-BC0_SO4;...
%       ul(4)-BC0_Fe;...
%       ul(5)-BC0_H2S;...
%       BC0_FeS/.05961;...
%       BC0_FeOH3/.05961/100;...
%       BC0_FeOH3/.05961/100];
pl = [ul(1)-BC0_O2;...
      BC0_FeOH3/.05961;...
      ul(3)-BC0_SO4;...
      ul(4)-BC0_Fe;...
      ul(5)-BC0_H2S;...
      BC0_FeS/.05961;...
      0;...
      0];

%1 for solid, 0 for solute
ql = [0;1;0;0;0;1;1;1];

pr = [w*ur(1);...
      w*ur(2);...
      w*ur(3);...
      0;...
      w*ur(5);...
      w*ur(6);...
      w*ur(7);...
      w*ur(8)];
  
% qr = [1;1;1;1;1;1;1;1;1];
qr=ones(1,12);
%----------------------------------------------------------------