
function pdex15
tic;
clear all;
clc;
global u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17...
       u18 u19 u20 u21 u22 u23 u24 u25 u26 u27 u28 u29 u30 u31 u32 u33 u34 u35...
       fD fIrr fT x t 
       
m = 0;
x = linspace(0,50,250);
t = linspace(0,50,250);%300 for fT calculations 100 for baseline

VarNames = {'OM(s)', ... %u1
            'O2(aq)',... %u2
            'Fe(OH)3(s)', ... %u3
            'SO4(2-)(aq)', ... %u4
            'Fe(2+)(aq)', ... %u5
            'S(-II)(aq)',... %u6
            'FeS(s)',... %u7
            'FeS-As(s)',... %u8
            'As(III,V)(aq)',... %u9
            'sulfidized OM',... %u10
            'FeS2(s)',... %u11
            'FeCO3(s)',... %u12
            'HCO3(aq)',... %u13
            'Zn(2+)(aq)',... %u14
            'ZnS(s)',... %u15
            'Fe(OH)3-Zn(s)',... %u16
            'T-Zn(s)',... %u17
            'Cu(2+)(aq)',... %u18
            'CuS(s)',... %u19
            'OM-Cu(s)',... %u20
            'Fe(OH)3-As(s)',... %u21
            'Co(2+)(aq)',... %u22
            'CoS(s)',... %u23
            'T-Co(s)',...%u24
            'Ni(2+)(aq)',... %u25
            'NiS(s)',... %u26
            'T-Ni(s)',... %u27
            'OM2(s)',... %u28
            'Al(3+)(aq)',... %u29
            'Al(OH)3(s)',... %u30
            'Mn(2+)(aq)',... %u31
            'MnO2(s)',... %u32
            'MnCO3(s)',... %u33
            'MnO2B(s)',... %u34
            'Fe(OH)3B(s)'}; %u35
            

NumVars = int16(length(VarNames)) ;

ql = [1;0;1;0;0;0;1;1;0;1;1;1;0;0;1;1;1;0;1;1;1;0;1;1;0;1;1;1;0;1;0;1;1;1;1];
NumPhases = int16(size(ql, 1));

%Checking input
if(NumPhases ~= NumVars)
    disp('size input does not fit');
    stop;
end

SimValues = cell(NumVars, 1);

disp('solving PDE ');
sol = pdepe(m,@pdex14pde,@pdex1ic,@pdex1bc,x,t);
% % Extract the first solution component as u.

for j=1:NumVars,
     MatValues = sol(:,:,j);
     [m, n]= size(MatValues);
     SimValues{j} = MatValues;
end

save ('Results.mat', 'SimValues', 'sol');

toc;

tic;

% Write results in excel
WriteSimValuesInExcel('simulation_results.xls', VarNames, SimValues, ...
                        NumVars, t, x );

                    
%per temps x1
timeStep = 250 ;
%DataValuesX1 = GetFieldVariables('C:\Users\Ester Torres\Desktop\Georgia Institute of Technology\Conceptual_model_RA\codes\working_codes\DEFINITIU\RA_2012_33\dades_camp_RA.xls', VarNames, NumVars );
%PlotSimVsData(VarNames, DataValuesX1, SimValues, NumVars, ...
%
%per temps x2
% time = 142.0 ;
% DataValuesX2 = GetFieldVariables('dades_campX2.xls', VarNames, NumVars );
% plotSimVsData(DataValuesX2, SimValues, NumVars, time );

toc;


% --------------------------------------------------------------
function [c,f,s] = pdex14pde(x,t,u,DuDx)

global BC0_OM BC0_FeOH3 BC0_FeOH3B BC0_O2 BC0_SO4 BC0_Fe BC0_H2S BC0_FeS BC0_FeS2...
       BC0_AsFeS BC0_AsO4 BC0_OMS BCO_FeCO3 BCO_HCO3...
       BC0_Zn BC0_ZnS BC0_FeOH3Zn BC0_TZn BC0_Cu BC0_Cu2S BC0_OMCu BC0_FeOH3As... 
       BC0_Co BC0_CoS BC0_TCo BC0_Ni BC0_NiS BC0_TNi BC0_OM2 BC0_Al BC0_AlOH3 BC0_Mn BC0_MnO2 BC0_MnCO3  BC0_MnO2B...
       D_bio D_O2 D_SO4 D_Fe D_H2S D_AsO4  F w phi Dphi v


%Boundary conditions 
BC0_OM=900; %900%580 calculat a partir del Corg a la primera llesca
BC0_OM2=200;
BC0_O2=.32; % microsensors
BC0_FeOH3=140;%240+120*sin(1.6*t);%240;
BC0_FeOH3B=100;
BC0_SO4=1.4;% resultats de l'aigua de fons de RA_oct_2010
BC0_Fe=4.28e-3; % dades del RA columna aigua, setembre 2009
BC0_H2S=0; % 
BC0_FeS=0;
BC0_AsO4=1.87e-5; % dades del RA columna aigua, setembre 2009
BC0_AsFeS=0; 
BC0_OMS=0;
BC0_FeS2=0;
BCO_FeCO3=0;
BCO_HCO3=0;
BC0_Zn=2.55e-2;%3.11e-2; % dades del RA columna aigua, setembre 2009
BC0_ZnS=0;
BC0_FeOH3Zn=0;
BC0_TZn=0;
BC0_Cu=8.25e-3; % dades del RA columna aigua, setembre 2009
BC0_Cu2S=0;
BC0_OMCu=0;
BC0_FeOH3As=0; % BC0_FeOH3*1.1e-3; 
BC0_Co=1.13e-3;%1.88e-3; % dades del RA columna aigua, setembre 2009
BC0_CoS=0;
BC0_TCo=0;
BC0_Ni=7.99e-4; % dades del RA columna aigua, setembre 2009
BC0_TNi=0.1;
BC0_NiS=0;
BC0_Al=5.67e-2;%1.65e-1; % dades del RA columna aigua, setembre 2009
BC0_AlOH3=0;
BC0_Mn=2.87e-2;  % dades del RA columna aigua, setembre 2009
BC0_MnO2=0; 
BC0_MnO2B=7;  
BC0_MnCO3=0;

% ALL UNITS of concentrations are uMOL and of volumes are cm-3 (cubic centimeter porewater !) (uMOL g-1 (gram solid phase) 

%Porosity
phi=-.008*x+.976; %depth dependant phi
Dphi=-.008; % phi derivada
%Porosity correction
porcoef=phi; %coefficient to modify the porous effects on diffusion coeffs

% Dry density of the sediment
rhob=2; 

% Conversion from solid to aqueous through density
F=rhob.*(1-phi)./phi;  %F=rhob*(1-phi)/phi;

% pH function
h_plus=(0.004*x*x-0.007*x+0.005)*(x<=1)+(6.0e-5*x+0.003)*(x>1); % pH is currently imposed in the model through concentration of H+ in uMOL cm-3 ! Can be step function or defined as a function of depth (x)

% All the diffusion coeffcients are modified for the porosity
D_O2=562*porcoef;  
D_SO4=261*porcoef;
D_Fe=176*porcoef;
D_Mn=176*porcoef;
D_Al=109*porcoef;
D_H2S=490*porcoef;
D_AsO4=500*porcoef; %inventat
D_Zn=200*porcoef; %inventat
D_Cu=180*porcoef; %inventat
D_Co=180*porcoef; %inventat
D_Ni=180*porcoef; %inventat
D_HCO3=2*D_Fe; %Li and Gregory no el donen a aquesta temperatura, és el valor de 25ºC

%Half saturation for SO4 reduction
KSO4=0.5;      % igual a retraso
%Half saturation for MnO2 reduction
KMnO2=16;    % umol/g valor promig de Wang and Cappellen, 1996. Rang: 4-32.
%Half saturation for FeOH3 reduction
KFeOH3=200;    % 200umol/g valor de Canavan
%Half saturation for oxic respiration
KO2=0.010;      % calculat i igual a retraso 
% Inhibition of oxic respiration
kinO2=1e-3;   % valor de Retraso i literatura
%Inhibition of Fe(III) reduction 
kinFeOH3=200;   %200 umol/g de Raoul
% Rate constant for sulfide oxidation by O2 
ktsox=3.15e4;      %Range 1e2 - 1.6e6 validated 3.15e4, as Retraso
% Rate constant for Fe(II) oxidation by O2 
kfeox=4e4;      %4e4
% Rate constant for sulfide oxidation by FeOH 
ktsfe=2.5;      %Range 2.5-95 validated 2.5

% FeS precipitation and disolution
KFeS=10^3.04;%9.6e3;   
kfedis=1e-3;    
kfepre=17;
Sat_FeS=u(5)*u(6)/(KFeS*h_plus^2); % saturation index

% FeS2 precipitation FeS2+ 2h+ = Fe2+ + 2H2S + S(0)
KFeS2=10^-9.4;
kpyr_pre=2.2e-12;
Sat_FeS2=u(5)*u(6)/(KFeS2*h_plus^2);

% MnCO3 precipitation MnCO3(s) + CO2 + H20 = Mn2+ + 2HCO3-
KMnCO3=10^-8.5;
kMnCO3dis=1e-3; %fitted
kMnCO3pre=1.3e-10; %fitted
Sat_MnCO3=u(31)*u(13)*u(13)/(KMnCO3);

% As sorption to Fe(OH3) 
kAsO4_ads=1; %Range 2.9-4.7 validated 1.35

% As sorption to FeS
kAs_FeS=0.3;      % 95 Canavan, 10 Raoul el nostre mes baix perquè tenim molts metalls competint?

% Rate constant for sulfidization of OM
kSorg=40;%4e-2;

% Siderite precipitation and disolution
KFeCO3=10^-8.4; %(Van Cappellen and Wang, 1996) 
kFeCO3dis=1e-3; %fitted
kFeCO3pre=2e-12; %fitted
Sat_FeCO3=u(13)*u(5)/(KFeCO3*h_plus); % saturation index

%Sphalerite precipitation and disolution
KZnS=10^0.5; % Retraso ZnS depth (amb H2S)
kZnSdis=1e-3; %fitted
kZnSpre=2.0; %fitted
Sat_ZnS=u(14)*u(6)/(KZnS*h_plus^2); % saturation index

%Zn adsorption to Fe(OH)3(s)
kZnFeOH=0.02; 

%Zn adsorption to OM(s). Actually is a sink of Zn to calculate the total Zn
kZnOM=0.1; % fitted

%CuS precipitation and disolution
KCuS=10^-4.3; % from Retraso
kCuSdis=1e-3; %fitted
kCuSpre=8e-5;%0.25e-5; %fitted
Sat_CuS=u(18)*u(6)/(KCuS*h_plus^2); % saturation index

%Cu adsorption to OM(s)
kOMCu=0.1; % fitted

%CoS precipitation and disolution
KCoS=10^-1.4; % from Retraso
kCoSdis=1e-3; %fitted
kCoSpre=5.0e-2; %fitted
Sat_CoS=u(22)*u(6)/(KCoS*h_plus^2); % saturation index

%NiS precipitation and disolution
KNiS=10^-7; % from Retraso
kNiSdis=1e-3; %fitted
kNiSpre=4.0e-8; %fitted
Sat_NiS=u(25)*u(6)/(KNiS*h_plus^2); % saturation index

%Al(OH)3 precipitation and disolution
KAlOH3=10^-10.8; % from Retraso
kAlOH3dis=1e-3; %fitted
kAlOH3pre=5e-15; %fitted
Sat_AlOH3=u(29)/(KAlOH3*h_plus^3); % saturation index

%Ni adsorption to other phases to fit Ni(aq)
kNisink=0.5; % inventat

% sedimentation rate
w=1.6;
% advection
v=w;

% bioirrigatoin coefficient
alfa0=2.0; 	
alfax=alfa0*exp(-.25*x);% formula for biorrigatoin decrease
% bioturbation coefficient
D_bio=5; %0.37; fitted

% Fraction of OM degradated by O2
fO2=u(2)/(KO2+u(2));

% Fraction of OM degradated by MnO2
fMnO2=u(32)/(KMnO2+u(32))*kinO2/(kinO2+u(2));
fMnO2B=u(34)/(KMnO2+u(34))*kinO2/(kinO2+u(2));


% Fraction of OM degradated by FeOH3
fFeOH3=u(3)/(KFeOH3+u(3))*kinO2/(kinO2+u(2));
fFeOH3B=u(35)/(KFeOH3+u(35))*kinO2/(kinO2+u(2));

% Fraction of OM degradated by SO4
fSO4=u(4)/(KSO4+u(4))*(kinO2/(kinO2+u(2))*kinFeOH3/(kinFeOH3+u(3)));

% ration of trace-metal to carrier phase, such as Fe oxides
% As sorbed to Fe(OH)3
fAsFe=1.0e-2;%1.0e-2; 

% ration of trace-metal to carrier phase, such as OM. Cu sorbed to OM
fCu_OM=0; %4.4e-3; trec la font de Cu de la OM


% ORGANIC MATTER RATE yr-1
kOM1=1.5;
kOM2=0.5;
% kOM3=0;

% R1: CH2O + O2 = CO2 + H2O
R1=u(1)*fO2*kOM1*2.5;% *35; % Range 17-34
R_1=u(27)*kOM2*fO2*2.5;

% R23: CH2O+2MnO2(s)+3CO2+H2O=2Mn(2+) +2HCO3- 
R23=u(1)*fMnO2*kOM1;
R_23=u(27)*kOM2*fMnO2;

R23_1=u(1)*fMnO2B*0.05*kOM1;
R_23_1=u(27)*kOM2*0.05*fMnO2B;

% R2: CH2O+4fe(OH)3(s)+7CO2=4Fe(II)+8HCO3- + 3H2O
R2=2.5*u(1)*fFeOH3*kOM1;
R_2=2.5*u(27)*kOM2*fFeOH3;

R27=u(1)*fFeOH3B*0.1*kOM1;
R_27=u(27)*fFeOH3B*0.1*kOM2;

% Acceleration factor applied to SO4 reduction due to the bacteria 
SO4accel=20; % fitted
% R3: 2CH2O+SO4=H2S+2HCO3
R3=u(1)*fSO4*SO4accel*kOM1; 
R_3=u(27)*kOM2*fSO4*SO4accel;

% R4: H2S+2O2+2HCO3=SO4+2CO2+2H2O
R4=ktsox*u(6)*u(2);

% R5: Fe(II)+0.25O2+2HCO3+1/2H2O=Fe(OH)3(s)+2CO2
R5=kfeox*u(2)*u(5);

% R25: Mn(II)+0.5O2+2HCO3+2H2O=MnO2(s)+2CO2+H2O 
kmnox=5e3; % Van Cappellen and Wang 1996. cm3/umol/yr
R25=kmnox*u(2)*u(31);

% R26: 2Fe(II)+ MnO2 + 2HCO3- + 2H2O = 2 Fe(OH)3 + Mn2+ + 2CO2
kfemno2=2.4e3;%3e3; % Van Cappellen and Wang 1996. cm3/umol/yr
R26=kfemno2*u(5)*(u(32)); %falta multiplicar o dividir per F!

%OXIDACIÓ DE SULFURS
% RFeS: FeS + 2O2 = SO42- + Fe2+
kfesox=2e3; % interval of Canavan et al., 2006
RFeS=kfesox*u(2)*u(7)*F;

% RZnS: ZnS + 2O2 = SO42- + Zn2+
kznsox=2e2;%2e-1; %Canavan et al., 2006
RZnS=kznsox*u(2)*u(15)*F;

% RCuS: CuS + 2O2 = SO42- + Cu2+
kcusox=1e3;%2e0; %Canavan et al., 2006
RCuS=kcusox*u(2)*u(19)*F;

% RCoS: CoS + 2O2 = SO42- + Co2+
kcosox=3e4;%2e0; %Canavan et al., 2006
RCoS=kcosox*u(2)*u(23)*F;

% RNiS: NiS + 2O2 = SO42- + Ni2+
knisox=2e2;%2e1; %Canavan et al., 2006
RNiS=knisox*u(2)*u(26)*F;


% R6: H2S+ 14CO2+8Fe(OH)3(s)=8Fe(II)+SO4+14HCO3+6H2O
R6=ktsfe*u(3)*F*u(6);

% FeS precipitation/dissolution: FeS + 2H+ = Fe(II)+ H2S
if (Sat_FeS>=1)
    R7=0;
    R_7=kfepre*(Sat_FeS-1);
else
    R7=kfedis*u(7)*(1-Sat_FeS);
    R_7=0;
end

% R8: As adsorption to FeOH3: As+Fe(OH)3(s)='Fe(OH)3(s)-As'
R8=kAsO4_ads*u(3)*F*u(9);

% R9: As adsorption to FeS: As+FeS='FeS-As'
R9=kAs_FeS*u(7)*u(9)*F;

% R10: OM sulfidization: OM+H2S='OM-H2S'
R10=(u(1)+u(27))*kSorg*u(6);%35000*kSorg*u(6);%(u(1)+u(27))*kSorg/F*u(6); %fitted. OM=10000

% R11: pyrite precipitation: FeS2 + 2H+ = Fe2+ + H2S + S(0)
% R11=kpyr_pre*u(6)*u(7);
if (Sat_FeS2>=1)
    R11=kpyr_pre*(Sat_FeS2-1);
else
    R11=0;
end

% FeCO3 precipitation/dissolution Fe(II)+HCO3-= FeCO3 + H+
if (Sat_FeCO3>=1)
    R12=0;
    R_12=kFeCO3pre*(Sat_FeCO3-1);
else
    R12=kFeCO3dis*u(13)*(1-Sat_FeCO3);
    R_12=0;
end

% ZnS precipitation/dissolution  ZnS + 2H+ = Zn(II)+ H2S
if (Sat_ZnS>=1)
    R13=0;
    R_13=kZnSpre*(Sat_ZnS-1);
else
    R13=kZnSdis*u(15)*(1-Sat_ZnS); %u(14) hauria de ser u(15)
    R_13=0;
end

% R14 Zn sorption to Fe(OH)3. Zn(2+)+ Fe(OH)3(s)='Fe(OH)3(s)-Zn+' + H+
R14=kZnFeOH*u(3)*u(14)/F;

% R15 TZn - ZnS - ZnFe(OH)3. Ho calculo a partir de la OM, però no és el
% lligat a la OM. És només per fer un sumidero on posar tota la resta de Zn
R15=u(14)/F*(u(1)+u(27))*kZnOM;

% R16 CuS precipitation/dissolution CuS + 2H+ = Cu2+ + H2S
if (Sat_CuS>=1)
    R16=0;
    R_16=kCuSpre*(Sat_CuS-1);
else
    R16=kCuSdis*u(19)*(1-Sat_CuS); 
    R_16=0;
end

% R17 Cu sorption to OM
R17=(u(1)+u(27))*u(18)/F*kOMCu; %


% R18 CoS precipitation/dissolution CoS + 2H+ = Co2+ + H2S
if (Sat_CoS>=1)
    R18=0;
    R_18=kCoSpre*(Sat_CoS-1);
else
    R18=kCoSdis*u(23)*(1-Sat_CoS); 
    R_18=0;
end


% R19: Co sink for all the phases except CoS:
R19=0.035*u(3)*u(22)/F; % to fit Co(aq)


% R20 NiS precipitation/dissolution NiS + 2H+ = Ni2+ + H2S
if (Sat_NiS>=1)
    R20=0;
    R_20=kNiSpre*(Sat_NiS-1);
else
    R20=kNiSdis*u(26)*(1-Sat_NiS); 
    R_20=0;
end

% R21: Ni sink
R21=kNisink*u(25)/F;

% R22 Al3+ + 3H2O = Al(OH)3 + 3H+
if (Sat_AlOH3>=1)
    R22=0;
    R_22=kAlOH3pre*(Sat_AlOH3-1);
else
    R22=kAlOH3dis*u(30)*(1-Sat_AlOH3); 
    R_22=0;
end

% R24 MnCO3 precipitation/dissolution: MnCO3(s) + CO2 + H20 = Mn2+ + 2HCO3-
if (Sat_MnCO3>=1)
    R24=0;
    R_24=kMnCO3pre*(Sat_MnCO3-1);
else
    R24=kMnCO3dis*u(33)*(1-Sat_MnCO3);
    R_24=0;
end
c = [1-phi;...
      phi;...
      1-phi;...
      phi;...
      phi;...
      phi;...
      1-phi;...
      1-phi;...
      phi;...
      1-phi;...
      1-phi;...
      1-phi;...
      phi;...
      phi;...
      1-phi;...
      1-phi;...
      1-phi;...
      phi;...
      1-phi;...
      1-phi;...
      1-phi;...
      phi;...
      1-phi;...
      1-phi;...
      phi;...
      1-phi;...
      1-phi;...
      1-phi;...
      phi;...
      1-phi;...
      phi;...
      1-phi;...
      1-phi;...
      1-phi;...
      1-phi];

% Transport
f = [-D_bio*Dphi*u(1)+D_bio*(1-phi)*DuDx(1)-w*(1-phi)*u(1);...
    D_O2*phi*DuDx(2)+D_bio*Dphi*u(2)+D_bio*phi*DuDx(2)-v*phi*u(2);...
    -D_bio*Dphi*u(3)+D_bio*(1-phi)*DuDx(3)-w*(1-phi)*u(3);...
    D_SO4*phi*DuDx(4)+D_bio*Dphi*u(4)+D_bio*phi*DuDx(4)-v*phi*u(4);...
    D_Fe*phi*DuDx(5)+D_bio*Dphi*u(5)+D_bio*phi*DuDx(5)-v*phi*u(5);...
    D_H2S*phi*DuDx(6)+D_bio*Dphi*u(6)+D_bio*phi*DuDx(6)-v*phi*u(6);...
    -D_bio*Dphi*u(7)+D_bio*(1-phi)*DuDx(7)-w*(1-phi)*u(7);...
    -D_bio*Dphi*u(8)+D_bio*(1-phi)*DuDx(8)-w*(1-phi)*u(8);...
    D_AsO4*phi*DuDx(9)+D_bio*Dphi*u(9)+D_bio*phi*DuDx(9)-v*phi*u(9);...
    -D_bio*Dphi*u(10)+D_bio*(1-phi)*DuDx(10)-w*(1-phi)*u(10);...
    -D_bio*Dphi*u(11)+D_bio*(1-phi)*DuDx(11)-w*(1-phi)*u(11);...
    -D_bio*Dphi*u(12)+D_bio*(1-phi)*DuDx(12)-w*(1-phi)*u(12);...
    D_HCO3*phi*DuDx(13)+D_bio*Dphi*u(13)+D_bio*phi*DuDx(13)-v*phi*u(13);...
    D_Zn*phi*DuDx(14)+D_bio*Dphi*u(14)+D_bio*phi*DuDx(14)-v*phi*u(14);...
    -D_bio*Dphi*u(15)+D_bio*(1-phi)*DuDx(15)-w*(1-phi)*u(15);...
    -D_bio*Dphi*u(16)+D_bio*(1-phi)*DuDx(16)-w*(1-phi)*u(16);...
    -D_bio*Dphi*u(17)+D_bio*(1-phi)*DuDx(17)-w*(1-phi)*u(17);...
    D_Cu*phi*DuDx(18)+D_bio*Dphi*u(18)+D_bio*phi*DuDx(18)-v*phi*u(18);...
    -D_bio*Dphi*u(19)+D_bio*(1-phi)*DuDx(19)-w*(1-phi)*u(19);...
    -D_bio*Dphi*u(20)+D_bio*(1-phi)*DuDx(20)-w*(1-phi)*u(20);...
    -D_bio*Dphi*u(21)+D_bio*(1-phi)*DuDx(21)-w*(1-phi)*u(21);...
    D_Co*phi*DuDx(22)+D_bio*Dphi*u(22)+D_bio*phi*DuDx(22)-v*phi*u(22);...
    -D_bio*Dphi*u(23)+D_bio*(1-phi)*DuDx(23)-w*(1-phi)*u(23);...
    -D_bio*Dphi*u(24)+D_bio*(1-phi)*DuDx(24)-w*(1-phi)*u(24);...
    D_Ni*phi*DuDx(25)+D_bio*Dphi*u(25)+D_bio*phi*DuDx(25)-v*phi*u(25);...
    -D_bio*Dphi*u(26)+D_bio*(1-phi)*DuDx(26)-w*(1-phi)*u(26);...
    -D_bio*Dphi*u(27)+D_bio*(1-phi)*DuDx(27)-w*(1-phi)*u(27);...
    -D_bio*Dphi*u(28)+D_bio*(1-phi)*DuDx(28)-w*(1-phi)*u(28);...
    D_Al*phi*DuDx(29)+D_bio*Dphi*u(29)+D_bio*phi*DuDx(29)-v*phi*u(29);...
    -D_bio*Dphi*u(30)+D_bio*(1-phi)*DuDx(30)-w*(1-phi)*u(30)
     D_Mn*phi*DuDx(31)+D_bio*Dphi*u(31)+D_bio*phi*DuDx(31)-v*phi*u(31);...
    -D_bio*Dphi*u(32)+D_bio*(1-phi)*DuDx(32)-w*(1-phi)*u(32);...
    -D_bio*Dphi*u(33)+D_bio*(1-phi)*DuDx(33)-w*(1-phi)*u(33);...
    -D_bio*Dphi*u(34)+D_bio*(1-phi)*DuDx(34)-w*(1-phi)*u(34);...
    -D_bio*Dphi*u(35)+D_bio*(1-phi)*DuDx(35)-w*(1-phi)*u(35)];
     
% Reaction 
s = [ (-R1-R2-R3)*(1-phi);...
    ((BC0_O2-u(2))*alfax-F*(R1+R_1)-.25*R5-0.5*R25-2*R4-2*RFeS-2*RZnS-2*RCuS-2*RCoS-2*RNiS)*phi;...
    (-4*(R2+R_2)+R5/F-8*R6/F+2*R26/F)*(1-phi);...
    ((BC0_SO4-u(4))*alfax+F*(-0.5*(R3+R_3))+R6+R4+RFeS+RZnS+RCuS+RCoS+RNiS)*phi;...
    ((BC0_Fe-u(5))*alfax+F*(4*R2+4*R_2+4*R27+4*R_27+R7-R_7-R11+R12-R_12)+8*R6-R5+RFeS-2*R26)*phi;...
    ((BC0_H2S-u(6))*alfax+F*(0.5*(R3+R_3)+R7-R_7-R11-R10)-R4-R6)*phi;...
    (-R7+R_7-RFeS/F)*(1-phi);...
    (R9/F)*(1-phi);...
    ((BC0_AsO4-u(9))*alfax+F*fAsFe*(4*R2+4*R_2+4*R27+4*R_27)+fAsFe*8*R6-R8-R9)*phi;...
    R10*(1-phi);...
    R11*(1-phi);...
    (-R12+R_12)*(1-phi);...
    ((BCO_HCO3-u(13))*alfax+F*(8*R2+8*R_2+8*R27+8*R_27+R3+R_3+4*R23+4*R23_1+4*R_23+4*R_23_1+2*R24+2*R_24+2*R12-2*R_12)+14*R6-2*R5-2*R4-2*R25-2*R26)*phi;...
    ((BC0_Zn-u(14))*alfax+F*(R13-R_13-R14-R15)+RZnS)*phi;...
    (-R13+R_13-RZnS/F)*(1-phi);...
    (RZnS/F)*(1-phi);...
    (R15+R14-R13+R_13)*(1-phi);...
    ((BC0_Cu-u(18))*alfax+F*fCu_OM*(R1+R_1+R2+R_2+R3+R_3)+F*(R16-R_16-R17)+RCuS)*phi;...
    (-R16+R_16-RCuS/F)*(1-phi);...
    R17*(1-phi);...
    R8*(1-phi);...
    ((BC0_Co-u(22))*alfax+F*(R18-R_18-R19)+RCoS)*phi;...
    (-R18+R_18-RCoS/F)*(1-phi);...
    (R19-R18+R_18)*(1-phi);...
    ((BC0_Ni-u(25))*alfax+F*(-R_20+R20-R21)+RNiS)*phi;...
    (-R20+R_20-RNiS/F)*(1-phi);...
    (R21-R20+R_20)*(1-phi);...
    (-R_1-R_2-R_3)*(1-phi);...
    ((BC0_Al-u(29))*alfax+F*(-R_22+R22))*phi;...
    (-R22+R_22)*(1-phi)
    ((BC0_Mn-u(31))*alfax+F*(2*R23+2*R_23+2*R23_1+2*R_23_1+R24-R_24)-R25+R26)*phi;...
    (-2*(R23+R_23)-R26/F+R25/F)*(1-phi);...
    (-R24+R_24)*(1-phi);...
    (-2*(R23_1+R_23_1))*(1-phi);...
    (-4*(R27+R_27))*(1-phi)];
       
                        
    
% --------------------------------------------------------------
function u0 = pdex1ic(x) 

u0=zeros(1,35);
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
  
global BC0_OM BC0_FeOH3 BC0_FeOH3B BC0_O2 BC0_SO4 BC0_Fe BC0_H2S BC0_FeS BC0_FeS2...
       BCO_FeCO3 BCO_HCO3 BC0_Zn BC0_ZnS BC0_FeOH3Zn BC0_TZn...
       BC0_AsFeS BC0_AsO4 BC0_OMS BC0_Cu BC0_Cu2S BC0_OMCu BC0_FeOH3As...
       BC0_Co BC0_CoS BC0_TCo BC0_Ni BC0_NiS BC0_TNi BC0_OM2 BC0_Al BC0_AlOH3...
       BC0_Mn BC0_MnO2 BC0_MnCO3 BC0_MnO2B w phi Dphi v  D_bio  F         
rhob=2;
pl = [BC0_OM/rhob;...
      ul(2)-BC0_O2;...
      BC0_FeOH3/rhob;...
      ul(4)-BC0_SO4;...
      ul(5)-BC0_Fe;...
      ul(6)-BC0_H2S;...
      BC0_FeS/rhob;...
      BC0_AsFeS/rhob;...
      ul(9)-BC0_AsO4;...
      BC0_OMS/rhob;...
      BC0_FeS2/rhob;...
      BCO_FeCO3/rhob;...
      ul(13)-BCO_HCO3;...
      ul(14)-BC0_Zn;...
      BC0_ZnS/rhob;...
      BC0_FeOH3Zn/rhob;...
      BC0_TZn/rhob;...
      ul(18)-BC0_Cu;...
      BC0_Cu2S/rhob;...
      BC0_OMCu/rhob;...
      BC0_FeOH3As/rhob;...
      ul(22)-BC0_Co;...
      BC0_CoS/rhob;...
      BC0_TCo/rhob;...
      ul(25)-BC0_Ni;...
      BC0_NiS/rhob;...
      BC0_TNi/rhob;... 
      BC0_OM2/rhob;...
      ul(29)-BC0_Al;...
      BC0_AlOH3/rhob
      ul(31)-BC0_Mn;...
      BC0_MnO2/rhob;...
      BC0_MnCO3/rhob;...
      BC0_MnO2B/rhob;...
      BC0_FeOH3B/rhob];
       

%1 for solid, 0 for solute
ql = [1;0;1;0;0;0;1;1;0;1;1;1;0;0;1;1;1;0;1;1;1;0;1;1;0;1;1;1;0;1;0;1;1;1;1];

pr = [w*(1-phi)*ur(1)+D_bio*Dphi*ur(1);...
      v*phi*ur(2)-D_bio*Dphi*ur(2);...
      w*(1-phi)*ur(3)+D_bio*Dphi*ur(3);...
      v*phi*ur(4)-D_bio*Dphi*ur(4);...
      v*phi*ur(5)-D_bio*Dphi*ur(5);...
      v*phi*ur(6)-D_bio*Dphi*ur(6);...
      w*(1-phi)*ur(7)+D_bio*Dphi*ur(7);...
      w*(1-phi)*ur(8)+D_bio*Dphi*ur(8);...
      v*phi*ur(9)-D_bio*Dphi*ur(9);...
      w*(1-phi)*ur(10)+D_bio*Dphi*ur(10);...
      w*(1-phi)*ur(11)+D_bio*Dphi*ur(11);...
      w*(1-phi)*ur(12)+D_bio*Dphi*ur(12);...
      v*phi*ur(13)-D_bio*Dphi*ur(13);...
      v*phi*ur(14)-D_bio*Dphi*ur(14);...
      w*(1-phi)*ur(15)+D_bio*Dphi*ur(15);...
      w*(1-phi)*ur(16)+D_bio*Dphi*ur(16);...
      w*(1-phi)*ur(17)+D_bio*Dphi*ur(17);...
      v*phi*ur(18)-D_bio*Dphi*ur(18);...
      w*(1-phi)*ur(19)+D_bio*Dphi*ur(19);...
      w*(1-phi)*ur(20)+D_bio*Dphi*ur(20);...
      w*(1-phi)*ur(21)+D_bio*Dphi*ur(21);...
      v*phi*ur(22)-D_bio*Dphi*ur(22);...
      w*(1-phi)*ur(23)+D_bio*Dphi*ur(23);...
      w*(1-phi)*ur(24)+D_bio*Dphi*ur(24);...
      v*phi*ur(25)-D_bio*Dphi*ur(25);...
      w*(1-phi)*ur(26)+D_bio*Dphi*ur(26);...
      w*(1-phi)*ur(27)+D_bio*Dphi*ur(27);...
      w*(1-phi)*ur(28)+D_bio*Dphi*ur(28);...
      v*phi*ur(29)-D_bio*Dphi*ur(29);...
      w*(1-phi)*ur(30)+D_bio*Dphi*ur(30);...
      v*phi*ur(31)-D_bio*Dphi*ur(31);...
      w*(1-phi)*ur(32)+D_bio*Dphi*ur(32);...
      w*(1-phi)*ur(33)+D_bio*Dphi*ur(33);...
      w*(1-phi)*ur(34)+D_bio*Dphi*ur(34);...
      w*(1-phi)*ur(35)+D_bio*Dphi*ur(35)];
      
      
qr=ones(1,35);
%----------------------------------------------------------------
% Last updated : 5/04/2013 

 

  
   
