%rate_calculate calculates the rate of reactions as a function of x and t
%
%   The rate is calculated directly from the concentration at each x and t,
%   and the rate law for each reaction. The 'sol' matrix must be present
%   for this function to work. This can be done either by running the
%   model, or by importing a MATLAB .mat file that contains a 'sol'
%   matrix
%
%   To access the rate matrix, type 'global SimRates'

%%Rate_caltulate_external RMC June 30 2013

    global sol SimRates F;
    SimRates = cell(6,1);          %declare the global cell matrix to store the rate
    
%     if isempty(para)==1             %if the para is not already defined, read it from file
%         fid = fopen('para.in','r'); %reading the parameters from para.in
%         C=textscan(fid,'%s%f');     %names of parameters, values of parameters
%         para = C{2};                %we only need the values of parameters
%         fclose(fid);                %close file
%     end
        
 % Load the information saved in Result.mat
    load results.mat

    %half saturation coefs-----------------------------------------------------

%Half saturation for SO4 reduction
KSO4=0.5;      % igual a retraso
%Half saturation for MnO2 reduction
KMnO2=16;    % umol/g valor promig de Wang and Cappellen, 1996. Rang: 4-32.
%Half saturation for FeOH3 reduction
KFeOH3=200;    % 200umol/g valor de Canavan
%Half saturation for oxic respiration
KO2=0.010;      % calculat i igual a retraso 


%inhibition coefs----------------------------------------------------------

% Inhibition of oxic respiration
kinO2=1e-3;   % valor de Retraso i literatura
%Inhibition of Fe(III) reduction 
kinFeOH3=200;   %200 umol/g de Raoul

%Secondary reaction constants----------------------------------------------

kOM1=1.5;                           %degredation rate constant for OM1, Rc1 Canavan et al
kOM2=0.5;                           %degredation rate constant for OM2, Rc2 Canavan et al
ktsox = 3.15e4;                     %S(II) oxidation by O2, R5
kfeox = 4e4;                        %Fe(II) oxidation by O2, R6
ktsfe=2.5;                          %Rate of H2S oxidation by Fe(OH)3
kfedis=1e-3;                        %dissolution rate of FeS, R7
kfepre=17;                          %precipitation rate of FeS, R_7
kSorg=40;                           %Sulfidization of OM, R8
kpyr_pre=2.2e-12;                   %Pyrite precipitation, R10
kFeCO3dis=1e-3;                     %dissolution rate of FeCO3
kFeCO3pre=2e-12;                    %precipitation rate of FeCO3
kZnSdis=1e-3;                       %dissolution rate of ZnS
kZnSpre=2.0;                        %precipitation rate of ZnS
kCuSdis=1e-3;                       %dissolution rate of CuS
kCuSpre=8e-5;                       %precipitation rate of CuS
kCoSdis=1e-3;                       %dissolution rate of CoS
kCoSpre=5.0e-2;                     %precipitation rate of CoS
kNiSdis=1e-3;                       %dissolution rate of NiS
kNiSpre=4.0e-8;                     %precipitation rate of NiS

 %saturation constants-----
    KFeS = 10^3.04;
    KFeS2=10^-9.4;
    KFeCO3=10^-8.4;
    KZnS=10^0.5;
    KCuS=10^-4.3;
    KCoS=10^-1.4;
    KNiS=10^-7;
    
  %other parameters----------------------------------------------------------
    x = linspace(0,50,250);
    h_plus=0.0033;%(0.004*x*x-0.007*x+0.005)*(x<=1)+(6.0e-5*x+0.003)*(x>1);  
    %Porosity
    phi=-.008*x+.976; 
    Dphi=-.008; % phi derivada
    %Porosity correction
    porcoef=phi; %coefficient to modify the porous effects on diffusion coeffs
    %Dry density of the sediment
    rhob=2; 
    %Conversion from solid to aqueous through density
    F=rhob.*(1-phi)./phi; 

  %indices-------------------------------------------------------------------
  %note that u(n) has been replaced by sol(end,:,n), and dot operation has
  %been used. e.g. 'u(1)*u(2)' has been replaced by sol(end,:,1).*sol(end,:,2)
    
    fO2 = sol(end,:,2)./(KO2+sol(end,:,2));     
    fMnO2B = sol(end,:,34)./(KMnO2+sol(end,:,34)).*kinO2./(kinO2+sol(end,:,2));
    fFeOH3 = sol(end,:,3)./(KFeOH3+sol(end,:,3)).*kinO2./(kinO2+sol(end,:,2));
    fFeOH3B = sol(end,:,35)./(KFeOH3+sol(end,:,35)).*kinO2./(kinO2+sol(end,:,2));
    fSO4 = sol(end,:,4)./(KSO4+sol(end,:,4)).*kinO2./(kinO2+sol(end,:,2)).*...
        kinFeOH3./(kinFeOH3+sol(end,:,3));
    
    Sat_FeS = sol(end,:,5).*sol(end,:,6)./(KFeS*h_plus^2);
    Sat_FeS2 = sol(end,:,5).*sol(end,:,6)./(KFeS2*h_plus^2);
    Sat_FeCO3 = sol(end,:,5).*sol(end,:,13)./(KFeCO3*h_plus);
    Sat_ZnS = sol(end,:,14).*sol(end,:,6)./(KZnS*h_plus^2);
    Sat_CuS = sol(end,:,18).*sol(end,:,6)./(KCuS*h_plus^2);
    Sat_CoS = sol(end,:,22).*sol(end,:,6)./(KCoS*h_plus^2);
    Sat_NiS = sol(end,:,25).*sol(end,:,6)./(KNiS*h_plus^2);


        %reaction rates------------------------------------------------------------
    
R1=sol(end,:,1).*F.*fO2*kOM1*2.5;                       % R1: OM1 + O2 = CO2 + H2O
R_1=sol(end,:,27).*F.*fO2*kOM2*2.5;                     % R_1: OM2 + O2 = CO2 + H2O

R23_1=sol(end,:,1).*F.*fMnO2B*0.05*kOM1;                % R23: CH2O+2MnO2B(s)+3CO2+H2O=2Mn(2+) +2HCO3- 
R_23_1=sol(end,:,27).*F.*fMnO2B*kOM2*0.05;              % R_23: CH2O+2MnO2B(s)+3CO2+H2O=2Mn(2+) +2HCO3- 

R2=sol(end,:,1).*F.*fFeOH3*2.5*kOM1;                    % R2: OM1+4fe(OH)3(s)+7CO2=4Fe(II)+8HCO3- + 3H2O
R_2=sol(end,:,27).*F.*fFeOH3*2.5*kOM2;                  % R_2: OM2+4fe(OH)3(s)+7CO2=4Fe(II)+8HCO3- + 3H2O

R27=sol(end,:,1).*F.*fFeOH3B*0.1*kOM1;                  % R27: OM1+4fe(OH)3B(s)+7CO2=4Fe(II)+8HCO3- + 3H2O
R_27=sol(end,:,27).*F.*fFeOH3B*0.1*kOM2;                % R_27: OM2+4fe(OH)3B(s)+7CO2=4Fe(II)+8HCO3- + 3H2O

SO4accel=20; % Acceleration factor applied to SO4 reduction due to the bacteria 
R3=sol(end,:,1).*F.*fSO4*SO4accel*kOM1;                 % R3: OM1+0.5SO4=0.5H2S+HCO3
R_3=sol(end,:,27).*F.*fSO4*kOM2*SO4accel;               % R_3: OM2+0.5SO4=0.5H2S+HCO3

R4=ktsox.*sol(end,:,6).*sol(end,:,2);                    % R4: H2S+2O2+2HCO3=SO4+2CO2+2H2O

R5=kfeox.*sol(end,:,2).*sol(end,:,5);                    % R5: Fe(II)+0.25O2+2HCO3+1/2H2O=Fe(OH)3(s)+2CO2

R6=ktsfe.*sol(end,:,3).*F.*sol(end,:,6);                 % R6: H2S+ 14CO2+8Fe(OH)3(s)=8Fe(II)+SO4+14HCO3+6H2O

kmnox=5e3;                                           % R25: Mn(II)+0.5O2+2HCO3+2H2O=MnO2(s)+2CO2+H2O 
R25=kmnox.*sol(end,:,2).*sol(end,:,31);

kfemno2=2.4e3;                                       % R26: 2Fe(II)+ MnO2 + 2HCO3- + 2H2O = 2 Fe(OH)3 + Mn2+ + 2CO2
R26=kfemno2.*sol(end,:,5).*(sol(end,:,32)).*F; 

R10=(sol(end,:,1)+sol(end,:,27)).*F.*kSorg.*sol(end,:,6);  % R10: OM sulfidization

R7 = (Sat_FeS<1)*kfedis.*sol(end,:,7).*(1-Sat_FeS);    %dissolution rate of FeS
R_7 = (Sat_FeS>=1)*kfepre.*(Sat_FeS-1);              %precipitation rate of FeS

R11 = (Sat_FeS2>=1)*kpyr_pre.*(Sat_FeS2-1);          %precipitation rate of FeS2

R12 = (Sat_FeCO3<1)*kFeCO3dis.*sol(end,:,13).*(1-Sat_FeCO3);    %dissolution rate of FeCO3
R_12 = (Sat_FeCO3>=1)*kFeCO3pre.*(Sat_FeCO3-1);              %precipitation rate of FeCO3
    
% TOTAL OF OM DEGRADATED
R_OM=(R1+R_1+R23_1+R_23_1+R2+R_2+R27+R_27+R3+R_3);

% TOTAL OF O2 CONSUMED
R_O2=(R1+R_1).*F+2*R4+0.25*R5+0.5*R25;

% TOTAL OF Fe(OH)3 DEGRADATED
R_FeOH=4*(R2+R_2+R27+R_27).*F+R5+8*R6;

% TOTAL OF Fe2+ GENERATED
R_Fe=(R_7-R7+R11+R_12-R12).*F+R5+R26;

% TOTAL OF SO4 CONSUMED
R_SO4=(0.5*(R3+R_3)+R7).*F;

% TOTAL OF H2S GENERATED
R_H2S=(R_7-R7+2*R11).*F+R4+R6+R10;
    
    SimRates{1} = R_OM;        %storing the rates in the cell matrix                
    SimRates{2} = R_O2;
    SimRates{3} = R_FeOH;
    SimRates{4} = R_Fe;
    SimRates{5} = R_SO4;
    SimRates{6} = R_H2S;
     
   
    IntRates = zeros(6,1);
    x = linspace(0,50,250);
    for i = [1:6]
        IntRates(i) = trapz(x,SimRates{i}(end,:));
    end
  
    save SimRates
    fprintf(('R_OM total OM degradated %0.2f\n'), IntRates(1));
    fprintf(('R_O2 total of O2 consumed %0.2f\n'), IntRates(2));
    fprintf(('R_FeOH total of Fe(OH)3 consumed %0.2f\n'), IntRates(3));
    fprintf(('R_Fe Fe distribution %0.2f\n'), IntRates(4));
    fprintf(('R_SO4 total of SO4 consumed %0.2f\n'), IntRates(5));
    fprintf(('R_H2S H2S distribution %0.2f\n'), IntRates(6));
    
          