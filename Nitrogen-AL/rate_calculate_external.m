%rate_calculate calculates the rate of reactions as a function of x and t
%
%   The rate is calculated directly from the concentration at each x and t,
%   and the rate law for each reaction. The 'sol' matrix must be present
%   for this function to work. This can be done either by running the
%   model, or by importing a MATLAB .mat file that contains a 'sol'
%   matrix
%
%   To access the rate matrix, type 'global SimRates'

%%Rate_caltulate_external RMC Jne 30 2013

    global sol SimRates para F;
    SimRates = cell(3,1);          %declare the global cell matrix to store the rate
    
    if isempty(para)==1             %if the para is not already defined, read it from file
        fid = fopen('para.in','r'); %reading the parameters from para.in
        C=textscan(fid,'%s%f');     %names of parameters, values of parameters
        para = C{2};                %we only need the values of parameters
        fclose(fid);                %close file
    end
        
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

%Secondary reaction constants----------------------------------------------
kOM1 = 1;                           %degredation rate constant for OM1, Rc1 Canavan et al
kOM2 = 0.01;                        %degredation rate constant for OM2, Rc2 Canavan et al
kfeox = para(19);                   %Fe(II) oxidation by O2, R7
knh4ox = 20; %max ammonium oxidation rate AML 
kno2ox = 20; %max nitrite oxidation rate AML 
kamx = 20; %max anammox rate AML 
 
    %indices-------------------------------------------------------------------
    %note that u(n) has been replaced by sol(:,:,n), and dot operation has
    %been used. e.g. 'u(1)*u(2)' has been replaced by sol(:,:,1).*sol(:,:,2)
    
fO2 = sol(:,:,1)./(KmO2+sol(:,:,1));
fNO3 = sol(:,:,12)./(KmNO3+sol(:,:,12)).*kinO2./(kinO2+sol(:,:,1));%nitrate reduction to nitrite
fNO2 = sol(:,:,14)./(KmNO2+sol(:,:,14)).*kinO2./(kinO2+sol(:,:,1)); %nitrite reduction to N2O AML
fNO2DNRA = sol(:,:,14)./(KmNO2DNRA+sol(:,:,14)).*kinO2./(kinO2+sol(:,:,1)); %nitrite reduction to ammonium AML
fN2O = sol(:,:,15)./(KmN2O+sol(:,:,15)).*kinO2./(kinO2+sol(:,:,1)) ; %nitrous oxide reduction to N2 AML
fFeOH3 = sol(:,:,2)./(KmFeOH3+sol(:,:,2)).*kinO2./(kinO2+sol(:,:,1)).*kinNO3./(kinNO3+sol(:,:,12));
%fSO4 = sol(:,:,4)./(KmSO4+sol(:,:,4)).*kinO2/(kinO2+sol(:,:,4)).*kinNO3./(kinNO3+sol(:,:,12)).*kinFeOH3./(kinFeOH3+sol(:,:,12));
%Sat_FeS = 0;% u(5)*u(6)/(KFeS*h_plus^2);    %saturation index of FeS    

        %reaction rates------------------------------------------------------------
    Rc1 = kOM1.*sol(:,:,9);        %first pool rate
    Rc2 = kOM2.*sol(:,:,10);         %second pool rate
     
    R2a = Rc1.*fNO3*25;             %OM oxiation by Fe(OH)3
    R2b = Rc2.*fNO3*25;             %OM oxiation by Fe(OH)3
   
    R6a = Rc1.*fNO2*25;
    R6b = Rc2.*fNO2*25;
    
    R7a = Rc1.*fN2O*25; 
    R7b = Rc2.*fN2O*25;   
    
    R9a = Rc1.*fNO2DNRA*25; 
    R9b = Rc2.*fNO2DNRA*25;
    
    R2 = R2a+R2b;
    R6 = R6a+R6b;
    R7 = R7a+R7b;
    R9 = R9a+R9b;
    
    R21 = knh4ox.*(sol(:,:,1)./(KmO2ao+sol(:,:,1))).*(sol(:,:,13)./(KmNH4ao+sol(:,:,13))); %ammonia oxidation with MM kinetics AML
    R22 = kno2ox.*(sol(:,:,1)./(KmO2no+sol(:,:,1))).*(sol(:,:,14)./(KmNO2no+sol(:,:,13))); %nitrite oxidation with MM kinetics AML
    R23 = kinO2./(kinO2+sol(:,:,1)).*kamx.*sol(:,:,13)./(KmNO2amx+sol(:,:,13)).*sol(:,:,14)./(KmNH4amx+sol(:,:,14)); %anammox with MM kinetics AML
    
    F = 0.26;
    
    SimRates{1} = R2 * F;        %storing the rates in the cell matrix                
    SimRates{2} = R6 * F;
    SimRates{3} = R7 * F;
    SimRates{4} = R9 * F;
    SimRates{5} = R21;
    SimRates{6} = R22;
    SimRates{7} = R23;
 
   
    IntRates = zeros(7,1);
    x = linspace(0,15,511);
    for i = [1:7]
        IntRates(i) = trapz(x,SimRates{i}(end,:));
    end
  
    save SimRates
    fprintf(('R2 NO3 red: %0.2f\n'), IntRates(1));
    fprintf(('R6 NO2-denit: %0.2f\n'), IntRates(2));
    fprintf(('R7 N2O: %0.2f\n'), IntRates(3));
    fprintf(('R9 DNRA red: %0.2f\n'), IntRates(4));
    fprintf(('R21 NH4_ox: %0.2f\n'), IntRates(5));
    fprintf(('R22 No2_ox: %0.2f\n'), IntRates(6));
    fprintf(('R23 annamox: %0.2f\n'), IntRates(7));
          