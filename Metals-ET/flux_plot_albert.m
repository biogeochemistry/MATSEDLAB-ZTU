function flux_plot
%plots the surface flux of solute species as a function of time
%
%flux_plot
%   This function calculate the surface flux of solute species using the
%   'pdeval' command, from the solution obtained by the pdepe solver. The
%   'sol' matrix must be present for this function to work. This can be
%   done either by running the model, or by importing a MATLAB .mat file
%   that contains a 'sol' matrix.
%
%   The calculated flux will be stored in the global cell matrix 'SimFlux',
%   which is accessible from the command window by typing 'global SimFlux'.
%   The function will then produce the plots of flux vs. time for all the
%   species.

    global SimFlux sol;
    SimFlux = cell(11,1);                %defining the cell matrix
    x = linspace(0,50,250);             %defining spatial domain, must be the same as for sol
    t = linspace(0,50,250);            %defining the years for which flux is wanted
    pd = zeros(11,250);                  %initializing the matrix for partial derivative
    index = 1;
    
    % Load the information saved in Result.mat
    load results.mat
    
    varPosition = [2,4,5,6,9,14,18,22,25,29,31];   %O2, SO4, Fe2+, S2-, As, Zn2+, Cu2+, Co2+, Ni2+, Al3+, Mn2+
    NumVars = size(varPosition, 2);
    VarNames = cell(1, NumVars);
    SimValues = cell(1, NumVars);

    
    for i1 = varPosition                %O2, SO4, Fe2+, S2-, As, Zn2+, Cu2+, Co2+, Ni2+, Al3+, Mn2+
        for t1 = 1:250
            [~,dudx] = pdeval(0,x,sol(t1,:,i1),0);  %calculating dudx at each t1
            pd(index,t1) = dudx;
        end
        switch i1
            case 2                                  %O2
                D_coef = 562;
                S_name = 'O2(aq)';
            case 4                                  %SO4
                D_coef = 261;
                S_name = 'SO4(2-)(aq)';
            case 5                                  %Fe2+
                D_coef = 176;
                S_name = 'Fe(2+)(aq)';
            case 6                                  %S2-
                D_coef = 490;
                S_name = 'S(-II)(aq)';
            case 9                                  %As
                D_coef = 500;
                S_name = 'As(III,V)(aq)';
            case 14                                 %Zn
                D_coef = 200;
                S_name = 'Zn(2+)(aq)';
            case 18                                 %Cu
                D_coef = 180;
                S_name = 'Cu(2+)(aq)';
            case 22                                 %Co
                D_coef = 180;
                S_name = 'Co(2+)(aq)';
            case 25                                 %Ni
                D_coef = 180;
                S_name = 'Ni(2+)(aq)';
            case 29                                 %Al
                D_coef = 109;
                S_name = 'Al(3+)(aq)';
            case 31                                 %Mn
                D_coef = 176;
                S_name = 'Mn(2+)(aq)';
        end
        phi=-0.008*x+0.976; %depth dependant POROSITY
        SimFlux{index} = pd(index,:)*(-0.976)*D_coef; %calculating the flux
        subplot(3,4,index);
        plot(t,SimFlux{index});
        hold on;
        set(gcf,'Position',get(0,'ScreenSize'));    %maximize the figure
        set(gcf,'PaperPosition',[3,5,24,11]);       %set the dimensions on paper when the plots are printed
        set(gca,'FontSize',12);                     %set the font size
        
        title(S_name, 'FontSize',12,'FontWeight','bold');
        xlabel('Time (yr)','FontSize',10);
        ylabel('Flux (umol/cm2/yr)','FontSize',10);

        
        
%         nameFile = strcat(S_name, '_fluxes.dat');
%         caca =  SimFlux{index};
%         save(nameFile, 'caca');
        
        VarNames{index} = S_name;
        SimValues{index} = SimFlux{index}';
        
        index = index+1;
    end
    
    %postprocess to excel


    
    % Write results in excel
    space = zeros(1,1);
    WriteSimValuesInExcel('fluxes_results.xls', VarNames, SimValues, ...
                        NumVars, t, space );

end