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
    SimFlux = cell(5,1);                %defining the cell matrix
    x = linspace(0,15,511);             %defining spatial domain, must be the same as for sol
    t = linspace(0,200,201);            %defining the years for which flux is wanted
    pd = zeros(5,201);                  %initializing the matrix for partial derivative
    index = 1;
    
    for i1 = [12,14,15,13,16]                %NO3, NO2, N2O, NH4, N2
        for t1 = 1:201 		% set t1 = 1:201 for the whole time range
            [~,dudx] = pdeval(0,x,sol(t1,:,i1),0);  %calculating dudx at each t1
            pd(index,t1) = dudx;
        end
        switch i1
            case 12                                  %NO3
                D_coef = 539;
                S_name = 'NO3(aq)';
            case 14                                  %NO2
                D_coef = 520;
                S_name = 'NO2(aq)';
            case 15                                  %N20
                D_coef = 666;
                S_name = 'N20(aq)';
            case 13                                  %NH4
                D_coef = 555;
                S_name = 'NH4(aq)';
            case 16                                  %N2
                D_coef = 550;
                S_name = 'N2';
        end
        SimFlux{index} = pd(index,:)*(-0.89)*D_coef; %calculating the flux
        subplot(2,3,index);
        plot(t,SimFlux{index});
        hold on;
        set(gcf,'Position',get(0,'ScreenSize'));    %maximize the figure
        set(gcf,'PaperPosition',[3,5,24,11]);       %set the dimensions on paper when the plots are printed
        set(gca,'FontSize',12);                     %set the font size
        
        title(S_name, 'FontSize',12,'FontWeight','bold');
        xlabel('Time (yr)','FontSize',10);
        ylabel('Flux (umol/cm2/yr)','FontSize',10);
        index = index+1;
    end
end
