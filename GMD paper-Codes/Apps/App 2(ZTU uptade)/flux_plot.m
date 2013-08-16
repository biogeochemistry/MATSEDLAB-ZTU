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

    close all;
    global SimFlux SimValues;
    SimFlux = cell(2,1);                %defining the cell matrix
    x = linspace(0,15,511);             %defining spatial domain, must be the same as for sol
    t = [0:99,100:0.5:200];            %defining the years for which flux is wanted
    pd = zeros(2,301);                  %initializing the matrix for partial derivative
    c_index = 1;
    
%    for i1 = [1,4,5,6,8,18,19,22]                %O2, SO4, Fe2+, S2-, S(0)
    for i1 = [13,14]                %O2, SO4, Fe2+, S2-, S(0)
        for t1 = 1:101
            [~,dudx] = pdeval(0,x,SimValues{i1,1},0);  %calculating dudx at each t1
            pd(c_index,t1) = dudx;
        end
        
        s_index = 2;
        
        for t2 = 102:301
            [~,dudx] = pdeval(0,x,SimValues{i1,s_index},0.1);  %calculating dudx at each t1
            pd(c_index,t2) = dudx;
            s_index = s_index+1;
        end
        
        switch i1
            case 1                                  %O2
                D_coef = 375;
                S_name = 'O2(aq)';
            case 4                                  %SO4
                D_coef = 175;
                S_name = 'SO4(2-)(aq)';
            case 5                                  %Fe2+
                D_coef = 118;
                S_name = 'Fe(2+)(aq)';
            case 6                                  %S2-
                D_coef = 284;
                S_name = 'S(-II)(aq)';
            case 8                                  %S(0)
                D_coef = 100;
                S_name = 'S(0)(aq)';
            case 13                                 %AsO4
                D_coef = 160;
                S_name = 'AsO4(aq)';
            case 14                                 %PO4
                D_coef = 232;
                S_name = 'PO4(aq)';
            case 17                                 %Ca2
                D_coef = 212;
                S_name = 'Ca2(aq)';
        end
        SimFlux{c_index} = pd(c_index,:)*(-0.9)*D_coef; %calculating the flux
        subplot(2,1,c_index);
%        yl = ylim;
        
%        ylim([0,yl(2)]);          
        plot(120:0.5:200,-SimFlux{c_index}(141:end),'LineWidth',1.5);
        hold on;
        set(gcf,'Position',get(0,'ScreenSize'));    %maximize the figure
        set(gcf,'PaperPosition',[3,5,24,11]);       %set the dimensions on paper when the plots are printed
        set(gca,'FontSize',12);                     %set the font size
        
        title(S_name, 'FontSize',16,'FontWeight','bold');
        xlabel('Time (yr)','FontSize',12);
        ylabel('Flux (umol/cm2/yr)','FontSize',12);
        c_index = c_index+1;
    end
end