function PlotSimVsData(VarNames, DataValuesX1, SimValues, ...
                                    NumVars, timeStep, x, phaseInfo )
%PLOTSIMVSDATA Summary of this function goes here
%   Detailed explanation goes here

    str = sprintf('plotting data for time step  %d ', timeStep);
    disp(str);

    for i=1:NumVars,
        
        % Plot simulation data
        h = figure; 
        y = cell2mat(SimValues(i));
        plot(y(timeStep,:), x);
        hold on;
        
        set(gca,'YDir','reverse')
         
        % Plot field data 
        fieldValues = cell2mat(DataValuesX1(i));
        
        % If we have field data draw it
        if (size(fieldValues , 2) == 2)
            scatter(fieldValues(:, 2), fieldValues(:, 1) );
            hleg = legend('simulated values', 'field data');
        else
            hleg = legend('simulated values');
        end

        %Legend settings
        set(hleg,'FontAngle','italic','TextColor',[.3 .2 .1])
        set(hleg, 'Location','NorthEastOutside');

        
        % Plot configuration
        if(phaseInfo(i)== 0)
            %aq phase
            xlabel('Concentration (\mumol/cm3)');
        else
            %s phase
            xlabel('Concentration (\mumol/g)');
        end
        
        
        xlim('auto')
        ylim([-3, 24]);
        
        %axis([xMin,Xmax,Ymin,Ymax]);
        %axis([-1.0,6.0,-1.2,5.55]);
        
        ylabel('Depth (cm)');
        
%         set(h ,'LineWidth',2,...
%                 'MarkerEdgeColor','k',...
%                 'MarkerFaceColor','g',...
%                 'MarkerSize',10)
        
        nameFile = strcat(char(VarNames(i)), '.ps');
        title(char(VarNames(i)), 'FontSize', 16, 'FontWeight','bold');


        % Save plot
        print(h,'-dps', nameFile );
        
        % Close plot (comment it if you want to keep them)
        % close(h);
    end

end

