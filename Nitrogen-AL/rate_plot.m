function rate_plot(varargin)
%plots the rate of reaction as a function of depth at the end of the
%simulation time
%
%rate_plot
%   plots all the rates of reaction, as a function of depth, at the end of
%   the simulation time (set to 400 years by default). The rate matrix
%   'SimRates' must be present in order for this function to work
%   correctly. This can be done by calling the function 'rate_calculate' to
%   calculate the concentration from the 'sol' matrix.
%
%rate_plot('ratename1','ratename2','ratename3',...)
%   plots the rates of reaction specified by the user. For example,
%       rate_plot('R1','R7','R_7')
%   will only plot these 3 rates.

    close all;
    global SimRates;
    
    x = linspace(0,15,511);
    RateIndex = {'R1';'R2';'R3';'R4';'R5';'R6';'R7';'R_7';...
        'R8';'R9';'R10';'R11';'R_11';'R12';'R13'};
    
    RateNames = {'Nitrate red',...                     %R1
                'Denit to N2O',...           %R2
                'Denit to Nitrogen gas',...             %R3
                'DNRA',...               %R4
                'NH4 ox to NO2',...             %R5
                'NO2 oxidation',...            %R6
                'annamox',...                   %R7
                'FeS precipitation',...                 %R_7
                'OM Sulfidization',...                  %R8
                'S8 formation',...                      %R9
                'Pyrite precipitation',...              %R10
                'Elemental sulfur dissolution',...      %R11
                'Elemental sulfur precipitation',...    %R_11
                'OMSO hydrolysis',...                   %R12
                'OMS Formation from SO4 directly'};     %R13

    RateUnits = [0;0;0;0;0;0;0;0;1;1;0;1;1;1;1];        %1 for umol/g/yr, 0 for umol/cm3/yr
    
    if nargin>0                                             %if user specifies the rate
        index = [];
        for i1 = 1:15
            for i2 = 1:nargin
                if size(RateIndex{i1})==size(varargin{i2})  %check whether the rate name matches
                    if RateIndex{i1}==varargin{i2}
                        index = [index,i1];                 %if a match is found, include that index
                    end                                     %in the index matrix
                end
            end
        end                 
    else
        index = 1:15;                                       %if no specification is given, the index
    end                                                     %matrix includes all the indices, 1 to 15
    
    p_index = 1;                                    %index of plot
    
    for i1 = index
        if rem(p_index,4)==1                        %create a new figure if the current figure
            figure;                                 %already contains 4 plots
        end
        y = SimRates{i1};
        subplot(1,4,rem(p_index+3,4)+1);
        plot(y(end,:),x);
        hold on;
        set(gca,'YDir','reverse');                  %reverse the y-axis, since depth should increase from top to bottom
        set(gcf,'Position',get(0,'ScreenSize'));    %maximize the figure
        set(gcf,'PaperPosition',[3,5,24,11]);       %set the dimensions on paper when the plots are printed
        set(gca,'FontSize',12);                     %set the font size

        xl = xlim;
        if xl(1)<0
            xlim([0,xl(2)]);                        %set the x-axis to start from 0 (if it's negative)
        end
        
        if RateUnits(i1)==1
            xlabel('\mumol/g/yr','FontSize',14);
        else
            xlabel('\mumol/cm3/yr','FontSize',14);
        end
        ylabel('Depth (cm)');
        title(char(RateNames(i1)), 'FontSize',16,'FontWeight','bold');  %title of the plot is the name of the rate
        
        p_index = p_index+1;                        %plot index increases by 1
    end
    
end

