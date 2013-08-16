function result_plot(varargin)
%result_plot plots the result of MATSEDLAB and saves them as .emf files
%
%By default, each figure contains the plots for 4 species, and the figures
%will be saved as 'plot1.emf', 'plot2.emf', etc.
%
%result_plot
%   plots the concentrations of all species at the end of the simulation
%   period (set to 400 years by default), as a function of depth. The 'sol'
%   matrix must be present for the plotting function to work. This can be
%   done either by running the model, or by importing a MATLAB .mat file
%   that contains an 'sol' matrix.
%
%result_plot(filename)
%   plots the data in the excel file 'filename' on top of those produced by
%   the model. e.g. result_plot('FIELD_DATA.xls')
%
%   Note that the data in the excel file has to be in the format similar to
%   the template provided, thus, a sheet for each species, while the first
%   column is the depth, and the second column is the corresponding
%   concentration at each depth. The names of the sheet must match those
%   used in the variable 'VarNames', and there has to be a sheet for each
%   species, even if there is no data for that species.
%
%result_plot(duration, interval)
%   plots the concentrations of all species as a function of depth, for the
%   last 'duration' years, with an interval of 'interval' years. e.g.
%               result_plot(40,10)
%   will plot the concentrations at the last 40 years, with an interval of
%   10 years. (So there will be five curves on the same plot for each species.)

    close all;
    global SimValues;                     %call the global variable storing the concentration
    
    VarNames = {'O2(aq)',...        %u1
            'Fe(OH)3(s)',...    %u2
            'FeOOH(s)', ...     %u3
            'SO4(2-)(aq)', ...  %u4
            'Fe(2+)(aq)', ...   %u5
            'S(-II)(aq)',...    %u6
            'FeS(s)',...        %u7
            'S(0)(aq)',...      %u8
            'OM1(s)',...        %u9
            'OM2(s)',...        %u10
			'AsFea(s)',...		%u11
			'AsFeb(s)',...		%u12
			'AsO4(aq)',...		%u13
			'PO4(aq)',...		%u14
			'PO4adsa(s)',...	%u15
			'PO4adsb(s)',...	%u16
			'Ca2(aq)',...		%u17
			'CaPO4(s)'};%u18
            
    phase = [0;1;1;0;0;0;1;0;1;1;1;1;0;0;1;1;0;1];  %1 for solid, 0 for solute
    Data = cell(18,1);                          %cell matrix for storing field data from excel sheet
    x = linspace(0,15,511);                     %spatial domain, same as the one used in model
    f_index = 1;                                %index for the plots (and for the corresponding .emf files)
    
    if nargin < 2                               %result_plot, or result_plot(filename)
        for i1 = 1:18
            if rem(i1,4)==1                     %create a new figure if the current figure
                figure;                         %already contains 4 plots
            end
            
            y = SimValues{i1,end}(end,:);
     
            subplot(1,4,rem(i1+3,4)+1);
            plot(y,x);
            hold on;
            set(gca,'YDir','reverse');              %reverse the y-axis, since depth should increase from top to bottom
            set(gcf,'Position',get(0,'ScreenSize'));%maximize the figure
            set(gcf,'PaperPosition',[3,5,24,11]);   %set the dimensions on paper when the plots are printed
            set(gca,'FontSize',8);                  %set the font size
            
            if nargin == 1                          %if filename is provided
                Data{i1} = xlsread(varargin{1}, char(VarNames(i1)));    %stores the data in excel sheet in the cell matrix Data
                field_data = cell2mat(Data(i1));                        %convert cell to matrix for easier access of data
                if (size(field_data,2) >= 2)
                    scatter(field_data(:,2), field_data(:,1) );         %plots a scatter diagram of conc. vs. depth
                end
            end
            
            xl = xlim;
            if xl(1)<0
                xlim([0,xl(2)]);                    %set the x-axis to start from 0 (if it's negative)
            end
            if phase(i1)==1
                xlabel('\mumol/g','FontSize',8);    %units for solid
            else
                xlabel('\mumol/cm3','FontSize',8);  %units for solute
            end
            ylabel('Depth (cm)');                   %name of y-axis
            title(char(VarNames(i1)), 'FontSize',16,'FontWeight','bold');   %title of the plots is the name of the species
        
            if rem(i1,4)==0
                print(gcf,'-dmeta',['plot',num2str(f_index)]);      %if the current figure contains 4 plots already,
                f_index = f_index + 1;                              %save it as a .emf file
            end
        end
    elseif nargin == 2                              %result_plot(duration, interval)
        duration = varargin{1};
        interval = varargin{2};
        n_curves = ceil((duration+1)/interval);     %no. of curves (needed for total FeOH and total s)
        Total_FeOH = zeros(n_curves,511);           %initializing total FeOH
        Total_S = zeros(n_curves,511);              %initializing total S
        for i1 = 1:16
            if rem(i1,4)==1                         %create a new figure if the current figure
                figure;                             %already contains 4 plots
            end
        
            if i1<15                                                    %if the species is stored in the sol matrix,
                y = sol((end-duration):interval:end,:,i1);              %let y be the concentration of that species at the time point specified
                if (i1==2)||(i1==3)
                    Total_FeOH = Total_FeOH + y;                        %if the species is FeOH or FeOOH
                elseif (i1==7)||(i1==9)||(i1==10)||(i1==11)||(i1==14)
                    Total_S = Total_S + y;                              %if the species is FeS, S8, FeS2, OMS or OMSO
                end
            elseif i1==15
                y = Total_FeOH;
            else
                y = Total_S;
            end
     
            subplot(1,4,rem(i1+3,4)+1);
            plot(y,x);
            hold on;
            set(gca,'YDir','reverse');              %reverse the y-axis, since depth should increase from top to bottom
            set(gcf,'Position',get(0,'ScreenSize'));%maximize the figure
            set(gcf,'PaperPosition',[3,5,24,11]);   %set the dimensions on paper when the plots are printed
            set(gca,'FontSize',8);                  %set the font size

            xl = xlim;
            if xl(1)<0
                xlim([0,xl(2)]);                    %set the x-axis to start from 0 (if it's negative)
            end
            if phase(i1)==1
                xlabel('\mumol/g','FontSize',8);    %units for solid
            else
                xlabel('\mumol/cm3','FontSize',8);  %units for solute
            end
            ylabel('Depth (cm)');                   %name of y-axis
            title(char(VarNames(i1)), 'FontSize',16,'FontWeight','bold');   %title of the plots is the name of the species
        
            if rem(i1,4)==0
                print(gcf,'-dmeta',['plot',num2str(f_index)]);      %if the current figure contains 4 plots already,
                f_index = f_index + 1;                              %save it as a .emf file
            end
        end
    else
        error('Too many inputs.')                   %if no. of inputs is greater than 2, give an error message
    end
end
