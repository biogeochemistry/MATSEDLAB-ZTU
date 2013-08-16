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
    global SimValues1 SimValues2;                     %call the global variable storing the concentration
    
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
			'PO4(aq)',...		%u11
			'PO4adsa(s)',...	%u12
			'PO4adsb(s)',...	%u13
			'Ca2(aq)',...		%u14
			'CaPO4(s)',...		%u15
			'NO3(aq)',...		%u16
			'NH4(aq)'};%u17
            
    phase = [0;1;1;0;0;0;1;0;1;1;0;1;1;0;1;0;0];  %1 for solid, 0 for solute
    Data = cell(17,1);                          %cell matrix for storing field data from excel sheet
    x = linspace(0,15,511);                     %spatial domain, same as the one used in model
    
    for season = 1:2    
        f_index = 1;                                %index for the plots (and for the corresponding .emf files)
        for i1 = 11:13
            if rem(i1,4)==1                     %create a new figure if the current figure
                figure;                         %already contains 4 plots
            end
            
            switch season
                case 1
                    y = SimValues1{i1,1}(end,:);
                    f_name = 'Field_Data_s.xlsx';
                    p_name = 'Spring';
                case 2
                    y = SimValues2{i1,1}(end,:);
                    f_name = 'Field_Data_w.xlsx';
                    p_name = 'Winter';
            end
     
            subplot(1,4,rem(i1+3,4)+1);
            plot(y,x);
            hold on;
            set(gca,'YDir','reverse');              %reverse the y-axis, since depth should increase from top to bottom
            set(gcf,'Position',get(0,'ScreenSize'));%maximize the figure
            set(gcf,'PaperPosition',[3,5,24,11]);   %set the dimensions on paper when the plots are printed
            set(gca,'FontSize',8);                  %set the font size
            
            Data{i1} = xlsread(f_name, char(VarNames(i1)));    %stores the data in excel sheet in the cell matrix Data
            field_data = cell2mat(Data(i1));                        %convert cell to matrix for easier access of data
            if (size(field_data,2) >= 2)
                scatter(field_data(:,2), field_data(:,1) );         %plots a scatter diagram of conc. vs. depth
            end
       
            if phase(i1)==1
                xlabel('\mumol/g','FontSize',8);    %units for solid
            else
                xlabel('\mumol/cm3','FontSize',8);  %units for solute
            end
            ylabel('Depth (cm)');                   %name of y-axis
            title(char(VarNames(i1)), 'FontSize',16,'FontWeight','bold');   %title of the plots is the name of the species
        
            if (rem(i1,4)==0)||(i1==17)
                print(gcf,'-dmeta',[p_name,num2str(f_index)]);      %if the current figure contains 4 plots already,
                f_index = f_index + 1;                              %save it as a .emf file
            end
        end
    end
end
