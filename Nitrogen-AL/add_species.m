function add_species(filename, species, index, phase, boundary, d_coef)
%Add a species into the current model
%   produces a new m-file named 'out_filename.m', with the species added,
%   while keeping the old m-file unchanged. If the user is satisfied with
%   the new m-file, he or she can delete the old one and rename the new
%   one with a desirable name.
%
%   add_species(filename, species, index, phase, boundary, d_coef)
%
%   The input parameters:
%       'filename' is a string indicating the name of the m-file
%       'species' is a string that represent the name of the species (alpha-numeric)
%       'index' is an integer, the index number of the species
%       'phase' is an integer, 1 for solid and 0 for solute
%       'boundary' is a number, the boundary condition of the species added
%       'd_coef' is a number, the diffusion coefficient of a solute
%           species. Only needed when the species is a solute.
%
%   Example:
%       add_species('MATSEDLAB.m', 'PO4', 15, 0, 0.5, 100)
%
%           produces a new m-file out_MATSEDLAB.m, with a new species of
%           'PO4' added into the model, as a solute with index number 15,
%           surface concentration 0.5 umol/cm3 and a diffusion coefficient
%           of 100 cm2/yr
%
%       add_species('MATSEDLAB.m', 'P2O5', 16, 1, 3)
%
%           produces a new m-file out_MATSEDLAB.m, with a new species of
%           'P2O5' added into the model, as a solid with index number 16
%           and surface flux of 3 umol/cm2/yr

    if nargin<5
        error('Not enough inputs.');
    end
    
    if phase==1
        if nargin~=5
            error('Too many inputs. Only 5 required for solid species.')
        end
        c_phase = '(s)';                            %solid phase
    elseif phase==0
        if nargin~=6
            error('Wrong number of inputs. 6 is required for solute species.')
        end
        c_phase = '(aq)';                           %solute phase
        D = ['D_',species];                         %defining diffusion coefficient
    else
        error('Invalid phase. 1 for solid and 0 for solute.');
    end
    
    B = ['BC0_',species];                                                       %Declare boundary variable
    tab_2 = [char(9),char(9)];                                                  %String of two tabs
    tab_3 = [char(9),char(9),char(9)];                                          %String of three tabs
    tab_s = [char(10),char(9),char(32),char(32)];                               %String, new line, 6 spaces
    i0 = num2str(index-1);
    i1 = num2str(index);
    
    f_in = fopen(filename,'r');                                                 %read the original m-file
    f_out = fopen(['out_',filename],'w');                                       %write to out_filename
    
    while (~feof(f_in))                                                         %loop, executed until the end of file is reached
        s = fgetl(f_in);                                                        %Read a new line from file
        if phase==0
            s = strrep(s, '%half saturation', [D,'=',num2str(d_coef),';',...    %D Coefficient
                char(10),'%half saturation']);
            f = [';...',char(10),char(9),'(D_bio+',D,')*DuDx(',...
                i1,')-w*u(',i1,')];%f'];
            pl = [';...',tab_s,...
                'ul(',i1,')-',B,'];%pl'];
            r = '(''R in umol/cm3/yr'') + (''R in umol/g/yr'')*F';
        else
            f = [';...',char(10),char(9),'D_bio*DuDx(',...
                i1,')-w*u(',i1,')];%f'];
            pl = [';...',tab_s,...
                B,'/F];%pl'];
            r = '(''R in umol/cm3/yr'')/F + (''R in umol/g/yr'')';
        end
        s = strrep(s, ['};%u',i0], [',...',tab_2,'%u',...
            i0,char(10),tab_3,'''',species,c_phase,'''};%u',i1]);               %Adding VarNames
        s = strrep(s, 'BC0_O2 BC0_FeOH3', ['BC0_O2 BC0_FeOH3 ',B]);             %Adding Boundary
        s = strrep(s, '%bioturbatiuon', [B,'=',num2str(boundary),';',...
            char(10),'%bioturbatiuon']);                                        %Setting Boundary
        s = strrep(s, '];%f',f);                                                %Modifying f Matrix
        s = strrep(s, '];%pr', [';...',tab_s,'w*ur(',i1,')];%pr']);             %Modifying pr
        s = strrep(s, '];%s', [';...',char(10),char(9),r,'];%s']);              %Modifying s Matrix
        s = strrep(s, '];%pl', pl);                                             %Modifying pl
        s = strrep(s, '];%ql', [';',num2str(phase),'];%ql']);                   %Modifying ql
        fprintf(f_out,'%s\r\n',s);                                              %Write modified line to file
    end
    fclose(f_in);                               %close file
    fclose(f_out);
end