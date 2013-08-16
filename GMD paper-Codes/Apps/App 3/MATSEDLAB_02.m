function [DataMatrix ] = MATSEDLAB_02 ( fieldName, VarName, numVar )
% Summary of this function goes here
%   Detailed explanation goes here


    DataMatrix = cell(numVar,1);

    for j=1:numVar,
        XMLSheetName = char(VarName(j));%converts into a MATLAB® character array. 
        [DataMatrix{j},txt,raw] = xlsread(fieldName, XMLSheetName);
    end

end

