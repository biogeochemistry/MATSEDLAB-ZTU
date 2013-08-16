function [DataMatrix ] = GetFieldVariables ( fieldName, VarName, numVar )
% Summary of this function goes here
%   Detailed explanation goes here


    DataMatrix = cell(numVar,1);

    for j=1:numVar,
        XMLSheetName = char(VarName(j));
        [DataMatrix{j},txt,raw] = xlsread(fieldName, XMLSheetName);
    end

end

