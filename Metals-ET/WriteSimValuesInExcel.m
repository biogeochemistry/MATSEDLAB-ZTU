function WriteSimValuesInExcel(fileName, VarNames, ...
                                    SimValues, NumVars, time, depth)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



for i=1:NumVars,
    MatVals = cell2mat(SimValues(i));
    [m, n]= size(MatVals);
    MatValsWithTime = zeros(m+1, n+1);
    MatValsWithTime(2:m+1,1)= time';
    MatValsWithTime(2:m+1,2:n+1)= MatVals;
    MatValsWithTime(1,2:n+1)= depth;
    MatValsWithTime(1,1) = -1;
    xlswrite(fileName, MatValsWithTime, char(VarNames(i)), 'B2');
end

 xlswrite(fileName, time', 'time', 'A1' );
 xlswrite(fileName, depth', 'depth', 'A1' );

end

