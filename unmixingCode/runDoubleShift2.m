function[doubleShiftReps]=runDoubleShift2(dataToCluster01Norm)
pathDists = get_path_metric_squared(dataToCluster01Norm);
    t_depth=1; %1 step of euc shift.
    hc1 = 1;
    hc2 = 1;
    
    %or larger hc2
   % hc2 = 0.01;
    doubleShiftReps = (doubleshift(dataToCluster01Norm,hc2,t_depth,...
                hc1));
end
