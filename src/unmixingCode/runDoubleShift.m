function[numClus,reps]=runDoubleShift(data,h)
    reps = doubleshift(data,h,1,h);
    %assumes same neighborhood size for euc and non-euc shifts
    %assumes 1 step of euc. medoidshift, then as many kernelmshifts as are
    %needed
    numClus = length(unique(reps(end,:)));
end