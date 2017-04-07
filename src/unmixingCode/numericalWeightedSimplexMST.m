function V = numericalWeightedSimplexMST(V0,data,gamma,reweights,optset)


%should be more samples than dimension
if size(data,1) > size(data,2)
    data = data';
end

minplex = @(v) (simplexfitErrorMultMSTWeighted(v,data,gamma,reweights));
if nargin ==3
    optset = optimset('Display','iter','MaxIter',100);
end

V = fminsearch(minplex,V0,optset);
end