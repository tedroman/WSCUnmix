function[clusterDimensionVec,dimStarts,dimEnds,dimList] = estDimSliver(T2,data,maxdim,sigma)
for i = 1:max(T2)
   currClustData = data(T2==i,:);
   clusterDimensionVec(i) = sliver_dim_est(sigma,currClustData,maxdim);
   dimStarts(i) = 1;
   dimEnds(i) = clusterDimensionVec(i);
   dimList{i}=1:dimEnds(i);
end

end