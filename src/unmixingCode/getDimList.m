function[dimList,clusterDimVec]=getDimList(numRep,numPts,repIdxAllReplicates,CIInt,repsBase,theta,dataSCCoeffs,dataSCBase01Norm,h,dataSCBase,data,dataPCCoeff,gamma,c,alphaVal,dispFlag,figIdx)
   %there is only 1 'winner' for the number of clusters.
   %We need a module to compute the stable clusters.
   disp('Computing Stable Clusters');
   T = computeStableClusters(numRep,numPts,repIdxAllReplicates,CIInt);
   
   %what if we use basereps?
   clusIdx = unique(repsBase);
   T2 = zeros(numPts,1);
   for i =1:length(clusIdx)
      T2(repsBase==clusIdx(i))=i; 
   end
   %scatterConsistentPoints(dataSCBase,T2);
   %jaccardMeanMatrix = computeJaccard(numRep,numPts,repIdxAllReplicates);
   %Z = linkage(1-jaccardMeanMatrix);
   %T = cluster(Z,'maxclust',max(CIInt));
   %There's an issue with this module.
   
   %Now we need to perform dimensionality estimation on the clusters
   %generate Gaussian data against which we can compare.
   %Does this cause issues when the estimated number of clusters is 1?
 
      
    gaussMVarExpl = zeros(size(dataSCBase,2),numRep);
    for p = 1:numRep
       %gData = randn(numPts,theta(1)); %this should be depend. on the dimensionality of the data.
       %new model for data generation.
        for q = 1:size(dataSCBase,2)
            gData = exprnd(mean(dataSCBase(:,q)),size(dataSCBase,1),1);
        end
       [~,~,~,~,expl]=pca(gData,'Economy',true,'Centered',false);
       gaussMVarExpl(:,p)=expl;
    end
    [clusterDimVec,dimStarts,dimEnds,dimList] = estDimClusterPCA3(gaussMVarExpl,T2,data,dataSCCoeffs,theta(1),theta(3),c,alphaVal,dispFlag,figIdx);
    
   

end