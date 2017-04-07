function[clusterDimEst,dimStart,dimEnd,dimList]=estDimClusterPCA3(gaussMVarExpl,T,data,dataSCCoeffs,maxPCs,numSigmas,c,alphaVal,dispFlag,figIdx)
for clusterIdx = 1:max(T)
    dataInCluster = data(T==clusterIdx,:);
    dimEstNum = getDim(dataInCluster,dispFlag,c,alphaVal,figIdx);
   clusterDimEst(clusterIdx) = dimEstNum;
    %now, WHICH PCS?
    varAmt = [];
    for i = 1:maxPCs
    
        varAmt(i) = std(dataInCluster(:,i));
    end
    [~,varIdx]=sort(varAmt);
    dimStart(clusterIdx) = varIdx(1);
    dimEnd(clusterIdx) = varIdx(min(clusterDimEst(clusterIdx),numel(varIdx)));
    dimList{clusterIdx}=varIdx(1:dimEnd(clusterIdx));
end
end