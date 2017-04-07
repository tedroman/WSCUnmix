function[logOddsLoss,minEIdx]=batchLogOddsLoss(scenarioStruct,gamma)
%presumes workspace is loaded.
for i=1:7
    cData = scenarioStruct{i}.data; %current data
    cE = [scenarioStruct{i}.E_h scenarioStruct{i}.E_KNN];
    cGamma = gamma;
    [logOddsLoss(i),minEIdx(i)] = computeLogOdds(cE,cData,cGamma,i);
    nComp =0;
    if i==1
        nComp =4;
    elseif i==2
        nComp = 5;
    elseif i==3
        nComp = 7;
    elseif i==4
        nComp = 6;
    elseif i==5
        nComp = 5;
    elseif i==6
        nComp = 5;
    elseif i==7
        nComp = 6;
    end
    [dataz,datamu,dataSigma]=zscore(cData);
    [datazc,datazs]=pca(dataz,'Economy',true);
    gmObj = gmdistribution.fit(datazs(:,1:10),nComp);
    %need to compute 'energy ' of gmm.
    
    
end
end