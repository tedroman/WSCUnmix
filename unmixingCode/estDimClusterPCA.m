function[clusterDimEst]=estDimClusterPCA(gaussMVarExpl,T,dataSCBase)
clusterDimEst = zeros(1,max(T));
for q = 1:max(T)
       dataCluster =dataSCBase(T==q,:);
       [~,~,~,~,explained]=pca(dataCluster,'Economy',true,'Centered',false); %this won't say WHICH dimensions.
       for s = 1:length(explained)
           if (explained(s) < mean(gaussMVarExpl(s,:))+3*std(gaussMVarExpl(s,:))) %should be p-val of < 0.01 in this way.
               break;
           end
       end
       clusterDimEst(q)=s;
       
end
end