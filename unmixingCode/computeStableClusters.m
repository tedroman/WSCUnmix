function [T] =computeStableClusters(numRep,numPts,repIdxAllReplicates,CIInt)
%this iteration of the method will use the adjusted rand index to compute a
%'stable' clustering over bootstrapped replicates.

%need to implement the adjusted rand index for bs replicates.
%do we consider "not selected" its own partition? otherwise we are not
%evaluating the same data, so ARI is invalid.


%#TODO

%coClusterAllPts = zeros(numPts,numPts,numRep);
repIdxAllReplicates = transpose(repIdxAllReplicates);
coCluster = zeros(numPts,numPts);
nTimesSelected = zeros(numPts,numPts);
for repInd = 1:numRep
   for p1 = 1:numPts
       for p2 = 1:numPts
           if repIdxAllReplicates(repInd,p1)==0 || repIdxAllReplicates(repInd,p2)==0
               nTimesSelected(p1,p2)=nTimesSelected(p1,p2);%do not increment, both were not selected.
           elseif repIdxAllReplicates(repInd,p1)==repIdxAllReplicates(repInd,p2)
               coCluster(p1,p2)=coCluster(p1,p2)+1;
               nTimesSelected(p1,p2)=nTimesSelected(p1,p2)+1;
           end
       end
   end

end
coCluster = coCluster./numRep;
hiCoCluster = zeros(numPts,numPts);
%what are the 'highly' coClustered.
highIdx = coCluster >0; %arbitrary parameter value?

%this whole setup should probably be changed to adjusted rand index to
%measure relatively stable clusters.

for i = 1:numPts
    for j = 1:numPts
        if highIdx(i,j)
            hiCoCluster(i,j)=coCluster(i,j);
        end
    end
end

%now compute the stable clusters from this.
Z = linkage(1-coCluster);
T = cluster(Z,'maxclust',max(CIInt));

    
end