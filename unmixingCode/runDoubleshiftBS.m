function[numClusVec,repIdxAllReplicates]=runDoubleshiftBS(dataBase,h,numRep)
    numPts = size(dataBase,1);
    numClusVec = zeros(1,numRep);
    repIdxAllReplicates = zeros(numPts,numRep);
for i = 1:numRep
        %obtain a sampling (bootstrap)
        replicateIdx = randi(numPts,1,numPts);
        replicateIdx = unique(replicateIdx);
        dataSC = dataBase(replicateIdx,:);
        [numClusVec(i),reps]= runDoubleShift(dataSC,h);
        repIdxAllReplicates(replicateIdx,i) = replicateIdx(reps(end,:));
end
end