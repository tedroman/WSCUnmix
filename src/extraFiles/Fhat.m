function[fVal]=Fhat(cParam,numRep,numPts,numClusBase,numClus)
fVal =0;
for j = 1:numRep
   fVal = fVal+(sqrt(numPts)*(numClus(j)-numClusBase)<=cParam);
end
fVal = fVal/numRep;
end