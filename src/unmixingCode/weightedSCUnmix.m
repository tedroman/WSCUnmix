function[V,A,VProj,I,E,Ks,CIp,Ps,Fs,RW,penDf,penPf]=weightedSCUnmix(data,dataCoeff,M,dimVec,dReps,dimStarts,dimEnds,dimList,vertexMergeThresh,numNearestNeighbors,maxIter,dataSCBase,gamma)
%This function will perform the weighted unmixing of simplicies designated
%by dimVec, and return the corresponding vertices and adjacency matrices in
%V and A.  In the beginning, we will only prototype the vertices and use
%the prior method of 1 fit on the entire dataset.  We will aim to update
%once a probabilistic framework is derive as to the best method to obtain
%the adjacency matrix, and how bootstraps or some other method (jacknife?)
%may be used in order to leverage the uncertainty in the data.

%cpyright 2015 Ted Roman, advised by Russell Schwartz, Carnegie Mellon
%Univ.

%%Data is a numpts * numPCs matrix of data in the principal components
%%space with simplicial structure to be exploited via soft geometric
%%unmixing.  M is a numpts * numSimplicies matrix of membership in each of
%%the simplicies, as determined previously.  dimVec represents the
%%dimensionality of each of the simplicies to be unmixed.  The goal with
%%this code will be to provide the unmixed choice for vertices and
%%adjacency matrix.


%let theta(5) = probability threshold for vertex merger
%let theta(6) = number of nearest neighbors (K) for KNN graph for vertex
%merger
%let theta(7) = maximum number of iterations(line and non-line the same)




%To start, we only provide vertices
V = {};
VProj = {};
A = [];
%We will aim to use the previous work in robust mst unmix to begin.
if ~exist('gamma','var')
    gamma = 10; %this is a bug. need to pass this in.
end



%forcing dimVec for debugging.
%dimVec = [2 2];  %need to remove this for real data--perhaps this was
%causing issues?


repIdx = [];
%find the reps that are the cluster centers--may be a faster way, but this
%is easy.
for j = 1:size(dReps,2)
    repIdx = [repIdx find(dReps(:,j)==0)];
end
for i =1:length(dimVec)
   currDim = dimVec(i);
   %The data are already in the PC space., so take the correct number of
   %PCs from the estimate.
   %similarly, may need to take the correct number of Coeff values--do we?

   %dataOfDim = data(:,1:currDim); %but this might be the "other" coeffs.
   
   %We also need the coeffs.
   [encl,a,F,P,errors,C0p,C0,K,enclProj,reweights,penDf,penPf]=robustUnmixWeightedMST(data',... %modified 19 Aug.
       (currDim+1),gamma,M(:,i),dataCoeff,dataSCBase,dReps,repIdx,dimStarts(i),dimEnds(i),dimList{i},maxIter); %need to modify for consideration of M matrix.
   I{i} = C0;
   V{i} = encl;    
   VProj{i} = enclProj;
   E{i}=errors;
   Ks{i}=K;
   CIp{i}=C0p;
   Ps{i}=P;
   Fs{i}=F;
   RW{i}=reweights;
end
%ignoring the other outputs from the unmixing for now.
%This is where the other module (solving for A) will go in the future.
k = numNearestNeighbors; %number of nearest neighbors to use (theoretically difficult to
%solve but should be relatively insensitive in the result.)

thresh = vertexMergeThresh; %at what probability density (or cum) threshold do we want to merge points?
[mergeIndicatorMatrix,newVMtrx,adjMtrxNew] = findSimilarVertices(data,V,k,thresh,dataSCBase);

%draw to test.%ignore drawing for now:
%{
carray =['k','c','m','g','b','r','y'];
    figure(1);
    hold on;
    
    scatter3(data(:,1),data(:,2),data(:,3),45,[M zeros(size(data,1),1)]);
    for i = 1:size(adjMtrxNew,1)
        scatter3(newVMtrx(i,1),newVMtrx(i,2),newVMtrx(i,3),60,strcat(carray(i),'o'),'filled');
        for j = 1:size(adjMtrxNew,2)
        if adjMtrxNew(i,j)
           line([newVMtrx(i,1) newVMtrx(j,1)],[newVMtrx(i,2) newVMtrx(j,2)],[newVMtrx(i,3) newVMtrx(j,3)],'Color','k'); 
        end
        end
    end
  %}  
   A = adjMtrxNew; %added this in 12 june 2015 (helpful?)
   V = newVMtrx;
   %need to also get a new VProj
   VProj = V*dataCoeff(:,1:size(V,2))' + repmat(mean(data),size(V,1),1);


end
function[mergeIndMtrx,newVMtrx2,adjMtrxNew]=findSimilarVertices(data,Vcell,k,thresh,dataScores)

%here we will use the number of overlapping nearest neighbors of inferred
%vertices compared against the total number of data points and number of
%nearest neighbors chosen in order to formulate relative to a
%hypergeometric test whether or not inferred vertices ought to be merged in
%order to form a simplicial complex.  This should represent an improvement
%against previous approaches in that (a) it does not require bootstrapping,
%(b) it is not ad hoc, and (c) it provides a statistically-reasoned
%approach as to how the mergers are made.
grpCell={};
currIdx = 1;
VMtrx = [];
for i = 1:length(Vcell)
    VMtrx = [VMtrx; Vcell{i}]; %this fails in the real data; 
    %why?  One of the points is 1-dimensional, and has only "0" as a
    %points--why would this happen?
    grpCell{i}=currIdx:currIdx+size(Vcell{i},1)-1; %which vertices are
    %originally connected in a fully connected fashion.
    currIdx = currIdx+size(Vcell{i},1);
end
mergeIndMtrx = zeros(size(VMtrx,1),size(VMtrx,1));

%find the knn for the vertices. don't want to use projection because shapes
%are in Pc space.  should be using pc of data.
nnidx = knnsearch(dataScores,VMtrx,'K',k);

%assign parameters for the hypergeometric test
bigK = k;%number of possible successes.
bigN = size(data,1); %size of population is number of data points.
littlen = k; %number of attempts
for v1 = 1:size(VMtrx,1)
    for v2 = v1+1:size(VMtrx,1)
        littlek = length(intersect(nnidx(v1,:),nnidx(v2,:))); %how many in common
        %compute probability of the observation assuming a random
        %assignment of nearest neighbors--this needs a better name, it's
        %more like the probability that the overlap is random, and it is <
        %because of cdf vs. pdf. and density vs. cumulative.
        hcdf = hygecdf(littlek,bigN,bigK,littlen); 
        %pRand = nchoosek(bigK,littlek)*nchoosek(bigN-bigK,littlen-littlek)...
            %/nchoosek(bigN,littlen) %should this computation be done in log space?
        mergeIndMtrx(v1,v2)= littlek > bigK^2/bigN; %littlek>1; %just for kicks.    %(1-hcdf) > thresh; %try switching to greater?
        %if the probability that the observed overlap is due to random
        %chance is small ( < threshold), then merge.  Otherwise, the
        %randomness may have contributed to the observed overlap--do not
        %merge.
        
    end
end
mergeIndSym = mergeIndMtrx;
for i = 1:size(mergeIndSym,1)
  for j = 1:size(mergeIndSym,2)
   if mergeIndMtrx(i,j)
       mergeIndSym(i,j)=1;
       mergeIndSym(j,i)=1;
   end
  end
end
%now here is where we could compute the new adjacency matrix let's first
%check that the rest is working though
adjMtrxOld = zeros(size(VMtrx,1),size(VMtrx,1));
for i = 1:length(grpCell)
   adjMtrxOld(grpCell{i},grpCell{i})=1; 
end

newTmpVMtrx = zeros(2*size(VMtrx,1),size(VMtrx,2));
newTmpVMtrx(1:size(VMtrx,1),:)=VMtrx;
for i = 1:size(mergeIndSym,1)
    toMergeWith = mergeIndSym(i,:); %logical vector of which points to
    %merge with the current point.
    if any(toMergeWith)
        newPt = mean([VMtrx(i,:);VMtrx(logical(toMergeWith),:)]);
    else 
        newPt = VMtrx(i,:);
    end
    newTmpVMtrx(i+size(VMtrx,1),:) = newPt; %??? why is this not correct?
    
end
newVMtrx = (newTmpVMtrx(size(VMtrx,1)+1:end,:));
[newVMtrx2,~,idxNew] = unique(newVMtrx,'rows');

%now we need to compute the new adjacency matrix.
%newVMtrx(idxOrig(i),:)) = newVMtrx2(i,:);
%how do we know to which row / column of the old adjacency matrix these
%pertain?  The old ordering was preserved.  How big is newVMtrx2?
%Finally newVMtrx2 is the points that we want, now we just have to map the
%computation of the new adjacency matrix to that space.

%idxNew(i) says that the ith row in the old is now the idxNewth row in the
%new, so for j = 1:max(idxNew), if sum(idxNew==j) > 1, some merger
%occurred.  If a merger occurred, then the new adjacency should be a
%composite of all of the connections of the old adjacency.  How do we do
%this?
%we can start by looping through the old adjacency matrix.  if (i,j) in the
%old, then idxNew(i),idxNew(j) in the new.  What happens?

adjMtrxNew = zeros(size(newVMtrx2,1),size(newVMtrx2,1));
for i = 1:size(adjMtrxOld,1)
    for j =1:size(adjMtrxOld,2)
        if adjMtrxOld(i,j)
            adjMtrxNew(idxNew(i),idxNew(j))=1;
        end
    end
end
%does the simple approach work?
%seems yes.


end