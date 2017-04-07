function [V,A,data,VProj,M,I,E,K,CIp,Ps,Fs,RW,penDf,penPf,dimList] = weightedSCUnmixMainSliver(data,theta,gamma,c,sigma)
%cpright 2015 Ted Roman, advised by Russell Schwartz, Carnegie Mellon Univ.

%This code takes as input a numPts * numDimensions dataset of presumed
%tumor data, called data.  The data is assumed to have some experimental
%noise.  The data is presumed to be bulk measurements of some genomic data,
%and as such the model with which we attempt to unmix the data is using a
%simplicial complex.  The simplicial complex is assumed to come from low
%dimensional simplicies, for which the low dimensionality and the number of
%such simplicies is unknown.  It is further unknown the fashion in which
%the simplicial complex maintains sub-simplicies.  We have several aims to
%process the data.  The first aim will be to ascertain the most likely
%estimates for the number of simplicies in the data and the minimum
%dimenisonality of such simplicies, following a maximum likelihood
%framework.  We will then aim to determine a fuzzy clustering of the data
%representing a confidence that a data point is a member of a particular
%simplex within the complex.  Finally, we will use this information
%combined with constrained optimization techniques to infer the best
%fitting simplicial complex, subject to a soft penalty functino
%(regularization) on the minimum spanning tree of each simplex, as well as
%the confidence of membership in each simplex.  Here we use the minimum
%spanning tree as it has an interpretation as the evolutionary distance
%that the tumor has spanned.  Lastly, we merge the simplicies together to
%form a coherent simplicial complex structure.

%theta is a vector of parameters, used later in the method.  

%from the batch caller: 
%call the parameters all together theta.
%let theta(1) = maximum number of dimensions
%let theta(2) = number of bootstrapped replicates
%let theta(3) = threshold for dimensionality estimation
%let theta(4) = h (h1 and h2 being the same)
%let theta(5) = probability threshold for vertex merger
%let theta(6) = number of nearest neighbors (K) for KNN graph for vertex
%merger
%let theta(7) = maximum number of iterations(line and non-line the same)


%the first step of the algorithm is to determine the dimensionality and
%number of clusters of the data.
%We use the derivations of MLE and CI presented in Wasserman (06), as these
%outputs represent parameters to other parts of the model.

%First, we need to determine the optimal hyperparameter for the
%determination of the dimensionality / clustering (bandwidth for
%medoidshift)

%Because the clustering is to be done in the PC space, we first consdier
%principal components of the data as our base.

%check to see if we should generate data.
I = [];
E = [];
K = [];
CIp = [];
Ps = [];
Fs = [];
RW = [];

if ~exist('data','var')
   numDim = 10; %only use 10 D so as not to compress the shapes.
   
    v1 = zeros(3,numDim);
    v2 = zeros(3,numDim);
    v1(2,1) =1;
    v1(3,2) =1;
    v2(2,3) =1;
    v2(3,2) =1;
    
    %both of these share the origin.
    numPts = 300;
    nPts1 = floor(numPts/2);
    
    dataS1 = getSimulatedData(nPts1,numDim,'unif',v1);
    dataS2 = getSimulatedData(numPts-nPts1,numDim,'unif',v2);
    
    dataSCBase = [dataS1;dataS2];
    
    %optional noise
    dataSCBase = dataSCBase+randn(numPts,numDim)*0.01;
    data = dataSCBase;
    dataSCBase = [];
else
    numPts = size(data,1);
end
disp('Computing PCs');
[dataPCCoeff,dataPC] = pca(data,'Centered',false,'Economy',true);
dataSCBase = dataPC(:,1:theta(1)); %this is the upper limit parameter.
dataSCCoeffs = dataPCCoeff(:,1:theta(1));
disp('Computing bandwidth');
h = theta(4); %will put this into place at some point.scanForBandwidth(dataSCBase); 

%Now that we have determined the optimal bandwidth hyperparameter, (h),
%we can issue a medoidshift clustering on the dataset
%try on the base to get the mle estimate.
disp('Running Mshift');
%[numClusBase,~] = runMedoidshift(dataSCBase,h);

%we should modify this module to be a two-stage unmixing process.

%doubleshift module:

%normalize: 
dataSCBaseMax = max(dataSCBase);
dataSCBaseMin = min(dataSCBase); %should be correct dimensions.
dataSCBase01Norm = (dataSCBase-repmat(dataSCBaseMin,size(dataSCBase,1),1))./(repmat(dataSCBaseMax-dataSCBaseMin,size(dataSCBase,1),1));


%what about zscore?
%dataSCBase01Norm = zscore(dataSCBase);

%can we try the ratio thing again
%l2squaredpath / euc is small for points in the same manifold, large if
%distant.
%eucSquaredD  = dist(dataSCBase01Norm').^2; %check that it is correct size.
%pathD = get_path_metric_squared(dataSCBase01Norm);

%pathRatio = pathD./eucSquaredD;
%pathRatio(isnan(pathRatio))=0; %self distances.
%do these cluster cleanly?
%Z = linkage(pathRatio);
%T = cluster(Z,'maxclust',10);



%[numClusBase,repsBase] = runDoubleShift(dataSCBase,h); %requires h1 and h2, but 
%we presume these are equal--can change later. (new 13 Aug 2015) %may need
%to add norm. messing with this 20Aug.

%this should be on the normalized clutering

repsBase = runDoubleShift2(dataSCBase);
numClusBase = length(unique(repsBase)); 
%BS Module
%numRep = 300; %what if we crank up the number of bootstrapped replicates
numRep = theta(2);
%numRep = 10000;
%numRep = 500;
disp('Running Bootstrapped Replicates');
%[numClus,repIdxAllReplicates] = runMedoidshiftBS(dataSCBase,h,numRep);  


%switched to doubleshift 13 Aug 2015
 [numClus,repIdxAllReplicates] = runDoubleshiftBS(dataSCBase01Norm,h,numRep);   
 

 %{
 %to see if the problem is occurring with the actual clustering, we will
 %examine the consensus clustering (via mode)
 consensusClustering = mode(repIdxAllReplicates,2);
 for i = 1:length(consensusClustering)
    if consensusClustering(i)==0
        toSelect = repIdxAllReplicates(i,:)~=0;
       consensusClustering(i) = mode(repIdxAllReplicates(i,toSelect));
    end
 end
 
%plot this just to examine for now --seems to be doing something reasonable
%w/r/t cluster centers.  so why are there issues downstream with where the
%cluster centers are?  they are just row numbers so whether it is in 01
%form or "regular" pc form should be irrelevant.

uniqueReps = unique(consensusClustering); 
carray = 'cbmrk';
hold on;
for i = 1:length(uniqueReps)
    if uniqueReps(i)~=0
        scatter3(dataSCBase01Norm(consensusClustering==uniqueReps(i),1),...
            dataSCBase01Norm(consensusClustering==uniqueReps(i),2),...
            dataSCBase01Norm(consensusClustering==uniqueReps(i),3),...
            125,carray(i),'filled');
        
        scatter3(dataSCBase01Norm(uniqueReps(i),1),...
            dataSCBase01Norm(uniqueReps(i),2),...
            dataSCBase01Norm(uniqueReps(i),3),400,...
            strcat(carray(i),'*'));
    end
end
     %scatter3(dataSCBase01Norm(:,1),...
     %   dataSCBase01Norm(:,2),dataSCBase01Norm(:,3),...
     %   24,'k.');
    
hold off;
%}
%CI Module
disp('Determining CI');
alpha = 0.01;
FhatDec = Fhat(1-alpha/2,numRep,numPts,numClusBase,numClus);
FhatInc = Fhat(alpha/2,numRep,numPts,numClusBase,numClus);
tDec = 1/FhatDec;
tInc = 1/FhatInc;
CI = [numClusBase-tDec/sqrt(numPts) numClusBase+tInc/sqrt(numPts)]
CIInt = ceil(min(CI)):floor(max(CI));
disp('Removing infeasible solutions')
CIInt = CIInt(CIInt>0); 


%do not allow 0 clusters--illogical.



if ~exist('gamma','var')
gamma = 10; %may be needed later.
end
%Now that we have the CI, we need to compute that number of stable clusters
%in the data for each of the integer number of clusters.  We first assume
%there will only be 1 number of clusters in the CI (as empirical evidence
%has suggested at depth / noise levels analyzed)
%We will handle the multiple possibilities case later.

if length(CIInt) > 1
    %disp('Handling this case later, no handle at current time');
    %V=[];
    %A=[];
    %VProj = [];
    %M = [];
    %Init = [];
    %return;
    disp('Multiple integers in interval')
    
    for i = 1:length(CIInt)
       currInt  = CIInt(i);
       [V,A,VProj,I,E,K,CIp,Ps,Fs,RW,M,penDf,penPf,dimList] = ...
           unmix1NumClus(numRep,numPts,repIdxAllReplicates,...
           currInt,repsBase,theta,dataSCCoeffs,dataSCBase01Norm,h,...
           dataSCBase,data,dataPCCoeff,gamma,c,sigma);
        modelsToCompare{i}.V = V;
        modelsToCompare{i}.A = A;
        modelsToCompare{i}.VProj = VProj;
        modelsToCompare{i}.I = I;
        modelsToCompare{i}.E = E;
        modelsToCompare{i}.K = K;
        modelsToCompare{i}.CIp = CIp;
        modelsToCompare{i}.Ps = Ps;
        modelsToCompare{i}.Fs = Fs;
        modelsToCompare{i}.RW = RW; 
        modelsToCompare{i}.M = M;
        modelsToCompare{i}.dimList = dimList;
    end
    
    %now compare the energy functions (probabilities) of each of the models
    Evals = zeros(length(modelsToCompare),1);
    for i = 1:length(modelsToCompare)
       Evals(i) = compute_energy_prob(modelsToCompare{i}.V,modelsToCompare{i}.A,modelsToCompare{i}.E,gamma);
    end
    [~,minIdx] = min(Evals);
    V = modelsToCompare{minIdx}.V;
    A = modelsToCompare{minIdx}.A;
    VProj = modelsToCompare{minIdx}.VProj;
    I = modelsToCompare{minIdx}.I;
    E = modelsToCompare{minIdx}.E;
    K = modelsToCompare{minIdx}.K;
    CIp = modelsToCompare{minIdx}.CIp;
    Ps = modelsToCompare{minIdx}.Ps;
    Fs = modelsToCompare{minIdx}.Fs;
    RW = modelsToCompare{minIdx}.RW;
    M = modelsToCompare{minIdx}.M;
    dimList = modelsToCompare{minIdx}.dimList;
    return;
elseif CIInt ==1
    disp('Only 1 cluster estimated, conflict with model assumptions, exiting');
    %how can we handle this better--> more bootstraps?  Why is this
    %happening?
    
    V=[];
    A=[];
    VProj = [];
    M = [];
    Init = [];
    penDf =0;
    penPf =0;
    return;
else
[V,A,VProj,I,E,K,CIp,Ps,Fs,RW,M,penDf,penPf,dimList]=...
    unmix1NumClus(numRep,numPts,repIdxAllReplicates,CIInt,...
    repsBase,theta,dataSCCoeffs,dataSCBase01Norm,h,dataSCBase,...
    data,dataPCCoeff,gamma,c,sigma);
    
end
    
    
disp('Completed unmixing...')    
    
    
end
function[LSCV]=crossValidationFunction(h,data)
LSCV = fn(h) - 2*f(h,data);
end

function[kernel]=fn(~)
%kernel = exp(-D/h) - 1; %need to fix this integral is improper but constant
%for all values, so choosing an arbitrary constant is fair (??)
kernel = 1;
end

function[kernel]=f(h,D)
kernel = sum(kernelHelper(h,D));
end

function[kernelVec]=kernelHelper(h,D)
projectionD = exp(-D/h)-1;
projectedDist = dist(projectionD');
kernelVec = zeros(size(D,1),1);
for i =1:size(D,1)
    for j = 1:size(D,1)
        kernelVec(i)=kernelVec(i) +projectedDist(i,j)/size(D,1); 
    end
end
end

function[fVal]=Fhat(cParam,numRep,numPts,numClusBase,numClus)
fVal =0;
for j = 1:numRep
   fVal = fVal+(sqrt(numPts)*(numClus(j)-numClusBase)<=cParam);
end
fVal = fVal/numRep;
end

function[hMinVal]=scanForBandwidth(dataSCBase)
hHatVec = 1:0.01:100; %how to choose this (parameter scanning) may be an issue.
    
    hMin = crossValidationFunction(1,dataSCBase);
    hMinVal = 1;
    for i = 1:length(hHatVec)
        %attemptedValue(i) = crossValidationFunction(hHatVec(i),dataSCBase);
        if crossValidationFunction(hHatVec(i),dataSCBase) < hMin
            hMin = crossValidationFunction(hHatVec(i),dataSCBase);
            hMinVal = hHatVec(i);
        end
    end
end

function[numClus,reps]=runMedoidshift(data,h)
    pathMetricBase = get_path_metric(data);
    reps = medoidshift(pathMetricBase,h);
    numClus = length(unique(reps(end,:)));
end

function[numClus,reps]=runDoubleShift(data,h)
    reps = doubleshift(data,h,1,h);
    %assumes same neighborhood size for euc and non-euc shifts
    %assumes 1 step of euc. medoidshift, then as many kernelmshifts as are
    %needed
    numClus = length(unique(reps(end,:)));
end

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

function[numClusVec,repIdxAllReplicates]=runMedoidshiftBS(dataBase,h,numRep)
    numPts = size(dataBase,1);
    numClusVec = zeros(1,numRep);
    repIdxAllReplicates = zeros(numPts,numRep);
for i = 1:numRep
        %obtain a sampling (bootstrap)
        replicateIdx = randi(numPts,1,numPts);
        replicateIdx = unique(replicateIdx);
        dataSC = dataBase(replicateIdx,:);
        [numClusVec(i),reps]= runMedoidshift(dataSC,h);
        repIdxAllReplicates(replicateIdx,i) = replicateIdx(reps(end,:));
end
end

function[jaccardMeanMatrix]=computeJaccard(numRep,numPts,repIdxAllReplicates)
jaccardIdxAllReplicates = zeros(numPts,numPts,numRep);
for repInd = 1:numRep
       for p1 = 1:numPts
          for p2 = p1:numPts
             if repIdxAllReplicates(p1,repInd)==0 || repIdxAllReplicates(p2,repInd)==0
                jaccardIdxAllReplicates(p1,p2,repInd) = 0;
                jaccardIdxAllReplicates(p2,p1,repInd) = 0;
                %if one of these was not selected for this replicate
             else
                 %we know both of these were selected for this replicate
                 p1ClusterRep = repIdxAllReplicates(p1,repInd);
                 p2ClusterRep = repIdxAllReplicates(p2,repInd);
                 
                 p1ClusterIdx = find(repIdxAllReplicates(:,repInd)==p1ClusterRep);
                 p2ClusterIdx = find(repIdxAllReplicates(:,repInd)==p2ClusterRep); %can we do away with find and use logical indexing?
                 
                 if isempty(intersect(p1ClusterIdx,p2ClusterIdx))
                     jaccardIdxAllReplicates(p1,p2,repInd) =0;
                     jaccardIdxAllReplicates(p2,p1,repInd) =0;
                 else
                     jaccardIdxAllReplicates(p1,p2,repInd) = length(intersect(p1ClusterIdx,p2ClusterIdx))/length(union(p1ClusterIdx,p2ClusterIdx));
                     jaccardIdxAllReplicates(p2,p1,repInd) = length(intersect(p1ClusterIdx,p2ClusterIdx))/length(union(p1ClusterIdx,p2ClusterIdx));
                 end
             end
          end
       end
end
    jaccardMeanMatrix = mean(jaccardIdxAllReplicates,3);
end

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


function[] = scatterConsistentPoints(dataSCBase,T)
colAr = 'gbrkcmy';
hold on;
for i = 1:max(T)
    scatter3(dataSCBase(T==(i),1),dataSCBase(T==(i),2),dataSCBase(T==(i),3),50,strcat(colAr(i),'o'),'filled');
end
end


function[T] =computeStableClustersARI(numRep,numPts,repIdxAllReplicates,CIInt)
 %CIInt should be a scalar and inform us as to how many clusters should be
 %formed via hierarchy
 
 %numRep is the number of replicates.
 
 %repIdxAllReplicates is the index of the representative for each of the bs
 %replicates, and 0 if not selected.
 %numPts is the number of points.
 %what is the dimensionality of repIdxAllReplicate?
 %size is numPts * numRep
 
 %What is the ARI of each clustering to every other clustering?
 for i = 1:(numPts)
    for j = 1:(numPts)
       ARIs(i,j) = computeARI(repIdxAllReplicates(i,:),repIdxAllReplicates(j,:)); %can we flip this treating points as replicates over bootstraps?
       disp(['replicate ' num2str(i) ' and replicate ' num2str(j) ' complete']);
    end
 end
 %'distance' is in some sense inversely proportional to ARI.
 ARId = 1-ARIs;
 Z = linkage(ARId);
 T = cluster(Z,'maxclust',max(CIInt)); %this tells us the partitions that are
 %consistently similar (i.e the replicates--not exactly the points).  
 %How do we go from repliates that are consistently partitioned together to
 %points that are consistently partitioned together.  In which case--why
 %use the CIInt maxclust?
 
 
 
 
end
 
function[ARI]=computeARI(reps1,reps2)
%compute the pairwise ARI between the clusterings who have representatives
%in the vectors reps1 and reps2

if length(reps1)~=length(reps2)
    disp('Error: representative vectors must have same length');
    ARI = 0;
    return;
end

clusIdx1 = unique(reps1);
clusIdx2 = unique(reps2);

contingencyTable = zeros(length(clusIdx1),length(clusIdx2));
for i = 1:length(clusIdx1)
    pointsOfCluster1 = (reps1==clusIdx1(i)); %which points are in cluster 1
    for j = 1:length(clusIdx2)
        pointsOfCluster2 = (reps2==clusIdx2(j)); %which points are in cluster 2
        pointsInCommon = and(pointsOfCluster1',pointsOfCluster2');
        contingencyTable(i,j) = sum(pointsInCommon);
    end
    
end

%compute the sums 
avec = sum(contingencyTable,2);
bvec = sum(contingencyTable,1);

idx = 0;
for i = 1:length(clusIdx1)
    for j =1:length(clusIdx2)
        if contingencyTable(i,j) > 1 %o.w. 0.
            idx = idx + nchoosek(contingencyTable(i,j),2);
        end
    end
end

expectIdx1 = 0;
for i =1:length(clusIdx1)
    if avec(i) > 1 %o.w. 0
        expectIdx1 = expectIdx1 + nchoosek(avec(i),2);
    end
end
expectIdx2 = 0;
for j = 1:length(clusIdx2)
    if bvec(j) > 1 % o.w. 0
        expectIdx2 = expectIdx2 + nchoosek(bvec(j),2);
    end
end
expectIdx = expectIdx1*expectIdx2/nchoosek(length(reps1),2);

maxIdx = 0.5*(expectIdx1+expectIdx2);
ARI = (idx-expectIdx)/(maxIdx-expectIdx);
end



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

function[T]=computeSimilarPointsAcrossClusters(numRep,numPts,repIdxAllReplicates,CIInt)


%simple approach at first, cluster how many times they agree vs. disagree +
%agree (total number of replicates)
agreements = zeros(numPts,numPts);
disagreements = zeros(numPts,numPts);
for pt1 = 1:numPts
        pt1Reps = repIdxAllReplicates(pt1,:);
    for pt2 = 1:numPts
        pt2Reps = repIdxAllReplicates(pt2,:);
        %how many times do they agree?
        for i = 1:numRep
            if pt1Reps(i)==pt2Reps(i)
                agreements(pt1,pt2) = agreements(pt1,pt2)+1;
            else
                disagreements(pt1,pt2) = disagreements(pt1,pt2)+1;
            end
        end
    end
end


%how can we correct for chance?  chance agreements and disagreements?
%probability of chance agreement--both are selected for the same cluster.
%If the cluster sizes are fixed but points are randomly assigned (w/o
%replacement, this is a hypergeometric distribution, but choosing 2)

%following rand index form, we can do it like this: 
%sum_i(sizeclus1 choose 2)) * sum_j (sizeclus2 choose 2) / n choose 2

%^call this expected value
%max is then that we match each time.
%treat agreements as contingency table.
%compute rand index on that??
avec = sum(agreements,2);
bvec = sum(agreements,1);
idx = zeros(numPts,numPts);
for i=1:numPts
    for j =1:numPts
            if agreements(i,j) > 1
               idx(i,j) =  nchoosek(agreements(i,j),2); 
            end
    end
end

expect1 = zeros(1,numPts);
expect2 = zeros(1,numPts);
for i=1:numPts
    if avec(i) > 1
       expect1(i) = nchoosek(avec(i),2); 
    end
    if bvec(i) > 1
        expect2(i) = nchoosek(bvec(i),2);
    end
end

adjustedAgreements = (idx-((expect1.*expect2)./nchoosek(numPts,2)))./((0.5*expect1+expect2)-((expect1.*expect2)./nchoosek(numPts,2)));

similarityCounts = agreements./(agreements+disagreements);
Z = linkage(1-similarityCounts);
T = cluster(Z,'maxclust',max(CIInt));
end


function[doubleShiftReps]=runDoubleShift2(dataToCluster01Norm)
pathDists = get_path_metric_squared(dataToCluster01Norm);
    t_depth=1; %1 step of euc shift.
    hc1 = 1;
    hc2 = 1;
    
    %or larger hc2
   % hc2 = 0.01;
    doubleShiftReps = (doubleshift(dataToCluster01Norm,hc2,t_depth,...
                hc1));
end


function[clusterDimEst,dimStart,dimEnd,dimList]=estDimClusterPCA2(gaussMVarExpl,T,data,dataSCCoeffs,maxPCs,numSigmas)
%how do we project into each principal component?  Each column of Coeffs is
%one PC, so we multiply column-wise
for i =1:max(T)
   currClus = data(T==i,:);
   %dataSCBase is NOT in PC space necessarily.
   for j =1:maxPCs %because we consider at most 10--nope this is the straight dope data.
    %disp(size(dataSCCoeffs))
    %disp(j)
       projScores =(currClus - repmat(mean(currClus),size(currClus,1),1))*dataSCCoeffs(:,j);
       varTot(i,j) = var(projScores);
   end
   varTotAll(i) = sum(varTot(i,:));
   [varTotSorted,varTotSortedIdx]=sort(varTot(i,:),'descend');
   dimStart(i) = varTotSortedIdx(1);
   for k = 1:maxPCs
      %disp(k)
      %disp(size(gaussMVarExpl))
      %if no variance i.e. 1 point, then we need a check
      if ~isnan(varTotSorted(k)/varTotAll(i))
       if (100*(varTotSorted(k))/varTotAll(i)) < mean(gaussMVarExpl(k,:)) + numSigmas*std(gaussMVarExpl(k,:)) %took away the 2
         break; 
       end
      else
          %each dimension works equally as well.??
          break; %just pick the first one because it is better for other points.
      end
      
   end
   clusterDimEst(i) = k;
   dimEnd(i) = dimStart(i) + k -1; %this makes additional assumptions.
   dimList{i} = varTotSortedIdx(1:k);
   %need to code in an exception for single points.
end

end


function[V,A,VProj,I,E,K,CIp,Ps,Fs,RW,M,penDf,penPf,dimList]=unmix1NumClus(numRep,numPts,repIdxAllReplicates,CIInt,repsBase,theta,dataSCCoeffs,dataSCBase01Norm,h,dataSCBase,data,dataPCCoeff,gamma,c,sigma)
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
    
    %color points by consistent clusters
    %scatterConsistentPoints(dataSCBase,T); %This is visually problematic.
    
    %we should have points that are more consistent with one another.
    %why is this happening in the computation of T?
    
    %size(GaussMVarExpl) = maxD * numReplicates
    %T2 = numPts * 1; value = clusterVal.
    %size(data) = numpts * total dims
    %theta(1) = maxD
    %theta(3) = boundary condition.
    
    
    %[clusterDimensionVec,dimStarts,dimEnds,dimList] = estDimClusterPCA2(gaussMVarExpl,T2,data,dataSCCoeffs,theta(1),theta(3));
    %this says how many, but does it say which ones?
    
    %[clusterDimensionVec,dimStarts,dimEnds,dimList] = estDimClusterPCA3(gaussMVarExpl,T2,data,dataSCCoeffs,theta(1),theta(3),c);
    [clusterDimensionVec,dimStarts,dimEnds,dimList] = estDimSliver(T2,data,theta(1),sigma);
    %now that we have this, can we RE-estimate the cluster centers? T is
    %ultimately faulty (does not give visually logical clusters, yet the
    %results sort of make sense.--why would this be?
    
    
    
    
    %this module should come prior to the unmixing modules, then we should
    %do the unmixing leveraging the max number of dimensions for a
    %subsimplex...right?
    
    
    %now that we have an estimate of the dimensionality of clusters, stable
    %regions within those clusters, and the number of clusters, we need to
    %derive representative points for each cluster, and the membership as
    %being part of a cluster based on the kernel function distance to
    %those representative points in a non-bootstrapped sense.
    disp('Getting Fuzzy Clusters');
    [M,dReps]=getFuzzyClusterMembership(T2,dataSCBase01Norm,h); %may need to modify to norm, or modify 
    %the computation altogether.
    %very little M for all but a handful of points--is this
    %expected or not? 
    
    %now dReps,2 sum is negative--returning distance in the kernel space.
    %how should this change the unmixing (if at all?)
    
    
    
    %Now that we have a weight (stored in M), we need to be able to use
    %this weight to do the unmixing.
    %This module will require modification to the unmixing code.
    %It should be the case then that all points are part of all of the
    %simplices in the complex, but with regularization parameter = M(i,j)
    %for the jth cluster, ith point.
    
    %try changing to just plain data.
    [V,A,VProj,I,E,K,CIp,Ps,Fs,RW,penDf,penPf]=weightedSCUnmix(data,dataPCCoeff,M,clusterDimensionVec,dReps,dimStarts,dimEnds,dimList,theta(5),theta(6),theta(7),dataSCBase,gamma); %should we send it just the max of the clusterDimVec?
   
    
    disp('unmixing complete');
    %[V,A,VProj,Init] = weightedSCUnmix(dataSCBase,dataPCCoeff, M,clusterDimensionVec);
    %vertices and adjacency matrix.  This is the part that will require
    %significant changes.  First check that the previous work is doing its
    %job.
    %{
   carray =['k','c','m','g','b','r','y'];
    figure(1);
    hold on;
    if size(M,2)==2
        scatter3(dataSCBase(:,1),dataSCBase(:,2),dataSCBase(:,3),45,[M zeros(size(dataSCBase,1),1)],'filled');
    elseif size(M,2)==3
        scatter3(dataSCBase(:,1),dataSCBase(:,2),dataSCBase(:,3),45,M,'filled'); 
    end
        %for i = 1:length(V)
        
        for i = 1:size(A,1)
            
            for j = 1:size(A,2)
                if A(i,j)
                   line([V(i,1) V(j,1)],[V(i,2) V(j,2)],[V(i,3) V(j,3)]); %will modify color later. 
                end
            end
            
        end
        scatter3(V(:,1),V(:,2),V(:,3),100,strcat('k','o'),'filled');
    %end
    %}
    %what can we do for an adjacency matrix?
    %try looking at the overlap and weights of points in the k-nearest
    %neighborhood of each vertex.
    

end

function[clusterDimEst,dimStart,dimEnd,dimList]=estDimClusterPCA3(gaussMVarExpl,T,data,dataSCCoeffs,maxPCs,numSigmas,c)
for clusterIdx = 1:max(T)
    dataInCluster = data(T==clusterIdx,:);
   clusterIdx
   dataInCluster
    dimEstNum = getDim(dataInCluster,false,c)
   clusterDimEst(clusterIdx) = dimEstNum;
    %now, WHICH PCS?
    varAmt = [];
    for i = 1:maxPCs
    
        varAmt(i) = std(dataInCluster(:,i));
    end
    [~,varIdx]=sort(varAmt);
    dimStart(clusterIdx) = varIdx(1);
    dimEnd(clusterIdx) = varIdx(clusterDimEst(clusterIdx));
    dimList{clusterIdx}=varIdx(1:clusterDimEst(clusterIdx));
end
end


function[dim]=getDim(data,dispFlag,c)
%%return the dimensionality of some data, using the gaussian-derived method
if ~exist('dispFlag','var')
    dispFlag = false;
end
disp('generating replicates...')

disp('computing pcs...')

%normalize
%data = data - repmat(min(data),size(data,1),1);
minimax = max(data)-min(data);
%data = data./repmat(minimax,size(data,1),1);
[~,~,~,~,dExpl]=pca(data,'Economy',true);
if isempty(dExpl)
    dExpl = zeros(1,(size(data,2)));
    dExpl(1)=1;
end
%generate some gaussian data of the same dimensionality.
%use 1000 replicates.
%expl = zeros(1000,min(size(data)));
%cumExpl = zeros(1000,min(size(data)));
datan=zeros(size(data));
for i = 1:1000
    for j = 1:size(data,2)
        %datan(:,j) = normrnd(0,1,[size(data,1),1]);
        %datan(:,j) = exprnd(mean(data(:,j)),[size(data,1),1]); %should this be data or datascores?
        %for RNASeq, can we mimic complexity?
       datan(:,j) = exprnd(mean(data(:,j)),[size(data,1) ,1]);
        nmin = min(datan(:,j));
        nmax = max(datan(:,j));
        szn = size(datan(:,j),1);
        %datan(:,j) = (datan(:,j) - repmat(nmin,szn,1))./(repmat(nmax,szn,1)-repmat(nmin,szn,1));
        %datan(:,j) = randn(size(data,1),1) + mean(data(:,j));
        
        %what if this data *is* the PC data?
        %then we need some coefficient matrix to see what the "real" data
        %in the noise model would look like, so that we can look at the 
        %variance explained.
        %we have now passed in c, the coefficients matrix.
        
        
    end
    datan(isnan(datan))=0; %this is actually happening.
    %datan = datan+randn(size(datan))*0.1; %10pct. noise estimate?
    
    
    %reconstruct the 'real' datan
    datanRNASp = datan*c';
    %assuming it is data, do pca (o.w look at just var. in comp.)
    [~,~,~,~,vExpl]=pca(datanRNASp,'Economy',true);
    if ~isempty(vExpl)
        expl(i,:) = vExpl;
    else
        expl(i,1) = 1;
        expl(i,2:end) = zeros(1,size(expl,2)-1);
    end
    for j = 1:size(expl,2)
      cumExpl(i,j) = sum(expl(i,1:j));
    end
end


if dispFlag
   figure(1);
   hold on;
   title('Variance Explanation Plot');
   boxplot(expl);
   plot(dExpl,'ko-');
   hold off;
end

ztarget = abs(norminv(1-(0.01/size(data,2)),0,1));
dim=1; %need to initialize.
if size(expl,2)>1
    
for j = 2:(size(expl,2)-1)
   if ( dExpl(j-1) <= (mean(expl(:,j-1)) - std(expl(:,j-1))) && ...
           dExpl(j) >= (mean(expl(:,j)) - std(expl(:,j))) && ...
           dExpl(j+1) >= ((mean(expl(:,j+1)) - std(expl(:,j+1)))) ) %|| ...
           %( dExpl(j) <=(mean(expl(:,j)) + std(expl(:,j))) && ...
           %dExpl(j) >= (mean(expl(:,j)) - std(expl(:,j))) ) %may need tuning
       break;
   end
end
dim = j-1;
end
end



function[E]=compute_energy_prob(V,A,ErrorsCell,gamma)
%errors cell already has the distance-based errors, and these are weighted.
Edist =0;
for i = 1:length(ErrorsCell)
    currCell = ErrorsCell{i};
    Edist = Edist+sum(currCell);
end

%we need to compute the MST, and the -logprob resultant from that.
distV =dist(V');
pathWeights = sparse(distV.*A);
MST = graphminspantree(pathWeights);
EMST = gamma*sum(reshape(MST,1,numel(MST)));
E = Edist + EMST;
end

