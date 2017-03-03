function [dimList,clusterDimVec] = getDimListMain(data,theta,gamma,c,alphaVal,dispFlag,figIdx)
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

    disp('Multiple integers in interval')
    
    for i = 1:length(CIInt)
       currInt  = CIInt(i);
       [dimList,clusterDimVec] = getDimList(numRep,numPts,repIdxAllReplicates,currInt,repsBase,theta,dataSCCoeffs,dataSCBase01Norm,h,dataSCBase,data,dataPCCoeff,gamma,c,alphaVal,dispFlag,figIdx);
        modelsToCompare{i}.dimList = dimList;
        dimList=dimList; % can't compare models in this module.
    end
    
   
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
    dimList=0;
    return;
else
[dimList,clusterDimVec]=getDimList(numRep,numPts,repIdxAllReplicates,CIInt,repsBase,theta,dataSCCoeffs,dataSCBase01Norm,h,dataSCBase,data,dataPCCoeff,gamma,c,alphaVal,dispFlag,figIdx);
    
end

end



