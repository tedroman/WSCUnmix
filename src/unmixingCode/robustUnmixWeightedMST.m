function[C,a,F,P,errors,C0p,C0,K,CProj,reweights,penDf,penPf] = robustUnmixWeightedMST(data,numVert,gamma,weights,dataCoeff,dataScores,dReps, repIdx,startDim,endDim,dimList,maxIter)
%Inputs: data--data that may or may not be in PC space (numDim * numPts)
%numVert-number of vertices to find
%gamma-hyperparameter for soft unmixing
%dataCoeff-coefficients to transform data to pc space (optional)
%dataScores-scores of data in PC space (optional)
%weights- a numPts * 1 vector of the weight of considering each point in
%the vertex.

%we aim to minimize a cost function:
% min_V mst(V) + \gamma sum_i weights(i)*(|x_i - V*\alpha_i|_1)

if ~exist('dataCoeff','var') || ~exist('dataScores','var')
    [dataCoeff,dataScores]=pca(data','Centered',false,'Economy',true);
end
%will add tolerance later.
%tol = 1e-3; what was the purpose of tolerance?

k = numVert;
if k==2 %cannot use MVES with a line.
    [C,a,F,P,errors,C0p,C0,K,CProj,reweights] = unmix_line(dataScores,repIdx,weights,k,gamma,dataCoeff,startDim,endDim,dimList,maxIter);
  penPf = 0;
  penDf = 0; %ought to be corrected later (edge case)
  
    return;
end
   
  
n = size(dataScores,1);
V = dataCoeff;
S = dataScores; 
%This is to keep symmetry with previous methods.

%tracer.
if length(S)==2
    idV = V;
end

if size(V,2)==1
    %There is only 1 data point.
    P = V';
    C0p = V; %No compression is possible.
else
    P = S(1:n,1:k-1)'; %we should take from start to end instead.-should it even be start-end (supposing these are cts?)
    
    %which dimensions should we take?
    %length of dimList should be k-1
    
    %should just be taking the dim list, as these best explain the current
    %cluster.
    %try doing this the other way
    P = transpose(S(:,1:k-1));
    P = transpose(S(:,dimList)); %added 20Aug
    C0p = S(:,1:k-1);
    C0p = S(:,dimList); %added 20Aug
    %P is the transpose( scores of all data points, desired dimension PCs)
    %C0p is all data points, desired dimension PCs
    %Shoudl be P = C0p';
end
normalizer = max(sqrt(sum(P.^2,2)));
disp(['size(P): ' num2str(size(P)) 'size(normalizer): ' num2str(size(normalizer))]);
P = P/normalizer; %This is to scale in the scores, given they were previously unscaled.
%We should only take the P's for the fast unmix in which the weights are >
%0.5 (high confidence).--does this need to be modified?


weightCutoff = 1/(length(repIdx));

disp(strcat('Number of High-Confidence points: ',num2str(sum(weights>weightCutoff)))); %How do we know this is the correct number?



if sum(weights > weightCutoff) < 5 %this may need to be changed for > 2 clusters. (change made 13 Aug 2015)
    disp('Not enough HC points');
    
    npoints = size(dataScores,2);
   a =0;
   F = 0;
   P = 0;
   errors = 0;
   C0p = 0;
   C0 = 0;
   K = 0; 
   reweights = zeros(npoints,1);
   penPf=0;
   penDf = 0;
   %placeholders; will code the algorithm later.
   
   % is this why?
   
   %just do a hard unmix here b/c so few points?
   %{
   
   for i =1:size(dataScores,1)
    d2 = size(C0p,2); %num col of K'
   F0(:,i) = lsqlin(C0p,dataScores(i,:),-eye(d2),zeros(d2,1),ones(1,d2),1,[],[],[],opts);
   E0(i) = sum(abs(dataScores(i,:)'-C0p*F0(:,i)))*reweights(i);
end
pointErrors0 = sum(E0);
MST0 = graphminspantree(sparse(dist(C0p')));
objectError0 = gamma*sum(reshape(MST0,numel(MST0),1));
%log(gamma) - log(sum(reshape(MST0,numel(MST0),1)));
totalError0 = pointErrors0 + objectError0;
%now fit the penalized deviation
 K = numericalWeightedSimplexMST(C0p',dataScores',gamma,reweights,optset);


errors = zeros(1,size(P,2));
F = zeros(size(K,1),size(P,2));

   %}
   C = 0;
   CProj = [];
    return;
end


%determine the 'reweights'
%use a Gaussian model for the mean point of all of the representatives.

%this should be modularized.
reweights = reweightPoints(dataScores,repIdx,weights);




%May not be the number per se, but that high confidence points should be
%'overweight'--the gradient in confidence may be distorting simplices
%found.  Can we turn to theory?

PHC = transpose(dataScores(weights>weightCutoff,:));
%There may be issues with the gradient of points (not dense enough?)
%Because points are seeded such that some triangles are overly generous and
%others not so.


optset = optimset('Display','Final','MaxIter',maxIter);%(k-2)*1000); %why this many?

opts = optimset('Display','Final','MaxIter',maxIter);

[C0p,a,F,P,V0,ECs,PCs]=fast_unmix3(PHC,k); %only initial guess so okay.

for i =1:size(dataScores,1)
    d2 = size(C0p,2); %num col of K'
   F0(:,i) = lsqlin(C0p,dataScores(i,:),-eye(d2),zeros(d2,1),ones(1,d2),1,[],[],[],opts);
   E0(i) = sum(abs(dataScores(i,:)'-C0p*F0(:,i)))*reweights(i);
end
pointErrors0 = sum(E0);
MST0 = graphminspantree(sparse(dist(C0p')));
objectError0 = gamma*sum(reshape(MST0,numel(MST0),1));
%log(gamma) - log(sum(reshape(MST0,numel(MST0),1)));
totalError0 = pointErrors0 + objectError0;
%now fit the penalized deviation
 K = numericalWeightedSimplexMST(C0p',dataScores',gamma,reweights,optset);

opts = optimset('Display','iter','MaxIter',maxIter); %set same as other.
errors = zeros(1,size(P,2));
F = zeros(size(K,1),size(P,2));
warning off;

for i =1:size(dataScores,1)
    d = numel(dataScores(i,:))+1;
    d2 = size(K,1); %num col of K'
   % F(:,i) = lsqlin(K',P(:,i),-eye(d),zeros(d,1),ones(1,d),1,[],[],[],opts);
   F(:,i) = lsqlin(K',dataScores(i,:),-eye(d2),zeros(d2,1),ones(1,d2),1,[],[],[],opts);
   
   errors(i) = sum(abs(dataScores(i,:)'-K'*F(:,i)))*reweights(i);%sigmoid(weights(i)); %l1.
   
   %sqrt(sum((dataScores(i,:)'-K'*F(:,i)).^2));
end
pointErrors = sum(errors);
MSTFinal = graphminspantree(sparse(dist(K')));
%trying different prior.
objectErrorFinal = gamma*sum(reshape(MSTFinal,numel(MSTFinal),1));
%log(gamma)-log(sum(reshape(MSTFinal,numel(MSTFinal),1)));
penDf = pointErrors;
penPf = objectErrorFinal;
totalErrorFinal = pointErrors + objectErrorFinal;
%finalTotalErrors = pointErrors+gamma*sum(reshape(MSTFinal,numel(MSTFinal),1));

warning on;
F=F';


%What if I take this chunk out because of the new pca stuff?
%promote to the full space.
%C = dataCoeff(:,1:(k))*(K*normalizer);
%C=V(:,1:(k-1))*(K'*normalizer);
%C0 = dataCoeff(:,1:(k))*(C0p'*normalizer);
%C0=V(:,1:(k-1))*(C0p*normalizer);
%because these are already in PCA space.

C = K;
C0 = C0p';

%{
for i=1:k
 C(:,i)=C(:,i)+mean(data,2);
 C0(:,i)=C0(:,i)+mean(data,2); 
end
%}

CProj = transpose(dataCoeff(1:size(C,2),:))*C';


a = 0; %placeholder.


%simplex_vol(C); 

end

function[y]=sigmoid(x,mu)
y = (1+exp(-20*(x - mu)))^-1; %just see what happens with this--experiment.
end


function[density]=mvnpdfpath(dataScores,currIdx,currRepIdx,sigma,dToCenter)
l2SquareDists = dist(dataScores,dataScores').^2;

%can we change to a path measurement?
disp(['Computing path costs for point number' num2str(currIdx)])
l2SquarePathDists = graphallshortestpaths(sparse(l2SquareDists));
%maybe we should set the median distance to this cluster cetner in l2 apath
%as the sigma
%newSigma =sqrt(median(l2SquarePathDists(currRepIdx,:)));
npdf = (normpdf(l2SquarePathDists(currIdx,currRepIdx),0,dToCenter));
density = npdf;
end

function[C0p]=line_fast_unmix(PHC)
%just use the maxes and mins
initMin = min(PHC,[],2);
initMax = max(PHC,[],2);
%construct a set of points of the mins and maxes
C0p =[initMin';initMax'];

end


function[C,a,F,P,errors,C0p,C0,K,CProj,reweights]=unmix_line(dataScores,repIdx,weights,k,gamma,dataCoeff,dimStart,dimEnd,dimList,maxIter)
    weightCutoff = 1/(length(repIdx));
    PHC = transpose(dataScores(weights>weightCutoff,:));
   
    %determine the 'reweights'
    %use a Gaussian model for the mean point of all of the representatives.
    
    reweights = reweightPoints(dataScores,repIdx,weights);
    
    
    [C0p] = line_fast_unmix(PHC); % F is unused, %uses ALL the dimensions
    C0p=C0p';
    a = 0; %placeholder.
    PHCt = transpose(PHC);
    P = transpose((PHCt(:,1:k-1))); %takes the first k components (should be 2)--we really want the 2 dimList components
    P = transpose((PHCt(:,dimList))); %added 20Aug.
    %V0 is the low-d endpoints,
    %C0p is projected out.
    %PCs is the first however many pcs (scores or data in pc
    %space)--unused, dont' worry about them now.
    
    
    opts = optimset('Display','iter','MaxIter',maxIter); %why is this taking so long?
    optset = opts; %will fix later.
    
    
    for i =1:size(dataScores,1)
         d2 = size(C0p,2); %num col of K'
         F0(:,i) = lsqlin(C0p,dataScores(i,:),-eye(d2),zeros(d2,1),ones(1,d2),1,[],[],[],opts);
         E0(i) = sum(abs(dataScores(i,:)'-C0p*F0(:,i)))*reweights(i);
    end   
    pointErrors0 = sum(E0);
    MST0 = graphminspantree(sparse(dist(C0p')));
    objectError0 = gamma*sum(reshape(MST0,numel(MST0),1));
    totalError0 = pointErrors0 + objectError0;
    %now fit the penalized deviation
    K = numericalWeightedSimplexMST(C0p',dataScores',gamma,reweights,opts);

    
    F = zeros(size(K,1),size(P,2));
    
    
    for i =1:size(dataScores,1)
        d = numel(dataScores(i,:))+1;
        d2 = size(K,1); %num col of K'
        F(:,i) = lsqlin(K',dataScores(i,:),-eye(d2),zeros(d2,1),ones(1,d2),1,[],[],[],opts);
        errors(i) = sum(abs(dataScores(i,:)'-K'*F(:,i)))*reweights(i);%sigmoid(weights(i)); %l1.
    end
    
    %finish up last assignments.

    C = K;
    C0 = C0p';
    CProj = transpose(dataCoeff(1:size(C,2),:))*C';   
   

end

function[reweights]=reweightPoints(dataScores,repIdx,weights)
meanOfAllClusters = mean(dataScores(repIdx,:));
%may need to check to see if transpose is needed.
dToClusCenter = dist(dataScores(repIdx(1),:),meanOfAllClusters'); %this is the variance/std.
currentRepIdx = find(weights==1);

for q = 1:size(weights,1)
    for r = 1:size(weights,2)
        rt = sqrt(dToClusCenter);
        sigma = diag(repmat(rt,size(dataScores,2),1));
        reweights(q,r) = mvnpdfpath(dataScores,q,currentRepIdx,sigma,dToClusCenter);
    end
end
reweights = reweights - min(reweights);
reweights = reweights / max(reweights); %to normalize.-this is still too powerful.
end
