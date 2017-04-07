function [logodds,minEIdx]=computeLogOdds(E,data,gamma,scenario)


%energy function of model
%this is pre-computed.
%how?
[minE,minEIdx]=min(E); %which was the lowest energy of the
%candidate models

%energy function of ground truth
%need to compute the weighted distance-based errors.
At = [];
Vt = [];
numDim = size(data,2);
sdims=[];
SVert={};
if scenario==1
        sdims= [2 2];
        At = [1 1 1 1 1 0;...
            1 1 1 1 1 1;...
            1 1 1 1 1 0;...
            1 1 0 1 1 1;...
            1 1 1 1 1 1;...
            1 1 0 1 1 1];
        v1 = zeros(3,numDim);
     v2 = zeros(3,numDim);
     v1(2,1) =1;
     v1(3,2) =1;
     v2(2,1) = 1;
     v2(3,3) = 1;
        Vt = [v1;v2];
        SVert{1}=Vt(1:3,:);
        SVert{2}=Vt(4:6,:);
    elseif scenario ==2
             At = [1 1 1 1 0 0;...
              1 1 1 1 0 0;...
              1 1 1 1 0 0;...
              0 0 1 1 1 1;...
              0 0 1 1 1 1;...
              0 0 1 1 1 1];
          
          
           v1 = zeros(3,numDim);
        v2 = zeros(3,numDim);
        v1(2,1) =1;
        v1(3,2) =1;
        v2(2,3) =1;
        v2(3,2) =1;
          
        
        Vt = [v1;v2];
        sdims = [2 2];
        SVert{1}=Vt(1:3,:);
        SVert{2}=Vt(4:6,:);
    elseif scenario ==3
        At = [1 1 1 1 1 0 0 0;...
            1 1 1 1 1 0 0 0;...
            1 1 1 1 1 0 0 0;...
            1 1 1 1 1 0 0 0;...
            1 0 0 0 1 1 1 1;...
            1 0 0 0 1 1 1 1;...
            1 0 0 0 1 1 1 1;...
            1 0 0 0 1 1 1 1];
        
        v1 = zeros(4,numDim);
        v2 = zeros(4,numDim);
        v1(2,1) =1;
        v1(3,2) =1;
        v1(4,3) =1;
        v2(2,4) = 1;
        v2(3,5) = 1;
        v2(4,6) = 1;
        
        Vt = [v1;v2];
        sdims = [3 3];
        SVert{1}=Vt(1:4,:);
        SVert{2}=Vt(5:8,:);
    elseif scenario==4
        At=[1 1 1 1 1 1 0 0;...
          1 1 1 1 1 1 0 0;...
          1 1 1 1 1 1 0 0;...
          1 1 1 1 1 1 0 0;...
          1 1 0 0 1 1 1 1;...
          1 1 0 0 1 1 1 1;...
          1 1 0 0 1 1 1 1;...
          1 1 0 0 1 1 1 1];
        
      v1 = zeros(4,numDim);
        v2 = zeros(4,numDim);
        v1(2,1) =1;
        v1(3,2) =1;
        v1(4,3) =1;
        v2(2,1) = 1;
        v2(3,5) = 1;
        v2(4,6) = 1;
      
        Vt = [v1;v2];
        sdims = [3 3];
        SVert{1}=Vt(1:4,:);
        SVert{2}=Vt(5:8,:);
    elseif scenario==5
        
        At = [1 1 1 1 1 1 1 0;...
            1 1 1 1 1 1 1 0;...
            1 1 1 1 1 1 1 0;...
            1 1 1 1 1 1 1 0;...
            1 1 1 0 1 1 1 1;...
            1 1 1 0 1 1 1 1;...
            1 1 1 0 1 1 1 1;...
            1 1 1 0 1 1 1 1];
        v1 = zeros(4,numDim);
        v2 = zeros(4,numDim);
        v1(2,1) =1;
        v1(3,2) =1;
        v1(4,3) =1;
        v2(2,1) = 1;
        v2(3,2) = 1;
        v2(4,6) = 1;
        
        Vt = [v1;v2];
        sdims = [3 3];
        SVert{1}=Vt(1:4,:);
        SVert{2}=Vt(5:8,:);
    elseif scenario==6
        At = [1 1 1 1 1 0 0;...
            1 1 1 1 1 0 0;...
            1 1 1 1 1 0 0;...
            1 1 0 1 1 1 1;...
            1 1 0 1 1 1 1;...
            1 1 0 1 1 1 1;...
            1 1 0 1 1 1 1];
        
          v1 = zeros(3,numDim);
        v2 = zeros(4,numDim);
        v1(2,1) =1;
        v1(3,2) =1;
        
        v2(2,1) = 1;
        v2(3,5) = 1;
        v2(4,6) = 1;
        
        Vt = [v1;v2];
        sdims = [2 3];
        SVert{1}=Vt(1:3,:);
        SVert{2}=Vt(4:7,:);
else
        
        At = [1 1 1 1 1 0 0;...
            1 1 1 1 1 0 0;...
            1 1 1 1 1 0 0;...
            1 1 1 1 1 0 0;...
            1 0 0 0 1 1 1;...
            1 0 0 0 1 1 1;...
            1 0 0 0 1 1 1];
        
        v1 = zeros(4,numDim);
        v2 = zeros(3,numDim);
        v1(2,1) =1;
        v1(3,2) =1;
        v1(4,3) =1;
        v2(2,4) = 1;
        v2(3,5) = 1;
        
        Vt = [v1;v2];
          sdims = [3 2];
        SVert{1}=Vt(1:4,:);
        SVert{2}=Vt(5:7,:);
end


%first step: compute the weights
%the weights are a challenge because there is no equivalent to a 
%shifted representative.
numPts = size(data,1);
wt_t  = zeros(numPts,length(SVert)); %to fill in

%here is the projecttion into pc space???
[dataz,datamu,datasigma]=zscore(data);
[datac,datas]=pca(dataz,'Economy',true);


%transform non-pca to pca space (by multiplication with coefficients)


for j = 1:length(SVert)

V_t_extrema = SVert{j};%maybe?
meanOfExtrema = mean(V_t_extrema);
distToMean = dist(V_t_extrema,meanOfExtrema'); %may need addtl '
stdVal = mean(distToMean);
for q = 1:size(data,1)
    %current rep is the nearest extrema
    distToExtrema = dist(data(q,:),V_t_extrema');
    [minD,minIdx]=min(distToExtrema);
    xtr = V_t_extrema(minIdx,:);
    ds = dist(data,xtr');
    [minD2,minIdx2]=min(ds);
wts(q,j) = mvnpdfpath(data,q,minIdx2,stdVal,minD);
end
%
maxwts = max(wts(:,j));
minwts = min(wts(:,j));
maxwtsr = repmat(maxwts,size(wts,1),1);
minwtsr = repmat(minwts,size(wts,1),1);
tnum = wts(:,j)-minwtsr;
tdenom = maxwtsr'-minwtsr';
wt_t(:,j) = (wts(:,j)-minwtsr)./transpose((maxwtsr'-minwtsr'));
%these also need to be per-simplex.
end
%second step: compute the distances
d_t = zeros(numPts,length(SVert)); %to fill in.
%this should be formulated as a least-squares problem.

%need PCs of the data (is that what data is?)

numPoints = size(data,1);
numSimplex = length(SVert);
F = zeros(size(Vt,1),numPts); %?
currDim = 1;
for j = 1:numSimplex
    dists = zeros(numPoints,numSimplex);
for i = 1:size(data,1)
   nDim = sdims(j);
   K = SVert{j}; %how to do this? 
   %(this should also depend on the GT and adj)
   %assume it is done
   opts = optimset('Display','Final','MaxIter',1000*max(sdims));
   %computing nDim and K is much easier with a scenario variable passed in.
    d2 = nDim+1;

    F0 = 1/d2*ones(d2,1);%[];% to fill in.
    display(['scenario= ' num2str(scenario)])

    disp(num2str(size(K)))
    disp(num2str(size(data)))
    disp(num2str(i))
    disp(num2str(d2))
    disp(num2str(size(F)))
    F(currDim:currDim+nDim,i) = lsqlin(K',data(i,:)',-eye(d2),zeros(d2,1),ones(1,d2),1,[],[],[],F0,opts);
   dists(i,j) = sum(abs(data(i,:)-transpose(K'*F(currDim:currDim+nDim,i))));
    currDim = currDim+d2+1;
end
d_t = dists;
end

%{
copy-pasted code segment
for i =1:size(dataScores,1)
        d = numel(dataScores(i,:))+1;
        d2 = size(K,1); %num col of K'
        F(:,i) = lsqlin(K',dataScores(i,:),-eye(d2),zeros(d2,1),ones(1,d2),1,[],[],[],opts);
        errors(i) = sum(abs(dataScores(i,:)'-K'*F(:,i)))*reweights(i);%sigmoid(weights(i)); %l1.
    end
    
%}


%third step: weight the distances by the weights

Er_t = d_t.*wt_t;
for i = 1:size(Er_t,2)
    Er_t_cell{i}=Er_t(:,i);
end

E_gt = compute_energy_prob(Vt,At,Er_t_cell,gamma);
logodds = E_gt - minE;
end

%copy-pasted from batch unmixer:
function[E]=compute_energy_prob(V,A,ErrorsCell,gamma)
%errors cell already has the distance-based errors, and these are weighted.
if isempty(V)
    E = 1e10; %arbitrary very large energy for single simplex 
    return;
end


%do we need to do a part of this in PC space?

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


function[density]=mvnpdfpath(dataScores,currIdx,currRepIdx,sigma,dToCenter)
l2SquareDists = dist(dataScores,dataScores').^2;

%can we change to a path measurement?

l2SquarePathDists = graphallshortestpaths(sparse(l2SquareDists));
%maybe we should set the median distance to this cluster cetner in l2 apath
%as the sigma
%newSigma =sqrt(median(l2SquarePathDists(currRepIdx,:)));
npdf = (normpdf(l2SquarePathDists(currIdx,currRepIdx),0,dToCenter));
density = npdf;
end
