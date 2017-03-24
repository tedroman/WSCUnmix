function[Vnew,Anew]=determine_min_sc(V,A,s)
%try a visualization mode:
baseSz = 300;
baseLn = 5;
step = 1; %every step, decrease the size
%use highlighting to show which points are being merged.
hold on;
%scatter3(s(:,1),s(:,2),s(:,3),'go','filled');
%first, is it a simplicial complex?

isSC = is_SC_A(A);
Vtmp=[];
Atmp=[];
distM = dist(V'); %how is best to encode fit? (PseudoMST + fit to surface given merging?)
%what would a pseudomst look like?

%what if we look at cliques, then fit data to each.
%use the distance, but the min distance to each clique.
%for the mst
%take the mst
%then add min distance to connect all cliquees.

%is this viable?
%get all the combinations
combos = nchoosek(1:size(A,1),2); %get all pairs

%trim the combos that are not "legal"
combosTmp = [];
for i = 1:size(combos,1)
    guess = combos(i,:);
    if A(guess(1),guess(2))==0
        combosTmp=[combosTmp;guess];
    end
end
combos = combosTmp;
        
%getMinIdx = minIdx(distM,A);
%instead of a min idx, we need to try them all.

%need to generate a lot of candidates, then iterate over these.
candV = cell(1,size(combos,1));
candA = cell(1,size(combos,1));
likelihoodVec = zeros(1,size(combos,1));

%try each combination:
for i = 1:size(combos,1)
   getMinIdx = combos(i,:);
   [Vcand,Acand]=getNewBody(getMinIdx,V,A);
   %what is the likelihood?
   candV{i}=Vcand;
   candA{i}=Acand;
  likelihoodVal = getLikelihood(Vcand,Acand,s);
   likelihoodVec(i) = likelihoodVal;
end

%which of these is the most likely?
[~,sortIdx]=sort(likelihoodVec);
bestCandA = candA{sortIdx(1)};
bestCandV = candV{sortIdx(1)};
getMinIdx = combos(sortIdx(1),:);
Atmp = bestCandA;
Vtmp = bestCandV;

%graph:
%scatter3(Vtmp(:,1),Vtmp(:,2),Vtmp(:,3),baseSz/step,'ko','filled');
%scatter3(Vtmp(getMinIdx,1),Vtmp(getMinIdx,2),Vtmp(getMinIdx,3),baseSz/step,'ro','filled');
%for i = 1:size(Atmp,1)
%    for j = 1:size(Atmp,2)
%        if Atmp(i,j)
%            line([Vtmp(i,1) Vtmp(j,1)],...
%                [Vtmp(i,2) Vtmp(j,2)],...
%                [Vtmp(i,3) Vtmp(j,3)],'linewidth',...
%                baseLn/step);
%        end
%    end
%end
isSC = is_SC_A(Atmp); %update is_SC
while isSC > 1
    step = step+1;
    %find the minimum distance between 2 unique points
    %then generate a candidate. NO--we are changing this to likelihoods...
    
    %generate all legal candidates, then choose the best one.
    combos = nchoosek(1:size(Atmp,1),2); %get all pairs

    %trim the combos that are not "legal"
    combosTmp = [];
    for i = 1:size(combos,1)
        guess = combos(i,:);
        if A(guess(1),guess(2))==0
            combosTmp=[combosTmp;guess];
        end
    end
    combos = combosTmp;

    %need to generate a lot of candidates, then iterate over these.
    candV = cell(1,size(combos,1));
    candA = cell(1,size(combos,1));
    likelihoodVec = zeros(1,size(combos,1));

    %try each combination:
    for i = 1:size(combos,1)
    getMinIdx = combos(i,:);
    [Vcand,Acand]=getNewBody(getMinIdx,Vtmp,Atmp);
    %what is the likelihood?
    candV{i}=Vcand;
    candA{i}=Acand;
    likelihoodVal = getLikelihood(Vcand,Acand,s);
    likelihoodVec(i) = likelihoodVal;
    end
    
    disp(['Likelihood of candidates: ' num2str(likelihoodVec)]);
    for i = 1:length(candA)
        disp(['Cand A: \n']);
        disp(num2str(candA{i}));
    end
    
    %never gets to this step.
    
    
    %which of these is the most likely?
    likelihoodVec(isinf(likelihoodVec))=max(likelihoodVec)+1e5; %overflow prevention.
    [~,sortIdx]=sort(likelihoodVec);
    bestCandA = candA{sortIdx(1)};
    bestCandV = candV{sortIdx(1)};
    getMinIdx = combos(sortIdx(1),:);
    Atmp = bestCandA;
    Vtmp = bestCandV;
    
    %[Vtmp,Atmp] = getNewBody(getMinIdx,Vtmp,Atmp);
    %update SC status
    isSC = is_SC_A(Atmp);
    distM = dist(Vtmp');
    
    %graph:
    scatter3(Vtmp(:,1),Vtmp(:,2),Vtmp(:,3),baseSz/step,'ko','filled');
    %getMinIdx = minIdx(distM,Atmp);
   % scatter3(Vtmp(getMinIdx,1),Vtmp(getMinIdx,2),Vtmp(getMinIdx,3),baseSz/step,'ro','filled');
    for i = 1:size(Atmp,1)
        for j = 1:size(Atmp,2)
            if Atmp(i,j)
             line([Vtmp(i,1) Vtmp(j,1)],...
                   [Vtmp(i,2) Vtmp(j,2)],...
                  [Vtmp(i,3) Vtmp(j,3)],'linewidth',...
                  baseLn/step);
            end
        end
    end

end
Vnew = Vtmp;
Anew = Atmp;
end

%function to test if A is representative of a simplicial complex.
%what does it mean to be a simplicial complex?
function[isSC]=is_SC_A(A)
%if there is 1 connected component, it is a simplicial complex.
isSC  = graphconncomp(sparse(A));
end

function[Vcand,Acand]=getNewBody(idx,V,A)
%goal: merge the two points in V specified by idx,
%and return the new vertex and adjacency matrix.
meanPoint = mean([V(idx(1),:);V(idx(2),:)]); %may need the transpose
%determine which points the new point should be connected to (oldidx)
newAVec = A(idx(1),:)+A(idx(2),:) > 0;
%make a map from old to new.
oldToNew = zeros(1,size(A,1));
mapIdxCounter = 1;
for i=1:size(A,1)
    if i~=idx(1)
        
        if i~=idx(2)
        %i passes through
        oldToNew(i) = mapIdxCounter;
        mapIdxCounter = mapIdxCounter + 1;
        else
            %i==idx(2)
            oldToNew(i) =0;
        end
    else
       %i==idx(1) 
       oldToNew(i) = 0;
    end
end
Acandtmp = zeros(size(A,1)+1,size(A,2)+1);
for i = 1:length(oldToNew)
   if oldToNew(i)~=0
       Vcand(oldToNew(i),:) = V(i,:);
   end
end
Vcand(size(A,1)-1,:) = meanPoint;
for i = 1:length(oldToNew)
    for j = 1:length(oldToNew)
        if oldToNew(i)~=0 && oldToNew(j)~=0
            Acandtmp(i,j) = A(i,j);
        end
        
    end
end
Acandtmp(size(A,1)+1,1:size(A,2)) = newAVec;
Acandtmp(1:size(A,1),size(A,2)+1) = newAVec;
Acandtmp(end,end) = 1;
Acandtmp(idx,end) =0;
Acandtmp(end,idx) =0;

%now get rid of rows and cols that are just 0's.
Atmp2 = [];%zeros(size(A,1)+1,size(A,2)-1);
for i = 1:size(Acandtmp,2)
   if sum(Acandtmp(:,i))>0
       Atmp2 =[Atmp2 Acandtmp(:,i)];
   end
end
Atmp3 = [];
for i = 1:size(Atmp2,1)
    if sum(Atmp2(i,:)) > 0
        Atmp3 = [Atmp3;Atmp2(i,:)];
    end
end
Acand = Atmp3;

end

function[LL]=getLikelihood(Vcand,Acand,s)
%return the value proportional to the negative log likelihood function
%for each subsimplex, for each data point, compute the distance, then the
%mst.
LL = 0;
gamma = 1.1;
opts = optimset('Display','Final');
%standard gamma value.

Acandtmp = Acand;
for i = 1:size(Acandtmp,1)
    Acandtmp(i,i) = 0;
end
%prerequisite forMC 
MC = maximalCliques(Acandtmp);
for i = 1:size(MC,2)
   subSimplexV = Vcand(MC(:,i)==1,:);
   for j = 1:size(s,1)
       d = size(subSimplexV,1);
       d2 = d;
       V = subSimplexV;
       datum = s(j,:);
      F = lsqlin(subSimplexV',s(j,:),-eye(d),zeros(d,1),ones(1,d),1,[],[],[],opts);
      %alpha =  lsqlin(V', datum, -eye(d2), zeros(d2,1), ones(1,d2),1,[],[],[], opts);
%     [~,penDf] = recompute_penalty_function(Vcand,Acand,ones(size(s,1),size(MC,2)),s);
      distL1(j,i) = sum(abs(s(j,:)' - subSimplexV'*F));
   end
   minSpan = graphminspantree(sparse(dist(subSimplexV')));
   minSpanWt = sum(reshape(minSpan,numel(minSpan),1));
   priorError = gamma*log(minSpanWt);
   LL = LL+priorError;
end

ratioMtrx = distL1./repmat(sum(distL1')',1,size(MC,2));
weightedDists = distL1.*ratioMtrx;
for i = 1:size(weightedDists,1)
    for j = 1:size(weightedDists,2)
        LL = LL + weightedDists(i,j); %the conditional part.
    end
end

end
