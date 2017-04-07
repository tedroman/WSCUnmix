function[theta,thetasKNN,elapse,VProjCell,Acell,E_KNN,E_h,penPf_vec,penDf_vec,datasets,Vcell,Fcell]=batchParameterSensitivitySCUnmixv2(theta0,maxNumDim,scenarioNum,gamma,datasets)
%call the parameters all together theta.
%let theta(1) = maximum number of dimensions
%let theta(2) = number of bootstrapped replicates
%let theta(3) = threshold for dimensionality estimation
%let theta(4) = h (h1 and h2 being the same)
%let theta(5) = probability threshold for vertex merger
%let theta(6) = number of nearest neighbors (K) for KNN graph for vertex
%merger
%let theta(7) = maximum number of iterations(line and non-line the same)

%initialize variables
et = tic;
VProjCell = {};
Acell = {};

%initialize Theta:
theta(1,:) = theta0;

%generate the data:
if ~exist('datasets','var')
    if ~exist('scenarioNum','var')
        [datasets] = generate_data(); %need the gt out of these also
    else
        [datasets] = generate_data(scenarioNum);
    end
end

if ~exist('maxNumDim','var')
    maxNumDim = 10;
end
%initialize Gamma
if ~exist('gamma','var')
    gamma = 1e-3;
end
%gamma = 10;

%what if we messed up and gamma is supposed to be 0.1?
%gamma =0.1;
%gamma=1;
%gamma = 0.01; %very small--should prefer more complex models.


%compute the values of the distance matrix.
numPlaces = 0; %changed to 0 because of too many neighborhoods.
distVec = roundDistanceMatrix(datasets,numPlaces);
disp(['There are ' num2str(length(unique(distVec))) ' unique possible neighborhood sizes']);

%re-initialize theta.
theta0(4) = distVec(1);
theta0(1) = maxNumDim;

%create an appropriate size of theta.
theta = repmat(theta0,length(distVec),1);
for i = 1:length(distVec)
   theta(i,4) = distVec(i); %modify to the correct value. 
end
endOfFirstLoop = size(theta,1);
%make the computation as to the energy function for altering the h
%parameter.
parfor i = 1:endOfFirstLoop
%for i = 1:endOfFirstLoop
    %because these are unique, we must compute the energy for each of them.
    %compute the energy for the given value of h.  We can then create a
    %mapping of the energy for a given value of h (assuming default values
    %of the other parameters to the model)
    
        cT = theta(i,:);
        [V,A,data,VProj,M,I,Errors,K,CIp,Ps,Fs,RW,penPf_vec(i),penDf_vec(i)] = weightedSCUnmixMain(datasets,cT,gamma); %assume at first only 1 replicate of each scenario.
        E_h(i) = compute_energy_prob(V,A,Errors,gamma);
        disp(['Computed Energy number ' num2str(i)]);
        VProjCell{i}=VProj;
        Acell{i}=A;
        Vcell{i}=V;
        Fcell{i}=Fs;
end

%compute the max based on the heuristic for KNN
maxNumDataPoints = computeMaxNN(datasets);

initialThetaKNN = theta(end,:);
thetasKNN = repmat(initialThetaKNN,maxNumDataPoints,1);
thetasKNN(:,6) = reshape(1:maxNumDataPoints,maxNumDataPoints,1);


parfor i = 1:maxNumDataPoints
%for i = 1:maxNumDataPoints
    [V,A,data,VProj,M,I,Errors,K,CIp,Ps,Fs,RW,penPf_vec(i),penDf_vec(i)] = weightedSCUnmixMain(datasets,thetasKNN(i,:),gamma);
    %record the energy of the solution
    E_KNN(i) = compute_energy_prob(V,A,Errors,gamma);
    VProjCell{endOfFirstLoop+i}=VProj;
    Acell{endOfFirstLoop+i}=A;
    Vcell{endOfFirstLoop+i}=V;
    Fcell{endOfFirstLoop+i}=Fs;
end
elapse = toc(et);
end



function[E] = compute_energy(V,A,gtV,gtA)
%let the energy be defined as the sum of the error from the vertex
%resolution and the error from the misses in the adjacency matrix.

%for the purposes of computation, we want to find the minimum matching from
%the inferred vertices to the ground truth vertices, then apply the sorting
%required to make the transformation to the adjacency matrix, to have as
%apples-to-appls of a comparison as possible.

%to enumerate all options is exponentially expensive (combinatorics
%problem)

%we want to consider both whether the number of vertices inferred is
%correct, and, additionally, the closeness of those vertices to the ground
%truth vertices.

%a la the previous work, the error per dimension per ground truth model
%component gives a percentage error like measurement.  however, we could
%also compute the error per dimension *  ( 1 + |inferred components - gt components|)
%to further penalize for missing the number of components.  Then we would
%only need to look at what the nearest gt component is, rather than trying
%to solve an exponential problem.

%we would then apply the same sorting to the adjacency matrix to compute
%the number of FP and FN in the adjacency matrix, adding 0's if needed to
%ensure the dimensions align.

%compute the (euclidean) distance between the gt and inferred vertices.

%gtA needs to be corrected?? no...
%V needs to be projected.
%this should be Vproj.

distanceMtrx = dist(V,gtV'); %make sure the dimensionality is correct.
[minDs,minDIdx] = min(distanceMtrx,[],2); %make sure the sorting is done in the correct way.
%minDs will tell us the errors., size of V tells us the 'true' dims,
%minDIdx is used later.

totalD = sum(minDs);
numDims = size(gtV,2); %number of true dimensions (if not, then this coudl be passed as a paramter)
errorInNumComponents = abs(size(V,1)-size(gtV,1));

normalizedVError = totalD/numDims *(1+errorInNumComponents);

%minDIdx tells us which ones we tried to mirror, so how did we do on the
%connectivity of those--because we mayn't have mapped to some gt vertices,
%we are punished in the form of errorin num components, so we mayn't need a
%second punishment here.

sortedA = A(minDIdx,minDIdx); 
if size(sortedA,1)~=size(gtA,1)
    %error in number inferred.
    if size(sortedA,1) < size(gtA,1)
       %sortedA is smaller
       tmp = zeros(size(gtA,1),size(gtA,2));
       tmp(1:size(sortedA,1),1:size(sortedA,2)) = sortedA;
       sortedA = tmp;
       clear tmp;
    else
        tmp = zeros(size(sortedA,1),size(sortedA,2));
        tmp(1:size(gtA,1),1:size(gtA,2)) = gtA;
        gtA = tmp;
        clear tmp;
    end
end
    mismatched = abs(sortedA-gtA);
    numMissedA = sum(sum(mismatched));
    
    E = numMissedA + normalizedVError;
end

function[step] = compute_step()
%let the step size be defined as [1 100 1 0.1 0.01 1 10].
%stepSize = [1 100 1 0.1 0.01 1 10];
stepSize = 0.1;
randomMove = rand(1,length(stepSize));
%we have 2 possibilities at each of the parameters as adjacent neighbor
%states--forward or backward, so we can use 0.5 as a cutoff.

%TRY ONLY VARYING THE NUMBER OF h (neighborhood size, theta(4), new
%modification 3 Sep 2015
while (sum(randomMove==0.5)~=0)
    randomMove = rand(1,length(stepSize));
    randomMove = rand(1,1);
end %don't want a tie.
for i = 1:length(stepSize)
   if randomMove(i) > 0.5
      step = stepSize(i); 
   else
       step = -1*stepSize(i); %changed 3 september
   end
end
end


function[notComputed]=energyNotComputed(theta,thetaMap)
notComputed = true;
for i = 1:length(thetaMap)
   currentTheta = thetaMap{i};
        disp(['theta =      ' num2str(theta)])
        disp(['c. theta =   ' num2str(currentTheta)])

        if length(currentTheta)~=length(theta)
            %well, the can't be equal.. but why is the current theta not
            %correct in some sense?
           disp('Error in computing whether or not energy is computed');
           
        end
        if currentTheta==theta
            notComputed = false;
            break;
        end
end
end

function[E] = lookupEnergy(theta,thetaMap,Eall)
%if we are here, we know the energy has been computed.
%which row of E should we look at?
for i = 1:length(thetaMap)
   currentTheta = thetaMap{i};
   if currentTheta==theta
      E = Eall(i); 
      return;
   end
   
end
E = -1;
disp('Error Computing E');
end
function[E] = lookupEnergyAll(allThetas,theta,Eall)
idx = 0;
for i = 1:size(allThetas,1)
   if allThetas(i,:)==theta
       idx = i;
       break;
   end
end
if idx ==0
    disp('Error Computing Energy');
    E = max(Eall) * 10;
else
   E = Eall(idx); 
end
end
function[notComputed] = energyNotComputedAll(allThetas,theta)
allThetaUnique = unique(allThetas,'rows');
notComputed = true;
for i = 1:size(allThetaUnique,1)
    if allThetaUnique(i,:)==theta
        notComputed = false;
        break;
    end
end
end
function[isValid]=isValidState(theta)
%want to make sure that we are not disobeying any edge cases
%let theta(1) = maximum number of dimensions
%let theta(2) = number of bootstrapped replicates
%let theta(3) = threshold for dimensionality estimation
%let theta(4) = h (h1 and h2 being the same)
%we can set it up so that if theta(4) > max dist, theta(4) = max dist.
%let theta(5) = probability threshold for vertex merger
%let theta(6) = number of nearest neighbors (K) for KNN graph for vertex
%merger
%let theta(7) = maximum number of iterations(line and non-line the same)
isValid = true;
if floor(theta(1))~=theta(1)
    isValid = false;
elseif theta(1) <=0
    isValid = false;
elseif floor(theta(2))~=theta(2)
    isValid = false;
elseif theta(2)<=0
    isValid = false;
elseif theta(3) < 0
    isValid = false;
elseif theta(4) <= 0
    isValid = false;
elseif theta(5) < 0
    isValid = false;
elseif floor(theta(6))~=theta(6)
    isValid = false;
elseif theta(6) <= 0 
    isValid = false;
elseif floor(theta(7))~=theta(7)
    isValid = false;
elseif theta(7)<=0
    isValid = false;
end

%we could also add a bound for neighborhood sizes (theta(4)) > max dist
%should be collapsed to the max distance, as they will yield the same
%result.
end

function[E]=compute_energy_prob(V,A,ErrorsCell,gamma)
%errors cell already has the distance-based errors, and these are weighted.
if isempty(V)
    E = 1e10; %arbitrary very large energy for single simplex 
    return;
end

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


function[distVec]=roundDistanceMatrix(datasets,numPlaces)
%numPlaces is the number of places after the decimal point.  e.g. 1 =
%tenths place.
distMtrx = dist(datasets');
distVec = reshape(distMtrx,1,numel(distMtrx));
distVec = unique(distVec); %these are the values of the neighborhoods that make sense.
distVec = round(distVec,numPlaces); %round to the nearest tenth. (doint whole number for speed right now)
distVec = unique(distVec);
distVec = distVec(distVec > 0);
end

function[num]=computeMaxNN(datasets)
num = floor(sqrt(size(datasets,1))) + floor(log(size(datasets,1))); %This is a rule of thumb from KNN
end