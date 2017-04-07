
function[theta,elapse,VProjCell,Acell,Echosen]=batchParameterSensitivitySCUnmix()
%call the parameters all together theta.
%let theta(1) = maximum number of dimensions
%let theta(2) = number of bootstrapped replicates
%let theta(3) = threshold for dimensionality estimation
%let theta(4) = h (h1 and h2 being the same)
%let theta(5) = probability threshold for vertex merger
%let theta(6) = number of nearest neighbors (K) for KNN graph for vertex
%merger
%let theta(7) = maximum number of iterations(line and non-line the same)

%let the step size be defined as [1 100 1 0.1 0.01 1 10].  Let the quality
%of result be the mean over 100 datasets for each of the 7 scenarios for
%vertex and adjacency matrix resolution.

%we can then say that the energy is in some sense the sum of the errors
%(RMSD / RMSE per component).  Then we cna apply the principles of
%metropolis sampling for some number of maximum moves, and we consider the
%number of times we are in each 'state' of the system.
allThetas = [];
allEs = [];
theta(1,:) = [10 500 2 1.2 0.01 15 100];

%try starting with only a few iterations to better debug.

%theta(1,:) = [10 500 2 1 0.01 15 10];

%this is the start state

%generate the data --- assume a noise level nearest to literature for real
%data.
[datasets,ground_truth_vertex,ground_truth_adjacency] = generate_data(); %need the gt out of these also
gamma = 10;
%for now, again, just assume 1 scenario--say two triangles at a point.
Eall = [];
thetaMap = {};
for i  = 1:1000%3 %100 "moves", so it's an even percentage %start with 5 to debug.
    %get the results of performing the unmixing on this
    %[results] = WSCUM(data,theta(i,:));
     %for scenario = 1:7 %7 scenarios--lets just assume 1 at first though.
     tstart = tic;
     
     %first, see if we have already computed this energy.
     %if energyNotComputed(theta(i,:),thetaMap)
      if energyNotComputedAll(allThetas,theta(i,:))
          [V,A,data,VProj,M,I,Errors,K,CIp,Ps,Fs,RW] = weightedSCUnmixMain(datasets,theta(i,:)); %assume at first only 1 replicate of each scenario.
         %E = compute_energy(VProj,A,ground_truth_vertex,ground_truth_adjacency); %where does the ground_truth come from?
         E = compute_energy_prob(V,A,Errors,gamma); %should be sufficient.
         %need to re-work the energy function, using the Errors (E) from
         %weighted SCUnmix main.
         thetaMap{i} = theta(i,:);
         Eall(i) = E;
         allThetas = [allThetas;theta(i,:)];
         allEs = [allEs E];
         
         disp(['Theta:      ' num2str(theta(i,:))])
         disp(['Energy:     ' num2str(E)])
     else
         %E = lookupEnergy(theta(i,:),thetaMap,Eall);
         E = lookupEnergyAll(allThetas,theta(i,:),allEs);
         disp(['Theta:      ' num2str(theta(i,:))])
         disp(['Energy:     ' num2str(E)])
         thetaMap{i} = theta(i,:);
     end
    %get the proposed next direction
    step = compute_step();
    newProposedT = theta(i,:) + [0 0 0 step 0 0 0];
    while ~isValidState(newProposedT)
       step = compute_step();
       disp('Invalid step; recomputing step');
       newProposedT = theta(i,:) + [0 0 0 step 0 0 0];
    end
    %what is the "new" E?
       
    %if energyNotComputed(newProposedT,thetaMap)
     if energyNotComputedAll(allThetas,newProposedT)
         [V2,A2,data2,VProj2,M2,I2,Errors2,K2,CIp2,Ps2,Fs2,RW2] = weightedSCUnmixMain(datasets,newProposedT); %assume at first only 1 replicate of each scenario.
       
        %need to recompute this.
        %E2 = compute_energy(VProj2,A2,ground_truth_vertex,ground_truth_adjacency);
        E2 = compute_energy_prob(V2,A2,Errors2,gamma);
        disp(['Theta:       ' num2str(newProposedT)])
        disp(['Energy:      ' num2str(E2)])
        allThetas = [allThetas;newProposedT];
        allEs = [allEs E2];
    else
        %E2 = lookupEnergy(newProposedT,thetaMap,Eall);
        E2 = lookupEnergyAll(allThetas,newProposedT,allEs);
        disp(['Theta:       ' num2str(newProposedT)])
        disp(['Energy:      ' num2str(E2)])
    end
        
    %if the new E < old E, move there with probability 1.
    if E2 < E
        theta(i+1,:) = newProposedT;
        VProjCell{i}=VProj2;
        Acell{i} = A2;
        Echosen(i) = E2;
    else
       randNum = rand(1);
       if randNum < exp(-1*(E2-E)) 
           theta(i+1,:) = newProposedT;
           VProjCell{i} = VProj2;
           Acell{i} = A2;
           Echosen(i) = E2;
       else %assumes kT = 1-->is this wrong?
           theta(i+1,:) = theta(i,:);
           VProjCell{i} = VProj;
           Acell{i} = A;
           Echosen(i) = E;
       end
    end
    %o.w. move there with some lesser probability.
    %theta(i+1,:) = theta(i,:)+step;
    disp(['Completed step number ' num2str(i)])
    elapse(i) = toc(tstart);
    disp(['Time taken for step: ' num2str(elapse(i))])
    save(strcat('tmp',num2str(i)));
end
end

function[data,gtV,gtA]=generate_data()

%for now, we just generate two triangles at a point, but we will modify
%this in the future to allow for multiple scenarios and repliates, to
%better evaluate.  We also assume 10 % noise, which should be changable in
%the future.

    numDim = 50; 
    v1 = zeros(3,numDim);%V needs corrected in order to 
    v2 = zeros(3,numDim);
    v1(2,1) =1;
    v1(3,2) =1;
    v2(2,3) =1;
    v2(3,2) =1;
    
    %both of these share the origin.
    numPts = 500;
    nPts1 = floor(numPts/2);
    
    dataS1 = getSimulatedData(nPts1,numDim,'unif',v1);
    dataS2 = getSimulatedData(numPts-nPts1,numDim,'unif',v2);
    
    dataSCBase = [dataS1;dataS2];
    
    %optional noise
    data = dataSCBase+randn(numPts,numDim)*0.1;
    gtV = [v1;v2];
    gtA = [1 1 1 1 0 0;...
           1 1 1 1 0 0;...
           1 1 1 1 0 0;...
           0 0 1 1 1 1;...
           0 0 1 1 1 1;...
           0 0 1 1 1 1];
   
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