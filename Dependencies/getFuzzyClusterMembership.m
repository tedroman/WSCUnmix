function[M,distToReps]=getFuzzyClusterMembership(T,data,h)

%given the cluster membership (T), the data, and the number of Clusters,
%this code will return a numPoints * numClus matrix M of membership of each
%point in each of the clusters, as a mixture fraction (normalized to 1
%across all clusters).  The algorithm will work by first identifying the
%representative point for each cluster as the point farthest on average
%from all other clusters, using the kernel function that was first used for
%hte medoidshift clustering, which led to the stable clusters input as T
%here.  Then, the distance of each point to each of the representative
%points is computed, using that same kernel function.  Finally, the
%normalization of these distances is computed and returned as the output
%parameter M, which provides a way to weight membership in each of the
%clusters.  %h is the bandwidth used previously.

%cpyright 2015 Ted Roman, advised by Russell Schwartz, Carnegie Mellon
%Univ.
numPts = size(data,1);
numClus = max(T);
reps = zeros(1,numClus);
M = zeros(numPts,numClus);
D = dist(data'); %get distance.
K = getKernelDistance(D,h); %get kernel distance.
for i = 1:numClus
   pointsInCluster = find(T==i);
   pointsOutsideCluster = find(T~=i);
   kdInToOut = K(pointsInCluster,pointsOutsideCluster);
   meanToOut = mean(kdInToOut,2); %may need to check dimension.
   [~,rep]=min(meanToOut);
   repIdx = pointsInCluster(rep);
   reps(i) = repIdx;
end
distToReps = K(:,reps);
M = (1./distToReps)./repmat(sum(1./distToReps,2),1,numClus); %may need to check dimension.
M(isnan(M))=1; %correct for Nan


end

function[K]=getKernelDistance(D,h) %get the kernel function of the dist.
K = exp(-D/h)-1;
end
