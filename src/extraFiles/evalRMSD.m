function[RMSD]=evalRMSD(V,scenario)
%evaluate the RMSD of each scenario
numDim = 50; 
if (scenario==1)
     v1 = zeros(3,numDim);
     v2 = zeros(3,numDim);
     v1(2,1) =1;
     v1(3,2) =1;
     v2(2,1) = 1;
     v2(3,3) = 1;
     
elseif (scenario==2)
        v1 = zeros(3,numDim);
        v2 = zeros(3,numDim);
        v1(2,1) =1;
        v1(3,2) =1;
        v2(2,3) =1;
        v2(3,2) =1;
elseif (scenario==3)
        v1 = zeros(4,numDim);
        v2 = zeros(4,numDim);
        v1(2,1) =1;
        v1(3,2) =1;
        v1(4,3) =1;
        v2(2,4) = 1;
        v2(3,5) = 1;
        v2(4,6) = 1;
elseif (scenario==4)
    v1 = zeros(4,numDim);
        v2 = zeros(4,numDim);
        v1(2,1) =1;
        v1(3,2) =1;
        v1(4,3) =1;
        v2(2,1) = 1;
        v2(3,5) = 1;
        v2(4,6) = 1;
elseif (scenario==5)
    v1 = zeros(4,numDim);
        v2 = zeros(4,numDim);
        v1(2,1) =1;
        v1(3,2) =1;
        v1(4,3) =1;
        v2(2,1) = 1;
        v2(3,2) = 1;
        v2(4,6) = 1;
elseif (scenario==6)
    v1 = zeros(3,numDim);
        v2 = zeros(4,numDim);
        v1(2,1) =1;
        v1(3,2) =1;
        
        v2(2,1) = 1;
        v2(3,5) = 1;
        v2(4,6) = 1;
elseif (scenario==7)
    v1 = zeros(4,numDim);
        v2 = zeros(3,numDim);
        v1(2,1) =1;
        v1(3,2) =1;
        v1(4,3) =1;
        v2(2,4) = 1;
        v2(3,5) = 1;
end
gtV = [v1;v2];
dists = dist(gtV,V);
minDists = min(dists);
RMSD = sum(minDists);

end