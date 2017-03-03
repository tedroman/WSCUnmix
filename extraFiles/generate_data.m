function[data,gtV,gtA]=generate_data(scenario)

%for now, we just generate two triangles at a point, but we will modify
%this in the future to allow for multiple scenarios and repliates, to
%better evaluate.  We also assume 10 % noise, which should be changable in
%the future.
if ~exist('scenario','var')
    scenario=1;
end
    numDim = 50; 
    numPts = 250; %changed for temp.
    nPts1 = floor(numPts/2);
    
    
    %build the datasets in correspondence with the 7 scenarios from the
    %Apbc paper (2015-2016)
    
    v1=[];
    v2 =[];
    gtA=[];
    gtV=[];
    if(scenario==1) %two triangles at an edge
     v1 = zeros(3,numDim);
     v2 = zeros(3,numDim);
     v1(2,1) =1;
     v1(3,2) =1;
     v2(2,1) = 1;
     v2(3,3) = 1;
     %shares both of the first components (edge).
     dataS1 = getSimulatedData(nPts1,...,
         numDim,'unif',v1);
     dataS2 = getSimulatedData(numPts-nPts1,...,
         numDim,'unif',v2);
     dataSCBase = [dataS1;dataS2];
     data = dataSCBase+randn(numPts,numDim)*0.1;
     gtA = [1 1 1 1 1 0;...
            1 1 1 1 1 1;...
            1 1 1 1 1 0;...
            1 1 0 1 1 1;...
            1 1 1 1 1 1;...
            1 1 0 1 1 1];
        
    elseif (scenario==2) %two triangles at a point
        v1 = zeros(3,numDim);
        v2 = zeros(3,numDim);
        v1(2,1) =1;
        v1(3,2) =1;
        v2(2,3) =1;
        v2(3,2) =1;
    
        dataS1 = getSimulatedData(nPts1,numDim,'unif',v1);
        dataS2 = getSimulatedData(numPts-nPts1,numDim,'unif',v2);
    
        dataSCBase = [dataS1;dataS2];
    
        %optional noise
        data = dataSCBase+randn(numPts,numDim)*0.1;
      
       gtA = [1 1 1 1 0 0;...
              1 1 1 1 0 0;...
              1 1 1 1 0 0;...
              0 0 1 1 1 1;...
              0 0 1 1 1 1;...
              0 0 1 1 1 1];
       %ground truths are not used right now.
    elseif (scenario==3) % two tetrahedra at a point
        v1 = zeros(4,numDim);
        v2 = zeros(4,numDim);
        v1(2,1) =1;
        v1(3,2) =1;
        v1(4,3) =1;
        v2(2,4) = 1;
        v2(3,5) = 1;
        v2(4,6) = 1;
     %shares origin
     gtA = [1 1 1 1 1 0 0 0;...
            1 1 1 1 1 0 0 0;...
            1 1 1 1 1 0 0 0;...
            1 1 1 1 1 0 0 0;...
            1 0 0 0 1 1 1 1;...
            1 0 0 0 1 1 1 1;...
            1 0 0 0 1 1 1 1;...
            1 0 0 0 1 1 1 1];
        
     dataS1 = getSimulatedData(nPts1,...,
         numDim,'unif',v1);
     dataS2 = getSimulatedData(numPts-nPts1,...,
         numDim,'unif',v2);
     dataSCBase = [dataS1;dataS2];
     data = dataSCBase+randn(numPts,numDim)*0.1;
        
        
    elseif (scenario==4) % two tetrahedra at an edge
        v1 = zeros(4,numDim);
        v2 = zeros(4,numDim);
        v1(2,1) =1;
        v1(3,2) =1;
        v1(4,3) =1;
        v2(2,1) = 1;
        v2(3,5) = 1;
        v2(4,6) = 1;
     %shares origin->2,1
     
     gtA=[1 1 1 1 1 1 0 0;...
          1 1 1 1 1 1 0 0;...
          1 1 1 1 1 1 0 0;...
          1 1 1 1 1 1 0 0;...
          1 1 0 0 1 1 1 1;...
          1 1 0 0 1 1 1 1;...
          1 1 0 0 1 1 1 1;...
          1 1 0 0 1 1 1 1];
      
     dataS1 = getSimulatedData(nPts1,...,
         numDim,'unif',v1);
     dataS2 = getSimulatedData(numPts-nPts1,...,
         numDim,'unif',v2);
     dataSCBase = [dataS1;dataS2];
     data = dataSCBase+randn(numPts,numDim)*0.1;
    elseif (scenario==5) % two tetrahedra at a face
        v1 = zeros(4,numDim);
        v2 = zeros(4,numDim);
        v1(2,1) =1;
        v1(3,2) =1;
        v1(4,3) =1;
        v2(2,1) = 1;
        v2(3,2) = 1;
        v2(4,6) = 1;
     %shares origin->2,1->3,2
     
     gtA = [1 1 1 1 1 1 1 0;...
            1 1 1 1 1 1 1 0;...
            1 1 1 1 1 1 1 0;...
            1 1 1 1 1 1 1 0;...
            1 1 1 0 1 1 1 1;...
            1 1 1 0 1 1 1 1;...
            1 1 1 0 1 1 1 1;...
            1 1 1 0 1 1 1 1];
     dataS1 = getSimulatedData(nPts1,...,
         numDim,'unif',v1);
     dataS2 = getSimulatedData(numPts-nPts1,...,
         numDim,'unif',v2);
     dataSCBase = [dataS1;dataS2];
     data = dataSCBase+randn(numPts,numDim)*0.1;
    elseif (scenario==6) % a triangle and a tetrahedron at an edge
        v1 = zeros(3,numDim);
        v2 = zeros(4,numDim);
        v1(2,1) =1;
        v1(3,2) =1;
        
        v2(2,1) = 1;
        v2(3,5) = 1;
        v2(4,6) = 1;
     %shares origin->2,1
     
     gtA = [1 1 1 1 1 0 0;...
            1 1 1 1 1 0 0;...
            1 1 1 1 1 0 0;...
            1 1 0 1 1 1 1;...
            1 1 0 1 1 1 1;...
            1 1 0 1 1 1 1;...
            1 1 0 1 1 1 1];
        
     dataS1 = getSimulatedData(nPts1,...,
         numDim,'unif',v1);
     dataS2 = getSimulatedData(numPts-nPts1,...,
         numDim,'unif',v2);
     dataSCBase = [dataS1;dataS2];
     data = dataSCBase+randn(numPts,numDim)*0.1;
    else % a triangle and a tetrahedron at a point.
        v1 = zeros(4,numDim);
        v2 = zeros(3,numDim);
        v1(2,1) =1;
        v1(3,2) =1;
        v1(4,3) =1;
        v2(2,4) = 1;
        v2(3,5) = 1;
     %shares origin
     gtA = [1 1 1 1 1 0 0;...
            1 1 1 1 1 0 0;...
            1 1 1 1 1 0 0;...
            1 1 1 1 1 0 0;...
            1 0 0 0 1 1 1;...
            1 0 0 0 1 1 1;...
            1 0 0 0 1 1 1];
        
     dataS1 = getSimulatedData(nPts1,...,
         numDim,'unif',v1);
     dataS2 = getSimulatedData(numPts-nPts1,...,
         numDim,'unif',v2);
     dataSCBase = [dataS1;dataS2];
     data = dataSCBase+randn(numPts,numDim)*0.1;
        
 
    end
gtV = [v1;v2];
    
end