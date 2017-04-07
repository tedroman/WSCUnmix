function[gtV]=getGtV(scenario,numDim)

gtV = [];
if scenario==1
        v1 = zeros(3,numDim);
     v2 = zeros(3,numDim);
     v1(2,1) =1;
     v1(3,2) =1;
     v2(2,1) = 1;
     v2(3,3) = 1;
        gtV = [v1;v2];
        %two triangles at a edge
elseif scenario==2
    v1 = zeros(3,numDim);
        v2 = zeros(3,numDim);
        v1(2,1) =1;
        v1(3,2) =1;
        v2(2,3) =1;
        v2(3,2) =1;
          %two triangles at piont
        
        gtV = [v1;v2];
elseif scenario==3
    v1 = zeros(4,numDim);
        v2 = zeros(4,numDim);
        v1(2,1) =1;
        v1(3,2) =1;
        v1(4,3) =1;
        v2(2,4) = 1;
        v2(3,5) = 1;
        v2(4,6) = 1;
        %two tet at point
        gtV = [v1;v2];
elseif scenario==4
     v1 = zeros(4,numDim);
        v2 = zeros(4,numDim);
        v1(2,1) =1;
        v1(3,2) =1;
        v1(4,3) =1;
        v2(2,1) = 1;
        v2(3,5) = 1;
        v2(4,6) = 1;
      %two tet edge.
        gtV = [v1;v2];
elseif scenario==5
    v1 = zeros(4,numDim);
        v2 = zeros(4,numDim);
        v1(2,1) =1;
        v1(3,2) =1;
        v1(4,3) =1;
        v2(2,1) = 1;
        v2(3,2) = 1;
        v2(4,6) = 1;
        %two tet face
        gtV = [v1;v2];
elseif scenario==6
     v1 = zeros(3,numDim);
        v2 = zeros(4,numDim);
        v1(2,1) =1;
        v1(3,2) =1;
        
        v2(2,1) = 1;
        v2(3,5) = 1;
        v2(4,6) = 1;
        %tri tet edge
        gtV = [v1;v2];
elseif scenario==7
    v1 = zeros(4,numDim);
        v2 = zeros(3,numDim);
        v1(2,1) =1;
        v1(3,2) =1;
        v1(4,3) =1;
        v2(2,4) = 1;
        v2(3,5) = 1;
        %tet tri point
        gtV = [v1;v2];
end
end