repeat = 10;
noise = 0;
t_depth = 1; %trace depth of the 1st medoidshift.
hc1 = 1:10; %bandwidth coefficients for the two medoidshifts;
hc2 = 0; %actual bandwidth = coefficient * mean distance.
score = zeros(repeat,5,numel(hc1));

dataPts1 = 100;
dataPts2 = 100; %500 dp for each simplex.

%set up the vertices.
v1 = sparse(3,200);
v2 = sparse(3,200);
v1(2,1)=1;
v1(3,2)=1;
v2(2,3)=1;
v2(3,4)=1;

for replicate = 1:50
   dataForRep = getSimulatedData(dataPts1,200,'unif',v1);
   dataForRep2 = getSimulatedData(dataPts2,200,'unif',v2);
   noise = randn(dataPts1+dataPts2,200)*0.3;
   dataRep = [dataForRep;dataForRep2]+noise;
   disp('Generated Data');
   %randProjMatrix = randn(20000,200);
   %'project' the data.
  % projData = dataRep*randProjMatrix';
   for pcs = 1:100
       [datac,datas]=pca(projData);
       
       repsVec = doubleshift(dataRep(:,1:pcs),1,1,1);
       disp('Computed mshift');
       numClus(replicate,pcs)=length(unique(repsVec));
   
   end
end