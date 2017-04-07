function[]=runWSCUnmix(nonPCdataDir,flag1,thetaStr,...
    flag2,gammaStr,flag3,sigmaStr,flag4,...
    k_upperStr,flag5,filename,mode)
%how exactly is the data being passed in? -- needs to be strings
 tic;


%for now, assume these are required to be in order
if ~strcmp(flag1,'-theta')
    disp(...
        ['wrong flag order:'...
        'usage runWSCUnmix inputDir'...
        '-theta thetaVal -gamma gammaVal'...
        '-sigma sigmaVal -kupper kupperVal -output outputDir mode']);
    exit(1);
elseif ~strcmp(flag2,'-gamma')
    disp(...
        ['wrong flag order:'...
        'usage runWSCUnmix inputDir'...
        '-theta thetaVal -gamma gammaVal'...
        '-sigma sigmaVal -kupper kupperVal -output outputDir mode']);
    exit(1);
elseif ~strcmp(flag3,'-sigma')
        disp(...
        ['wrong flag order:'...
        'usage runWSCUnmix inputDir'...
        '-theta thetaVal -gamma gammaVal'...
        '-sigma sigmaVal -kupper kupperVal -output outputDir mode ']);
    exit(1);
elseif ~strcmp(flag4,'-kupper')
        disp(...
        ['wrong flag order:'...
        'usage runWSCUnmix inputDir'...
        '-theta thetaVal -gamma gammaVal'...
        '-sigma sigmaVal -kupper kupperVal -output outputDir mode ']);
    exit(1);
elseif ~strcmp(flag5,'-output')
        disp(...
        ['wrong flag order:'...
        'usage runWSCUnmix inputDir'...
        '-theta thetaVal -gamma gammaVal'...
        '-sigma sigmaVal -kupper kupperVal -output outputDir mode']);
    exit(1);
end

%if we get here we know the flags are all correct
theta = str2num(thetaStr); %assumes could be non-scalar 
%should double-check this.
%seems to work if comma-separated.

sigma = str2double(sigmaStr);

gamma = str2double(gammaStr);

k_upper = str2double(k_upperStr);

%might need some sort of parsing script or to leverage the previous parsing
%script.  probably should be distinct from this program.

%now get the nonPCData
if strcmp(mode,'-d')
[nonPCdata,ids,geneList] = parseTCGAv4Data(nonPCdataDir);
elseif strcmp(mode,'-r')
   [nonPCdata,ids,geneList] = parseTCGARNAv3(nonPCdataDir); 
elseif strcmp(mode,'-dr') || strcmp(mode,'-rd')
    disp('WARNING: Heterogeneous mode untested!');
    [nonPCdata,ids,geneList] = parseTCGADNARNA(nonPCdataDir);
    %exit(1); % this is only temporary
else
    disp(['mode flag not properly set']);
    exit(1);
end

%need to toss out the 'bad' data marked by 0's and o's (sentinels)
%create a vector of those to keep, and only keep the ones that are marked.
keepRows = ones(size(nonPCdata,1),1);
for i = 1:size(nonPCdata,1)
   if ids(i,1)=='o'
       keepRows(i)=0;
   elseif sum(any(nonPCdata(i,:)))==0
       keepRows(i)=0;
   end
end
nonPCdataTmp = nonPCdata(logical(...
    keepRows),:);
idstmp = ids(logical(keepRows),:);
nonPCdata = nonPCdataTmp;
ids = idstmp;
clearvars nonPCdataTmp idstmp;

if strcmp(mode,'-d')
    [dataCoeff,dataScore]=pca(10.^nonPCdata,'Economy',true);
elseif strcmp(mode,'-r')
    [dataCoeff,dataScore]=pca(2*nonPCdata,'Economy',true); %test to see if this is the correct normalization.
elseif strcmp(mode,'-dr') || strcmp(mode,'-rd')
    [dataCoeff,dataScore]=pca(nonPCdata,'Economy',true);
else
    disp('ERROR: mode flag improperly set');
    exit(1);    
end

%need to run a check for orphaned data points (i.e. orphaned clusters)
%here.

repsBase = runDoubleShift2(dataScore(:,1:k_upper));

uniqueReps = unique(repsBase);
keepIdx = ones(length(repsBase),1);

for i = 1:size(uniqueReps,2)
   if sum(repsBase==uniqueReps(i))<=k_upper %we cannot unmix this
       keepIdx(repsBase==uniqueReps(i))=0;
   end
   
end

dataScore = dataScore(logical(keepIdx),:);
minSum = size(dataScore,1);

repsBase = runDoubleShift2(dataScore(:,1:k_upper));
uniqueReps = unique(repsBase);
keepIdx = ones(length(repsBase),1);
for i = 1:size(uniqueReps,2)
    if sum(repsBase==uniqueReps(i))<=k_upper
        keepIdx(repsBase==uniqueReps(i))=0;
    end
    if sum(repsBase==uniqueReps(i))<=minSum
        minSum=sum(repsBase==uniqueReps(i));
    end
end
dataScore=dataScore(logical(keepIdx),:);
while minSum <=k_upper
   minSum = size(dataScore,1);
   repsBase = runDoubleShift2(dataScore(:,1:k_upper));
   uniqueReps = unique(repsBase);
   keepIdx = ones(length(repsBase),1);
   for i= 1:size(uniqueReps,2)
       if sum(repsBase==uniqueReps(i))<=k_upper
            keepIdx(repsBase==uniqueReps(i))=0;
       end
       if sum(repsBase==uniqueReps(i))<=minSum
           minSum=sum(repsBase==uniqueReps(i));
       end
       if length(uniqueReps)==length(repsBase)
          disp('Cycles in pre-clustering; exiting');
          return;
       end
   end
   %now we need to do the trimming
   dataScore = dataScore(logical(keepIdx),:);
 
end

[V,A,data,VProj,...
    M,I,E,K,CIp,...
    Ps,Fs,RW,penDf,penPf,...
    dimList] = ...
    weightedSCUnmixMainSliver(dataScore(:,1:k_upper),...
    theta,gamma,dataCoeff,sigma);

%timesaver for debuging.
%save('TmpCheckspace');

%check to see if A is one connected component
numConnComp = graphconncomp(sparse(A),'DIRECTED',false);

if numConnComp > 1
   [Vnew,Anew] = determine_min_sc(VProj(:,1:1+min(k_upper,max(sum(A)))),...
       A,dataScore(:,1:1+min(k_upper,max(sum(A))))); %dunno if this is to be max(sum(A)) or sum(sum(A))
   M = recompute_M(Vnew,Anew,dataScore(:,1:size(Vnew,2)));
   [penPf,penDf]=recompute_penalty_function(Vnew,...
       Anew,M,dataScore(:,1:size(Vnew,2)),gamma); %looks like an edge-case bug here for lines. in the automated version.
   
   Atmp = Anew;
   for i = 1:size(Anew,1)
       Atmp(i,i)=0;
   end
   MC = maximalCliques(Atmp); 
   [FNew,errors,errorRatio]=recomputeFs(Vnew,MC,dataScore(:,1:size(Vnew,2)));
else
    Vnew = V;
    Anew = A;
    FNew = Fs;
    errors = E;
    errorRatio = 1;
end


%output the result to a file

csvwrite(strcat(filename,'_V.csv'),Vnew);


csvwrite(strcat(filename,'_A.csv'),Anew);

%think about how the new VProj will be computed
VProjNew = ...
    Vnew*transpose(dataCoeff(:,1:size(Vnew,2)));

csvwrite(strcat(filename,'_VProj.csv'),VProjNew);

csvwrite(strcat(filename,'_coeff.csv'),dataCoeff);

csvwrite(strcat(filename,'_F.csv'),FNew);

fid = fopen(strcat(filename,'_ids.txt'),'w+');
ids = char(ids);
for i = 1:size(ids,1)
    fprintf(fid,'%s\n',ids(i,:));
end
fclose(fid);
%csvwrite(strcat(filename,'_ids.csv'),char(ids)); %need to be converted to char 
%equivalents of the ascii codes
fid = fopen(strcat(filename,'_gene_list.txt'),'w+');
geneList = char(geneList);
for i = 1:size(geneList,1)
      fprintf(fid,'%s\n',geneList(i,:));
end
fclose(fid);
%dlmwrite(strcat(filename,'_gene_list.csv'),char(geneList),'newline','unix'); % this does not work


disp(['Neg. Log Likelihood: ' num2str(penDf+penPf)]);
analysisTime = toc;
disp(['Time for Analysis in s ' num2str(analysisTime)]);
disp(['Time for Analysis in min ' num2str(analysisTime/60)]);
disp(['Time for Analysis in hours ' num2str(analysisTime/60/60)]);


end