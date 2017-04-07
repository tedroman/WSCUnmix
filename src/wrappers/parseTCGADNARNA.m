function[nonPCData,dataIds,geneList]=parseTCGADNARNA(nonPCDir)
nonPCData=[];
nonPCDataRNA = [];
nonPCDataDNA = [];
dnaCt =1;
rnaCt =1;
geneList = {};
geneListRNA = {};
geneListDNA = {};
dataIds=[];

%we will need to maintain separate 
%data ids for dna and rna then merge

dataIdsDNA = [];
dataIdsRNA = [];

%an alternative way to encode this is 
%a vector to indicate DNA or RNA

%first -- we need to strip away the directory structure,
%as there may be multiple layers to the subdirectories under the root
[status,fileNamesCell]=...
    system(['ls -R ' nonPCDir ' | grep "\.txt$"']);
%this assumes there are not additional txt files, and only the
%files used for importation

if status~=0
    disp('Error: could not execute file name extraction');
    return;
end
fileNamesCell = strsplit(fileNamesCell,'\n');


for i = 1:length(fileNamesCell)
   currentFileName = fileNamesCell{i};
   [status,currentFileNameFull]=...
       system(...
       ['find '...
       nonPCDir ' -type f -not -path '...
       '''*/\.*'''...
       '| grep "' currentFileName '"']);
   if status~=0
       disp(['Error: could not determine full file '...
           'path for ' currentFileName]);
       return;
   end
   
  tempFNCell = strsplit(currentFileNameFull,'\n');
  currentFileNameFull = tempFNCell{1};
   
   [status,currentDir]=system('pwd');
   if status ~=0
       disp('Error: could not print current directory');
   return;
   end;
   %next step: 
   %read in the file contents, which are 
   %gene name \t locus id (not used) \t cytoband (not used) ...
   %\t sampleid/cnv ratio (sample id is on the header line)
   %disp(['length of file: '...
   %    num2str(length([currentDir(1:end-1)...
   %    '/' currentFileNameFull(1:end-1)]))]);
   
   %disp([currentDir(1:end-1)...
   %    '/' currentFileNameFull(1:end-1)]);
   
   fullFS = [currentDir(1:end-1) currentFileNameFull];
   fid = fopen(currentFileNameFull,'r');
   
   %to debug:
   %disp([currentDir(1:end-1)...
   %    '/' currentFileNameFull(1:end-1)]);
   
   if fid~=-1
   outStruct = textscan(fid,...
       '%s%f','Delimiter','\t','Headerlines',1); %this has to be 
   %customized for RNA and for DNA DOH!
   
   fclose(fid);
   
   %go through, and see if any of the 'gene' column has '.' characters
   %these are indicative that this was actually a DNA file and not an RNA
   %file -- the '.' represents the cytoband.
   
   atLeastOneNonBar = false;
   for j = 1:length(outStruct{1})
      if isempty(strfind(outStruct{1}{j},'|'))
          atLeastOneNonBar = true;
          disp(outStruct{1}{j})
          break;
      end
   end
   
  
   if atLeastOneNonBar
      fid = fopen(currentFileNameFull,'r');
      if fid~=-1
         outStruct = textscan(fid,...
             '%s%s%s%f','Delimiter','\t',...
             'Headerlines',1);
      end
      fclose(fid);
       
   end
   end
   
   
   %how big is the outstruct? -- this tells us if it is DNA or RNA
   %in the DNA Case: outstruct length = 4
   %in the RNA Case: outstruct length = 2
   disp(currentFileName)
   if length(currentFileName)>28
   dataIds(i,:) = currentFileName(1:28);
   else
   dataIds(i,:) = repmat('o',1,28);
   end
   
   %have a checkpoint here.
   if length(outStruct)==2
       %RNA case
       if length(currentFileName) > 28
            dataIdsRNA(rnaCt,:) = currentFileName(1:28);
       else
           dataIdsRNA(rnaCt,:)= repmat('o',1,28);
       end
       if length(geneListRNA) < 10
      %store the gene List if it has the sentinel value
            geneListRNA = outStruct{1};
       end
       
       %now process:
       if rnaCt > 1 && length(outStruct{2})~=size(nonPCDataRNA,2)
        outStruct = [{'XYZtemp'} zeros(size(nonPCDataRNA,2),1)];
       end
        nonPCDataRNA(rnaCt,:) = outStruct{2}';
        disp(['Completed parsing ' ...
            num2str(i) ' of '...
            num2str(length(fileNamesCell))]); 
        rnaCt = rnaCt + 1;
      
    else
       %DNA case
       if length(geneListDNA) < 10
           %store the gene List if it has the sentinel value
           geneListDNA = outStruct{1}; %need separate DNA 
       end
       
       if length(currentFileName) > 28
           dataIdsDNA(dnaCt,:) = currentFileName(1:28);
       else
           dataIdsDNA(dnaCt,:) = repmat('o',1,28);
       end
       
       
       if dnaCt > 1 && length(outStruct{4})~=size(nonPCDataDNA,2)
        outStruct = [{'XYZtemp'} {'XYZtemp'} {'XYZtemp'} zeros(size(nonPCDataDNA,2),1)];
       end
   nonPCDataDNA(dnaCt,:) = outStruct{4}';
   disp(['Completed parsing ' ...
       num2str(i) ' of '...
       num2str(length(fileNamesCell))]);
   dnaCt = dnaCt + 1;
       
    end
   
   
  
   
   %store the data values
   
   
  
end

dataIdsDNAChar= char(dataIdsDNA);
dataIdsRNAChar = char(dataIdsRNA);



rnaZ = zscore(2*nonPCDataRNA);
dnaZ = zscore(10.^nonPCDataDNA);

%if the sample is from the same person, then the first 
%12+3 = 15 characters will be the samej

% we need to create a map of all the samples to keep
%and how the DNA samples are indexed in RNA and vica versa

dnaToRNAMap = zeros(size(dataIdsDNAChar,1),1);
rnaToDNAMap = zeros(size(dataIdsRNAChar,1),1);

for r = 1:size(dataIdsDNAChar,1)
    dnaPatient = dataIdsDNAChar(r,1:15);
    for s = 1:size(dataIdsRNAChar,1) % this is inefficient but these are small so it's just easy to code
        rnaPatient = dataIdsRNAChar(s,1:15);
        if strcmp(dnaPatient,rnaPatient)
           dnaToRNAMap(r)=s;
           rnaToDNAMap(s)=r;
        end
    end
    
end

dataIdsFinal = [];
dataCt = 1;
for t = 1:length(dnaToRNAMap)
   if dnaToRNAMap(t)>0
       dataIdsFinal(dataCt,:) = ...
           dataIdsDNA(t,:);
       nonPCData(dataCt,:) = ...
           [dnaZ(t,:) rnaZ(dnaToRNAMap(t,:))];
       dataCt = dataCt + 1;
       
   end
end
%test this module
%nonPCData = [dnaZ rnaZ];

%forgotten point -- can only merge data points
%that have both RNA and DNA data.

dataIds = dataIdsFinal;
geneList = [geneListDNA' geneListRNA'];




end