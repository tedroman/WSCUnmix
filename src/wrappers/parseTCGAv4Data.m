function[nonPCData,dataIds,geneList]=parseTCGAv4Data(nonPCDir)


%input: a string describing the root directory, under which
%level 4 TCGA CNV data is stored


%output: a matrix form of the data, along with
%a cell array of strings that contain the ids.
nonPCData=[];
geneList = {};
dataIds=[];
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

%otherwise, begin to read the files.
%use readlink to find the full path name for each file

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
       '%s%s%s%f','Delimiter','\t','Headerlines',1);
   
   fclose(fid);
   end
   %outStruct has all fields-- we need to shrink
   %use only the first and fourth fields.
   if exist('outStruct','var')
    outStructSmall = [outStruct(1) outStruct(4)];
    clearvars outStruct;
   else
    outStructSmall = [{'XYZtemp'} zeros(size(nonPCData,2),1)];
   end
   %first 28 chars of currentFileName are the 
   %TCGA ID in full
   disp(currentFileName)
   if length(currentFileName)>28
   dataIds(i,:) = currentFileName(1:28);
   else
   dataIds(i,:) = repmat('o',1,28);
   end
   
   
   %now we process the outStruct.  
   %if the gene list is not stored yet, we should store it
   %but not overwrite it
   
   if length(geneList) < 10
      %store the gene List if it has the sentinel value
      
      geneList = outStructSmall{1};
   end
   
   %store the data values
   if i > 1 && length(outStructSmall{2})~=size(nonPCData,2)
       outStructSmall = [{'XYZtemp'} zeros(size(nonPCData,2),1)];
   end
   nonPCData(i,:) = outStructSmall{2}';
   disp(['Completed parsing ' ...
       num2str(i) ' of '...
       num2str(length(fileNamesCell))]); 
   
  
end



end