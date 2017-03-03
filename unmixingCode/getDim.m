function[dim]=getDim(data,dispFlag,c,alphaVal,figIdx)
%%return the dimensionality of some data, using the gaussian-derived method
if ~exist('dispFlag','var')
    dispFlag = false;
end
disp('generating replicates...')

disp('computing pcs...')

%normalize
%data = data - repmat(min(data),size(data,1),1);
minimax = max(data)-min(data);
%data = data./repmat(minimax,size(data,1),1);
[~,~,~,~,dExpl]=pca(data,'Economy',true);
if isempty(dExpl)
    dExpl = zeros(1,(size(data,2)));
    dExpl(1)=1;
end
%generate some gaussian data of the same dimensionality.
%use 1000 replicates.
%expl = zeros(1000,min(size(data)));
%cumExpl = zeros(1000,min(size(data)));
datan=zeros(size(data));
for i = 1:1000
    for j = 1:size(data,2)
        %datan(:,j) = normrnd(0,1,[size(data,1),1]);
        %datan(:,j) = exprnd(mean(data(:,j)),[size(data,1),1]); %should this be data or datascores?
        %for RNASeq, can we mimic complexity?
       datan(:,j) = exprnd(mean(data(:,j)),[size(data,1) ,1]);
        nmin = min(datan(:,j));
        nmax = max(datan(:,j));
        szn = size(datan(:,j),1);
        %datan(:,j) = (datan(:,j) - repmat(nmin,szn,1))./(repmat(nmax,szn,1)-repmat(nmin,szn,1));
        %datan(:,j) = randn(size(data,1),1) + mean(data(:,j));
        
        %what if this data *is* the PC data?
        %then we need some coefficient matrix to see what the "real" data
        %in the noise model would look like, so that we can look at the 
        %variance explained.
        %we have now passed in c, the coefficients matrix.
        
        
    end
    datan(isnan(datan))=0; %this is actually happening.
    %datan = datan+randn(size(datan))*0.1; %10pct. noise estimate?
    
    
    %reconstruct the 'real' datan
    datanRNASp = datan*c';
    %assuming it is data, do pca (o.w look at just var. in comp.)
    [~,~,~,~,vExpl]=pca(datanRNASp,'Economy',true);
    if ~isempty(vExpl)
        expl(i,:) = vExpl;
    else
        expl(i,1) = 1;
        expl(i,2:end) = zeros(1,size(expl,2)-1);
    end
    for j = 1:size(expl,2)
      cumExpl(i,j) = sum(expl(i,1:j));
    end
end


if dispFlag
   figure(figIdx);
   hold on;
   title('Variance Explanation Plot');
   boxplot(expl);
   plot(dExpl,'ko-');
   hold off;
end

ztarget = abs(norminv(1-(0.01/size(data,2)),0,1));
dim=1; %need to initialize.
%does some parameter go here???

if size(expl,2)>1
    
for j = 2:(size(expl,2)-1)
  % if ( dExpl(j-1) <= (mean(expl(:,j-1)) - std(expl(:,j-1))) && ...
   %        dExpl(j) >= (mean(expl(:,j)) - std(expl(:,j))) && ...
   %        dExpl(j+1) >= ((mean(expl(:,j+1)) - std(expl(:,j+1)))) ) %|| ...
           %( dExpl(j) <=(mean(expl(:,j)) + std(expl(:,j))) && ...
           %dExpl(j) >= (mean(expl(:,j)) - std(expl(:,j))) ) %may need tuning
           
      %make a probability distribution object
      pdPast = makedist('Normal','mu',mean(expl(:,j-1)),'sigma',std(expl(:,j-1)));
      pdCurr = makedist('Normal','mu',mean(expl(:,j)),'sigma',std(expl(:,j)));
      pdNext = makedist('Normal','mu',mean(expl(:,j+1)),'sigma',std(expl(:,j)));
      %ciPast = paramci(pdPast,'Alpha',alphaVal);
      %ciCurr = paramci(pdCurr,'Alpha',alphaVal);
      %ciNext = paramci(pdNext,'Alpha',alphaVal);
      disp('Computing confidence intervals...');
      ciPast(1) = norminv(alphaVal,mean(expl(:,j-1)),std(expl(:,j-1)));
      ciPast(2) = norminv(1-alphaVal,mean(expl(:,j-1)),std(expl(:,j-1)));
      ciCurr(1) = norminv(alphaVal,mean(expl(:,j)),std(expl(:,j)));
      ciCurr(2) = norminv(1-alphaVal,mean(expl(:,j)),std(expl(:,j)));
      ciNext(1) = norminv(alphaVal,mean(expl(:,j+1)),std(expl(:,j+1)));
      ciNext(2) = norminv(1-alphaVal,mean(expl(:,j+1)),std(expl(:,j+1)));
      disp('Ci past')
      disp(num2str(ciPast))
      disp(['dExpl_past:' num2str(dExpl(j-1))])
      disp('Ci curr')
      disp(num2str(ciCurr))
      disp(['dExpl_curr:' num2str(dExpl(j))])
      disp('Ci next')
      disp(num2str(ciNext))
      disp(['dExpl_next:' num2str(dExpl(j+1))])
      
      figure(figIdx);
      line([j-1 j-1], [ciPast(1) ciPast(2)],'linewidth',2,'color','k');
      line([j j],[ciCurr(1) ciCurr(2)],'linewidth',2,'color','k');
      line([j+1 j+1],[ciNext(1) ciNext(2)],'linewidth',2,'color','k');
      
      if ( dExpl(j-1) <=ciPast(2) && ...
              ( dExpl(j) >=ciCurr(1) )  && ...
              ( dExpl(j+1) >=ciNext(1) ) ) 
       break;
      end
end
disp('Dim computed...');
dim = j-1;
end
end