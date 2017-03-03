function[FNew,errors,errorRatio]=recomputeFs(V,MC,dataScores)
opts = optimset('Display','Final');
FNew = cell(1,size(MC,2));
for q = 1:size(MC,2)
    currVert = find(MC(:,q)==1);

for i =1:size(dataScores,1)
    d = numel(dataScores(i,:))+1;
    d2 = size(V,2); %num col of K'
    i2 = i
    cv = currVert
    
   % F(:,i) = lsqlin(K',P(:,i),-eye(d),zeros(d,1),ones(1,d),1,[],[],[],opts);
   FNew{q}(:,i) = lsqlin(transpose(V(currVert,:)),transpose(dataScores(i,1:d2)),-eye(length(currVert)),zeros(length(currVert),1),ones(1,length(currVert)),1,[],[],[],opts);
   
   errors{q}(i) = sum(...
       abs(...
       dataScores(i,:)'-transpose(V(currVert,:))*FNew{q}(:,i)));%*reweights(i);%sigmoid(weights(i)); %l1.
   
   %sqrt(sum((dataScores(i,:)'-K'*F(:,i)).^2));
end

end

errorMtrx = zeros(1080,3);
for i = 1:size(MC,2)
    errorMtrx(:,i)= errors{i}';
end
%errorInv = ones(size(errorMtrx))./errorMtrx;
errorTot = sum(errorMtrx'); %smaller dist -> larger membership
errorRatio = errorMtrx./repmat(errorTot',1,size(MC,2));

