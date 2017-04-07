%iterate over the TCGA labels from our data
%find a map from our data to the data we have imported.
mapToPubd = zeros(472,1);
for i = 1:length(TCGAIDPub)
    
       tcgaIdPubCurr = char(TCGAIDPub{i});
    for j = 1:size(tcgaToRowIDChar,1)
       if strcmp(tcgaIdPubCurr(1:16),tcgaToRowIDChar(j,1:16))
           mapToPubd(j)=i;
       end
    end
end

ourPurityEst = zeros(sum(mapToPubd>0),1);
theirPurityEst = zeros(sum(mapToPubd>0),1);
counter = 1;
for i = 1:size(tcgaToRowIDChar,1)
    if mapToPubd(i)>0 && ~isnan(CPE(mapToPubd(i)))
       ourPurityEst(counter)=M(i,1)*purityEstTri(i)+M(i,2)*purityEstTet(i);
       
       if ~isnan( CPE(mapToPubd(i)))
          theirPurityEst(counter) = CPE(mapToPubd(i));
       else
           theirPurityEst(counter) =IHC(mapToPubd(i));
       end
       
       counter=counter+1;
    end
end

scatter(ourPurityEst,theirPurityEst,'ko');
axis([0 1 0 1]);
xlabel('Our Purity');
ylabel('Their Purity');

%{
ybar = mean(ourPurityEst);
ssTot=sum((theirPurityEst-repmat(ybar,length(theirPurityEst),1)).^2);
ssReg=sum((ourPurityEst-repmat(ybar,length(ourPurityEst),1)).^2);
ssRes=sum((theirPurityEst-ourPurityEst).^2);
r2=1-ssRes/ssTot;
%}
