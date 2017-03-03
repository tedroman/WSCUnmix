alphaVal9=zeros(1,20);
dimListAll9=cell(1,20);
i=1;
for i=1:length(alphaVal9)
    alphaVal9(i) = i/20;
   [dimListAll9{i},clusterDimVec9{i}]=getDimListMain(s,theta9,gamma9,c,alphaVal9(i),true,i); 
end