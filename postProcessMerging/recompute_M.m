function[Mnew]=recompute_M(Vnew,Anew,s)
%recompute the ratios for each of the 
%subsimplices in the simplicial complex denoted by 
%Vnew,Anew as a pair
Atst = Anew;
opts = optimset('Display','Final');
for i = 1:size(Anew,1)
    Atst(i,i)=0;
end
MC = maximalCliques(Atst);
Mnew = zeros(size(s,1),size(MC,2));
distMtrx = zeros(size(s,1),size(MC,2));
for i = 1:size(MC,2)
    VSub = Vnew(MC(:,i)==1,:);
    d = size(VSub,1);
    for j = 1:size(s,1)
        point = s(j,:);
        alpha =  lsqlin(VSub', point, -eye(d), zeros(d,1),...
            ones(1,d),1,[],[],[], opts);
        
        distMtrx(j,i) = sum(abs(VSub'*alpha - point'));
      
    end
end
Mnew = distMtrx./...
    repmat(sum(...
    transpose(distMtrx))',1,size(MC,2));


end