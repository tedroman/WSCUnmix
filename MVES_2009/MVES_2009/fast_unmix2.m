function [C,a,F,P,lowDimEnds, ECs, PCs] = fast_unmix2(M, k); 
% fast_unmix : calls Chan's convex program for M.V.E.S.
%
% [C,a,F,P] = fast_unmix(M, k); 
% 
% uses MVES, requires SEDUMI 
% 
% David Tolliver (c) Carnegie Mellon University 2009 
% 


if(0==exist('sedumi'))
    error('fast_unmix requires that the SEDUMI pkg be in the path');
end

[V,S]=princomp(M','econ');
P=S(:,1:k-1)';

% rescale the points
nrps = max(sqrt(sum(P.^2,1))) 
P = P./nrps; 
%tic
%disp('computing minimum volume simplex'); 
[C, F] = MVES(P, k,0); 
C = nrps*C;


lowDimEnds = C; 
C=V(:,1:(k-1))*C;
for i=1:k
 C(:,i)=C(:,i)+mean(M,2);
end

%to add the constraint that one of the C's must be fully diploid, consider
%stretching the C most close to fully diploidy.

%minResid=1e10;
%minResidInd=0;
%for i=1:size(C,2)
%    resid=norm(C(:,i)-ones(size(C,1),1),2);
%    if resid <minResid
%        minResid=resid;
%        minResidInd=i;
%    end
%end

%C(:,minResidInd)=ones(size(C,1),1);



F = F'; 
% promote to full dimension

P = nrps*P; 
%keyboard 
a = simplex_vol(C); 
%disp([' complete in ' num2str(toc) ' seconds']); 
%save(filename,'C','a','F','P','Cc');

ECs = V(:,1:(k-1)); 
PCs = P(:,1:(k-1)); 
