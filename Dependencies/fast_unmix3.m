function [C,a,F,P,lowDimEnds, ECs, PCs] = fast_unmix3(M, k)
% fast_unmix : calls Chan's convex program for M.V.E.S.
%
% [C,a,F,P] = fast_unmix(M, k); 
% 
% uses MVES, requires SEDUMI 
% 
% David Tolliver (c) Carnegie Mellon University 2009 
% 

M = M'; % to put into poitns * dims.

if(0==exist('sedumi'))
    error('fast_unmix requires that the SEDUMI pkg be in the path');
end
[coeffs,scores]=pca(M,'Economy',true,'Centered',false);

%[V,S]=princomp(M','econ');
P = transpose(scores(:,1:k-1));
%take the first k-1 pcs.

%P=S(:,1:k-1)';

% rescale the points --no longer needed

%nrps = max(sqrt(sum(P.^2,1))) 
%P = P./nrps; 
%tic
%disp('computing minimum volume simplex'); 

[C, F] = MVES(P, k,0); 
%P = C * F 
%M = ~ P'*coeffs(:,1:k-1)'


%C = nrps*C; 

lowDimEnds = C; 
C=coeffs(:,1:(k-1))*C;
%for i=1:k
% C(:,i)=C(:,i)+mean(M,2);
%end
F = F'; 
% promote to full dimension

%P = nrps*P; 
%keyboard 

%a = simplex_vol(C); 
a = 0; %Placeholder.
%disp([' complete in ' num2str(toc) ' seconds']); 
%save(filename,'C','a','F','P','Cc');

ECs = scores(:,1:(k-1)); 
PCs = P(:,1:(k-1));
end