function [C,a,F,P,lowDimEnds, ECs, PCs] = fast_unmix4(M, k,pcFlag)
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

if ~pcFlag
[coeffs,scores]=pca(M,'Economy',true,'Centered',false);
else
 scores = M;
 coeffs = eye(size(M));
end

P = transpose(scores(:,1:k-1));
%take the first k-1 pcs.

[C, F] = MVES(P, k,0); 
lowDimEnds = C; 
C=coeffs(:,1:(k-1))*C;
F = F'; 

% promote to full dimension
a = 0; %no longer using vol.
ECs = scores(:,1:(k-1)); 
PCs = P(:,1:(k-1));

end