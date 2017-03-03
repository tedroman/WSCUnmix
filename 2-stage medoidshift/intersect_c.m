function [ I, CA, CB, nI, nCA, nCB ] = intersect_c ( A, B )

%Input: two sets A, B
%Output:
%  I: intersect of A and B
%  CA: A - I
%  CB: B - I

   nA = numel(A);
   nB = numel(B);
   CA = zeros(1,nA);
   CB = zeros(1,nB);
   I = zeros(1,min(nA,nB));
   nCA = 0;
   nCB = 0;
   nI = 0;
   
   %for i = 1 : nA-1, if A(i) > A(i+1), A = sort(A); break, end, end
   
   %for i = 1 : nB-1, if B(i) > B(i+1), B = sort(B); break, end, end
   
   %no need to sort if A and B are indices given by "find".
   
   i = 1;
   j = 1;
   
   while i <= nA && j <= nB
       
       if A(i) == B(j)
           
           nI = nI + 1;
           I(nI) = A(i);
           i = i + 1;
           j = j + 1;
           
       elseif A(i) < B(j)
           
           nCA = nCA + 1;
           CA(nCA) = A(i);
           i = i + 1;
           
       else
           
           nCB = nCB + 1;
           CB(nCB) = B(j);
           j = j + 1;
           
       end
       
   end
   
   I = I(1:nI);
   CA = [CA(1:nCA),A(i:end)];
   CB = [CB(1:nCB),B(j:end)];

end