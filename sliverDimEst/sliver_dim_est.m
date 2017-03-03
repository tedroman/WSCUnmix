function[dimEst]=sliver_dim_est(sigma,s,maxDim)
dimEst=-1;
for j = 2:maxDim
   %get the enclosure
   [C,a,F,P,lowDimEnds, ECs, PCs] = fast_unmix3(s(:,1:j), j+1,true);

   %what about this simplex?
   
   %find the longest edge.
   %the longest edge can be computed
   %with a distance matrix
   dMtrx = dist(C);
   longestEdge = max(max(dMtrx));
   C_to_compute_area = transpose(C(1:j,:));
   [~,area]=convhulln(C_to_compute_area);
   sliverFactor = (longestEdge^j)/(factorial(j));
   if area < (sigma^j)*sliverFactor
      disp(['Found a sliver in dimension: ' num2str(j)]);
      dimEst = j-1;
      break;
   end
   disp('finished a step');

    
end
   if dimEst==-1
       dimEst = maxDim;
   end
end