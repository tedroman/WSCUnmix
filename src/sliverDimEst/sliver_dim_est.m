function[dimEst]=sliver_dim_est(sigma,s,maxDim)
dimEst=-1;
for j = 2:maxDim
   %get the enclosure
   disp(['number of points: ' num2str(size(s,1))]);
   if size(s,1)==1
       dimEst = 0;
       return;
   elseif size(s,1)==2
       dimEst = 1;
       return;
   else
   [C,a,F,P,lowDimEnds, ECs, PCs] = fast_unmix4(s(:,1:j), j+1,true);
   end
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