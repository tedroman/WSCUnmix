sigmaVals = 0.05:0.05:1;
dimEsts = zeros(length(sigmaVals));
for i = 1:length(sigmaVals)
   dimEsts(i) = sliver_dim_est(sigmaVals(i),snorm,12); 
   disp(['Estimate complete for sigma = ' num2str(sigmaVals(i))]);
end