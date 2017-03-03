function[errVec]=batch_compute_rmsd_vertices(scenarioStruct,ss2)
for scenario= 1:7
   gtV = getGtV(scenario,size(ss2{scenario}.data,2));
    %allEs = [scenarioStruct{scenario}.E_h scenarioStruct{scenario}.E_KNN];
    %[minE,minEIdx]=min(allEs);
    %minE_vert = scenarioStruct{scenario}.VProjCell{minEIdx};
        %minE_vert = scenarioStruct{scenario}.sc_C';%scenarioStruct{scenario}.VProjCell{minEIdx};
    minE_vert = scenarioStruct{scenario}.mu;
   %minE_vert = scenarioStruct{scenario};   
   errVec(scenario)=compute_rmsd_vertices(minE_vert,gtV);
    
end

end
