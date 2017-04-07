function[errVecGM]=getGMError(scenarioStruct3,GM3)
for scenario = 1:7
    coeff = scenarioStruct3{scenario}.datacoeff;
    currMu = GM3{scenario}.mu;
    currMuProj = currMu*coeff(:,1:size(currMu,2))';
    gtV = getGtV(scenario,size(scenarioStruct3{scenario}.data,2));
    errVecGM(scenario)=compute_rmsd_vertices(currMuProj,gtV);
end