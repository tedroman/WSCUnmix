function[GM,GMProj]=computeGMCorrectNumComponents(scenarioStruct)

for scenario= 1:7
  numComponents = 0;
    if scenario==1
       numComponents = 4;
   elseif scenario==2
           numComponents = 5;
    elseif scenario==3
        numComponents = 7;
   elseif scenario==4
            numComponents = 6;
   elseif scenario==5
            numComponents = 5;
   elseif scenario==6
            numComponents = 5;
   elseif scenario==7
            numComponents = 6;
    end
    datascore = scenarioStruct{scenario}.datascore;
    coeff = scenarioStruct{scenario}.datacoeff;
   GM{scenario} = fitgmdist(datascore(:,1:7),...
      numComponents,'RegularizationValue',1e-3);
  GMProj{scenario}=GM{scenario}.mu*coeff(:,1:7)';
end