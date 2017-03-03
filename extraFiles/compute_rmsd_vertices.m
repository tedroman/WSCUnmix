%do the old way of measuring error:
function[normErr]=compute_rmsd_vertices(minSolutionVert,gtVert)
%find the distance between the min Solution vertices and gt for each
distMtrx = dist(minSolutionVert,gtVert'); %this is in gene space
%do we need to convert to pc space first?
distMtrx =distMtrx';
totaldist = 0;
for i = 1:size(distMtrx,1)
    totaldist = totaldist + min(distMtrx(i,:));
end
normErr = totaldist/(size(gtVert,1)*size(minSolutionVert,2));
normErr = normErr*(1+abs(size(distMtrx,1)-size(minSolutionVert,1))); %compensating for error in estimation of number of components.
end