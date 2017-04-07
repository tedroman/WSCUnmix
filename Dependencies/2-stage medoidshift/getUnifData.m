function [ data ] = getUnifData ( n_point, vertex )

% generate uniformly distributed data points within the convex hull defined
% by vertex.
% n_point: number of data points
% vertex: n-by-p matrix that defines n vertices in p-dimension
% data: n_point-by-p matrix that contains n_point data points

    n_vertex = size(vertex,1);
    raw_data = gamrnd(1,1,n_point,n_vertex);
    S = sum(raw_data,2);
    data = raw_data./repmat(S,1,n_vertex) * vertex;

end