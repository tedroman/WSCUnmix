function [ path_metric ] = get_path_metric ( coordinates )

    path_metric = graphallshortestpaths(...
        sparse(tril(squareform(pdist(coordinates)))),'Directed',false);

end