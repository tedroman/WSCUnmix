function [ path_metric ] = get_path_metric_squared ( coordinates )

    path_metric = graphallshortestpaths(...
        sparse(tril(squareform(pdist(coordinates).^2))),'Directed',false);

end