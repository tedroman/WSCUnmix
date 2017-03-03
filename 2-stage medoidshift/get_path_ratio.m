function [ path_ratio ] = get_path_ratio ( coordinates )

    l1 = squareform(pdist(coordinates,'cityblock'));
    %l2 = squareform(pdist(coordinates));
    %l2 = squareform(pdist(coordinates).^2);
    linf = squareform(pdist(coordinates,'chebychev'));
    %l2p = graphallshortestpaths(sparse(tril(l2)),'Directed',false);
    path_ratio = l1.^2./linf.^0.5;
    path_ratio(linf==0) = 1;
    %path_ratio=linf.^2;
        
end