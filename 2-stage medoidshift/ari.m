function [ similarity ] = ari ( rep, n_point )

%compute ari similarity between inferred class (rep) and ground truth
%(n_point). Rep is a 1*n vector for the label of n data points, and n_point
%is a 1*m vector for the size of m classes.

    n_class = numel(n_point);
    infer_label = unique(rep);
    n_infer_class = numel(infer_label);
    true_set = cell(1,n_class);
    infer_set = cell(1,n_infer_class);
    n =  sum(n_point);
    N = zeros(n_class,n_infer_class);
    
    for i = 1 : n_class
        
        true_set{i} = sum(n_point(1:i-1))+1 : sum(n_point(1:i));
        
    end
        
    for j = 1 : n_infer_class
        
        infer_set{j} = find(rep==infer_label(j));
        
    end
    
    for i = 1 : n_class
        
        for j = 1 : n_infer_class
            
            %N(i,j) = numel(intersect(true_set{i},infer_set{j}));
            [~,true_set{i},infer_set{j},N(i,j)] = ...
                intersect_c(true_set{i},infer_set{j});
            
        end
        
    end
            
    Ci = sum(N,2);
    Cj = sum(N,1);
    t0 = sum(sum(N.*(N-1)/2));
    t1 = sum(Ci.*(Ci-1)/2);
    t2 = sum(Cj.*(Cj-1)/2);
    t3 = 2*t1*t2/n/(n-1);
    similarity = (t0 - t3) / (t1/2 + t2/2 - t3);

end