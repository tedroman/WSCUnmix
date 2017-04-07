function [ rep2 ] = doubleshift ( data, hc2, t_depth, hc1 )

% implement clustering with two steps of medoidshift
% first one with regular kernal exp(-D/h)
% second one with improper kernel exp(-D/h)-1
% data: nxm matrix with n points in m dimension
% hc2: bandwidth coefficient of second medoidshift
%      bandwidth = coefficient * mean distance
% t_depth: steps of back-trace for the first medoidshift
% hc1: bandwidth coefficient of first medoidshift
% rep1 & rep2: representative vector of data points after first and second
%              medoidshift

    if nargin < 4, hc1 = 1; end
    
    if nargin < 3, t_depth = 1; end
    
    if nargin < 2, hc2 = 1; end
    
    if t_depth > 0 || hc2 == 0
            
        d1 = squareform(pdist(data).^2);
        %d1 = get_path_metric(data);
        h1 = hc1 * mean(d1(:));
        trace1 = medoidshift(d1,h1,0);
            
        if hc2 > 0
        
            rep1 = trace1(min(t_depth,size(trace1,1)),:);
                
        else
                
            rep2 = trace1(end,:);
                
        end
            
    else
            
        rep1 = 1 : size(data,1);
            
    end
        
    if hc2 > 0
            
        u_rep1 = unique(rep1);
        n_rep1 = numel(u_rep1);
            
        if n_rep1 == 1
            
            rep2 = rep1;
            
        else
            
            d2 = get_path_metric_squared(data(u_rep1,:)); % L2-squared path
            %d2 = squareform(pdist(data(u_rep1,:)).^2); % L2-squared
            %d2 = get_path_ratio(data(u_rep1,:)); % L2^2 / L2^2 path
            %d2 = squareform(pdist(data(u_rep1,:)).^2);
            h2 = hc2 * mean(d2(:));
            %h2 = hc2;
            trace2 = medoidshift(d2,h2,1);
            hash = containers.Map(u_rep1,u_rep1(trace2(end,:)));
      
            rep2 = cell2mat(hash.values(num2cell(rep1)));
                
        end
            
    end
        
end