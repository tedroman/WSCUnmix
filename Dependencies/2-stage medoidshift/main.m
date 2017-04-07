repeat = 10;
noise = 0;
t_depth = 1; %trace depth of the 1st medoidshift.
hc1 = 1:10; %bandwidth coefficients for the two medoidshifts;
hc2 = 0; %actual bandwidth = coefficient * mean distance.
score = zeros(repeat,5,numel(hc1));

for i = 1 : 5

    [vertex,n_point] = get_vertex(i);
    dim = size(vertex{1},2); %ambient dimension.
    data = zeros(sum(n_point),dim);
    
    for j = 1 : repeat
        
        for k = 1 : numel(n_point);
            
            data(sum(n_point(1:k-1))+1:sum(n_point(1:k)),:) = ...
                getUnifData(n_point(k),vertex{k});
        
        end
        
        noise_data = data + randn(sum(n_point),dim) * noise;
        
        for h = 1 : numel(hc1)
            
            fprintf(['Structure %d out of 5, repeat %d out of %d, '...
                    'hc1 = %d\n'],i,j,repeat,hc1(h));
            score(j,i,h) = ari(doubleshift(noise_data,hc2,t_depth,...
                hc1(h)),n_point);
    
        end
    
    end
    
end

figure;boxplot(score);ylim([-0.2,1.2]);