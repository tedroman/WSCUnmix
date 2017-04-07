function [ mode_path ] = medoidshift( D, h, a )

% implement medoidshift
% D: distance metric; h: kernel bandwidth; a = 0 or 1.
% a = 0: regular medoidshift with kernel 1-D/h for D<=h
% a = 1: medoidshift with improper kernel exp(-D/h)-1

    if nargin == 2 || (a ~= 0 && a ~= 1)
        
        error(['Insufficient or incorrect input arguments.\n',...
               'Usage: mp = medoidshift(D,h,a)\n',...
               'D: distance metric; h: kernel bandwidth; a = 0 or 1.'],...
               nargin);
           
    end
    
    if a == 0
        
        K = 1 - D/h;
        K(K<0) = 0;
        
    else
        
        K = 1 - exp(D/h);
        %K = -(D/h).^2;
        %K = exp(-D/h) - 1;
        %K = -D/h;
        
    end
    if (size(D,2)~=size(K,1))
        disp('Error in size--mismatch')
    end
    S = D * K;
    [~,mapping] = min(S);
    mode_path = zeros(size(D));
    mode_path(1,:) = mapping;
    old_label = 0;
    count = 1;
    
    while any(old_label~=mode_path(count,:))
        
        old_label = mode_path(count,:);
        count = count + 1;
        mode_path(count,:) = mapping(old_label);
        
    end
    
    mode_path(count:end,:) = [];

end