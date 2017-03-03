function totalError = simplexfitErrorMultMSTWeighted(V,data,gamma,reweights)

if ~exist('gamma','var')
    gamma = 1;
end
if ~exist('reweights','var')
    reweights = ones(size(data,2),1);
end

opts = optimset('Display','off');
fit_error = 0;
warning off;

max_point_error = 0;

for i =1:size(data,2)
    %visualize--initial fits are not bad.
    %hold on;
    %line(V(:,1)',V(:,2)',V(:,3)','color','c');
    %line([V(end,1) V(1,1)],[V(end,2) V(1,2)],[V(end,3) V(1,3)],'color','c');
    %scatter3(data(1,i),data(2,i),data(3,i),25,'k','filled');
    disp(num2str(i))
    point_error = l1_dist2simplex(V,data(:,i),reweights(i));
    if point_error > max_point_error
        max_point_error = point_error;
    end
    fit_error = fit_error + point_error;
end
warning on;
%hold on;
%line(V(:,1)',V(:,2)',V(:,3)','color','c');
%line([V(end,1) V(1,1)],[V(end,2) V(1,2)],[V(end,3) V(1,3)],'color','c');

mst = graphminspantree(sparse(dist(V')));
mst_cost = sum(reshape(mst,numel(mst),1));
prior_penalty = gamma*(log(mst_cost)+size(data,2)*size(V,1));
totalError = fit_error + prior_penalty;%gamma*(log(mst_cost)+size(data,2)*(size(V,1))); %modified 23 Oct 2015
%trying new prior
%totalError = fit_error + log(gamma) - log(mst_cost);
% iso-tropic \ell_1 distance to the simplex of the sample point datum
    function error = l1_dist2simplex(V, datum,reweight)
        d = numel(datum)+1; 
        d2 = size(V,1);
        alpha =  lsqlin(V', datum, -eye(d2), zeros(d2,1), ones(1,d2),1,[],[],[], opts);
      
        %argument order: C, d, A, b, Aeq, beq, lb, ub, x0, options
        %finds x such that least squares of C*x-d s.t. Ax<=b, Aeqx=beq,
        %lb<=x<=ub
        
        %this means we want our program to find the minimal vector such
        %that norm( V'*alpha - datum) is minimized, subject to
        %-eye(d)*alpha < = zeros(d,1) (that is, alpha > = 0) for all
        %dimensions, ones(1,d)*alpha=1 (alpha has desirable geometric behavior)
        %if weight > 0 
            %to get l1 we need to take sum of abs value of points.
            error = sum(abs(datum-V'*alpha))*reweight;%(sigmoid(reweight));
            if error < 0
                disp('problem with error');
            end
       % error = sqrt(sum((datum-V'*alpha).^2))*log(weight);
        %else
         %   error = sum(abs(datum-V'*alpha))*((1e-2));
            %error = sqrt(sum((datum-V'*alpha).^2))*log(1e-2); %small substitute to avoid -inf
        %end
        
        %seems to be working for nearby points.
        %what happens when the weight is low?
        if (reweight < 0.3)
          %  disp('Low Weight point');
        end
    end

end
function[y]=sigmoid(x)
y = (1+exp(-20*(x - 0.5)))^-1; %just see what happens with this--experiment.
end