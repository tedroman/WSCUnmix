function[penPf,penDf]=recompute_penalty_function(V,A,M,data,gamma)

%need to find the cliques--use maximal cliques (assume it is on the path)
Aprime = A;
for i = 1:size(Aprime,1)
    Aprime(i,i)=0;
end
%precondition for clique finding
MCA = maximalCliques(Aprime);
numCliques = size(MCA,2);
%need to fit ALL data as before.
%weight by M.

numPoints = size(data,1);
fit_error = 0;
mst_cost =0;
for i = 1:numCliques
   for j = 1:numPoints
       point_error = l1_dist2simplex(V(logical(MCA(:,i)),:),data(j,:),M(j,i));
       fit_error = fit_error+point_error;
   end
   mst = graphminspantree(sparse(dist(V(logical(MCA(:,i)),:)')));
   mst_cost = mst_cost + sum(reshape(mst,numel(mst),1));
end

prior_penalty = gamma*(log(mst_cost)+size(data,1)*size(V,1));
penPf = prior_penalty;
penDf = fit_error;
end

 function error = l1_dist2simplex(V, datum,reweight)
        d = numel(datum)+1; 
        d2 = size(V,1);
        opts = optimset('Display','off');
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
            error = sum(abs(transpose(datum)-V'*alpha))*reweight;%(sigmoid(reweight));
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