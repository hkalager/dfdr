% return p-val grid for one-sided permutation test
get_p_grid_perm_one_sided=@(nPerm) (1:nPerm)/nPerm;

% if two support are numerically tied, prefer the left/smaller one to break the tie consistently
tol=1e-16;
function outarg=which_min_left(x, tol)
[~,index]=min(x);
if index==1
    outarg=1;
elseif (x(index-1)-x(index))<tol
    outarg=index-1;
else
    outarg=index;
end
end

% DRB method
%input: 
% p_count: the observed p-value counts
% p_grid: p-value grid (from small to large)
% n: number of bins
% lambda_max: largest lambda considered    
function pi_opt=est_pi0_disc(p_count, p_grid, n, lambda_max)
m=sum(p_count);
cum_count=cumsum(p_count);
target=1/n*(1:(n-1));
candidate=p_grid(and(p_grid>0,p_grid<1));
% for each target, find the closest candidate
lambda_vec=unique(candidate(which_min_left(abs(candidate-target))));
if length(lambda_vec)==1 
    pi0=(m-cum_count(lambda_vec==p_grid))/(1-lambda_vec)/m; 
    lambda=lambda_vec;
    pi_opt={pi0;lambda};
    return;
end
if lambda_max<1
    lambda_vec=lambda_vec(lambda_vec<=lambda_max);
end
gap = [lambda_vec(1), diff(lambda_vec')];
index_lambda= find(p_grid==lambda_vec);
cum_bin_count=cum_count(index_lambda);
% all but the last bin, length B-1
bin_counts=[cum_bin_count(1), diff(cum_bin_count)];
% number of bins
B=length(lambda_vec)+1;
R = cumsum(bin_counts);
tail_m0 = (m-R)./(1-lambda_vec');
temp = bin_counts./gap - tail_m0;
if sum(temp <= 0)>0
            index =find(temp<=0,1,'first');
else
            index = B-1;
        
end
pi0=tail_m0(index)/m; 
lambda=lambda_vec(index);
pi_opt={pi0,lambda};
end