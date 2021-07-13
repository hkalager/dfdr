function AA_adjusted=adjustperf_SR(AA_zero, SR)
% input: SR, YEARLY SR!!!

% convert to daily SR:
SR = SR / sqrt(252);

nbdays = size(AA_zero, 1);

sigma = std(AA_zero);

mu = SR * sigma;

AA_adjusted = AA_zero + kron(ones(nbdays, 1), mu);











