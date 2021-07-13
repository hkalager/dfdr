function pvalues = compute_pvalues(Perfs, Perfs_B, sigma)
% Perfs need to be centred around the benchmark
% Perfs_B need to be centred around the Perfs

[nbstrats, B] = size(Perfs_B);

pvalues = zeros(nbstrats, 1);

for i = 1:nbstrats
    pvalues(i) = 2 * min(sum((Perfs_B(i,: ) - Perfs(i)) > Perfs(i)),...
        sum((Perfs_B(i,:) - Perfs(i)) < Perfs(i))) / B;
end

% assume that strategy with volatility=0 from the null hypothesis:
pvalues(sigma == 0) = rand(sum(sigma == 0), 1);
