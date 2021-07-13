function [AA_final, nboutperf, nbunderperf, nbnull]=preparadata_auxfun_mu(AA, outperf_mu, ...
    underperf_mu, prop_outperf, prop_underperf)

% average daily log-return:
mu = mean(AA);

nbdays = size(AA, 1);
nbstrats = size(AA, 2);

nboutperf = round(nbstrats * prop_outperf);
nbunderperf = round(nbstrats * prop_underperf);
nbnull = nbstrats - nboutperf - nbunderperf;

% center all trajectories:
AA_zero = AA - kron(ones(nbdays, 1), mu);

thresh = 5E-4;
% AA_zero = AA_zero/1.5;
AA_zero_notnull = AA_zero(:, std(AA_zero) > thresh);

nbstrats_notnull = size(AA_zero_notnull, 2);

mu_notnull = mean(AA(:, std(AA_zero) > thresh));

BB = [mu_notnull; 1:nbstrats_notnull]';
BB = sortrows(BB, 1);
BB = [BB (nbstrats_notnull:-1:1)'];
BB = sortrows(BB, 2);
rankings = BB(:, 3);

% select subset of outperforming strategies:
averageranking = nbstrats;
outperf_startidx = 1;
for i = 1:(nbstrats_notnull - nboutperf + 1)
    averageranking_tmp = mean(rankings(i:(i + nboutperf - 1)));
    if averageranking_tmp < averageranking
        outperf_startidx = i;
        averageranking = averageranking_tmp;
    end
end

idx = 1:nbstrats_notnull;
idx_outperf = (idx >= outperf_startidx) & (idx < outperf_startidx + nboutperf);

% select subset of underperforming strategies:
averageranking = 1;
underperf_startidx = 1;
enough_succesive_underperf = 0;
% skip block of selected outperforming strategies:
if (outperf_startidx - nbunderperf) > 0
    enough_succesive_underperf = 1;
    for i = 1:(outperf_startidx-nbunderperf)
        averageranking_tmp = mean(rankings(i:(i + nbunderperf - 1)));
        if averageranking_tmp > averageranking
            underperf_startidx = i;
            averageranking = averageranking_tmp;
        end
    end
end
if (outperf_startidx + nboutperf - 1 + nbunderperf) <= nbstrats_notnull
    enough_succesive_underperf = 1;
    for i = (outperf_startidx + nboutperf):(nbstrats_notnull - nbunderperf + 1)
        averageranking_tmp = mean(rankings(i:(i + nbunderperf - 1)));
        if averageranking_tmp > averageranking
            underperf_startidx = i;
            averageranking = averageranking_tmp;
        end
    end
end

nbmissingunderperf = 0;
if enough_succesive_underperf == 1
    idx_underperf = (idx >= underperf_startidx) & (idx < underperf_startidx + nbunderperf);
else
    if nbstrats_notnull < (nboutperf + nbunderperf)
        idx_underperf = ~ idx_outperf;
        % take some strategies from outperforming block to make up for missing underperforming strategies:
        % (Problem arises because we impose that the volatility of strategies we select as outperforming 
        % or underperforming is > thresh. With 20% outperforming, and 30% underperforming, we have to
        % make up for only 228 strategies, which is negligible.)
        nbmissingunderperf = nbunderperf - sum(idx_underperf);
        if nbmissingunderperf > nboutperf
            fprintf('Problem!\n')
        end
        idx_underperf = idx_underperf | (idx >= outperf_startidx) & (idx < outperf_startidx + nbmissingunderperf);
    else
        idx_underperf = ~ idx_outperf;
        tmp_a = cumsum(idx_underperf);
        tmp_b = find(tmp_a > nbunderperf, 1, 'first');
        idx_underperf(tmp_b:end) = 0;
    end
end

% select null strategies: 
idx_null = ~(idx_outperf | idx_underperf);
AA_sigmazero = AA_zero(:, std(AA_zero) <= thresh);
nbexcessnull = nbmissingunderperf;
AA_null = [AA_zero_notnull(:, idx_null) AA_sigmazero(:, 1:(size(AA_sigmazero, 2) - nbexcessnull))];

% shift outperforming and underperforming trajectories:
AA_outperf = adjustperf_mu(AA_zero_notnull(:, idx_outperf), outperf_mu);
AA_underperf =  adjustperf_mu(AA_zero_notnull(:, idx_underperf), underperf_mu);

AA_final = [AA_outperf AA_underperf AA_null];





