close all
clc

% INPUT PARAMETERS:
outperf_mu = .3;
underperf_mu = -.1;
pi_aplus = .2;
pi_aminus = .3;
%%%%%%%%%%%%%%%%% 

load('MCSAMPLE.mat')
AA = MCSAMPLE;

nbstrats = 7846;
nbdays = 126;
nboutperf = round(nbstrats * pi_aplus);
nbunderperf = round(nbstrats * pi_aminus);

idx = 1:497;
[AA_fr, nb_outperf_fr, nb_underperf_fr, nb_null_fr] = preparadata_auxfun_mu(AA(:, idx), ...
    outperf_mu, underperf_mu, pi_aplus, pi_aminus);

idx = 498:2546;
[AA_ma, nb_outperf_ma, nb_underperf_ma, nb_null_ma] = preparadata_auxfun_mu(AA(:, idx), ...
    outperf_mu, underperf_mu, pi_aplus, pi_aminus);

idx = 2547:3766;
[AA_sr, nb_outperf_sr, nb_underperf_sr, nb_null_sr] = preparadata_auxfun_mu(AA(:, idx), ...
    outperf_mu, underperf_mu, pi_aplus, pi_aminus);

idx = 3767:5806;
[AA_cb, nb_outperf_cb, nb_underperf_cb, nb_null_cb] = preparadata_auxfun_mu(AA(:, idx), ...
    outperf_mu, underperf_mu, pi_aplus, pi_aminus);

idx = 5807:7846;
[AA_obv, nb_outperf_obv, nb_underperf_obv, nb_null_obv] = preparadata_auxfun_mu(AA(:, idx), ...
    outperf_mu, underperf_mu, pi_aplus, pi_aminus);

AA_final = [AA_fr(:, 1:nb_outperf_fr)...
    AA_ma(:, 1:nb_outperf_ma)...
    AA_sr(:, 1:nb_outperf_sr) ...
    AA_cb(:, 1:nb_outperf_cb)...
    AA_obv(:, 1:nb_outperf_obv) ...
    AA_fr(:, (nb_outperf_fr + 1):(nb_outperf_fr + nb_underperf_fr))...
    AA_ma(:, (nb_outperf_ma + 1):(nb_outperf_ma + nb_underperf_ma))...
    AA_sr(:, (nb_outperf_sr + 1):(nb_outperf_sr + nb_underperf_sr)) ...
    AA_cb(:, (nb_outperf_cb + 1):(nb_outperf_cb + nb_underperf_cb))...
    AA_obv(:, (nb_outperf_obv + 1):(nb_outperf_obv + nb_underperf_obv)) ...
    AA_fr(:, (nb_outperf_fr + nb_underperf_fr + 1):end)...
    AA_ma(:, (nb_outperf_ma + nb_underperf_ma + 1):end)...
    AA_sr(:, (nb_outperf_sr + nb_underperf_sr + 1):end)...
    AA_cb(:, (nb_outperf_cb + nb_underperf_cb + 1):end)...
    AA_obv(:, (nb_outperf_obv + nb_underperf_obv + 1):end)];

save(['AA_final_' num2str(100*pi_aplus) '_' num2str(100*pi_aminus) '_mu_'...
        num2str(100*outperf_mu) '_' num2str(100*(-underperf_mu)) '.mat'], 'AA_final')




