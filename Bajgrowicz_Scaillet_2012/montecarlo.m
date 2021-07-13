% execute the following lines once to setup parallel computing environment:

% matlabpool close force local
% matlabpool close force
% matlabpool close
% matlabpool open 4
% matlabpool open 2

% add files for java progress bar (from matlab central):
% addpath C:\UNIGE\projet1\matlab\mc
% pctRunOnAll javaaddpath C:\UNIGE\projet1\matlab\mc
% addpath D:\UNIGE\matlab\mc
% pctRunOnAll javaaddpath D:\UNIGE\matlab\mc

clc
clear all
close all

tic

% INPUT PARAMETERS:
MC = 1000;
B = 1000;

% perf='_SR';
perf=''; % mu

% outperf_SR = 4;
% underperf_SR = -3;
outperf_mu = .3;
underperf_mu = -.1;

pi_aplus = .2;
pi_aminus = .3;
pi_0 = 1 - pi_aplus - pi_aminus;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q = 1/10;
sigma_thresh = 1e-3;

if strcmp(perf,'_SR')
    load(['AA_final_' num2str(100*pi_aplus) '_' num2str(100*pi_aminus) '_SR_'...
        num2str(100*outperf_SR) '_' num2str(100*(-underperf_SR)) '.mat'])
else
    load(['AA_final_' num2str(100*pi_aplus) '_' num2str(100*pi_aminus) '_mu_'...
        num2str(100*outperf_mu) '_' num2str(100*(-underperf_mu)) '.mat'])
end

nbstrats = 7846;
nbdays = 126;

nboutperf = round(nbstrats * pi_aplus);
nbunderperf = round(nbstrats * pi_aminus);
nbnull = nbstrats - nboutperf - nbunderperf;

RETS=zeros(nbdays, nbstrats);

% ***************************************************************
% compute p-values by bootstrap
IDX=zeros(nbdays,B);
for b=1:B
    IDX(1,b)=unidrnd(nbdays);
    U=unifrnd(0,1,nbdays-1,1);
    ind1=find(U<q)+1;
    ind2=find(U>=q)+1;
    IDX(ind1,b)=unidrnd(nbdays,length(ind1),1);
    for i=1:length(ind2)
        t=ind2(i);
        IDX(t,b)=IDX(t-1,b)+1;
        if IDX(t,b)>nbdays
            IDX(t,b)=1;
        end
    end
end

IDX_MC=zeros(nbdays,MC);
for mc=1:MC
    IDX_MC(1,mc)=unidrnd(nbdays);
    U=unifrnd(0,1,nbdays-1,1);
    ind1=find(U<q)+1;
    ind2=find(U>=q)+1;
    IDX_MC(ind1,mc)=unidrnd(nbdays,length(ind1),1);
    for i=1:length(ind2)
        t=ind2(i);
        IDX_MC(t,mc)=IDX_MC(t-1,mc)+1;
        if IDX_MC(t,mc)>nbdays
            IDX_MC(t,mc)=1;
        end
    end
end

pi_0hat=zeros(MC,1);
pi_aplushat=zeros(MC,1);
pi_aminushat=zeros(MC,1);
FDRhat20=zeros(MC,1);
FDRhatb=zeros(MC,1);
FDRrealFDRb=zeros(MC,1);
FDRrealFDR20=zeros(MC,1);
FDRrealRW5=zeros(MC,1);
FDRrealRW20=zeros(MC,1);
powerFDRb=zeros(MC,1);
powerFDR20=zeros(MC,1);
powerRW5=zeros(MC,1);
powerRW20=zeros(MC,1);
portsizeFDRb=zeros(MC,1);
portsizeFDR20=zeros(MC,1);
portsizeRW5=zeros(MC,1);
portsizeRW20=zeros(MC,1);

CC_outperf = zeros(3, MC);
CC_underperf = zeros(3, MC);

rankingbestlucky=zeros(MC,1);

% for i=1:MC
%ppm = ParforProgMon('Bob', MC);
parfor i=1:MC
    disp(i)
    %     simul by bootstrap:
    RETS = AA_final(IDX_MC(:, i), :);
    
    if strcmp(perf,'_SR')
        mu_final = mean(RETS);
        BB_outperf = 252 * mu_final(1:nboutperf);
        BB_underperf = 252 * mu_final((nboutperf + 1):(nboutperf + nbunderperf));
        BB_null = 252 * mu_final((nboutperf + nbunderperf + 1):end);
    else
        SR_final = mean(RETS) ./ std(RETS);
        BB_outperf = sqrt(252) * SR_final(1:nboutperf);
        BB_underperf = sqrt(252) * SR_final((nboutperf + 1):(nboutperf + nbunderperf));
        BB_null = sqrt(252) * SR_final((nboutperf + nbunderperf + 1):end);
    end
    CC_outperf(:, i) = [quantile(BB_outperf, .25); quantile(BB_outperf, .5); quantile(BB_outperf, .75)];
    CC_underperf(:, i) = [quantile(BB_underperf, .75); quantile(BB_underperf, .5); quantile(BB_underperf, .25)];
    
    Perfs_B=zeros(nbstrats,B);
    for b=1:B
        RETS_B = RETS(IDX(:, b), :);
        if strcmp(perf,'_SR')
            sigma = std(RETS_B)';
            mu = mean(RETS_B)';
            Perfs_B(:, b) = mu ./ sigma;
            
            % if sigma == 0 => Perfs_B = 0:
            Perfs_B(sigma == 0, b) = 0;
            
            % avoid situation with sigma too small:
            idx_tmp = sigma > 0 & sigma < sigma_thresh;
            Perfs_B(idx_tmp, b) = mu(idx_tmp) / sigma_thresh;
        else
            Perfs_B(:, b) = mean(RETS_B)';
        end
    end
    
    if strcmp(perf,'_SR')
        sigma = std(RETS)';
        mu = mean(RETS)';
        Perfs = mu ./ sigma;
        
        % if sigma == 0 => Perfs = 0:
        Perfs(sigma == 0) = 0;
        
        % avoid situation with sigma too small:
        idx_tmp = sigma > 0 & sigma < sigma_thresh;
        Perfs(idx_tmp) = mu(idx_tmp) / sigma_thresh;
    else
        Perfs = mean(RETS)';
    end
    
    pvalues = compute_pvalues(Perfs, Perfs_B, std(RETS));
    
    lambda = .6;
    gamma = .4;
    pi_0hat(i) = compute_pi_0hat(pvalues, lambda);
    [pi_aplushat(i), pi_aminushat(i)] = compute_pi_ahat(pvalues, Perfs, pi_0hat(i), gamma);
    
    [PORTFDR20, FDRhat20(i)] = portfolio_FDR(.2, Perfs, pvalues, pi_0hat(i));
    [PORTFDRb, FDRhatb(i)] = portfolio_FDR(.1, Perfs, pvalues, pi_0hat(i));

    PORTRW20 = portfolio_RW(.2, Perfs, Perfs_B);
    PORTRW5 = portfolio_RW(.05, Perfs, Perfs_B);
    
    % ranking of best lucky rule:
    idx = (1:nbstrats)';
    A = [Perfs idx];
    A = sortrows(A, -1);
    rankingbestlucky(i) = find(A(:, 2)  > nboutperf, 1, 'first');
    
    [FDRrealFDRb(i), portsizeFDRb(i)] = ComputeRealFDR(PORTFDRb, nboutperf);
    [FDRrealFDR20(i), portsizeFDR20(i)] = ComputeRealFDR(PORTFDR20, nboutperf);
    [FDRrealRW5(i), portsizeRW5(i)] = ComputeRealFDR(PORTRW5, nboutperf);
    [FDRrealRW20(i), portsizeRW20(i)] = ComputeRealFDR(PORTRW20, nboutperf);

    powerFDRb(i)=sum(PORTFDRb(1:nboutperf))/nboutperf;
    powerFDR20(i)=sum(PORTFDR20(1:nboutperf))/nboutperf;
    powerRW5(i)=sum(PORTRW5(1:nboutperf))/nboutperf;
    powerRW20(i)=sum(PORTRW20(1:nboutperf))/nboutperf;
    
    % java advance bar:
    %ppm.increment();
end

AA=[mean(FDRhatb) std(FDRhatb); mean(FDRhat20) std(FDRhat20)]';

BB=[mean(FDRrealFDRb) std(FDRrealFDRb);...
    mean(FDRrealFDR20) std(FDRrealFDR20);...
    mean(FDRrealRW5) std(FDRrealRW5);...
    mean(FDRrealRW20) std(FDRrealRW20)];

CC=[mean(powerFDRb) std(powerFDRb);...
    mean(powerFDR20) std(powerFDR20);...
    mean(powerRW5) std(powerRW5);...
    mean(powerRW20) std(powerRW20)];

DD=[mean(portsizeFDRb) std(portsizeFDRb);...
    mean(portsizeFDR20) std(portsizeFDR20);...
    mean(portsizeRW5) std(portsizeRW5);...
    mean(portsizeRW20) std(portsizeRW20)];

EE=[mean(pi_0hat) std(pi_0hat);...
    mean(pi_aplushat) std(pi_aplushat);...
    mean(pi_aminushat) std(pi_aminushat)];

ratio_pi_a = ((mean(pi_aplushat) + mean(pi_aminushat)) / (1 - mean(pi_0hat)))^-1;

FF=[mean(pi_0hat) std(pi_0hat);...
    mean(pi_aplushat) * ratio_pi_a std(pi_aplushat) * ratio_pi_a;...
    mean(pi_aminushat) * ratio_pi_a std(pi_aminushat) * ratio_pi_a];

GG=[median(rankingbestlucky) mean(rankingbestlucky) std(rankingbestlucky)];

HH=[BB CC DD];

II = [mean(CC_outperf, 2) std(CC_outperf, 0, 2) mean(CC_underperf, 2) std(CC_underperf, 0, 2)];

if strcmp(perf,'_SR')
    save(['MCresults_' num2str(100*pi_aplus) '_' num2str(100*pi_aminus) '_SR_'...
        num2str(100*outperf_SR) '_' num2str(100*(-underperf_SR)) '.mat'])
else
    save(['MCresults_' num2str(100*pi_aplus) '_' num2str(100*pi_aminus) '_mu_'...
        num2str(100*outperf_mu) '_' num2str(100*(-underperf_mu)) '.mat'])
end

toc


