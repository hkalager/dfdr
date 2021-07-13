%  Monte Carlo simulation script for the new data-snooping paper
%  Script developed by Arman Hassanniakalager for the 3rd chapter of PhD,
%  Based on the BS, 2012 paper
%  Last modified 01 Feb 2018 12:13 GMT.
addpath(genpath(cd));
clc
clear
close all
delete(gcp('nocreate'))
%% Options
% BS+BS= Original BS= BB
% BS+LIANG= BL
pi0typ={'BL','BB'}; 
%% Initiating parallel pool
tic
%poolobj=parpool('EC2EU',18,'AttachedFiles','Arman_montecarlo4.m');
poolobj=parpool('local', 4);
%% INPUT PARAMETERS:
MC =1000; % Number of simulations
disp(['Running ',num2str(MC),' Monte Carlo simulation...']);
B = 1000; % Number of bootstrap
%% Parameters
gamma=0.4;
% BB setting
lambda=0.55;
% BL setting
lambda_max_BL=0.99;
N_BL=10;
% Sharpe Ratio or Mu
perf='_SR';
%% Loops
for outer=4:-1:2
    for under=-4:1:-2
        
        outperf_SR = outer;
        underperf_SR = under;
        %  perf=''; % mu
        %  outperf_mu = .3;
        %  underperf_mu = -.1;
        disp(['Running performance measure ',perf,' outperformers ',...
            num2str(outer),' underperformer ' num2str(under)]);
        pi_aplus = .2;
        pi_aminus = .3;
        pi_0 = 1 - pi_aplus - pi_aminus;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        q = 1/10;
        sigma_thresh = 5e-4;
        
        if strcmp(perf,'_SR')
            nameSR=['AA_final_' num2str(100*pi_aplus) '_' num2str(100*pi_aminus) '_SR_'...
                num2str(100*outperf_SR) '_' num2str(100*(-underperf_SR)) '.mat'];
            if ~exist(nameSR,'file')
                prepare_SR_Arman(outperf_SR,underperf_SR,pi_aplus,pi_aminus);
            end
            load(nameSR)
        else
            namemu=['AA_final_' num2str(100*pi_aplus) '_' num2str(100*pi_aminus) '_mu_'...
                num2str(100*outperf_mu) '_' num2str(100*(-underperf_mu)) '.mat'];
            if ~exist(namemu,'file')
                prepare_mu_Arman(outperf_mu,underperf_mu,pi_aplus,pi_aminus);
            end
            load(namemu)
        end
        
        nbstrats = 21195;
        nbdays = 155;
        
        nboutperf = round(nbstrats * pi_aplus);
        nbunderperf = round(nbstrats * pi_aminus);
        nbnull = nbstrats - nboutperf - nbunderperf;
        
        RETS=zeros(nbdays, nbstrats);
        
        finalret=sum(AA_final);
        %% Bootstrap generation
        %  Firstly for the main index
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
        % Secondly for the whole simulation
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
        %% Declaring the variables
        pi_0hatBB=zeros(MC,1);
        pi_0hatBL=zeros(MC,1);
        lambdaBL=zeros(MC,1);
        pi_aplushatBB=zeros(MC,1);
        pi_aminushatBB=zeros(MC,1);
        pi_aplushatBL=zeros(MC,1);
        pi_aminushatBL=zeros(MC,1);
        FDRhat20BB=zeros(MC,1);
        FDRhatbBB=zeros(MC,1);
        FDRhat20BL=zeros(MC,1);
        FDRhatbBL=zeros(MC,1);
        FDRrealFDRbBB=zeros(MC,1);
        FDRrealFDR20BB=zeros(MC,1);
        FDRrealFDRbBL=zeros(MC,1);
        FDRrealFDR20BL=zeros(MC,1);
        FDRrealRW5=zeros(MC,1);
        FDRrealRW20=zeros(MC,1);
        PORTFDR20BB=nan(nbstrats,MC);
        PORTFDRbBB=nan(nbstrats,MC);
        PORTFDR20BL=nan(nbstrats,MC);
        PORTFDRbBL=nan(nbstrats,MC);
        powerFDRbBB=zeros(MC,1);
        powerFDR20BB=zeros(MC,1);
        powerFDRbBL=zeros(MC,1);
        powerFDR20BL=zeros(MC,1);
        powerRW5=zeros(MC,1);
        powerRW20=zeros(MC,1);
        portsizeFDRbBB=zeros(MC,1);
        portsizeFDR20BB=zeros(MC,1);
        portsizeFDRbBL=zeros(MC,1);
        portsizeFDR20BL=zeros(MC,1);
        portsizeRW5=zeros(MC,1);
        portsizeRW20=zeros(MC,1);
        pvalues=nan(nbstrats,MC);
        ISFDRbBB=nan(1,MC);
        ISFDR20BB=nan(1,MC);
        ISFDRbBL=nan(1,MC);
        ISFDR20BL=nan(1,MC);
        RSRW5=nan(1,MC);
        RSRW20=nan(1,MC);
        CC_outperf = zeros(3, MC);
        CC_underperf = zeros(3, MC);
        
        rankingbestlucky=zeros(MC,1);
        
        %% Performing the simulation
        parfor i=1:MC
            
            disp(i)
            tic;
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
            
            pvalues(:,i) = compute_pvalues(Perfs, Perfs_B, std(RETS));
            
            % BS Procedure for calculation of pi0
            lambdaBB = lambda;
            pi_0hatBB(i) = compute_pi_0hat(pvalues(:,i), lambdaBB);
            gammaBB = gamma;
            [pi_aplushatBB(i), pi_aminushatBB(i)] = compute_pi_ahat(pvalues(:,i), Perfs, pi_0hatBB(i), gammaBB);
            % Liang Procedure
            [pi_0hatBL(i),lambdaBL(i)]=est_pi0_disc(pvalues(:,i), N_BL,lambda_max_BL);
            gammaBL = gammaBB;
            [pi_aplushatBL(i), pi_aminushatBL(i)] = compute_pi_ahat(pvalues(:,i), Perfs, pi_0hatBL(i), gammaBL);
            
            % Portfolio construction based on the FDR+        
            [PORTFDR20BB(:,i), FDRhat20BB(i)] = my_portfolio_FDR(.2, Perfs, pvalues(:,i), pi_0hatBB(i));
            [PORTFDRbBB(:,i), FDRhatbBB(i)] = my_portfolio_FDR(.1, Perfs, pvalues(:,i), pi_0hatBB(i));
            [PORTFDR20BL(:,i), FDRhat20BL(i)] = my_portfolio_FDR(.2, Perfs, pvalues(:,i), pi_0hatBL(i));
            [PORTFDRbBL(:,i), FDRhatbBL(i)] =my_portfolio_FDR(.1, Perfs, pvalues(:,i), pi_0hatBL(i));
            
            PORTRW20(:,i) = portfolio_RW(.2, Perfs, Perfs_B);
            PORTRW5(:,i) = portfolio_RW(.05, Perfs, Perfs_B);
            
            % ranking of best lucky rule:
            idx = (1:nbstrats)';
            A = [Perfs idx];
            A = sortrows(A, -1);
            rankingbestlucky(i) = find(A(:, 2)  > nboutperf, 1, 'first');
            
            [FDRrealFDRbBB(i), portsizeFDRbBB(i)] = ComputeRealFDR(PORTFDRbBB(:,i), nboutperf);
            [FDRrealFDR20BB(i), portsizeFDR20BB(i)] = ComputeRealFDR(PORTFDR20BB(:,i), nboutperf);
            
            [FDRrealFDRbBL(i), portsizeFDRbBL(i)] = ComputeRealFDR(PORTFDRbBL(:,i), nboutperf);
            [FDRrealFDR20BL(i), portsizeFDR20BL(i)] = ComputeRealFDR(PORTFDR20BL(:,i), nboutperf);
            
            [FDRrealRW5(i), portsizeRW5(i)] = ComputeRealFDR(PORTRW5(:,i), nboutperf);
            [FDRrealRW20(i), portsizeRW20(i)] = ComputeRealFDR(PORTRW20(:,i), nboutperf);
            
            
%             % In-sample return
%             ISFDRbBB(i)=finalret*PORTFDRbBB(:,i)/portsizeFDRbBB(i);
%             ISFDR20BB(i)=finalret*PORTFDR20BB(:,i)/portsizeFDR20BB(i);
%             ISFDRbBL(i)=finalret*PORTFDRbBL(:,i)/portsizeFDRbBL(i);
%             ISFDR20BL(i)=finalret*PORTFDR20BL(:,i)/portsizeFDR20BL(i);
%             ISRW5(i)=finalret*PORTRW5/max(portsizeRW5(i),1);
%             ISRW20(i)=finalret*PORTRW20/max(portsizeRW20(i),1);
            toc;
        end
        
        for i=1:MC
            powerFDRbBB(i)=sum(PORTFDRbBB(1:nboutperf,i))/nboutperf;
            powerFDR20BB(i)=sum(PORTFDR20BB(1:nboutperf,i))/nboutperf;
            powerFDRbBL(i)=sum(PORTFDRbBL(1:nboutperf,i))/nboutperf;
            powerFDR20BL(i)=sum(PORTFDR20BL(1:nboutperf,i))/nboutperf;
            powerRW5(i)=sum(PORTRW5(1:nboutperf,i))/nboutperf;
            powerRW20(i)=sum(PORTRW20(1:nboutperf,i))/nboutperf;
        end
        %% Outputs
        AA_BB=[mean(FDRhatbBB) std(FDRhatbBB); mean(FDRhat20BB) std(FDRhat20BB)]';
        AA_BL=[mean(FDRhatbBL) std(FDRhatbBL); mean(FDRhat20BL) std(FDRhat20BL)]';
        
        BB_BB=[mean(FDRrealFDRbBB) std(FDRrealFDRbBB);...
            mean(FDRrealFDR20BB) std(FDRrealFDR20BB);...
            mean(FDRrealRW5) std(FDRrealRW5);...
            mean(FDRrealRW20) std(FDRrealRW20)];
        BB_BL=[mean(FDRrealFDRbBL) std(FDRrealFDRbBL);...
            mean(FDRrealFDR20BL) std(FDRrealFDR20BL);...
            mean(FDRrealRW5) std(FDRrealRW5);...
            mean(FDRrealRW20) std(FDRrealRW20)];
        
        CC_BB=[mean(powerFDRbBB) std(powerFDRbBB);...
            mean(powerFDR20BB) std(powerFDR20BB);...
            mean(powerRW5) std(powerRW5);...
            mean(powerRW20) std(powerRW20)];
        CC_BL=[mean(powerFDRbBL) std(powerFDRbBL);...
            mean(powerFDR20BL) std(powerFDR20BL);...
            mean(powerRW5) std(powerRW5);...
            mean(powerRW20) std(powerRW20)];
        
        DD_BB=[mean(portsizeFDRbBB) std(portsizeFDRbBB);...
            mean(portsizeFDR20BB) std(portsizeFDR20BB);...
            mean(portsizeRW5) std(portsizeRW5);...
            mean(portsizeRW20) std(portsizeRW20)];
        DD_BL=[mean(portsizeFDRbBL) std(portsizeFDRbBL);...
            mean(portsizeFDR20BL) std(portsizeFDR20BL);...
            mean(portsizeRW5) std(portsizeRW5);...
            mean(portsizeRW20) std(portsizeRW20)];
        
        EE_BB=[mean(pi_0hatBB) std(pi_0hatBB);...
            mean(pi_aplushatBB) std(pi_aplushatBB);...
            mean(pi_aminushatBB) std(pi_aminushatBB)];
        EE_BL=[mean(pi_0hatBL) std(pi_0hatBL);...
            mean(pi_aplushatBL) std(pi_aplushatBL);...
            mean(pi_aminushatBL) std(pi_aminushatBL)];
        
        ratio_pi_a_BB = ((mean(pi_aplushatBB) + mean(pi_aminushatBB)) / (1 - mean(pi_0hatBB)))^-1;
        ratio_pi_a_BL = ((mean(pi_aplushatBL) + mean(pi_aminushatBL)) / (1 - mean(pi_0hatBL)))^-1;

        FF_BB=[mean(pi_0hatBB) std(pi_0hatBB);...
            mean(pi_aplushatBB) * ratio_pi_a_BB std(pi_aplushatBB) * ratio_pi_a_BB;...
            mean(pi_aminushatBB) * ratio_pi_a_BB std(pi_aminushatBB) * ratio_pi_a_BB];
        FF_BL=[mean(pi_0hatBL) std(pi_0hatBL);...
            mean(pi_aplushatBL) * ratio_pi_a_BL std(pi_aplushatBL) * ratio_pi_a_BL;...
            mean(pi_aminushatBL) * ratio_pi_a_BL std(pi_aminushatBL) * ratio_pi_a_BL];
        
        GG=[median(rankingbestlucky) mean(rankingbestlucky) std(rankingbestlucky)];
        
        HH_BB=[BB_BB CC_BB DD_BB];
        HH_BL=[BB_BL CC_BL DD_BL];
        
        II = [mean(CC_outperf, 2) std(CC_outperf, 0, 2) mean(CC_underperf, 2) std(CC_underperf, 0, 2)];
        KK=mean(lambdaBL);
        
%         LL=[mean(ISFDRbBL),std(ISFDRbBL);...
%             mean(ISFDRbBB),std(ISFDRbBB);...
%             mean(ISFDR20BL),std(ISFDR20BL);...
%             mean(ISFDR20BB),std(ISFDR20BB);...
%             mean(ISRW5),std(ISRW5);...
%             mean(ISRW20),std(ISRW20)];
        
        %% Write to file
        if exist('adj','var')
            switch adj
                case 1
                    adjlbl='UNADJ';
                case 2
                    adjlbl='ADJ';
            end
        else
            adjlbl='RAW';
        end
        suff=['-BB_BL-' num2str(N_BL)];
        if strcmp(perf,'_SR')
            save(['MCresults_' num2str(MC) '_' num2str(100*pi_aplus) '_' num2str(100*pi_aminus) '_SR_'...
                num2str(100*outperf_SR) '_' num2str(100*(-underperf_SR)) suff '.mat'])
        else
            save(['MCresults_' num2str(MC) '_' num2str(100*pi_aplus) '_' num2str(100*pi_aminus) '_mu_'...
                num2str(100*outperf_mu) '_' num2str(100*(-underperf_mu)) suff '.mat'])
        end
    end
end
toc
delete(poolobj);