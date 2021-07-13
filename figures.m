figure,
load('MCresults_1000_10_20_SR_400_400-BS-BS-RAW.mat');
subplot(1,2,1)
histogram(pi_0hat);title('BS-BS $$\hat{\pi}$$0 Histogram ($${\pi}$$0=0.7) 1000 MC',...
    'Interpreter','Latex')
load('MCresults_1000_10_20_SR_400_400-BS-LIANG-RAW.mat')
subplot(1,2,2)
histogram(pi_0hat);title('BS-Liang $$\hat{\pi}$$0 Histogram ($${\pi}$$0=0.7) 1000 MC',...
    'Interpreter','Latex')

figure,
histogram(lambdaLIANG);title('$$\hat{\lambda}_{BL}$$ Histogram','Interpreter','Latex')

