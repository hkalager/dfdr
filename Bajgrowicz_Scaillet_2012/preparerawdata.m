clear all
close all
clc

nbstrats = 7846;
AA = csvread('retsstrats_mu_new.txt');
nbdays_tmp = length(AA) / nbstrats;
AA = reshape(AA, nbdays_tmp, nbstrats);

nbdays = 126;
tt = 450;
MCSAMPLE = AA((tt-nbdays+1):tt, :);

save('MCSAMPLE', 'MCSAMPLE')

