% Created 10 Oct 2017
% Last Modified 19 Oct 2017 12:37
clear;
clc;
labels={'NDDUUS'	'NDDUUK'	'NDDUJN'	'NDEUSRU'	'NDEUCHF'	'NDUEBRAF'...
    'MSEIESUN'	'M1MA'	'MIMUJORN'};
rftable=readtable('FEDL01.csv');
rftable{:,2}=rftable{:,2}/100/260;
transser=[25*ones(1,3),50*ones(1,6)];
delete(gcp('nocreate'))
disp('Checking the datasets');
start=2002;
finish=2016;
oosperiod=1; % Out of sample holding period
% Type of the performance measure (Simple mean or Sharpe ratio (studentized))
Studentized = true;
sigmathresh=1e-3;
% Bootstrap specification
Bsize=1000;
Bwindow=10;
%% Liang specification
Max_lambda=0.99;
N_bins=20;
gamma = 0.4;
%% FDR setting
fdrtarget=0.1;
%% switch function between STUDENTIZED and GENERAL test statistic
whichstd=@(x,y) y*std(x)+(1-y)*ones(1,size(x,2));
yearser=2006:2015;
tic
% parfor i=1:numel(labels)
for i=1:numel(labels)
    lbl=labels{i};
    transbasis=transser(i);
    fllbl=[lbl,'_',num2str(start),'_',num2str(finish),'-insample-',num2str(transbasis),'.mat'];
    if ~exist(fllbl,'file')
        disp(['Dataset is not available and will be generated for ',lbl,' ...']);
        rulegenerator21kR(lbl,transbasis,1,start,finish); % 1 stands in-sample
    else
        disp(['Dataset ready for ',lbl,' and being loaded...']);
        
    end
end