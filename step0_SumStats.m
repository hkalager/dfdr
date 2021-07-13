%% Summary statistics for securities under study
% Script developed by Arman Hassannia Kalager for completing 3rd chapter,
% Created on 02 Mar 2018,
% Last modified 02 Mar 2018 16:47 GMT.

%% Basics 
clear;
clc;
labels={'MXWO'  'NDDUUS'    'NDDUUK'	'NDDUJN'    ...
    'MXEF'  'NDUEBRAF'  'NDEUCHF'  'NDEUSRU'	...
    'MSEIESUN'	'M1MA'	'MIMUJORN'  'MXFEM'};
rftable=readtable('FEDL01.csv');
rftable{:,2}=rftable{:,2}/100/260;
start=2004;
finish=2016;
tbl1=table();
%% loading the MSCI indexes for loops
for i=1:numel(labels)
    lbl1=labels{i};
    sym=[lbl1,'.csv'];
    master=importdata(sym,',');
    TF1 = contains(master.textdata(:,1),mat2str(start));
    TF1=find(TF1==1,1,'first');
    startdt=master.textdata(TF1,1);
    TF2 = contains(master.textdata(:,1),mat2str(finish));
    TF2=find(TF2==1,1,'last');
    finishdt=master.textdata(TF2,1);
    priceser=master.data(find(strcmp(master.textdata,startdt)):find(strcmp(master.textdata,finishdt)),max(1,end-1));
    data=price2ret(priceser);
    tbl1(i,'Market')={lbl1}; 
    tbl1{i,'Mean'}=mean(data)*1e2;
    tbl1{i,'Max'}=max(data)*1e2;
    tbl1{i,'Min'}=min(data)*1e2;
    tbl1{i,'Std'}=std(data)*1e2;
    tbl1{i,'Kurt'}=kurtosis(data);
    tbl1{i,'Skew'}=skewness(data);
    mdl=arima(1,0,0);
    [estmdl,a,b,c]=estimate(mdl,data);
    cell1=estmdl.AR;
    tbl1{i,'FirstAC'}=cell1{:};   
    [h,pval] = archtest(data);
    tbl1{i,'pval'}=pval;
end
%% Effective Federal Funds rate
i=i+1;
rfind1=find(rftable{:,1}==char(startdt),1,'first');
rfind2=find(rftable{:,1}==char(finishdt),1,'first');
effr=rftable{rfind1:rfind2,2};
tbl1{i,'Market'}={'FFR'};
tbl1{i,'Mean'}=mean(effr)*1e2;
tbl1{i,'Max'}=max(effr)*1e2;
tbl1{i,'Min'}=min(effr)*1e2;
tbl1{i,'Std'}=std(effr)*1e2;
tbl1{i,'Kurt'}=kurtosis(effr);
tbl1{i,'Skew'}=skewness(effr);
mdl=arima(1,0,0);
estmdl=estimate(mdl,effr);
cell1=estmdl.AR;
tbl1{i,'FirstAC'}=cell1{:};
[h,pval] = archtest(data);
tbl1{i,'pval'}=pval;
%% Recording to the table
writetable(tbl1,['SummaryStats_',num2str(start),'-',num2str(finish),'.xlsx'])