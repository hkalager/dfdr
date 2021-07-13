% Based on a code created 10 Oct 2017
% Last Modified 28 August 2018 16:18 IRST
clear;
clc;
labels={'MXWO' 'NDDUUS'	'NDDUUK'	'NDDUJN'	'MXEF' 'NDEUSRU'	'NDEUCHF'	'NDUEBRAF'...
    'MXFEM' 'MSEIESUN'	'M1MA'	'MIMUJORN'};
regions={'Advanced' 'US' 'UK' 'Japan' 'Emerging' 'Russia' 'China' 'Brazil' ...
    'Frontier' 'Estonia' 'Morocco' 'Jordan'};

rftable=readtable('FEDL01.csv');
rftable{:,2}=rftable{:,2}/100/260;
transser=[25*ones(1,4),50*ones(1,8)];
delete(gcp('nocreate'))
%poolobj=parpool('EC2',16,'AttachedFiles','ModifiedRunner.m');
%poolobj=parpool('EC2compatibility',18,'AttachedFiles','ModifiedRunner.m');
poolobj=parpool('local',8);
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
%% FDR setting
fdrtarget=0.1;
%% Years range
yearser=2006:2015;
tic
for i=5%1:numel(labels)
    lbl=labels{i};
    transbasis=transser(i);
    fllbl=[lbl,'_',num2str(start),'_',num2str(finish),'-insample-',num2str(transbasis),'.mat'];
    if ~exist(fllbl,'file')
        disp(['Dataset is not available and will be generated for ',lbl,' ...']);
        rulegenerator21kR(lbl,transbasis,1,start,finish); % 1 stands in-sample
    else
        disp(['Dataset ready for ',lbl,' and being loaded...']);
        
    end
    load(fllbl,'master','ret','TF1');
    % Calculating the returns without transaction costs
    disp(['Dataset loaded for ',lbl]);
    iter=0;
    tic;
    for yearOOS=yearser(1):yearser(end)%2006:2016
        for month=1:12
            for IS=2
                
                year=yearOOS-IS;
                begindtstrIS=[num2str(month),'/',num2str(year)];
                TF1SMP = contains(master.textdata(:,1),begindtstrIS);
                TF1SMP=find(TF1SMP==1,1,'first');
                startdtIS=master.textdata(TF1SMP,1);
                finishdtstrIS=[num2str(month),'/',num2str(year+IS)];
                TF2SMP = contains(master.textdata(:,1),finishdtstrIS);
                TF2SMP=find(TF2SMP==1,1,'first')-1;
                finishdtIS=master.textdata(TF2SMP,1);
                startdtOOS=master.textdata(TF2SMP+1,1);
                %% 1M OOS
                finishdtstrOOS1M=[num2str(rem(month+oosperiod,12+1e-10)),'/',num2str(year+IS+(month+oosperiod>12))];
                TF3SMP1M = contains(master.textdata(:,1),finishdtstrOOS1M);
                TF3SMP1M=find(TF3SMP1M==1,1,'first')-1;
                finishdtOOS1M=master.textdata(TF3SMP1M,1);
                %% 3M OOS
                finishdtstrOOS3M=[num2str(rem(month+3*oosperiod,12+1e-10)),'/',num2str(year+IS+(month+3*oosperiod>12))];
                TF3SMP3M = contains(master.textdata(:,1),finishdtstrOOS3M);
                TF3SMP3M=find(TF3SMP3M==1,1,'first')-1;
                finishdtOOS3M=master.textdata(TF3SMP3M,1);
                %% 6M OOS
                finishdtstrOOS6M=[num2str(rem(month+6*oosperiod,12+1e-10)),'/',num2str(year+IS+(month+6*oosperiod>12))];
                TF3SMP6M = contains(master.textdata(:,1),finishdtstrOOS6M);
                TF3SMP6M=find(TF3SMP6M==1,1,'first')-1;
                finishdtOOS6M=master.textdata(TF3SMP6M,1);
                % IS return
                retIS=ret(TF1SMP-TF1+1:TF2SMP-TF1+1,:);
                %% OOS return
                % 1M
                retOOS1M=ret(TF2SMP+1-TF1+1:TF3SMP1M-TF1+1,:);
                % 3M
                retOOS3M=ret(TF2SMP+1-TF1+1:TF3SMP3M-TF1+1,:);
                % 6M
                retOOS6M=ret(TF2SMP+1-TF1+1:TF3SMP6M-TF1+1,:);
                %% Excess return calculation
                % IS
                rfind1=find(rftable{:,1}==char(startdtIS),1,'first');
                rfind2=find(rftable{:,1}==char(finishdtIS),1,'first');
                % OOS
                rfind31M=find(rftable{:,1}==char(finishdtOOS1M),1,'first');
                rfind33M=find(rftable{:,1}==char(finishdtOOS3M),1,'first');
                rfind36M=find(rftable{:,1}==char(finishdtOOS6M),1,'first');
                excessretIS=retIS-rftable{rfind1:rfind2,2};
                excessretOOS1M=retOOS1M-rftable{rfind2+1:rfind31M,2};
                excessretOOS3M=retOOS3M-rftable{rfind2+1:rfind33M,2};
                excessretOOS6M=retOOS6M-rftable{rfind2+1:rfind36M,2};
                
                %% Generating the bootstraps
                indices=stationary_bootstrap((1:size(retIS,1))',Bsize,Bwindow);
                
                %% Generating performance index
                
                modelscount=size(excessretIS,2);
                Perf_B=zeros(Bsize,modelscount);
                Perf_B_Sort=zeros(Bsize,modelscount);
                Perf_B_Mu=zeros(Bsize,modelscount);
                disp(['Bootstrap generation initiated for IS ', ...
                    num2str(month),'/',num2str(year),'-',num2str(month),'/',num2str(yearOOS)]);

                
                mu=mean(excessretIS);
                sigma=std(excessretIS);                
                
                SortSig=zeros(1,modelscount);
                for s=1:modelscount
                    neg_ind=excessretIS(:,s)<0;
                    SortSig(s)=std(excessretIS(neg_ind,s));
                end
                
                Perf=mu./sigma;
                Perf_Sort=mu./SortSig;
                Perf_Mu=mu;
                parfor b=1:Bsize
                    bsdata=excessretIS(indices(:,b),:);
                    mu=mean(bsdata);
                    sigma=std(bsdata);
                    Perf_B(b,:)=mu./sigma;
                    Perf_B_Mu(b,:)=mu;
                    
                    for s=1:modelscount
                        neg_ind=bsdata(:,s)<0;
                        SortSig_B=std(excessretIS(neg_ind,s));
                        Perf_B_Sort(b,s)=mu(s)/SortSig_B;
                    end
                    
                end
                pvalues = compute_pvalues(Perf', Perf_B', std(excessretIS));
                pvalues_Mu = compute_pvalues(Perf_Mu', Perf_B_Mu', std(excessretIS));
                pvalues_Sort = compute_pvalues(Perf_Sort', Perf_B_Sort', std(excessretIS));
                
                % Liang procedure
                [pi_0hat,lambda]=est_pi0_disc(pvalues, N_bins,Max_lambda);
                [pi_0hat_Mu,lambda_Mu]=est_pi0_disc(pvalues_Mu, N_bins,Max_lambda);
                [pi_0hat_Sort,lambda_Sort]=est_pi0_disc(pvalues_Sort, N_bins,Max_lambda);
                % Barras et al 2010 estimation of FDR+/-
                [PORTFDR, FDRhat] = my_portfolio_FDR(fdrtarget, Perf', pvalues, pi_0hat);
                [PORTFDR_Mu, FDRhat_Mu] = my_portfolio_FDR(fdrtarget, Perf_Mu', pvalues_Mu, pi_0hat_Mu);
                [PORTFDR_Sort, FDRhat_Sort] = my_portfolio_FDR(fdrtarget, Perf_Sort', pvalues_Sort, pi_0hat_Sort);
                
                %% IS and OOS return for outperformers
                % IS portfolio return
                finalretIS=sum(excessretIS);
                PortRetIS=finalretIS*PORTFDR/max(sum(PORTFDR),1);
                PortRetIS_Mu=finalretIS*PORTFDR_Mu/max(sum(PORTFDR_Mu),1);
                PortRetIS_Sort=finalretIS*PORTFDR_Sort/max(sum(PORTFDR_Sort),1);
                % OOS portfolio return
                % 1M
                finalretOOS1M=sum(excessretOOS1M);
                PortRetOOS1M=finalretOOS1M*PORTFDR/max(sum(PORTFDR),1);
                PortRetOOS1M_Mu=finalretOOS1M*PORTFDR_Mu/max(sum(PORTFDR_Mu),1);
                PortRetOOS1M_Sort=finalretOOS1M*PORTFDR_Sort/max(sum(PORTFDR_Sort),1);
                
                %% Recording the results
                %iter=IS+(month-1)*2+(yearOOS-2006)*24;
                iter=iter+1;
                resvar{iter,1}=PORTFDR;
                resvar{iter,2}=PORTFDR_Mu;
                resvar{iter,3}=PORTFDR_Sort;
                resvar{iter,4}=pvalues;
                resvar{iter,5}=pvalues_Mu;
                resvar{iter,6}=pvalues_Sort;
                resvar{iter,7}=Perf;
                resvar{iter,8}=Perf_Mu;
                resvar{iter,9}=Perf_Sort;
                % Written to table
                resvar{iter,1+9}=yearOOS;
                resvar{iter,2+9}=sum(PORTFDR);
                resvar{iter,3+9}=sum(PORTFDR_Mu);
                resvar{iter,4+9}=sum(PORTFDR_Sort);
                resvar{iter,5+9}=(1+PortRetIS)^(1/IS)-1;
                resvar{iter,6+9}=(1+PortRetIS_Mu)^(1/IS)-1;
                resvar{iter,7+9}=(1+PortRetIS_Sort)^(1/IS)-1;
                resvar{iter,8+9}=(1+mean(PortRetOOS1M))^12-1;
                resvar{iter,9+9}=(1+mean(PortRetOOS1M_Mu))^12-1;
                resvar{iter,10+9}=(1+mean(PortRetOOS1M_Sort))^12-1;
                resvar{iter,11+9}=regions{i};
                resvar{iter,12+9}=month;
                toc
                %xlswrite([lbl,'xlsx'],result_table);
            end
        end
        
    end
    result_table=cell2table(resvar(:,10:end));
    result_table.Properties.VariableNames={'Year' ...
        'PortSizeSharpe' 'PortSizeMu' 'PortSizeSort' ...
        'FDRRetISSharpe' 'FDRRetISMu' 'FDRRetISSort' ...
        'FDRRetOOS1MSharpe' 'FDRRetOOS1MMu' 'FDRRetOOS1MSort' 'Country' 'Month'};

    %drpadd='C:\Dropbox\Dropbox\Shared_Folder\Chapter3\';
    %drpadd='C:\Users\Hassannia\Dropbox\Shared_Folder\Chapter3\';
    drpadd='C:\Users\Arman\Dropbox\Shared_Folder\Chapter3\';
    writetable(result_table,[drpadd,'Metric1','_',lbl,'_',num2str(yearser(1)),'_',num2str(yearser(end)),'.xlsx'])
    save(['Metric1','_',lbl,'_',num2str(yearser(1)),'_',num2str(yearser(end)),'.mat'],'resvar');
    
    toc
    clear resvar result_table;
end
delete(poolobj);