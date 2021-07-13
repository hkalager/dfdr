% Created 10 Oct 2017 by Arman Hassanniakalager
% Last Modified 13 July 2021 

clear;
clc;
warning off;
addpath(genpath(cd));
labels={'NDDUUS'	'NDDUUK'	'NDDUJN'	'NDEUSRU'	'NDEUCHF'	'NDUEBRAF'...
    'MSEIESUN'	'M1MA'	'MIMUJORN' 'MXWO' 'MXEF' 'MXFEM'};
rftable=readtable('FEDL01.csv');
rftable{:,2}=rftable{:,2}/100/260;
transser=[25*ones(1,3),50*ones(1,6),25,50,50];
delete(gcp('nocreate'))
%poolobj=parpool('local',4);
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
%% Years range
yearser=2006:2015;
tic
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
    load(fllbl);
    % Calculating the returns without transaction costs
    retnocost=TransCost(inputser,dataret,0);
    disp(['Dataset loaded for ',lbl]);
    iter=0;
    tic;
    for yearOOS=yearser(1):yearser(end)%2006:2016
        for month=1:12
            for IS=1:2
                
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
                retISnocost=retnocost(TF1SMP-TF1+1:TF2SMP-TF1+1,:);
                %% OOS return
                % 1M
                retOOS1M=ret(TF2SMP+1-TF1+1:TF3SMP1M-TF1+1,:);
                retOOS1Mnocost=retnocost(TF2SMP+1-TF1+1:TF3SMP1M-TF1+1,:);
                % 3M
                retOOS3M=ret(TF2SMP+1-TF1+1:TF3SMP3M-TF1+1,:);
                retOOS3Mnocost=retnocost(TF2SMP+1-TF1+1:TF3SMP3M-TF1+1,:);
                % 6M
                retOOS6M=ret(TF2SMP+1-TF1+1:TF3SMP6M-TF1+1,:);
                retOOS6Mnocost=retnocost(TF2SMP+1-TF1+1:TF3SMP6M-TF1+1,:);
                %% Excess return calculation
                % IS
                rfind1=find(rftable{:,1}==char(startdtIS),1,'first');
                rfind2=find(rftable{:,1}==char(finishdtIS),1,'first');
                excessretISnocost=retISnocost-rftable{rfind1:rfind2,2};
                % OOS
                rfind31M=find(rftable{:,1}==char(finishdtOOS1M),1,'first');
                rfind33M=find(rftable{:,1}==char(finishdtOOS3M),1,'first');
                rfind36M=find(rftable{:,1}==char(finishdtOOS6M),1,'first');
                excessretIS=retIS-rftable{rfind1:rfind2,2};
                excessretOOS1M=retOOS1M-rftable{rfind2+1:rfind31M,2};
                excessretOOS3M=retOOS3M-rftable{rfind2+1:rfind33M,2};
                excessretOOS6M=retOOS6M-rftable{rfind2+1:rfind36M,2};
                excessretOOS1Mnocost=retOOS1Mnocost-rftable{rfind2+1:rfind31M,2};
                excessretOOS3Mnocost=retOOS3Mnocost-rftable{rfind2+1:rfind33M,2};
                excessretOOS6Mnocost=retOOS6Mnocost-rftable{rfind2+1:rfind36M,2};
                %% Generating the bootstraps
                indices=stationary_bootstrap((1:size(retIS,1))',Bsize,Bwindow);
                
                %% Generating performance index
                
                modelscount=size(excessretIS,2);
                Perf_B=zeros(Bsize,modelscount);
                Perf_Bnocost=zeros(Bsize,modelscount);
                disp(['Bootstrap generation initiated for IS ', ...
                    num2str(month),'/',num2str(year),'-',num2str(month),'/',num2str(yearOOS)]);

                
                mu=mean(excessretIS);
                sigma=std(excessretIS);
                munocost=mean(excessretISnocost);
                sigmanocost=std(excessretISnocost);
                
                
                Perf=mu./sigma;
                Perfnocost=munocost./sigmanocost;
                
%                 Perf(sigma == 0) = 0;
%                 Perfnocost(sigmanocost == 0) = 0;
%                 % avoid situation with sigma too small:
%                 idx_tmp = sigma > 0 & sigma < sigmathresh;
%                 Perf(idx_tmp) = mu(idx_tmp) / sigmathresh;
%                 
%                 idx_tmpnocost = sigmanocost > 0 & sigmanocost < sigmathresh;
%                 Perfnocost(idx_tmpnocost) = mu(idx_tmpnocost) / sigmathresh;
%                 
                
                parfor b=1:Bsize
                    bsdata=excessretIS(indices(:,b),:);
                    bsdatanocost=excessretISnocost(indices(:,b),:);
                    mu=mean(bsdata);
                    sigma=std(bsdata);
                    Perf_B(b,:)=mu./sigma;
%                     Perf_B(sigma==0,b)=0;
%                     idx_tmp = sigma > 0 & sigma < sigmathresh;
%                     Perf_B(b,idx_tmp) = mu(idx_tmp) / sigmathresh;
                    munocost=mean(bsdata);
                    sigmanocost=std(bsdata);
                    Perf_Bnocost(b,:)=munocost./sigmanocost;
%                     Perf_Bnocost(sigmanocost==0,b)=0;
%                     idx_tmpnocost = sigmanocost > 0 & sigmanocost < sigmathresh;
%                     Perf_Bnocost(b,idx_tmpnocost) = mu(idx_tmpnocost) / sigmathresh;
                    
                end
                
                pvalues = compute_pvalues(Perf', Perf_B', std(excessretIS));
                pvaluesnocost = compute_pvalues(Perfnocost', Perf_Bnocost', std(excessretISnocost));
                % Liang procedure
                [pi_0hat,lambda]=est_pi0_disc(pvalues, N_bins,Max_lambda);
                [pi_0hatnocost,lambdanocost]=est_pi0_disc(pvaluesnocost, N_bins,Max_lambda);
                % Barras et al 2010 estimation of FDR+/-
                [pi_aplushat, pi_aminushat] = compute_pi_ahat(pvalues, Perf', pi_0hat, gamma);
                [pi_aplushatnocost, pi_aminushatnocost] = compute_pi_ahat(pvaluesnocost, Perf', pi_0hat, gamma);
                [PORTFDR, FDRhat] = my_portfolio_FDR(fdrtarget, Perf', pvalues, pi_0hat);
                [PORTFDRnocost, FDRhatnocost] = my_portfolio_FDR(fdrtarget, Perfnocost', pvaluesnocost, pi_0hat);
                
                
                
                %% Break even transaction costs
                
                [max_perf,topmodel]=max(Perf);
                newcost=transbasis;
                totalrettop=sum(excessretIS(:,topmodel));
                while totalrettop>0 && newcost<=1e4
                    newcost=newcost+1;
                    [retser,trdcount]=TransCost(inputser(TF1SMP-TF1+1:TF2SMP-TF1+1,topmodel)...
                        ,dataret(TF1SMP-TF1+1:TF2SMP-TF1+1),newcost*1e-4);
                    excessretser=retser-rftable{rfind1:rfind2,2};
                    totalrettop=sum(excessretser);
                    newcost=newcost+(trdcount==0)*1e4;
                end
                
                %% IS and OOS return for outperformers
                % IS portfolio return
                finalretIS=sum(excessretIS);
                finalretISnocost=sum(excessretISnocost);
                PortRetIS=finalretIS*PORTFDR/max(sum(PORTFDR),1);
                PortRetISnocost=finalretISnocost*PORTFDRnocost/max(sum(PORTFDRnocost),1);
                % OOS portfolio return
                % 1M
                finalretOOS1M=sum(excessretOOS1M);
                finalretOOS1Mnocost=sum(excessretOOS1Mnocost);
                PortRetOOS1M=finalretOOS1M*PORTFDR/max(sum(PORTFDR),1);
                PortRetOOS1Mnocost=finalretOOS1Mnocost*PORTFDRnocost/max(sum(PORTFDRnocost),1);
                % 3M
                finalretOOS3M=sum(excessretOOS3M);
                finalretOOS3Mnocost=sum(excessretOOS3Mnocost);
                PortRetOOS3M=finalretOOS3M*PORTFDR/max(sum(PORTFDR),1);
                PortRetOOS3Mnocost=finalretOOS3Mnocost*PORTFDRnocost/max(sum(PORTFDRnocost),1);
                % 6M
                finalretOOS6M=sum(excessretOOS6M);
                finalretOOS6Mnocost=sum(excessretOOS6Mnocost);
                PortRetOOS6M=finalretOOS6M*PORTFDR/max(sum(PORTFDR),1);
                PortRetOOS6Mnocost=finalretOOS6Mnocost*PORTFDRnocost/max(sum(PORTFDRnocost),1);
                
                
                %% Recording the results
                %iter=IS+(month-1)*2+(yearOOS-2006)*24;
                iter=iter+1;
                resvar{iter,1}=PORTFDR;
                resvar{iter,2}=PORTFDRnocost;
                resvar{iter,3}=pvalues;
                resvar{iter,4}=pvaluesnocost;
                resvar{iter,5}=Perf;
                resvar{iter,6}=Perfnocost;
                resvar{iter,7}=pi_0hat;
                resvar{iter,8}=pi_0hatnocost;
                resvar{iter,9}=pi_aplushat;
                resvar{iter,10}=pi_aplushatnocost;
                resvar{iter,11}=pi_aminushat;
                resvar{iter,12}=pi_aminushatnocost;
                resvar{iter,13}=lambda;
                resvar{iter,14}=lambdanocost;
                resvar{iter,15}=yearOOS;
                resvar{iter,16}=month;
                resvar{iter,17}=IS;
                resvar{iter,1+17}=yearOOS;
                resvar{iter,2+17}=month;
                resvar{iter,3+17}=IS;
                resvar{iter,4+17}=sum(PORTFDR);
                resvar{iter,5+17}=sum(PORTFDRnocost);
                resvar{iter,6+17}=trdcount;
                resvar{iter,7+17}=newcost;
                resvar{iter,8+17}=pi_0hat;
                resvar{iter,9+17}=pi_0hatnocost;
                resvar{iter,10+17}=PortRetIS;
                resvar{iter,11+17}=PortRetISnocost;
                resvar{iter,12+17}=PortRetOOS1M;
                resvar{iter,13+17}=PortRetOOS1Mnocost;
                resvar{iter,14+17}=PortRetOOS3M;
                resvar{iter,15+17}=PortRetOOS3Mnocost;
                resvar{iter,16+17}=PortRetOOS6M;
                resvar{iter,17+17}=PortRetOOS6Mnocost;
                
                toc
                %xlswrite([lbl,'xlsx'],result_table);
            end
        end
        
    end
    result_table=cell2table(resvar(:,18:end));
    result_table.Properties.VariableNames={'Year' 'Month' 'IS_Period' 'FDRPortfolioSizeAT' ...
        'FDRPortSizeBT' 'TradeCount' 'BreakEvenCost' 'pi0hatAT' 'pi0hatBT'...
        'FDRRetISAT' 'FDRRetISBT' 'FDRRetOOS1MAT' 'FDRRetOOS1MBT'...
        'FDRRetOOS3MAT' 'FDRRetOOS3MBT' 'FDRRetOOS6MAT' 'FDRRetOOS6MBT'};

    writetable(result_table,[lbl,'_',num2str(yearser(1)),'_',num2str(yearser(end)),'.xlsx'])
    save([lbl,'_',num2str(yearser(1)),'_',num2str(yearser(end)),'.mat'],'resvar');
    
    toc
    clear resvar result_table;
end
delete(poolobj);