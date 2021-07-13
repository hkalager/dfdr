% Based on a code created 10 Oct 2017
% Created 15 Nov 2017
% Last Modified 13 July 2021 
clear;
clc;
labels={'NDDUUS'	'NDDUUK'	'NDDUJN'	'NDEUSRU'	'NDEUCHF'	'NDUEBRAF'...
    'MSEIESUN'	'M1MA'	'MIMUJORN' 'MXWO' 'MXEF' 'MXFEM'};
rftable=readtable('FEDL01.csv');
rftable{:,2}=rftable{:,2}/100/260;
transser=[25*ones(1,3),50*ones(1,6),25,50,50];
delete(gcp('nocreate'));
poolobj=parpool('local',4);
disp('Checking the datasets');
start=2002;
finish=2016;
oosperiod=1; % Out of sample holding period (1 Month)
sigmathresh=1e-3;
% Bootstrap specification
Bsize=1000;
Bwindow=10;
%% Liang specification
Max_lambda=0.99;
N_bins=20;
gamma = 0.4;
%% FDR setting
fdrtarget=0.2;
%% Year range
yearser=2006:2015;
tic
%% Main loop
for i=11%:numel(labels)
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
            for IS=2
                yr=yearOOS-IS;
                begindtstrIS=[num2str(month),'/',num2str(yr)];
                TF1SMP = contains(master.textdata(:,1),begindtstrIS);
                TF1SMP=find(TF1SMP==1,1,'first');
                startdtIS=master.textdata(TF1SMP,1);
                %% IS+1M OOS
                finishdtstrIOS1M=[num2str(rem(month+oosperiod,12+1e-10)),'/',num2str(yr+IS+(month+oosperiod>12))];
                TF3SMP1M = contains(master.textdata(:,1),finishdtstrIOS1M);
                TF3SMP1M=find(TF3SMP1M==1,1,'first')-1;
                finishdtIOS1M=master.textdata(TF3SMP1M,1);
                %% IS+3M OOS
                finishdtstrIOS3M=[num2str(rem(month+3*oosperiod,12+1e-10)),'/',num2str(yr+IS+(month+3*oosperiod>12))];
                TF3SMP3M = contains(master.textdata(:,1),finishdtstrIOS3M);
                TF3SMP3M=find(TF3SMP3M==1,1,'first')-1;
                finishdtIOS3M=master.textdata(TF3SMP3M,1);
                %% IS+6M OOS
                finishdtstrIOS6M=[num2str(rem(month+6*oosperiod,12+1e-10)),'/',num2str(yr+IS+(month+6*oosperiod>12))];
                TF3SMP6M = contains(master.textdata(:,1),finishdtstrIOS6M);
                TF3SMP6M=find(TF3SMP6M==1,1,'first')-1;
                finishdtIOS6M=master.textdata(TF3SMP6M,1);
                %% IOS return
                % 1M
                retIOS1M=ret(TF1SMP-TF1+1:TF3SMP1M-TF1+1,:);
                % 3M
                retIOS3M=ret(TF1SMP-TF1+1:TF3SMP3M-TF1+1,:);
                % 6M
                retIOS6M=ret(TF1SMP-TF1+1:TF3SMP6M-TF1+1,:);
                %% Excess return calculation
                % IS
                rfind1=find(rftable{:,1}==char(startdtIS),1,'first');
                % IOS
                rfind31M=find(rftable{:,1}==char(finishdtIOS1M),1,'first');
                rfind33M=find(rftable{:,1}==char(finishdtIOS3M),1,'first');
                rfind36M=find(rftable{:,1}==char(finishdtIOS6M),1,'first');
                excessretIOS1M=retIOS1M-rftable{rfind1:rfind31M,2};
                excessretIOS3M=retIOS3M-rftable{rfind1:rfind33M,2};
                excessretIOS6M=retIOS6M-rftable{rfind1:rfind36M,2};
                
                %% Generating the bootstraps
                indicesIOS1M=stationary_bootstrap((1:size(retIOS1M,1))',Bsize,Bwindow);
                indicesIOS3M=stationary_bootstrap((1:size(retIOS3M,1))',Bsize,Bwindow);
                indicesIOS6M=stationary_bootstrap((1:size(retIOS6M,1))',Bsize,Bwindow);
                %% Generating performance index
                
                modelscount=size(excessretIOS1M,2);
                Perf_BIOS1M=zeros(Bsize,modelscount);
                Perf_BIOS3M=zeros(Bsize,modelscount);
                Perf_BIOS6M=zeros(Bsize,modelscount);
                disp(['Bootstrap generation initiated for IOS ', ...
                    num2str(month),'/',num2str(yr),'-',num2str(month),'/',num2str(yearOOS)]);

                %% Performance measure 
                % IOS 1M
                muIOS1M=mean(excessretIOS1M);
                sigmaIOS1M=std(excessretIOS1M);
                PerfIOS1M=muIOS1M./sigmaIOS1M;
                % IOS 3M
                muIOS3M=mean(excessretIOS3M);
                sigmaIOS3M=std(excessretIOS3M);
                PerfIOS3M=muIOS3M./sigmaIOS3M;
                % IOS 6M
                muIOS6M=mean(excessretIOS6M);
                sigmaIOS6M=std(excessretIOS6M);
                PerfIOS6M=muIOS6M./sigmaIOS6M;
                
                %% Bootstrapping 
                parfor b=1:Bsize
                    bsdataIOS1M=excessretIOS1M(indicesIOS1M(:,b),:);
                    muIOS1M=mean(bsdataIOS1M);
                    sigmaIOS1M=std(bsdataIOS1M);
                    Perf_BIOS1M(b,:)=muIOS1M./sigmaIOS1M;
                    
                    bsdataIOS3M=excessretIOS3M(indicesIOS3M(:,b),:);
                    muIOS3M=mean(bsdataIOS3M);
                    sigmaIOS3M=std(bsdataIOS3M);
                    Perf_BIOS3M(b,:)=muIOS3M./sigmaIOS3M;
                    
                    bsdataIOS6M=excessretIOS6M(indicesIOS6M(:,b),:);
                    muIOS6M=mean(bsdataIOS6M);
                    sigmaIOS6M=std(bsdataIOS6M);
                    Perf_BIOS6M(b,:)=muIOS6M./sigmaIOS6M;

                end
                % P value calculation
                pvaluesIOS1M = compute_pvalues(PerfIOS1M', Perf_BIOS1M', std(excessretIOS1M));
                pvaluesIOS3M = compute_pvalues(PerfIOS3M', Perf_BIOS3M', std(excessretIOS3M));
                pvaluesIOS6M = compute_pvalues(PerfIOS6M', Perf_BIOS6M', std(excessretIOS6M));
                % Liang procedure
                [pi_0hatIOS1M,lambdaIOS1M]=est_pi0_disc(pvaluesIOS1M, N_bins,Max_lambda);
                [pi_0hatIOS3M,lambdaIOS3M]=est_pi0_disc(pvaluesIOS3M, N_bins,Max_lambda);
                [pi_0hatIOS6M,lambdaIOS6M]=est_pi0_disc(pvaluesIOS6M, N_bins,Max_lambda);
                % Barras et al 2010 estimation of FDR+/-
                [pi_aplushatIOS1M, pi_aminushatIOS1M] = compute_pi_ahat(pvaluesIOS1M, PerfIOS1M', pi_0hatIOS1M, gamma);
                [PORTFDRIOS1M, FDRhatIOS1M] = my_portfolio_FDR(fdrtarget, PerfIOS1M', pvaluesIOS1M, pi_0hatIOS1M);
                
                [pi_aplushatIOS3M, pi_aminushatIOS3M] = compute_pi_ahat(pvaluesIOS3M, PerfIOS3M', pi_0hatIOS3M, gamma);
                [PORTFDRIOS3M, FDRhatIOS3M] = my_portfolio_FDR(fdrtarget, PerfIOS3M', pvaluesIOS3M, pi_0hatIOS3M);
                
                [pi_aplushatIOS6M, pi_aminushatIOS6M] = compute_pi_ahat(pvaluesIOS6M, PerfIOS6M', pi_0hatIOS6M, gamma);
                [PORTFDRIOS6M, FDRhatIOS6M] = my_portfolio_FDR(fdrtarget, PerfIOS6M', pvaluesIOS6M, pi_0hatIOS6M);
                                
                %% Number of trades for the best strategy
                
                [~,topmodelIOS1M]=max(PerfIOS1M);
                [~,trdcountIOS1M]=TransCost(inputser(TF1SMP-TF1+1:TF3SMP1M-TF1+1,topmodelIOS1M)...
                        ,dataret(TF1SMP-TF1+1:TF3SMP1M-TF1+1),transbasis*1e-4);
                
                [~,topmodelIOS3M]=max(PerfIOS3M);
                [~,trdcountIOS3M]=TransCost(inputser(TF1SMP-TF1+1:TF3SMP3M-TF1+1,topmodelIOS3M)...
                        ,dataret(TF1SMP-TF1+1:TF3SMP3M-TF1+1),transbasis*1e-4);

                [~,topmodelIOS6M]=max(PerfIOS6M);
                [~,trdcountIOS6M]=TransCost(inputser(TF1SMP-TF1+1:TF3SMP6M-TF1+1,topmodelIOS6M)...
                        ,dataret(TF1SMP-TF1+1:TF3SMP6M-TF1+1),transbasis*1e-4);

                %% IS and OOS return for outperformers
                % IS portfolio return
                finalretIOS1M=sum(excessretIOS1M);
                PortRetIOS1M=finalretIOS1M*PORTFDRIOS1M/max(sum(PORTFDRIOS1M),1);
                
                finalretIOS3M=sum(excessretIOS3M);
                PortRetIOS3M=finalretIOS3M*PORTFDRIOS3M/max(sum(PORTFDRIOS3M),1);
                
                finalretIOS6M=sum(excessretIOS6M);
                PortRetIOS6M=finalretIOS6M*PORTFDRIOS6M/max(sum(PORTFDRIOS6M),1);
                
                %% Recording the results
                %iter=IS+(month-1)*2+(yearOOS-2006)*24;
                iter=iter+1;
                dim=1;
                resvar{iter,dim}=PORTFDRIOS1M;
                dim=dim+1;
                resvar{iter,dim}=PORTFDRIOS3M;
                dim=dim+1;
                resvar{iter,dim}=PORTFDRIOS6M;
                dim=dim+1;
                resvar{iter,dim}=pvaluesIOS1M;
                dim=dim+1;
                resvar{iter,dim}=pvaluesIOS3M;
                dim=dim+1;
                resvar{iter,dim}=pvaluesIOS6M;
                dim=dim+1;
                resvar{iter,dim}=PerfIOS1M;
                dim=dim+1;
                resvar{iter,dim}=PerfIOS3M;
                dim=dim+1;
                resvar{iter,dim}=PerfIOS6M;
                dim=dim+1;
                resvar{iter,dim}=pi_0hatIOS1M;
                dim=dim+1;
                resvar{iter,dim}=pi_0hatIOS3M;
                dim=dim+1;
                resvar{iter,dim}=pi_0hatIOS6M;
                dim=dim+1;
                resvar{iter,dim}=pi_aplushatIOS1M;
                dim=dim+1;
                resvar{iter,dim}=pi_aplushatIOS3M;
                dim=dim+1;
                resvar{iter,dim}=pi_aplushatIOS6M;
                dim=dim+1;
                resvar{iter,dim}=pi_aminushatIOS1M;
                dim=dim+1;
                resvar{iter,dim}=pi_aminushatIOS3M;
                dim=dim+1;
                resvar{iter,dim}=pi_aminushatIOS6M;
                dim=dim+1;
                resvar{iter,dim}=lambdaIOS1M;
                dim=dim+1;
                resvar{iter,dim}=lambdaIOS3M;
                dim=dim+1;
                resvar{iter,dim}=lambdaIOS6M;
                dim=dim+1;
                resvar{iter,dim}=yearOOS;
                dim=dim+1;
                resvar{iter,dim}=month;
                %% Outcome table
                resvar{iter,1+dim}=yearOOS;
                resvar{iter,2+dim}=month;
                resvar{iter,3+dim}=sum(PORTFDRIOS1M);
                resvar{iter,4+dim}=sum(PORTFDRIOS3M);
                resvar{iter,5+dim}=sum(PORTFDRIOS6M);
                resvar{iter,6+dim}=trdcountIOS1M;
                resvar{iter,7+dim}=trdcountIOS3M;
                resvar{iter,8+dim}=trdcountIOS6M;
                resvar{iter,9+dim}=pi_0hatIOS1M;
                resvar{iter,10+dim}=pi_0hatIOS3M;
                resvar{iter,11+dim}=pi_0hatIOS6M;
                resvar{iter,12+dim}=PortRetIOS1M;
                resvar{iter,13+dim}=PortRetIOS3M;
                resvar{iter,14+dim}=PortRetIOS6M;

                toc
                %xlswrite([lbl,'xlsx'],result_table);
            end
        end
        
    end
    result_table=cell2table(resvar(:,end-13:end));
    result_table.Properties.VariableNames={'Year' 'Month' 'FDRPortSizeIS_1M' ...
        'FDRPortSizeIS_3M' 'FDRPortSizeIS_6M' 'TradeCountIS_1M' 'TradeCountIS_3M' ...
        'TradeCountIS_6M' 'pi0hatIS_1M' 'pi0hatIS_3M' 'pi0hatIS_6M'...
        'FDRRetIS_1M' 'FDRRetIS_3M' 'FDRRetIS_6M'};
    
    writetable(result_table,['CV_',lbl,'_',num2str(yearser(1)),'_',num2str(yearser(end)),'.xlsx'])
    save(['CV_',lbl,'_',num2str(yearser(1)),'_',num2str(yearser(end)),'.mat'],'resvar');
    
    toc
    clear resvar result_table;
end
delete(poolobj);