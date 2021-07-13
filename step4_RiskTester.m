% Script by Arman Hassanniakalager
% Last Modified 25 May 2018 16:00 BST
clear;
clc;
labels={'NDDUUS' 'MXWO' 'MXEF'};
regions={ 'US' 'Advanced' 'Emerging'};
rftable=readtable('FEDL01.csv');
rsk_table=readtable('fsi1.csv');
t=datetime(rsk_table.Date,'Format','dd/MM/uuuu');
rsk_table.Date=t;
rftable{:,2}=rftable{:,2}/100/260;
transser=[25,25,50];
%% Year range
yearser=2006:2015;
startyr=2002;
finishyr=2017;
IS=2; % In year
oosperiod=1;
tic
iter=0;
%% Main loop
for i=1:numel(labels)
    lbl=labels{i};
    transbasis=transser(i);
    fllbl=[lbl,'_',num2str(startyr),'_',num2str(finishyr),'-insample-',num2str(transbasis),'.mat'];
    load(fllbl,'master','ret','TF1');
    S1=load([lbl,'_',num2str(yearser(1)),'_',num2str(yearser(end)),'.mat'],'resvar');
    resvarIS=S1.resvar;
    disp(['Dataset loaded for ',lbl]);
    iterplus=0;
    tic;
    PortRetOOS1M_Y=[];
    PortRetOOS3M_Y=[];
    PortRetOOS6M_Y=[];
    PortRetOOS1MHighRisk_Y=[];
    PortRetOOS3MHighRisk_Y=[];
    PortRetOOS6MHighRisk_Y=[];
    PortRetOOS1MLowRisk_Y=[];
    PortRetOOS3MLowRisk_Y=[];
    PortRetOOS6MLowRisk_Y=[];
    
    for yearOOS=yearser(1):yearser(end)%2006:2016
        iter=iter+1;
        PortRetIS=[];
        PortRetOOS1M=[];
        PortRetOOS3M=[];
        PortRetOOS6M=[];
        PortRetOOS1MHighRisk=[];
        PortRetOOS3MHighRisk=[];
        PortRetOOS6MHighRisk=[];
        PortRetOOS1MLowRisk=[];
        PortRetOOS3MLowRisk=[];
        PortRetOOS6MLowRisk=[];
        ISportFDRsize=[];
        ISportFDRser=[];
        OOS1MportFDRser=[];
        OOS1MportFDRserHighRisk=[];
        OOS1MportFDRserLowRisk=[];
        OOS3MportFDRser=[];
        OOS3MportFDRserHighRisk=[];
        OOS3MportFDRserLowRisk=[];
        OOS6MportFDRser=[];
        OOS6MportFDRserHighRisk=[];
        OOS6MportFDRserLowRisk=[];
        
        risk1MHsum=0;
        risk1MLsum=0;
        risk3MHsum=0;
        risk3MLsum=0;
        risk6MHsum=0;
        risk6MLsum=0;
        for month=1:12
            iterplus=iterplus+1;
            yr=yearOOS-IS;
            begindtstrIS=[num2str(month),'/',num2str(yr)];
            TF1SMP = contains(master.textdata(:,1),begindtstrIS);
            TF1SMP=find(TF1SMP==1,1,'first');
            startdtIS=master.textdata(TF1SMP,1);
            
            %% IS
            finishdtstrIS=[num2str(month),'/',num2str(yr+IS)];
            TF2SMP = contains(master.textdata(:,1),finishdtstrIS);
            TF2SMP=find(TF2SMP==1,1,'first')-1;
            finishdtIS=master.textdata(TF2SMP,1);
            
            %% IS+1M OOS
            finishdtstrIOS1M=[num2str(rem(month+oosperiod,12+1e-10)),'/',num2str(yr+IS+(month+oosperiod>12))];
            TF3SMP1M = contains(master.textdata(:,1),finishdtstrIOS1M);
            TF3SMP1M=find(TF3SMP1M==1,1,'first')-1;
            finishdtOOS1M=master.textdata(TF3SMP1M,1);
            %% IS+3M OOS
            finishdtstrIOS3M=[num2str(rem(month+3*oosperiod,12+1e-10)),'/',num2str(yr+IS+(month+3*oosperiod>12))];
            TF3SMP3M = contains(master.textdata(:,1),finishdtstrIOS3M);
            TF3SMP3M=find(TF3SMP3M==1,1,'first')-1;
            finishdtOOS3M=master.textdata(TF3SMP3M,1);
            %% IS+6M OOS
            finishdtstrIOS6M=[num2str(rem(month+6*oosperiod,12+1e-10)),'/',num2str(yr+IS+(month+6*oosperiod>12))];
            TF3SMP6M = contains(master.textdata(:,1),finishdtstrIOS6M);
            TF3SMP6M=find(TF3SMP6M==1,1,'first')-1;
            finishdtOOS6M=master.textdata(TF3SMP6M,1);
            
            %% OOS return
            % 1M
            retOOS1M=ret(TF2SMP+1-TF1+1:TF3SMP1M-TF1+1,:);
            % 3M
            retOOS3M=ret(TF2SMP+1-TF1+1:TF3SMP3M-TF1+1,:);
            % 6M
            retOOS6M=ret(TF2SMP+1-TF1+1:TF3SMP6M-TF1+1,:);
            
            %% Excess return
            rfind1=find(rftable{:,1}==char(startdtIS),1,'first');
            rfind2=find(rftable{:,1}==char(finishdtIS),1,'first');
            
            rkind1=find(rsk_table{:,1}==startdtIS,1,'first');
            rkind2=find(rsk_table{:,1}==finishdtIS,1,'first');
            
            rfind31M=find(rftable{:,1}==char(finishdtOOS1M),1,'first');
            rfind33M=find(rftable{:,1}==char(finishdtOOS3M),1,'first');
            rfind36M=find(rftable{:,1}==char(finishdtOOS6M),1,'first');
            
            dtoos1=char(finishdtOOS1M);
            dtoos3=char(finishdtOOS3M);
            dtoos6=char(finishdtOOS6M);
            
            rkfind31M=find(rsk_table{:,1}==dtoos1,1,'last');
            rkfind33M=find(rsk_table{:,1}==dtoos3,1,'last');
            rkfind36M=find(rsk_table{:,1}==dtoos6,1,'last');
            
            excessretOOS1M=retOOS1M-rftable{rfind2+1:rfind31M,2};
            excessretOOS3M=retOOS3M-rftable{rfind2+1:rfind33M,2};
            excessretOOS6M=retOOS6M-rftable{rfind2+1:rfind36M,2};
            
            %% High risk low risk excess return
            riskser1M=sign(rsk_table{rkind2+1:rkfind31M,7+i});
            riskser3M=sign(rsk_table{rkind2+1:rkfind33M,7+i});
            riskser6M=sign(rsk_table{rkind2+1:rkfind36M,7+i});
            
            %% Loading the FDR portfolios
            ISportFDR=resvarIS{(iterplus-1)*2+IS,1};
            
            %% Portfolio generation
            % OOS portfolio return
            % 1M
            finalretOOS1M=sum(excessretOOS1M);
            PortRetOOS1M(end+1)=finalretOOS1M*ISportFDR/max(sum(ISportFDR),1);
            OOS1MportFDRser(end+1:end+size(excessretOOS1M,1))=mean(excessretOOS1M(:,ISportFDR==1),2);
            add_count1MH=(sum(riskser1M>0)/numel(riskser1M)*22);
            risk1MHsum=risk1MHsum+add_count1MH;
            if sum(riskser1M>0)
                OOS1MH=sum(excessretOOS1M(riskser1M>0,:),1)*add_count1MH;
                OOS1MportFDRserHighRisk(end+1:end+size(excessretOOS1M(riskser1M>0,:),1))=mean(excessretOOS1M(riskser1M>0,ISportFDR==1),2);
            else
                OOS1MH=nan(1,size(excessretOOS1M,2));
         
            end
            PortRetOOS1MHighRisk(end+1)=OOS1MH*ISportFDR/max(sum(ISportFDR),1);
            
            add_count1ML=(sum(riskser1M<=0)/numel(riskser1M)*22);
            risk1MLsum=risk1MLsum+add_count1ML;
            if sum(riskser1M<=0)
                OOS1ML=sum(excessretOOS1M(riskser1M<=0,:),1)*add_count1ML;
                OOS1MportFDRserLowRisk(end+1:end+size(excessretOOS1M(riskser1M<=0,:),1))=mean(excessretOOS1M(riskser1M<=0,ISportFDR==1),2);
            else
                OOS1ML=nan(1,size(excessretOOS1M,2));
            end
            
            PortRetOOS1MLowRisk(end+1)=OOS1ML*ISportFDR/max(sum(ISportFDR),1);

            % 3M
            finalretOOS3M=sum(excessretOOS3M);
            PortRetOOS3M(end+1)=finalretOOS3M*ISportFDR/max(sum(ISportFDR),1);
            OOS3MportFDRser(end+1:end+size(excessretOOS3M,1))=mean(excessretOOS3M(:,ISportFDR==1),2);
            
            add_count3MH=(sum(riskser3M>0)/numel(riskser3M)*65);
            risk3MHsum=risk3MHsum+add_count3MH;
            
            if sum(riskser3M>0)
                OOS3MH=sum(excessretOOS3M(riskser3M>0,:),1)*add_count3MH;
                OOS3MportFDRserHighRisk(end+1:end+size(excessretOOS3M(riskser3M>0,:),1))=mean(excessretOOS3M(riskser3M>0,ISportFDR==1),2);
            else
                OOS3MH=nan(1,size(excessretOOS3M,2));
            end
            PortRetOOS3MHighRisk(end+1)=OOS3MH*ISportFDR/max(sum(ISportFDR),1);
            
            add_count3ML=(sum(riskser3M<=0)/numel(riskser3M)*65);
            risk3MLsum=risk3MLsum+add_count3ML;
            
            if sum(riskser3M<=0)
                OOS3ML=sum(excessretOOS3M(riskser3M<=0,:),1)*add_count3ML;
                OOS3MportFDRserLowRisk(end+1:end+size(excessretOOS3M(riskser3M<=0,:),1))=mean(excessretOOS3M(riskser3M<=0,ISportFDR==1),2);
            else
                OOS3ML=nan(1,size(excessretOOS3M,2));
            end
            
            PortRetOOS3MLowRisk(end+1)=OOS3ML*ISportFDR/max(sum(ISportFDR),1);

            
            % 6M
            finalretOOS6M=sum(excessretOOS6M);
            PortRetOOS6M(end+1)=finalretOOS6M*ISportFDR/max(sum(ISportFDR),1);
            OOS6MportFDRser(end+1:end+size(excessretOOS6M,1))=mean(excessretOOS6M(:,ISportFDR==1),2);
            
            add_count6MH=(sum(riskser6M>0)/numel(riskser6M)*130);
            risk6MHsum=risk6MHsum+add_count6MH;            
            if sum(riskser6M>0)
                OOS6MH=sum(excessretOOS6M(riskser6M>0,:),1)*add_count6MH;
                OOS6MportFDRserHighRisk(end+1:end+size(excessretOOS6M(riskser6M>0,:),1))=mean(excessretOOS6M(riskser6M>0,ISportFDR==1),2);
            else
                OOS6MH=nan(1,size(excessretOOS6M,2));
            end
            PortRetOOS6MHighRisk(end+1)=OOS6MH*ISportFDR/max(sum(ISportFDR),1);
            
            add_count6ML=(sum(riskser6M<=0)/numel(riskser6M)*130);
            risk6MLsum=risk6MLsum+add_count6ML;
            
            if sum(riskser6M<=0)
                OOS6ML=sum(excessretOOS6M(riskser6M<=0,:),1)*add_count6ML;
                OOS6MportFDRserLowRisk(end+1:end+size(excessretOOS6M(riskser6M<=0,:),1))=mean(excessretOOS6M(riskser6M<=0,ISportFDR==1),2);
            else
                OOS6ML=nan(1,size(excessretOOS6M,2));
            end
            
            PortRetOOS6MLowRisk(end+1)=OOS6ML*ISportFDR/max(sum(ISportFDR),1);

            
        end
        %% Recording annual results
        PortRetOOS1M_Y((iter-1)*12+1:iter*12)=PortRetOOS1M;
        PortRetOOS1MHighRisk_Y((iter-1)*12+1:iter*12)=PortRetOOS1MHighRisk/22;
        PortRetOOS1MLowRisk_Y((iter-1)*12+1:iter*12)=PortRetOOS1MLowRisk/22;
        
        PortRetOOS3M_Y((iter-1)*12+1:iter*12)=PortRetOOS3M;
        PortRetOOS3MHighRisk_Y((iter-1)*12+1:iter*12)=PortRetOOS3MHighRisk/65;
        PortRetOOS3MLowRisk_Y((iter-1)*12+1:iter*12)=PortRetOOS3MLowRisk/65;
        
        PortRetOOS6M_Y((iter-1)*12+1:iter*12)=PortRetOOS6M;
        PortRetOOS6MHighRisk_Y((iter-1)*12+1:iter*12)=PortRetOOS6MHighRisk/130;
        PortRetOOS6MLowRisk_Y((iter-1)*12+1:iter*12)=PortRetOOS6MLowRisk/130;
        
        %% Sharpe ratio
        OOS1MSR=sharpe(OOS1MportFDRser,0)*260^.5;
        OOS1MSRH=(mean(OOS1MportFDRserHighRisk,'omitnan')/std(OOS1MportFDRserHighRisk,1,'omitnan'))*260^.5;
        OOS1MSRL=(mean(OOS1MportFDRserLowRisk,'omitnan')/std(OOS1MportFDRserLowRisk,1,'omitnan'))*260^.5;
        OOS3MSR=sharpe(OOS3MportFDRser,0)*260^.5;
        OOS3MSRH=(mean(OOS3MportFDRserHighRisk,'omitnan')/std(OOS3MportFDRserHighRisk,1,'omitnan'))*260^.5;
        OOS3MSRL=(mean(OOS3MportFDRserLowRisk,'omitnan')/std(OOS3MportFDRserLowRisk,1,'omitnan'))*260^.5;
        OOS6MSR=sharpe(OOS6MportFDRser,0)*260^.5;
        OOS6MSRH=(mean(OOS6MportFDRserHighRisk,'omitnan')/std(OOS6MportFDRserHighRisk,1,'omitnan'))*260^.5;
        OOS6MSRL=(mean(OOS6MportFDRserLowRisk,'omitnan')/std(OOS6MportFDRserLowRisk,1,'omitnan'))*260^.5;
        
        %% Recording the results
        % Period
        dim=1;
        resvar{iter,dim}=regions{i};
        dim=dim+1;
        resvar{iter,dim}=yr+IS;
        % OOS return (survivors)
        dim=dim+1;
        resvar{iter,dim}=(1+mean(PortRetOOS1M))^12-1;
        dim=dim+1;
        resvar{iter,dim}=(1+sum(PortRetOOS1MHighRisk,'omitnan')/risk1MHsum)^12-1;
        dim=dim+1;
        resvar{iter,dim}=(1+sum(PortRetOOS1MLowRisk,'omitnan')/risk1MLsum)^12-1;
        dim=dim+1;
        resvar{iter,dim}=(1+mean(PortRetOOS3M))^4-1;
        dim=dim+1;
        resvar{iter,dim}=(1+sum(PortRetOOS3MHighRisk,'omitnan')/risk3MHsum)^4-1;
        dim=dim+1;
        resvar{iter,dim}=(1+sum(PortRetOOS3MLowRisk,'omitnan')/risk3MLsum)^4-1;
        dim=dim+1;
        resvar{iter,dim}=(1+mean(PortRetOOS6M))^2-1;
        dim=dim+1;
        resvar{iter,dim}=(1+sum(PortRetOOS6MHighRisk,'omitnan')/risk6MHsum)^2-1;
        dim=dim+1;
        resvar{iter,dim}=(1+sum(PortRetOOS6MLowRisk,'omitnan')/risk6MLsum)^2-1;
        dim=dim+1;
        resvar{iter,dim}=OOS1MSR;
        dim=dim+1;
        resvar{iter,dim}=OOS1MSRH;
        dim=dim+1;
        resvar{iter,dim}=OOS1MSRL;
        dim=dim+1;
        resvar{iter,dim}=OOS3MSR;
        dim=dim+1;
        resvar{iter,dim}=OOS3MSRH;
        dim=dim+1;
        resvar{iter,dim}=OOS3MSRL;
        dim=dim+1;
        resvar{iter,dim}=OOS6MSR;
        dim=dim+1;
        resvar{iter,dim}=OOS6MSRH;
        dim=dim+1;
        resvar{iter,dim}=OOS6MSRL;
        disp(['Processing done for ',regions{i},' period ',num2str(yearOOS), ' completed ...']);
        toc;
    end
    [~,p1H(i)]=ttest(PortRetOOS1M_Y,PortRetOOS1MHighRisk_Y);
    [~,p1L(i)]=ttest(PortRetOOS1M_Y,PortRetOOS1MLowRisk_Y);
    [~,p3H(i)]=ttest(PortRetOOS3M_Y,PortRetOOS3MHighRisk_Y);
    [~,p3L(i)]=ttest(PortRetOOS3M_Y,PortRetOOS3MLowRisk_Y);
    [~,p6H(i)]=ttest(PortRetOOS6M_Y,PortRetOOS6MHighRisk_Y);
    [~,p6L(i)]=ttest(PortRetOOS6M_Y,PortRetOOS6MLowRisk_Y);
end
Rejections=table(p1H',p1L',p3H',p3L',p6H',p6L');
Rejections.Properties.VariableNames={'pval1M_High' 'pval1M_Low' ...
    'pval3M_High' 'pval3M_Low' 'pval6M_High' 'pval6M_Low'};
result_table=cell2table(resvar);
result_table.Properties.VariableNames={'Country' 'Year'  ...
    'OOS1_R' 'OOS1H_R' 'OOS1L_R' 'OOS3_R' 'OOS3H_R' 'OOS3L_R' ...
    'OOS6_R' 'OOS6H_R' 'OOS6L_R' 'OOS1_SR' 'OOS1H_SR' 'OOS1L_SR' 'OOS3_SR' 'OOS3H_SR' 'OOS3L_SR' ...
    'OOS6_SR' 'OOS6H_SR' 'OOS6L_SR'};
writetable(result_table,['RiskAnalysis_',num2str(yearser(1)),'_',num2str(yearser(end)),'_IS_',num2str(IS),'R9.xlsx'])
toc