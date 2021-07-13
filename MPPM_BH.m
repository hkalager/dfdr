% Script by Arman Hassanniakalager
% Last Modified 04 Jan 2019 14:31
clear;
clc;
labels={'MXWO' 'NDDUUS'	'NDDUUK'	'NDDUJN'	'MXEF' 'NDEUSRU'	'NDEUCHF'	'NDUEBRAF'...
    'MXFEM' 'MSEIESUN'	'M1MA'	'MIMUJORN'};
regions={'Advanced' 'US' 'UK' 'Japan' 'Emerging' 'Russia' 'China' 'Brazil' ...
    'Frontier' 'Estonia' 'Morocco' 'Jordan'};
rftable=readtable('FEDL01.csv');
rftable{:,2}=rftable{:,2}/100/260;
transser=[25*ones(1,4),50*ones(1,8)];
%% MPPM definition
rho=2;
d_t=260;
MPPM=@(r_t,r_f) (1/((1-rho)*(1/d_t))*log(mean(((1+r_t)./(1+r_f)).^(1-rho))));

%% Year range
yearser=2006:2015;
startyr=2002;
finishyr=2016;
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
    for yearOOS=yearser(1):yearser(end)%2006:2016
        iter=iter+1;
        PortRetIS=[];
        PortRetOOS1M=[];
        PortRetOOS3M=[];
        PortRetOOS6M=[];
        ISportFDRsize=[];
        ISportFDRser=[];
        OOS1MportFDRser=[];
        OOS3MportFDRser=[];
        OOS6MportFDRser=[];
        rf_IS=[];
        rf_1M=[];
        rf_3M=[];
        rf_6M=[];
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
            retIS=ret(TF1SMP-TF1+1:TF2SMP-TF1+1,:);
            
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
            excessretIS=retIS-rftable{rfind1:rfind2,2};
            rfind31M=find(rftable{:,1}==char(finishdtIOS1M),1,'first');
            rfind33M=find(rftable{:,1}==char(finishdtIOS3M),1,'first');
            rfind36M=find(rftable{:,1}==char(finishdtIOS6M),1,'first');
            rf_IS(end+1:end+size(retIS,1),1)=rftable{rfind1:rfind2,2};
            rf_1M(end+1:end+size(retOOS1M,1),1)=rftable{rfind2+1:rfind31M,2};
            rf_3M(end+1:end+size(retOOS3M,1),1)=rftable{rfind2+1:rfind33M,2};
            rf_6M(end+1:end+size(retOOS6M,1),1)=rftable{rfind2+1:rfind36M,2};
            
            %% Loading the FDR portfolios
            ISportFDR=resvarIS{(iterplus-1)*2+IS,1};
            
            %% Portfolio generation
            % IS portfolio return
            finalretIS=sum(excessretIS);
            PortRetIS(end+1)=finalretIS*ISportFDR/max(sum(ISportFDR),1);
            ISportFDRser(end+1:end+size(excessretIS,1))=mean(excessretIS(:,ISportFDR==1),2);
            % OOS portfolio return
            % 1M
            finalretOOS1M=sum(retOOS1M);
            PortRetOOS1M(end+1)=finalretOOS1M*ISportFDR/max(sum(ISportFDR),1);
            OOS1MportFDRser(end+1:end+size(retOOS1M,1))=mean(retOOS1M(:,ISportFDR==1),2);
            % 3M
            finalretOOS3M=sum(retOOS3M);
            PortRetOOS3M(end+1)=finalretOOS3M*ISportFDR/max(sum(ISportFDR),1);
            OOS3MportFDRser(end+1:end+size(retOOS3M,1))=mean(retOOS3M(:,ISportFDR==1),2);
            % 6M
            finalretOOS6M=sum(retOOS6M);
            PortRetOOS6M(end+1)=finalretOOS6M*ISportFDR/max(sum(ISportFDR),1);
            OOS6MportFDRser(end+1:end+size(retOOS6M,1))=mean(retOOS6M(:,ISportFDR==1),2);
            
            
        end
        %% B&H Ret & Sharpe
        T_S = find(contains(master.textdata(:,1),num2str(yearOOS)),1,'first')-1;
        T_L = find(contains(master.textdata(:,1),num2str(yearOOS)),1,'last')-1;
        t_ser=price2ret(master.data(T_S-1:T_L,end));
        t_excess=t_ser-rf_1M;
        B_H_ret=sum(t_excess);
        B_H_Sharpe=sharpe(t_excess,0)*260^.5;
        
        %% Manipulation ratio
        
        IS_MPPM=MPPM(ISportFDRser',rf_IS);
        OOS1M_SMPPM=MPPM(OOS1MportFDRser',rf_1M);
        OOS3M_SMPPM=MPPM(OOS3MportFDRser',rf_3M);
        OOS6M_SMPPM=MPPM(OOS6MportFDRser',rf_6M);
        
        %% Recording the results
        % Period
        dim=1;
        resvar{iter,dim}=yearOOS;
        % B&H RET
        dim=dim+1;
        resvar{iter,dim}=B_H_ret;
        
        % B&H Sharpe
        dim=dim+1;
        resvar{iter,dim}=B_H_Sharpe;
        
        % MPPM
        dim=dim+1;
        resvar{iter,dim}=IS_MPPM;
        
        dim=dim+1;
        resvar{iter,dim}=OOS1M_SMPPM;
        dim=dim+1;
        resvar{iter,dim}=OOS3M_SMPPM;
        dim=dim+1;
        resvar{iter,dim}=OOS6M_SMPPM;
        
        dim=dim+1;
        resvar{iter,dim}=regions{i};
        disp(['Processing done for ',regions{i},' period ',num2str(yearOOS), ' completed ...']);
        toc;
    end
    
end
result_table=cell2table(resvar);
result_table.Properties.VariableNames={'Year' 'BH_RET' ...
    'BH_SHARPE' 'IS_MPPM' 'OOS1M_MPPM' 'OOS3M_MPPM' 'OOS6M_MPPM' ...
    'Country'};
%drpadd='/Empirics3/';
drpadd='C:\Users\ah2493\Dropbox\Shared_Folder\Chapter3\';
writetable(result_table,[drpadd,'AnnualTable_',num2str(yearser(1)),'_',num2str(yearser(end)),'_IS_',num2str(IS),'_MPPM.xlsx'])
toc