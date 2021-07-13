% Script by Arman Hassanniakalager
% Last Modified 22 May 2018 12:12 BST
clear;
clc;
labels={'MXWO' 'NDDUUS'	'NDDUUK'	'NDDUJN'	'MXEF' 'NDEUSRU'	'NDEUCHF'	'NDUEBRAF'...
    'MXFEM' 'MSEIESUN'	'M1MA'	'MIMUJORN'};
regions={'Advanced' 'US' 'UK' 'Japan' 'Emerging' 'Russia' 'China' 'Brazil' ...
    'Frontier' 'Estonia' 'Morocco' 'Jordan'};
rftable=readtable('FEDL01.csv');
rftable{:,2}=rftable{:,2}/100/260;
transser=[25*ones(1,4),50*ones(1,8)];
%% Year range
yearser=2006:2015;
startyr=2002;
finishyr=2017;
oosperiod=1; 
tic
iter=0;
%% Main loop
for i=1:numel(labels)
    lbl=labels{i};
    transbasis=transser(i);
    fllbl=[lbl,'_',num2str(startyr),'_',num2str(finishyr),'-insample-',num2str(transbasis),'.mat'];
    load(fllbl,'master','ret','TF1','inputser','dataret');
    S1=load([lbl,'_',num2str(yearser(1)),'_',num2str(yearser(end)),'.mat'],'resvar');
    resvarIS=S1.resvar;
    S2=load(['CV_',lbl,'_',num2str(yearser(1)),'_',num2str(yearser(end)),'.mat'],'resvar');
    resvarIOS=S2.resvar;
    % Calculating the returns without transaction costs
    retnocost=TransCost(inputser,dataret,0);
    disp(['Dataset loaded for ',lbl]);
    iterplus=0;
    tic;
    for yearOOS=yearser(1):yearser(end)%2006:2016
        iter=iter+1;
        PortRetIS=[];
        PortRetOOS1M=[];
        PortRetOOS3M=[];
        PortRetOOS6M=[];
        InterPortRetOOS1M=[];
        InterPortRetOOS3M=[];
        InterPortRetOOS6M=[];
        ISportFDRsize=[];
        IOS1MportFDRsize=[];
        IOS3MportFDRsize=[];
        IOS6MportFDRsize=[];
        
        ISportFDRser=[];
        OOS1MportFDRser=[];
        OOS3MportFDRser=[];
        OOS6MportFDRser=[];
        IOS1MportFDRser=[];
        IOS3MportFDRser=[];
        IOS6MportFDRser=[];
        
        IOS1Mratio=[];
        IOS3Mratio=[];
        IOS6Mratio=[];
        mnthchk1Mser=[];
        mnthchk3Mser=[];
        mnthchk6Mser=[];
        Intermnthchk1Mser=[];
        Intermnthchk3Mser=[];
        Intermnthchk6Mser=[];
        
        for month=1:12
            for IS=2
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
                retIOS1M=ret(TF1SMP-TF1+1:TF3SMP1M-TF1+1,:);
                %% IS+3M OOS
                finishdtstrIOS3M=[num2str(rem(month+3*oosperiod,12+1e-10)),'/',num2str(yr+IS+(month+3*oosperiod>12))];
                TF3SMP3M = contains(master.textdata(:,1),finishdtstrIOS3M);
                TF3SMP3M=find(TF3SMP3M==1,1,'first')-1;
                finishdtIOS3M=master.textdata(TF3SMP3M,1);
                retIOS3M=ret(TF1SMP-TF1+1:TF3SMP3M-TF1+1,:);
                %% IS+6M OOS
                finishdtstrIOS6M=[num2str(rem(month+6*oosperiod,12+1e-10)),'/',num2str(yr+IS+(month+6*oosperiod>12))];
                TF3SMP6M = contains(master.textdata(:,1),finishdtstrIOS6M);
                TF3SMP6M=find(TF3SMP6M==1,1,'first')-1;
                finishdtIOS6M=master.textdata(TF3SMP6M,1);
                retIOS6M=ret(TF1SMP-TF1+1:TF3SMP6M-TF1+1,:);
                
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
                excessretIOS1M=retIOS1M-rftable{rfind1:rfind31M,2};
                excessretIOS3M=retIOS3M-rftable{rfind1:rfind33M,2};
                excessretIOS6M=retIOS6M-rftable{rfind1:rfind36M,2};
                excessretOOS1M=retOOS1M-rftable{rfind2+1:rfind31M,2};
                excessretOOS3M=retOOS3M-rftable{rfind2+1:rfind33M,2};
                excessretOOS6M=retOOS6M-rftable{rfind2+1:rfind36M,2};
                
                %% Loading the FDR portfolios
                ISportFDR=resvarIS{iterplus*2,1};
                IOS1MportFDR=resvarIOS{iterplus,1};
                IOS3MportFDR=resvarIOS{iterplus,2};
                IOS6MportFDR=resvarIOS{iterplus,3};
                %% Portfolio generation
                % IS portfolio return
                finalretIS=sum(excessretIS);
                PortRetIS(end+1)=finalretIS*ISportFDR/max(sum(ISportFDR),1);
                ISportFDRser(end+1:end+size(excessretIS,1))=mean(excessretIS(:,ISportFDR==1),2);
                % IOS+OOS portfolio return
                % 1M
                finalretOOS1M=sum(excessretOOS1M);
                PortRetOOS1M(end+1)=finalretOOS1M*ISportFDR/max(sum(ISportFDR),1);
                OOS1MportFDRser(end+1:end+size(excessretOOS1M,1))=mean(excessretOOS1M(:,ISportFDR==1),2);
                IOS1MportFDRser(end+1:end+size(excessretIOS1M,1))=mean(excessretIOS1M(:,IOS1MportFDR==1),2);
                % 3M
                finalretOOS3M=sum(excessretOOS3M);
                PortRetOOS3M(end+1)=finalretOOS3M*ISportFDR/max(sum(ISportFDR),1);
                OOS3MportFDRser(end+1:end+size(excessretOOS3M,1))=mean(excessretOOS3M(:,ISportFDR==1),2);
                IOS3MportFDRser(end+1:end+size(excessretIOS3M,1))=mean(excessretIOS3M(:,IOS3MportFDR==1),2);
                % 6M
                finalretOOS6M=sum(excessretOOS6M);
                PortRetOOS6M(end+1)=finalretOOS6M*ISportFDR/max(sum(ISportFDR),1);
                OOS6MportFDRser(end+1:end+size(excessretOOS6M,1))=mean(excessretOOS6M(:,ISportFDR==1),2);
                IOS6MportFDRser(end+1:end+size(excessretIOS6M,1))=mean(excessretIOS6M(:,IOS6MportFDR==1),2);
                
                %% Cross-validating the results
                intersect1=and(ISportFDR==1,and(finalretOOS1M'>0,ISportFDR==IOS1MportFDR));
                IOS1Mratio(end+1)=sum(intersect1)/21195;
                intersect3=and(ISportFDR==1,and(finalretOOS3M'>0,ISportFDR==IOS3MportFDR));
                IOS3Mratio(end+1)=sum(intersect3)/21195;
                intersect6=and(ISportFDR==1,and(finalretOOS6M'>0,ISportFDR==IOS6MportFDR));
                IOS6Mratio(end+1)=sum(intersect6)/21195;              
                
                %% Robustness check (survivors)
                % 1M 
                mnthchk1M=0;
                PortRetOOS=PortRetOOS1M(end);
                mnthchk1M=mnthchk1M+(PortRetOOS>0);
                while mnthchk1M<=6 && PortRetOOS>0
                    finishdtstrOOS=[num2str(rem(month+(mnthchk1M+1)*oosperiod,12+1e-10)),...
                        '/',num2str(yr+IS+(month+(mnthchk1M+1)*oosperiod>12))];
                    TF3SMPOOS = contains(master.textdata(:,1),finishdtstrOOS);
                    TF3SMPOOS=find(TF3SMPOOS==1,1,'first')-1;
                    finishdtOOS=master.textdata(TF3SMPOOS,1);
                    retOOS=ret(TF2SMP+1-TF1+1:TF3SMPOOS-TF1+1,:);
                    rfind3OOS=find(rftable{:,1}==char(finishdtOOS),1,'first');
                    excessretOOS=retOOS-rftable{rfind2+1:rfind3OOS,2};
                    finalretOOS=sum(excessretOOS);
                    PortRetOOSbuffer=finalretOOS*ISportFDR/max(sum(ISportFDR),1);
                    PortRetOOS=(PortRetOOSbuffer>PortRetOOS)*PortRetOOSbuffer+...
                        (1-(PortRetOOSbuffer>PortRetOOS))*(PortRetOOSbuffer-PortRetOOS);
                    mnthchk1M=mnthchk1M+(PortRetOOS>0);
                    
                end
                % 3M 
                mnthchk3M=0;
                PortRetOOS=PortRetOOS3M(end);
                mnthchk3M=mnthchk3M+3*(PortRetOOS>0);
                while mnthchk3M<=12 && PortRetOOS>0
                    finishdtstrOOS=[num2str(rem(month+(mnthchk3M+3)*oosperiod,12+1e-10)),...
                        '/',num2str(yr+IS+floor((month+mnthchk3M+3)*oosperiod/12))];
                    TF3SMPOOS = contains(master.textdata(:,1),finishdtstrOOS);
                    TF3SMPOOS=find(TF3SMPOOS==1,1,'first')-1;
                    finishdtOOS=master.textdata(TF3SMPOOS,1);
                    retOOS=ret(TF2SMP+1-TF1+1:TF3SMPOOS-TF1+1,:);
                    rfind3OOS=find(rftable{:,1}==char(finishdtOOS),1,'first');
                    excessretOOS=retOOS-rftable{rfind2+1:rfind3OOS,2};
                    finalretOOS=sum(excessretOOS);
                    PortRetOOSbuffer=finalretOOS*ISportFDR/max(sum(ISportFDR),1);
                    PortRetOOS=(PortRetOOSbuffer>PortRetOOS)*PortRetOOSbuffer+...
                        (1-(PortRetOOSbuffer>PortRetOOS))*(PortRetOOSbuffer-PortRetOOS);
                    mnthchk3M=mnthchk3M+3*(PortRetOOS>0);
                    
                end
                % 6M 
                mnthchk6M=0;
                PortRetOOS=PortRetOOS6M(end);
                mnthchk6M=mnthchk6M+6*(PortRetOOS>0);
                while mnthchk6M<=12 && PortRetOOS>0
                    finishdtstrOOS=[num2str(rem(month+(mnthchk6M+6)*oosperiod,12+1e-10)),...
                        '/',num2str(yr+IS+floor((month+mnthchk6M+6)*oosperiod/12))];
                    TF6SMPOOS = contains(master.textdata(:,1),finishdtstrOOS);
                    TF6SMPOOS=find(TF6SMPOOS==1,1,'first')-1;
                    finishdtOOS=master.textdata(TF6SMPOOS,1);
                    retOOS=ret(TF2SMP+1-TF1+1:TF6SMPOOS-TF1+1,:);
                    rfind6OOS=find(rftable{:,1}==char(finishdtOOS),1,'first');
                    excessretOOS=retOOS-rftable{rfind2+1:rfind6OOS,2};
                    finalretOOS=sum(excessretOOS);
                    PortRetOOSbuffer=finalretOOS*ISportFDR/max(sum(ISportFDR),1);
                    PortRetOOS=(PortRetOOSbuffer>PortRetOOS)*PortRetOOSbuffer+...
                        (1-(PortRetOOSbuffer>PortRetOOS))*(PortRetOOSbuffer-PortRetOOS);
                    mnthchk6M=mnthchk6M+6*(PortRetOOS>0);
                end
                mnthchk1Mser(end+1)=mnthchk1M;
                mnthchk3Mser(end+1)=mnthchk3M;
                mnthchk6Mser(end+1)=mnthchk6M;
                %% Intersection Portfolio Performance
                InterPortOOS1M=mean(excessretOOS1M(:,intersect1),2);
                InterPortOOS3M=mean(excessretOOS3M(:,intersect3),2);
                InterPortOOS6M=mean(excessretOOS6M(:,intersect6),2);
                InterPortRetOOS1M(end+1)=sum(InterPortOOS1M);
                InterPortRetOOS3M(end+1)=sum(InterPortOOS3M);
                InterPortRetOOS6M(end+1)=sum(InterPortOOS6M);
                %% Robustness check (Intersection Portfolio)
                % 1M 
                Intermnthchk1M=0;
                InterPortRetOOS=InterPortRetOOS1M(end);
                Intermnthchk1M=Intermnthchk1M+(InterPortRetOOS>0);
                while Intermnthchk1M<=6 && InterPortRetOOS>0
                    finishdtstrOOS=[num2str(rem(month+(Intermnthchk1M+1)*oosperiod,12+1e-10)),...
                        '/',num2str(yr+IS+(month+(Intermnthchk1M+1)*oosperiod>12))];
                    TF3SMPOOS = contains(master.textdata(:,1),finishdtstrOOS);
                    TF3SMPOOS=find(TF3SMPOOS==1,1,'first')-1;
                    finishdtOOS=master.textdata(TF3SMPOOS,1);
                    retOOS=ret(TF2SMP+1-TF1+1:TF3SMPOOS-TF1+1,:);
                    rfind3OOS=find(rftable{:,1}==char(finishdtOOS),1,'first');
                    excessretOOS=retOOS-rftable{rfind2+1:rfind3OOS,2};
                    finalretOOS=sum(excessretOOS);
                    InterPortRetOOSbuffer=finalretOOS*intersect1/max(sum(intersect1),1);
                    InterPortRetOOS=(InterPortRetOOSbuffer>InterPortRetOOS)*InterPortRetOOSbuffer+...
                        (1-(InterPortRetOOSbuffer>InterPortRetOOS))*(InterPortRetOOSbuffer-InterPortRetOOS);
                    Intermnthchk1M=Intermnthchk1M+(InterPortRetOOS>0);
                    
                end
                % 3M 
                Intermnthchk3M=0;
                InterPortRetOOS=InterPortRetOOS3M(end);
                Intermnthchk3M=Intermnthchk3M+3*(InterPortRetOOS>0);
                while Intermnthchk3M<=12 && InterPortRetOOS>0
                    finishdtstrOOS=[num2str(rem(month+(Intermnthchk3M+3)*oosperiod,12+1e-10)),...
                        '/',num2str(yr+IS+floor((month+Intermnthchk3M+3)*oosperiod/12))];
                    TF3SMPOOS = contains(master.textdata(:,1),finishdtstrOOS);
                    TF3SMPOOS=find(TF3SMPOOS==1,1,'first')-1;
                    finishdtOOS=master.textdata(TF3SMPOOS,1);
                    retOOS=ret(TF2SMP+1-TF1+1:TF3SMPOOS-TF1+1,:);
                    rfind3OOS=find(rftable{:,1}==char(finishdtOOS),1,'first');
                    excessretOOS=retOOS-rftable{rfind2+1:rfind3OOS,2};
                    finalretOOS=sum(excessretOOS);
                    InterPortRetOOSbuffer=finalretOOS*intersect3/max(sum(intersect3),1);
                    InterPortRetOOS=(InterPortRetOOSbuffer>InterPortRetOOS)*InterPortRetOOSbuffer+...
                        (1-(InterPortRetOOSbuffer>InterPortRetOOS))*(InterPortRetOOSbuffer-InterPortRetOOS);
                    Intermnthchk3M=Intermnthchk3M+3*(InterPortRetOOS>0);
                    
                end
                % 6M 
                Intermnthchk6M=0;
                InterPortRetOOS=InterPortRetOOS6M(end);
                Intermnthchk6M=Intermnthchk6M+6*(InterPortRetOOS>0);
                while Intermnthchk6M<=12 && InterPortRetOOS>0
                    finishdtstrOOS=[num2str(rem(month+(Intermnthchk6M+6)*oosperiod,12+1e-10)),...
                        '/',num2str(yr+IS+floor((month+Intermnthchk6M+6)*oosperiod/12))];
                    TF6SMPOOS = contains(master.textdata(:,1),finishdtstrOOS);
                    TF6SMPOOS=find(TF6SMPOOS==1,1,'first')-1;
                    finishdtOOS=master.textdata(TF6SMPOOS,1);
                    retOOS=ret(TF2SMP+1-TF1+1:TF6SMPOOS-TF1+1,:);
                    rfind6OOS=find(rftable{:,1}==char(finishdtOOS),1,'first');
                    excessretOOS=retOOS-rftable{rfind2+1:rfind6OOS,2};
                    finalretOOS=sum(excessretOOS);
                    InterPortRetOOSbuffer=finalretOOS*intersect6/max(sum(intersect6),1);
                    InterPortRetOOS=(InterPortRetOOSbuffer>InterPortRetOOS)*InterPortRetOOSbuffer+...
                        (1-(InterPortRetOOSbuffer>InterPortRetOOS))*(InterPortRetOOSbuffer-InterPortRetOOS);
                    Intermnthchk6M=Intermnthchk6M+6*(InterPortRetOOS>0);
                end
                Intermnthchk1Mser(end+1)=Intermnthchk1M;
                Intermnthchk3Mser(end+1)=Intermnthchk3M;
                Intermnthchk6Mser(end+1)=Intermnthchk6M;
                ISportFDRsize(end+1)=sum(ISportFDR);
                IOS1MportFDRsize(end+1)=sum(IOS1MportFDR);
                IOS3MportFDRsize(end+1)=sum(IOS3MportFDR);
                IOS6MportFDRsize(end+1)=sum(IOS6MportFDR);
            end           
        end
        %% Sharpe ratio
        ISSR=sharpe(ISportFDRser,0)*260^.5;
        OOS1MSR=sharpe(OOS1MportFDRser,0)*260^.5;
        OOS3MSR=sharpe(OOS3MportFDRser,0)*260^.5;
        OOS6MSR=sharpe(OOS6MportFDRser,0)*260^.5;
        IOS1MSR=sharpe(IOS1MportFDRser,0)*260^.5;
        IOS3MSR=sharpe(IOS3MportFDRser,0)*260^.5;
        IOS6MSR=sharpe(IOS6MportFDRser,0)*260^.5;
        
        InterPortRetOOS1M(isnan(InterPortRetOOS1M))=0;
        InterPortRetOOS3M(isnan(InterPortRetOOS3M))=0;
        InterPortRetOOS6M(isnan(InterPortRetOOS6M))=0;
        %% Recording the results
        % Period
        dim=1;
        resvar{iter,dim}=yr;
        % IS return (survivors)
        dim=dim+1;
        resvar{iter,dim}=(1+mean(PortRetIS))^(1/IS)-1;
        % IOS return (full dataset survivors)
        dim=dim+1;
        resvar{iter,dim}=(1+sum(IOS1MportFDRser)/12)^(12/(IS*12+1))-1;
        dim=dim+1;
        resvar{iter,dim}=(1+sum(IOS3MportFDRser)/12)^(12/(IS*12+3))-1;
        dim=dim+1;
        resvar{iter,dim}=(1+sum(IOS6MportFDRser)/12)^(12/(IS*12+6))-1;
        % OOS return (survivors)
        dim=dim+1;
        resvar{iter,dim}=(1+mean(PortRetOOS1M))^12-1;
        dim=dim+1;
        resvar{iter,dim}=(1+mean(PortRetOOS3M))^4-1;
        dim=dim+1;
        resvar{iter,dim}=(1+mean(PortRetOOS6M))^2-1;
        % OOS return (intersection portfolio)
        dim=dim+1;
        resvar{iter,dim}=(1+mean(InterPortRetOOS1M))^12-1;
        dim=dim+1;
        resvar{iter,dim}=(1+mean(InterPortRetOOS3M))^4-1;
        dim=dim+1;
        resvar{iter,dim}=(1+mean(InterPortRetOOS6M))^2-1;
        % Portfolio size
        dim=dim+1;
        resvar{iter,dim}=mean(ISportFDRsize);
        dim=dim+1;
        resvar{iter,dim}=mean(IOS1MportFDRsize);
        dim=dim+1;
        resvar{iter,dim}=mean(IOS3MportFDRsize);
        dim=dim+1;
        resvar{iter,dim}=mean(IOS6MportFDRsize);
        % Sharpe ratio
        dim=dim+1;
        resvar{iter,dim}=ISSR;
        dim=dim+1;
        resvar{iter,dim}=OOS1MSR;
        dim=dim+1;
        resvar{iter,dim}=OOS3MSR;
        dim=dim+1;
        resvar{iter,dim}=OOS6MSR;
        dim=dim+1;
        resvar{iter,dim}=IOS1MSR;
        dim=dim+1;
        resvar{iter,dim}=IOS3MSR;
        dim=dim+1;
        resvar{iter,dim}=IOS6MSR;
        dim=dim+1;
        % Intersection Ratio
        resvar{iter,dim}=mean(IOS1Mratio);
        dim=dim+1;
        resvar{iter,dim}=mean(IOS3Mratio);
        dim=dim+1;
        resvar{iter,dim}=mean(IOS6Mratio);
        % Robustness check (IS survivors)
        dim=dim+1;
        resvar{iter,dim}=mean(mnthchk1Mser);
        dim=dim+1;
        resvar{iter,dim}=mean(mnthchk3Mser)/3;
        dim=dim+1;
        resvar{iter,dim}=mean(mnthchk6Mser)/6;
        % Robustness check (intersection)
        dim=dim+1;
        resvar{iter,dim}=mean(Intermnthchk1Mser);
        dim=dim+1;
        resvar{iter,dim}=mean(Intermnthchk3Mser);
        dim=dim+1;
        resvar{iter,dim}=mean(Intermnthchk6Mser);
        dim=dim+1;
        resvar{iter,dim}=regions{i};
        disp(['Processing done for ',regions{i},' period ',num2str(yearOOS), ' completed ...']);
        toc;
    end
    
end
result_table=cell2table(resvar);
result_table.Properties.VariableNames={'Year' 'FDRReturnIS' ...
    'FDRReturnIOS1M' 'FDRReturnIOS3M' 'FDRReturnIOS6M'...
    'FDRReturnOOS1M' 'FDRReturnOOS3M' 'FDRReturnOOS6M' ...
    'InterReturnOOS1M' 'InterReturnOOS3M' 'InterReturnOOS6M'...
    'FDRPortSizeIS' 'FDRPortSizeIOS1M' 'FDRPortSizeIOS3M' 'FDRPortSizeIOS6M' ...
    'SharpeRatioIS' 'SharpeRatioOOS1M' 'SharpeRatioOOS3M' 'SharpeRatioOOS6M' ...
    'SharpeRatioIOS1M' 'SharpeRatioIOS3M' 'SharpeRatioIOS6M' ...
    'MatchRatio1M' 'MatchRatio3M' 'MatchRatio6M' ...
    'Robustness_1MN' 'Robustness_3MN' 'Robustness_6MN'...
    'InterRobustness_1MN' 'InterRobustness_3MN' 'InterRobustness_6MN' 'Country'};
%drpadd='/Empirics3/';
drpadd='C:\Dropbox\Dropbox\Shared_Folder\Chapter3\';
writetable(result_table,[drpadd,'CV_AnnualTable_',num2str(yearser(1)),'_',num2str(yearser(end)),'FF.xlsx'])
toc