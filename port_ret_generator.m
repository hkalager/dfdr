% Created 10 Oct 2017
% Last Modified 22 Oct 2017 23:31
clear;
clc;
labels={'NDDUUS'	'NDDUUK'	'NDDUJN'	'NDEUSRU'	'NDEUCHF'	'NDUEBRAF'...
    'MSEIESUN'	'M1MA'	'MIMUJORN' 'MXWO' 'MXEF' 'MXFEM'};
rftable=readtable('/Users/arman/Dropbox (Personal)/Shared_Folder/Chapter3/FEDL01.csv');
rftable{:,2}=rftable{:,2}/100/260;
transser=[25*ones(1,3),50*ones(1,6),25,50,50];

disp('Checking the datasets');
start=2002;
finish=2016;
oosperiod=1; % Out of sample holding period

%% FDR setting
fdrtarget=0.1;
%% switch function between STUDENTIZED and GENERAL test statistic
whichstd=@(x,y) y*std(x)+(1-y)*ones(1,size(x,2));
yearser=2006:2015;
cd('/Volumes/My Passport/Empirics3');

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
    
    load([lbl,'_',num2str(yearser(1)),'_',num2str(yearser(end)),'.mat'],'resvar');
    tbl_main=table();
    % Calculating the returns without transaction costs
    disp(['Dataset loaded for ',lbl]);
    iter=0;
    for yearOOS=yearser(1):yearser(end)%2006:2016
        for month=1:12
            for IS=2
                iter=iter+1;
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
                finishdtstrOOS1M=[num2str(rem(month,12)+oosperiod),'/',num2str(year+IS+(month+oosperiod>12))];
                TF3SMP1M = contains(master.textdata(:,1),finishdtstrOOS1M);
                TF3SMP1M=find(TF3SMP1M==1,1,'first')-1;
                finishdtOOS1M=master.textdata(TF3SMP1M,1);
                %% 3M OOS
                finishdtstrOOS3M=[num2str(rem(month+3*oosperiod,12)),'/',num2str(year+IS+(month+3*oosperiod>12))];
                TF3SMP3M = contains(master.textdata(:,1),finishdtstrOOS3M);
                TF3SMP3M=find(TF3SMP3M==1,1,'first')-1;
                finishdtOOS3M=master.textdata(TF3SMP3M,1);
                %% 6M OOS
                finishdtstrOOS6M=[num2str(rem(month+6*oosperiod,12)),'/',num2str(year+IS+(month+6*oosperiod>12))];
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
                
                PORTFDR=resvar{iter,1};
                
                dt_ser=rftable{rfind2+1:rfind31M,1};
                PortRetOOS1M_ser= mean(retOOS1M(:,logical(PORTFDR)),2);
                PortRetOOS1M_excess_ser= mean(excessretOOS1M(:,logical(PORTFDR)),2);
                tb1=table(dt_ser,PortRetOOS1M_ser,PortRetOOS1M_excess_ser);
                
                tbl_main=[tbl_main;tb1];
                
                
                
                %xlswrite([lbl,'xlsx'],result_table);
            end
        end
        
    end
    tbl_main.Properties.VariableNames={'Date','Return','Excess_Return'};
    drpadd='/Users/arman/Dropbox (Personal)/Shared_Folder/Chapter3/';
    %drpadd='C:\Users\Hassannia\Dropbox\Shared_Folder\Chapter3\';
    writetable(tbl_main,[drpadd,lbl,'_',num2str(yearser(1)),'_',num2str(yearser(end)),'_port_returns.xlsx'])
    
end
cd(drpadd);
