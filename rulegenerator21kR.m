% Create Technical trading rules for input symbol,
% Script developed by Arman Hassannia Kalager for completing 3rd chapter,
% Created on 28 Oct 2015 11:15 GMT,
% Last modified 06 Dec 2017 15:30 GMT.
% ****************************************
% Input arguments:
% sym1: The desired ticker name available in directory e.g. 'EURUSD''1440'
% transbasis: Transaction cost in basis points (Optional, default=3)
% typ: Type of data, 1 for insample & 2 for out-of-sample
% start: Start date of dataset ('YYYY.MM.DD')
% finish: Last date of dataset ('YYYY.MM.DD')
% ****************************************
% Outputs:
% A file with all information required named as tickername-in/outsample

function rulegenerator21kR(sym1,transbasis,typ,start,finish)
tic
rng(0);
sym=[sym1,'.csv'];
if  exist(sym,'file')==0
    error('No file with such name available for daily period, retry!');
end
master=importdata(sym,',');
if ischar(start)
    TF1 = contains(master.textdata(:,1),start);
else 
    TF1 = contains(master.textdata(:,1),mat2str(start));
end
TF1=find(TF1==1,1,'first');
startdt=master.textdata(TF1,1);
if ischar(finish)
    TF2 = contains(master.textdata(:,1),finish);
else 
    TF2 = contains(master.textdata(:,1),mat2str(finish));
end
TF2=find(TF2==1,1,'last');
finishdt=master.textdata(TF2,1);
switch typ
    case 1
        sym3='insample';
        initstart=startdt;
    case 2
        sym3='outsample';
        fllbl=[sym1,'_',num2str(start),'_',num2str(finish),'-',sym3,'-',num2str(transbasis)];
        load(fllbl,'initstart');
        gap=find(strcmp(master.textdata,startdt))-find(strcmp(master.textdata,initstart));
    
end

data=master.data(find(strcmp(master.textdata,initstart)):find(strcmp(master.textdata,finishdt)),max(1,end-1));
%highvals=master.data(find(strcmp(master.textdata,initstart)):find(strcmp(master.textdata,finishdt)),2);
%lowvals=master.data(find(strcmp(master.textdata,initstart)):find(strcmp(master.textdata,finishdt)),3);
%volumes=master.data(find(strcmp(master.textdata,initstart)):find(strcmp(master.textdata,finishdt)),5);
data1=master.data(find(strcmp(master.textdata,initstart)):find(strcmp(master.textdata,finishdt))+1,max(1,end-1));
disp(' Data imported! Processing initiated');
datsiz=max(size(data));
dataret=price2ret(data1,[],'continuous');
datadir=sign(dataret);
paramstrct=struct('type',[],'case',[],'params',[]); % parameters of each trading rule
if isempty(transbasis)
    transbasis=25;
end
trncost=transbasis*1e-4;
% Alternative in case of CFD calculation:
%transcost=transbasis*1e-4*10^floor(log10(mean(data)));

%% Variable "iter1" keeps account of trading rule
%  Variable "iter2" keeps account of observations (datapoints)
% These variables appear on all trading strategies
%% Relative Strength Index (RSI)
disp(' RSI rule generation started...');
hserRSI=[5,10,15,20,25,50,100,150,200,250];
vserRSI=10:5:25;
dserRSI=[1,2,5];
kserRSI=[1,5,10,25];
rrh=numel(hserRSI);
rrv=numel(vserRSI);
rrd=numel(dserRSI);
rrk=numel(kserRSI);
%RSI case 1 (h,v,d)
iter1=1;
pos=0;
inputser(1,1:rrh*rrv*rrd)=pos;
ope=nan;
aaa=1;
for iter31=1:rrh
    for iter32=1:rrv
        for iter33=1:rrd
            RSIser=rsindex(data,hserRSI(iter31));
            for iter2=hserRSI(iter31)+dserRSI(iter33)+1:datsiz
                if pos==0
                    if min(RSIser(iter2-dserRSI(iter33)-1:iter2-1))>50+vserRSI(iter32) && RSIser(iter2)<50+vserRSI(iter32)
                        ope=data(iter2);
                        hld=iter2;
                        pos=-1;
                    elseif max(RSIser(iter2-dserRSI(iter33)-1:iter2-1))<50-vserRSI(iter32) && RSIser(iter2)>50-vserRSI(iter32)
                        ope=data(iter2);
                        hld=iter2;
                        pos=1;
                    end
                elseif pos==1
                    if RSIser(iter2)<50+vserRSI(iter32)
                        close=data(iter2);
                        pos=-1;
                        hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                        aaa=aaa+1;
                        ope=close;
                    end
                elseif pos==-1
                    if RSIser(iter2)>50-vserRSI(iter32)
                        close=data(iter2);
                        pos=1;
                        hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                        aaa=aaa+1;
                        ope=close;
                    end
                    
                end
                inputser(iter2,iter1)=pos;
            end
            paramstrct.type{iter1}='RSI';
            paramstrct.case(iter1)=1;
            paramstrct.params{iter1}=['h=',num2str(hserRSI(iter31)),',v=',num2str(vserRSI(iter32)),',d=',num2str(dserRSI(iter33))];
            %disp([' Iteration ',num2str(iter1),' finished...']);
            iter1=iter1+1;
            aaa=1;
            pos=0;
        end
    end
end

%RSI case 2 (h,v,d,k)
pos=0;
inputser(1,1:rrh*rrv*rrd)=pos;
ope=nan;
aaa=1;
for iter31=1:rrh
    for iter32=1:rrv
        for iter33=1:rrd
            for iter34=1:rrk
                RSIser=rsindex(data,hserRSI(iter31));
                for iter2=hserRSI(iter31)+dserRSI(iter33)+1:datsiz
                    if pos==0
                        if min(RSIser(iter2-dserRSI(iter33)-1:iter2-1))>50+vserRSI(iter32) && RSIser(iter2)<50+vserRSI(iter32)
                            ope=data(iter2);
                            hld=iter2;
                            pos=-1;
                        elseif max(RSIser(iter2-dserRSI(iter33)-1:iter2-1))<50-vserRSI(iter32) && RSIser(iter2)>50-vserRSI(iter32)
                            ope=data(iter2);
                            hld=iter2;
                            pos=1;
                        end
                    elseif pos==1
                        if iter2-hld==kserRSI(iter34)
                            close=data(iter2);
                            pos=0;
                            hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                            aaa=aaa+1;
                            ope=close;
                        end
                    elseif pos==-1
                        if iter2-hld==kserRSI(iter34)
                            close=data(iter2);
                            pos=0;
                            hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                            aaa=aaa+1;
                            ope=close;
                        end
                        
                    end
                    inputser(iter2,iter1)=pos;
                end
                paramstrct.type{iter1}='RSI';
                paramstrct.case(iter1)=2;
                paramstrct.params{iter1}=['h=',num2str(hserRSI(iter31)),',v=',num2str(vserRSI(iter32)),',d=',num2str(dserRSI(iter33)),',k=',num2str(kserRSI(iter34))];
                %disp([' Iteration ',num2str(iter1),' finished...']);
                iter1=iter1+1;
                aaa=1;
                pos=0;
            end
        end
    end
end
disp(' RSI rule generation completed successfully!');
toc

%% Filter rule generation (filt)
disp(' Filter rule generation started...');
xserfilt=[0.0005,0.001,0.005,0.01,0.05,0.10,0.2]; % alternative in %
yserfilt=[0.0005,0.001,0.005,0.01,0.05,0.10,0.2]; % alternative in %
dserfilt=0:5;
jserfilt=[1,2,5,10,20];
kserfilt=[5,10,15,20,25];
aa=numel(xserfilt);
%Filter case 1 (x,d,j)
pos=0;
aab=numel(jserfilt);
aaz=numel(dserfilt);
zzy=find(inputser(1,:)==0, 1, 'last' );
inputser(1,zzy+1:zzy+aa*aab*aaz)=pos;
ope=nan;
aaa=1;
for iter91=1:aaz % for all d
    delay=dserfilt(iter91);
    for iter3=1:aab % for all j
        hlperiod=jserfilt(iter3);
        for iter4=1:aa % for all x
            filter=xserfilt(iter4);
            for iter2=hlperiod+delay+1:datsiz
                if pos==0
                    if min(data(iter2-delay:iter2))/min(data(iter2-delay-hlperiod:iter2-delay-1))-1>filter % Upward trend
                        ope=data(iter2);hld=iter2;
                        pos=1;
                    elseif max(data(iter2-delay:iter2))/max(data(iter2-delay-hlperiod:iter2-delay-1))-1<-filter % Downward trend
                        ope=data(iter2);hld=iter2;
                        pos=-1;
                    end
                elseif pos==1
                    if max(data(iter2-delay:iter2))/max(data(iter2-delay-hlperiod:iter2-delay-1))-1<-filter % Downward trend
                        close=data(iter2);
                        pos=-1;
                        hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                        aaa=aaa+1;
                        ope=close;
                    end
                    
                elseif pos==-1
                    if min(data(iter2-delay:iter2))/min(data(iter2-delay-hlperiod:iter2-delay-1))-1>filter % Upward trend
                        close=data(iter2);
                        pos=1;
                        hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                        aaa=aaa+1;
                        ope=close;
                    end
                end
                inputser(iter2,iter1)=pos;
            end
            paramstrct.type{iter1}='Filter';
            paramstrct.case(iter1)=1;
            paramstrct.params{iter1}=['x=',num2str(xserfilt(iter4)),',j=',num2str(jserfilt(iter3)),',d=',num2str(dserfilt(iter91))];
            iter1=iter1+1;
            aaa=1;
            pos=0;
        end
    end
end

%Filter case 2 (x-y,dx-dy,j)
pos=0;
zzy=find(inputser(1,:)==0, 1, 'last' );
xycombcount=sum(sum(xserfilt(:)>yserfilt));
dxdycombcount=sum(sum(dserfilt(:)>dserfilt));
inputser(1,zzy+1:zzy+xycombcount*dxdycombcount*aab)=pos;
ope=nan;
aaa=1;
for iter911=2:aaz % for all dx
    for iter912=1:iter911-1 % for dy<dx
        delay1=dserfilt(iter911);
        delay2=dserfilt(iter912);
        for iter3=1:aab % for all j
            for iter41=2:aa % for all x
                for iter42=1:iter41-1 % for all y<x
                    filter1=xserfilt(iter41);
                    filter2=yserfilt(iter42);
                    hlperiod=jserfilt(iter3);
                    for iter2=hlperiod+delay1+1:datsiz
                        if pos==0
                            if min(data(iter2-delay1:iter2))/min(data(iter2-delay1-hlperiod:iter2-delay1-1))-1>filter1 % Opening Long position
                                ope=data(iter2);hld=iter2;
                                pos=1;
                            elseif max(data(iter2-delay1:iter2))/max(data(iter2-delay1-hlperiod:iter2-delay1-1))-1<-filter1 % Opening Short position
                                ope=data(iter2);hld=iter2;
                                pos=-1;
                                
                            end
                        elseif pos==1
                            if max(data(iter2-delay2:iter2))/max(data(iter2-delay2-hlperiod:iter2-delay2-1))-1<-filter2 % Downward trend
                                close=data(iter2);
                                pos=-1;
                                hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                                aaa=aaa+1;
                                ope=close;
                            end
                            
                        elseif pos==-1
                            if min(data(iter2-delay2:iter2))/min(data(iter2-delay2-hlperiod:iter2-delay2-1))-1>filter2 % Upward trend
                                close=data(iter2);
                                pos=1;
                                hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                                aaa=aaa+1;
                                ope=close;
                            end
                        end
                        inputser(iter2,iter1)=pos;
                        
                    end
                    paramstrct.type{iter1}='Filter';
                    paramstrct.case(iter1)=2;
                    paramstrct.params{iter1}=['x=',num2str(xserfilt(iter41)),'y=',num2str(xserfilt(iter42)),',j=',num2str(jserfilt(iter3)),',dx=',num2str(dserfilt(iter911)),',dy=',num2str(dserfilt(iter912))];
                    iter1=iter1+1;
                    aaa=1;
                    pos=0;
                end
            end
        end
    end
end

%Filter case 3 (x,j,d,k)
pos=0;
aac=numel(kserfilt);
zzy=find(inputser(1,:)==0, 1, 'last' );
inputser(1,zzy+1:zzy+aa*aab*aaz*aac)=pos;
ope=nan;
aaa=1;
for iter44=1:aac % for all k
    for iter91=1:aaz % for all d
        for iter3=1:aab % for all j
            for iter4=1:aa % for all x
                filter=xserfilt(iter4);
                delay=dserfilt(iter91);
                hlperiod=jserfilt(iter3);
                neutper=kserfilt(iter44); % Time to neutralize (fixed holding period)
                for iter2=hlperiod+delay+1:datsiz
                    if pos==0
                        if min(data(iter2-delay:iter2))/min(data(iter2-delay-hlperiod:iter2-delay-1))-1>filter % Upward trend
                            ope=data(iter2);hld=iter2;
                            pos=1;
                        elseif max(data(iter2-delay:iter2))/max(data(iter2-delay-hlperiod:iter2-delay-1))-1<-filter % Downward trend
                            ope=data(iter2);hld=iter2;
                            pos=-1;
                            
                        end
                    elseif pos==1
                        if iter2-hld==neutper % Neutralizing condition (fixed time)
                            close=data(iter2);
                            pos=0;
                            hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                            aaa=aaa+1;
                            ope=close;
                        end
                        
                    elseif pos==-1
                        if iter2-hld==neutper % Neutralizing condition (fixed time)
                            close=data(iter2);
                            pos=0;
                            hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                            aaa=aaa+1;
                            ope=close;
                        end
                    end
                    inputser(iter2,iter1)=pos;
                end
                paramstrct.type{iter1}='Filter';
                paramstrct.case(iter1)=3;
                paramstrct.params{iter1}=['x=',num2str(xserfilt(iter4)),',j=',num2str(jserfilt(iter3)),',d=',num2str(dserfilt(iter91)),',k=',num2str(kserfilt(iter44))];
                iter1=iter1+1;
                aaa=1;
                pos=0;
            end
        end
    end
end

disp(' Filter rule generation completed successfully!');
toc

%% Moving average crossover rules generation (cross)
disp(' MA Crossover rules generation started...');
qserMA=[2,5,10,15,20,25,50,100,150,200,250];
pqMAcmb=combnk(qserMA,2);
npqMAcmb=combnk(qserMA,3);
xserMA=[0,.05,.1,.5,1,5]/100;
dserMA=[0,2:5];
kserMA=[5,10,25];
% MA crosover case 1 (q,d,x)
pos=0;
zzy=find(inputser(1,:)==0, 1, 'last' );
inputser(1,zzy+1:zzy+numel(qserMA)*numel(xserMA)*numel(dserMA))=pos;
ope=nan;
aaa=1;
for iter3=1:numel(qserMA) % for all q
    lagper=qserMA(iter3);
    slowma=tsmovavg(data,'s',lagper,1);
    for iter91=1:numel(dserMA) % for all d
        for iter4=1:numel(xserMA) % for all x
            filter=xserMA(iter4);
            delay=dserMA(iter91);
            for iter2=lagper+delay+1:datsiz
                if pos==0
                    if min(data(iter2-delay:iter2))/max(slowma(iter2-delay:iter2-1))-1>filter % Upward trend
                        ope=data(iter2);hld=iter2;
                        pos=1;
                    elseif max(data(iter2-delay:iter2))/min(slowma(iter2-delay:iter2-1))-1<-filter % Downward trend
                        ope=data(iter2);hld=iter2;
                        pos=-1;
                        
                    end
                elseif pos==1
                    if max(data(iter2-delay:iter2))/min(slowma(iter2-delay:iter2-1))-1<-filter % Downward trend
                        close=data(iter2);
                        pos=-1;
                        hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                        aaa=aaa+1;
                        ope=close;
                    end
                    
                elseif pos==-1
                    if min(data(iter2-delay:iter2))/max(slowma(iter2-delay:iter2-1))-1>filter % Upward trend
                        close=data(iter2);
                        pos=1;
                        hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                        aaa=aaa+1;
                        ope=close;
                    end
                end
                inputser(iter2,iter1)=pos;
            end
            paramstrct.type{iter1}='MA Cross';
            paramstrct.case(iter1)=1;
            paramstrct.params{iter1}=['x=',num2str(filter),',q=',num2str(lagper),',d=',num2str(delay)];
            iter1=iter1+1;
            aaa=1;
            pos=0;
        end
    end
end

% MA crosover case 2 (q,x,d,k)
zzy=find(inputser(1,:)==0, 1, 'last' );
inputser(1,zzy+1:zzy++numel(qserMA)*numel(xserMA)*numel(dserMA)*numel(kserMA))=pos;
ope=nan;
aaa=1;
for iter3=1:numel(qserMA) % for all q
    lagper=qserMA(iter3);
    slowma=tsmovavg(data,'s',lagper,1);
    for iter44=1:numel(kserMA) % for all k
        neutper=kserMA(iter44); % Time to neutralize (fixed holding period)
        for iter91=1:numel(dserMA) % for all d
            delay=dserMA(iter91);
            for iter4=1:numel(xserMA) % for all x
                filter=xserMA(iter4);
                for iter2=lagper+delay+1:datsiz
                    if pos==0
                        if min(data(iter2-delay:iter2))/max(slowma(iter2-delay:iter2-1))-1>filter % Upward trend
                            ope=data(iter2);hld=iter2;
                            pos=1;
                        elseif max(data(iter2-delay:iter2))/min(slowma(iter2-delay:iter2-1))-1<-filter % Downward trend
                            ope=data(iter2);hld=iter2;
                            pos=-1;
                            
                        end
                    elseif pos==1
                        if iter2-hld==neutper % Neutralizing condition (fixed time)
                            close=data(iter2);
                            pos=0;
                            hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                            aaa=aaa+1;
                            ope=close;
                        end
                        
                    elseif pos==-1
                        if iter2-hld==neutper % Neutralizing condition (fixed time)
                            close=data(iter2);
                            pos=0;
                            hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                            aaa=aaa+1;
                            ope=close;
                        end
                    end
                    inputser(iter2,iter1)=pos;
                end
                paramstrct.type{iter1}='MA Cross';
                paramstrct.case(iter1)=2;
                paramstrct.params{iter1}=['x=',num2str(filter),',q=',num2str(lagper),',d=',num2str(delay),',k=',num2str(kserfilt(iter44))];
                iter1=iter1+1;
                aaa=1;
                pos=0;
            end
        end
    end
end
zzy=find(inputser(1,:)==0, 1, 'last' );
% MA crosover case 3 (p-q,d,x)
pqcmbMA=size(pqMAcmb,1);
inputser(1,zzy+1:zzy+pqcmbMA*numel(xserMA)*numel(dserMA))=pos;
ope=nan;
aaa=1;
for iter3=1:pqcmbMA % for all p-q
    lagper1=pqMAcmb(iter3,1);
    lagper2=pqMAcmb(iter3,2);
    fastma=tsmovavg(data,'s',lagper1,1);
    slowma=tsmovavg(data,'s',lagper2,1);
    for iter91=1:numel(dserMA) % for all d
        for iter4=1:numel(xserMA) % for all x
            filter=xserMA(iter4);
            delay=dserMA(iter91);
            for iter2=lagper+delay+1:datsiz
                if pos==0
                    if min(fastma(iter2-delay:iter2-1))/max(slowma(iter2-delay:iter2-1))-1>filter % Upward trend
                        ope=data(iter2);hld=iter2;
                        pos=1;
                    elseif max(fastma(iter2-delay:iter2-1))/min(slowma(iter2-delay:iter2-1))-1<-filter % Downward trend
                        ope=data(iter2);hld=iter2;
                        pos=-1;
                        
                    end
                elseif pos==1
                    if max(fastma(iter2-delay:iter2-1))/min(slowma(iter2-delay:iter2-1))-1<-filter % Downward trend
                        close=data(iter2);
                        pos=-1;
                        hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                        aaa=aaa+1;
                        ope=close;
                    end
                    
                elseif pos==-1
                    if min(fastma(iter2-delay:iter2-1))/max(slowma(iter2-delay:iter2-1))-1>filter % Upward trend
                        close=data(iter2);
                        pos=1;
                        hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                        aaa=aaa+1;
                        ope=close;
                    end
                end
                inputser(iter2,iter1)=pos;
            end
            paramstrct.type{iter1}='MA Cross';
            paramstrct.case(iter1)=3;
            paramstrct.params{iter1}=['x=',num2str(filter),',q=',num2str(lagper2),',p=',num2str(lagper1),',d=',num2str(delay)];
            iter1=iter1+1;
            aaa=1;
            pos=0;
        end
    end
end

% MA crosover case 4 (p-q,d,x,k)

inputser(1,zzy+1:zzy+pqcmbMA*numel(xserMA)*numel(dserMA)*numel(kserMA))=pos;
ope=nan;
aaa=1;
for iter3=1:pqcmbMA % for all p-q
    lagper1=pqMAcmb(iter3,1);
    lagper2=pqMAcmb(iter3,2);
    fastma=tsmovavg(data,'s',lagper1,1);
    slowma=tsmovavg(data,'s',lagper2,1);
    for iter44=1:numel(kserMA) % for all k
        neutper=kserMA(iter44); % Time to neutralize (fixed holding period)
        for iter91=1:numel(dserMA) % for all d
            delay=dserMA(iter91);
            for iter4=1:numel(xserMA) % for all x
                filter=xserMA(iter4);
                for iter2=lagper+delay+1:datsiz
                    if pos==0
                        if min(fastma(iter2-delay:iter2-1))/max(slowma(iter2-delay:iter2-1))-1>filter % Upward trend
                            ope=data(iter2);hld=iter2;
                            pos=1;
                        elseif max(fastma(iter2-delay:iter2-1))/min(slowma(iter2-delay:iter2-1))-1<-filter % Downward trend
                            ope=data(iter2);hld=iter2;
                            pos=-1;
                            
                        end
                    elseif pos==1
                        if iter2-hld==neutper % Neutralizing condition (fixed time)
                            close=data(iter2);
                            pos=0;
                            hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                            aaa=aaa+1;
                            ope=close;
                        end
                        
                    elseif pos==-1
                        if iter2-hld==neutper % Neutralizing condition (fixed time)
                            close=data(iter2);
                            pos=0;
                            hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                            aaa=aaa+1;
                            ope=close;
                        end
                    end
                    inputser(iter2,iter1)=pos;
                end
                paramstrct.type{iter1}='MA Cross';
                paramstrct.case(iter1)=4;
                paramstrct.params{iter1}=['x=',num2str(filter),',q=',num2str(lagper),',p=',num2str(lagper1),',d=',num2str(delay),',k=',num2str(kserfilt(iter44))];
                iter1=iter1+1;
                aaa=1;
                pos=0;
            end
        end
    end
end

zzy=find(inputser(1,:)==0, 1, 'last' );

% MA crosover case 5 (n-p-q,x,d)
npqcmbMA=size(npqMAcmb,1);
inputser(1,zzy+1:zzy+npqcmbMA*numel(xserMA)*numel(dserMA))=pos;
ope=nan;
aaa=1;
for iter3=1:npqcmbMA % for all n-p-q
    lagper1=npqMAcmb(iter3,1); %n
    lagper2=npqMAcmb(iter3,2); %p
    lagper3=npqMAcmb(iter3,3); %q
    fastma=tsmovavg(data,'s',lagper1,1);
    midma=tsmovavg(data,'s',lagper2,1);
    slowma=tsmovavg(data,'s',lagper3,1);
    for iter91=1:numel(dserMA) % for all d
        for iter4=1:numel(xserMA) % for all x
            filter=xserMA(iter4);
            delay=dserMA(iter91);
            for iter2=lagper3+delay+1:datsiz
                longcount=sum([min(data(iter2-delay:iter2))/max(slowma(iter2-delay:iter2-1))-1>filter,...
                    min(data(iter2-delay:iter2))/max(midma(iter2-delay:iter2-1))-1>filter,...
                    min(data(iter2-delay:iter2))/max(fastma(iter2-delay:iter2-1))-1>filter]);
                shortcount=sum([max(data(iter2-delay:iter2))/min(slowma(iter2-delay:iter2-1))-1<-filter,...
                    max(data(iter2-delay:iter2))/min(midma(iter2-delay:iter2-1))-1<-filter,...
                    max(data(iter2-delay:iter2))/min(fastma(iter2-delay:iter2-1))-1<-filter]);
                if pos==0
                    if  longcount>=2 % Upward trend
                        ope=data(iter2);hld=iter2;
                        pos=2/3*longcount-1;
                    elseif shortcount>=2 % Downward trend
                        ope=data(iter2);hld=iter2;
                        pos=-1*(2/3*longcount-1);
                        
                    end
                elseif pos==1/3
                    if longcount==3
                        pos=1;
                    end
                    if shortcount>=2 % Downward trend
                        close=data(iter2);
                        pos=-1*(2/3*longcount-1);
                        hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                        aaa=aaa+1;
                        ope=close;
                    end
                elseif pos==1
                    if shortcount>=2 % Downward trend
                        close=data(iter2);
                        
                        pos=-1*(2/3*longcount-1);
                        hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                        aaa=aaa+1;
                        ope=close;
                    end
                elseif pos==-1/3
                    if shortcount==3
                        pos=-1;
                    end
                    if longcount>=2 % Downward trend
                        close=data(iter2);
                        pos=(2/3*longcount-1);
                        hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                        aaa=aaa+1;
                        ope=close;
                    end
                elseif pos==-1
                    if longcount>=2 % Downward trend
                        close=data(iter2);
                        pos=(2/3*longcount-1);
                        hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                        aaa=aaa+1;
                        ope=close;
                    end
                end
                inputser(iter2,iter1)=pos;
            end
            paramstrct.type{iter1}='MA Cross';
            paramstrct.case(iter1)=5;
            paramstrct.params{iter1}=['x=',num2str(filter),',q=',num2str(lagper3),',p=',num2str(lagper2),',n=',num2str(lagper1),',d=',num2str(delay)];
            iter1=iter1+1;
            aaa=1;
            pos=0;
        end
    end
end
zzy=find(inputser(1,:)==0, 1, 'last' );
disp(' MA Crossover rules generation completed successfully!');
toc

%% Support and resilience trading rules generation (S_R)
disp(' Trading Ranges rules generation started...');
xserS_R=[0.0005,0.001,0.005,0.01,0.025,0.05,0.1]; % alternative in %
dserS_R=0:5;
jserS_R=[2,5,10,15,20,25,50,100,250];
kserS_R=[1,5,10,25];
aasr=numel(xserS_R);

% S&R case 1, (x,d,j)
aabsr=numel(jserS_R);
aazsr=numel(dserS_R);
aaysr=numel(kserS_R);
inputser(1,zzy+1:zzy+aasr*aabsr*aazsr)=pos;
ope=nan;
aaa=1;
for iter91=1:aazsr % for all d
    for iter3=1:aabsr % for all j
        for iter4=1:aasr % for all x
            filter=xserS_R(iter4);
            delay=dserS_R(iter91);
            hlperiod=jserS_R(iter3);
            for iter2=hlperiod+delay+1:datsiz
                if pos==0
                    if min(data(iter2-delay:iter2))/max(data(iter2-delay-hlperiod:iter2-delay-1))-1>filter % Upward trend
                        ope=data(iter2);hld=iter2;
                        pos=1;
                    elseif max(data(iter2-delay:iter2))/min(data(iter2-delay-hlperiod:iter2-delay-1))-1<-filter % Downward trend
                        ope=data(iter2);hld=iter2;
                        pos=-1;
                        
                    end
                elseif pos==1
                    if max(data(iter2-delay:iter2))/min(data(iter2-delay-hlperiod:iter2-delay-1))-1<-filter % Downward trend
                        close=data(iter2);
                        pos=-1;
                        hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                        aaa=aaa+1;
                        ope=close;
                    end
                    
                elseif pos==-1
                    if min(data(iter2-delay:iter2))/max(data(iter2-delay-hlperiod:iter2-delay-1))-1>filter % Upward trend
                        close=data(iter2);
                        pos=1;
                        hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                        aaa=aaa+1;
                        ope=close;
                    end
                end
                inputser(iter2,iter1)=pos;
            end
            paramstrct.type{iter1}='S&R';
            paramstrct.case(iter1)=1;
            paramstrct.params{iter1}=['x=',num2str(filter),',j=',num2str(hlperiod),',d=',num2str(delay)];
            iter1=iter1+1;
            aaa=1;
            pos=0;
        end
    end
end
zzy=find(inputser(1,:)==0, 1, 'last' );

% S&R case 2, (x,d,j,k)
inputser(1,zzy+1:zzy+aasr*aabsr*aazsr*aaysr)=pos;
ope=nan;
aaa=1;
for iter44=1:aaysr % for all k
    for iter91=1:aazsr % for all d
        for iter3=1:aabsr % for all j
            for iter4=1:aasr % for all x
                filter=xserS_R(iter4);
                delay=dserS_R(iter91);
                hlperiod=jserS_R(iter3);
                neutper=kserS_R(iter44); % Time to neutralize (fixed holding period)
                for iter2=hlperiod+delay+1:datsiz
                    if pos==0
                        if min(data(iter2-delay:iter2))/max(data(iter2-delay-hlperiod:iter2-delay-1))-1>filter % Upward trend
                            ope=data(iter2);hld=iter2;
                            pos=1;
                        elseif max(data(iter2-delay:iter2))/min(data(iter2-delay-hlperiod:iter2-delay-1))-1<-filter % Downward trend
                            ope=data(iter2);hld=iter2;
                            pos=-1;
                            
                        end
                    elseif pos==1
                        if iter2-hld==neutper % Neutralizing condition (fixed time)
                            close=data(iter2);
                            pos=0;
                            hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                            aaa=aaa+1;
                            ope=close;
                        end
                        
                    elseif pos==-1
                        if iter2-hld==neutper % Neutralizing condition (fixed time)
                            close=data(iter2);
                            pos=0;
                            hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                            aaa=aaa+1;
                            ope=close;
                        end
                    end
                    inputser(iter2,iter1)=pos;
                end
                paramstrct.type{iter1}='S&R';
                paramstrct.case(iter1)=3;
                paramstrct.params{iter1}=['x=',num2str(filter),',j=',num2str(hlperiod),',d=',num2str(delay),',k=',num2str(neutper)];
                iter1=iter1+1;
                aaa=1;
                pos=0;
            end
        end
    end
end
zzy=find(inputser(1,:)==0, 1, 'last' );

disp(' Trading Ranges rules generation completed successfully!');
toc

%% Channel breakout trading rules generation (CHB)
disp(' Channel breakout rules generation started...');
cserchb=[.1,.5,1,5,10]/100;
xserchb=[.05,.1,.5,1,5]/100;
jserchb=[5,10,15,20,25,50,100,200];
dserchb=0:2;
kserchb=[1,5,10,25];

% CHB case 1, (c,x,d,j)
inputser(1,zzy+1:zzy+numel(cserchb)*numel(xserchb)*numel(dserchb)*numel(jserchb))=pos;
ope=nan;
aaa=1;
for iter96=1:numel(cserchb) % For all c
    channel=cserchb(iter96);
    for iter91=1:numel(dserchb) % for all d
        delay=dserchb(iter91);
        for iter3=1:numel(jserchb) % for all j
            hlperiod=jserchb(iter3);
            for iter4=1:numel(xserchb) % for all x
                filter=xserchb(iter4);
                for iter2=hlperiod+delay+1:datsiz
                    upperbound=min(data(iter2-delay-hlperiod:iter2-delay-1))*(1+channel);
                    lowerbound=max(data(iter2-delay-hlperiod:iter2-delay-1))/(1+channel);
                    if pos==0
                        if min(data(iter2-delay:iter2))/upperbound-1>filter % Upward breakout
                            ope=data(iter2);hld=iter2;
                            pos=1;
                        elseif max(data(iter2-delay:iter2))/lowerbound-1<-filter % Downward breakout
                            ope=data(iter2);hld=iter2;
                            pos=-1;
                            
                        end
                    elseif pos==1
                        if  max(data(iter2-delay:iter2))/lowerbound-1<-filter % Downward breakout
                            close=data(iter2);
                            pos=-1;
                            hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                            aaa=aaa+1;
                            ope=close;
                        end
                        
                    elseif pos==-1
                        if min(data(iter2-delay:iter2))/upperbound-1>filter % Upward breakout
                            close=data(iter2);
                            pos=1;
                            hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                            aaa=aaa+1;
                            ope=close;
                        end
                    end
                    inputser(iter2,iter1)=pos;
                end
                paramstrct.type{iter1}='CHB';
                paramstrct.case(iter1)=1;
                paramstrct.params{iter1}=['x=',num2str(filter),',j=',num2str(hlperiod),',d=',num2str(delay),',c=',num2str(channel)];
                iter1=iter1+1;
                aaa=1;
                pos=0;
            end
        end
    end
end

zzy=find(inputser(1,:)==0, 1, 'last' );
% CHB case 2, (c,x,d,j,k)
inputser(1,zzy+1:zzy+numel(cserchb)*numel(xserchb)*numel(dserchb)*numel(jserchb)*numel(kserchb))=pos;
ope=nan;
aaa=1;
for iter96=1:numel(cserchb)
    channel=cserchb(iter96);
    for iter44=1:numel(kserchb) % for all k
        neutper=kserchb(iter44); % Time to neutralize (fixed holding period)
        for iter91=1:numel(dserchb) % for all d
            delay=dserchb(iter91);
            for iter3=1:numel(jserchb) % for all j
                hlperiod=jserchb(iter3);
                for iter4=1:numel(xserchb) % for all x
                    filter=xserchb(iter4);
                    for iter2=hlperiod+delay+1:datsiz
                        upperbound=min(data(iter2-delay-hlperiod:iter2-delay-1))*(1+channel);
                        lowerbound=max(data(iter2-delay-hlperiod:iter2-delay-1))/(1+channel);
                        if pos==0
                            if min(data(iter2-delay:iter2))/upperbound-1>filter % Upward breakout
                                ope=data(iter2);hld=iter2;
                                pos=1;
                            elseif max(data(iter2-delay:iter2))/lowerbound-1<-filter % Downward breakout
                                ope=data(iter2);hld=iter2;
                                pos=-1;
                                
                            end
                        elseif pos==1
                            if iter2-hld==neutper % Neutralizing condition (fixed time)
                                close=data(iter2);
                                pos=0;
                                hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                                aaa=aaa+1;
                                ope=close;
                            end
                            
                        elseif pos==-1
                            if iter2-hld==neutper % Neutralizing condition (fixed time)
                                close=data(iter2);
                                pos=0;
                                hldp(aaa,iter1)=iter2-hld;hldt(aaa,iter1)=iter2;hld=iter2;
                                aaa=aaa+1;
                                ope=close;
                            end
                        end
                        inputser(iter2,iter1)=pos;
                    end
                    paramstrct.type{iter1}='CHB';
                    paramstrct.case(iter1)=2;
                    paramstrct.params{iter1}=['x=',num2str(filter),',j=',num2str(hlperiod),',d=',num2str(delay),',k=',num2str(neutper),',c=',num2str(channel)];
                    iter1=iter1+1;
                    aaa=1;
                    pos=0;
                end
            end
        end
    end
end
zzy=find(inputser(1,:)==0, 1, 'last' );
disp(' Channel breakout rules generation completed successfully!');
toc

%% Output arguments & finalizing
% Applying transaction cost
ret=TransCost(inputser,dataret,trncost);
tts=zzy*numel(dataret); % Total trading signal
cumretadj=cumsum(ret); % Cumulative return with transaction cost for adjusted returns
finalret=cumretadj(end,:);
disp(' Processing completed!');
disp([' Total number of trading signals generated is ',mat2str(roundn(tts/1e6,-1)),' millions']);
strstart=char(startdt);
strfinish=char(finishdt);
fllbl=[sym1,'_',strstart(end-3:end),'_',strfinish(end-3:end),'-',sym3,'-',num2str(transbasis)];
switch typ
    case 1
        save(fllbl);
    case 2 % For out-of-sample one step before formal start is considered to avoid losing any observation
        inputser(1:gap-1,:)=[];
        data(1:gap-1)=[];
        datadir(1:gap-1)=[];
        ret(1:gap-1)=[];
        save(fllbl,'inputser','data','datadir','start','finish','initstart','gap','ret','transbasis');
        
end
toc