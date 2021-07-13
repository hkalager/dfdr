% Based on a code created 10 Oct 2017
% Created 10 May 2019
% Last Modified 10 May 2019 16:50 BST
clear;
clc;
labels={'MXWO' 'NDDUUS'	'NDDUUK'	'NDDUJN'	'MXEF' 'NDEUSRU'	'NDEUCHF'	'NDUEBRAF'...
    'MXFEM' 'MSEIESUN'	'M1MA'	'MIMUJORN'};
regions={'Advanced' 'US' 'UK' 'Japan' 'Emerging' 'Russia' 'China' 'Brazil' ...
    'Frontier' 'Estonia' 'Morocco' 'Jordan'};

transser=[25*ones(1,4),50*ones(1,8)];
%% Year range
yearser=2006:2015;
startyr=2002;
finishyr=2017;
oosperiod=1; 
tic
%% Main loop
for i=1:numel(labels)
    lbl=labels{i};
    transbasis=transser(i);
    fllbl=[lbl,'_',num2str(startyr),'_',num2str(finishyr),'-insample-',num2str(transbasis),'.mat'];
    load(fllbl);
    S1=load([lbl,'_',num2str(yearser(1)),'_',num2str(yearser(end)),'.mat'],'resvar');
    resvarIS=S1.resvar;

    % Calculating the returns without transaction costs
    disp(['Dataset loaded for ',lbl]);
    iter=0;
    tic;
    for yearOOS=yearser(1):yearser(end)%2006:2016
        for month=1:12
            for IS=2
                
                
                %% Loading the FDR portfolios
                iter=iter+1;
                ISportFDR=resvarIS{iter*2,1};
                
                %% Portfolio generation
                D_Index=find(ISportFDR);
                D_ind_count=sum(ISportFDR);
                % Relative Strength Index
                D_ind_rsi=sum(D_Index<=600);
                % Filter rule
                D_ind_fr=sum(and(D_Index>600,D_Index<=3435));
                % Moving average
                D_ind_ma=sum(and(D_Index>3435,D_Index<=16305));
                % Support and resilience
                D_ind_sr=sum(and(D_Index>16305,D_Index<=18195));
                % Channel breakout
                D_ind_cb=sum(and(D_Index>18195,D_Index<=21195));
                
                         
                %% Recording the results
                % Period
                dim=1;
                resvar{iter,dim}=yearOOS;
                dim=dim+1;
                resvar{iter,dim}=month;
                
                % Portfolio size
                dim=dim+1;
                resvar{iter,dim}=sum(ISportFDR);
                dim=dim+1;
                resvar{iter,dim}=D_ind_rsi;
                dim=dim+1;
                resvar{iter,dim}=D_ind_fr;
                dim=dim+1;
                resvar{iter,dim}=D_ind_ma;
                dim=dim+1;
                resvar{iter,dim}=D_ind_sr;
                dim=dim+1;
                resvar{iter,dim}=D_ind_cb;
                dim=dim+1;
                resvar{iter,dim}=regions{i};
                
            end
            
        end
    end
    
    result_table=cell2table(resvar);
    result_table.Properties.VariableNames={'Year' 'Month' ...
        'Discoveries' 'RSI' 'FR' ...
        'MA' 'SR' 'CB' 'Region'};
    drpadd='/Empirics3/';
    %drpadd='C:\Users\Hassannia\Dropbox\Shared_Folder\Chapter3\';
    writetable(result_table,[drpadd,'ID_',lbl,'_',num2str(yearser(1)),'_',num2str(yearser(end)),'.xlsx'])
    save(['ID_',lbl,'_',num2str(yearser(1)),'_',num2str(yearser(end)),'.mat'],'resvar');      
    toc
end
