% Test statistic generator based on stationary bootstrap of Potilis and
% Romano (1994)
% 
%  Script developed by Arman Hassanniakalager for the 3rd chapter of PhD,
%  Created on 04 Jul 2017 11:45 BST,
%  Last modified 06 Jul 2017 18:07 BST.
% ****************************************
%% Input arguments:
%   -dataset: the models time series
%   -bench: the Benchmark model
%   -typ: test statistic type 'student' for STUDENTIZED or anything else
%   for GENERAL
%   -Bsize: number of Bootstrap replications
%   -Bwindow: average numbr of block lengeth in the bootstrap replications
%% Output arguments:
%   -tststat: test statistic of the hyposethes
%   -tststatB: Bootstrap test statistic of the hyposethes

function [tststat,tststatB]=tststatboots(dataset,bench,typ,Bsize,Bwindow)
%% Input arguments
if nargin<3
    error('Please call the function with models and benchmark properly')
elseif nargin==3
    Bsize=1000; % Bootstrap replications
    Bwindow=10; % Bootstrap length
elseif nargin==4
    Bwindow=10; % Bootstrap length
end
if strcmpi(typ,'student')
    isStudentized = true;
    disp('Test statistic type is STUDENTIZED');
else
    isStudentized = false;
    disp('Test statistic type is GENERAL');
end

%% switch function between STUDENTIZED and GENERAL test statistic
whichstd=@(x,y) y*std(x)+(1-y)*ones(1,size(x,2));
%% Bootstrap indices
indices=stationary_bootstrap((1:size(dataset,1))',Bsize,Bwindow);
%%
modelscount=size(dataset,2);
tststatB=nan(Bsize,modelscount);
tststat=nan(1,modelscount);
tic;
disp('Bootstrap generation initiated');
for i=1:modelscount
    if mod(i,ceil(modelscount./20)) == 0
        disp([num2str(round(100*(i/modelscount))) '% completed'])
        toc;
    end
    candmodel=dataset(:,i);
    zeroind=find(cumsum(candmodel)==0,1,'last');
    if zeroind<=size(dataset,1)-Bwindow
        tststat(i)=mean(candmodel-bench)/whichstd(candmodel,isStudentized);
        for b=1:Bsize
            bsdata=candmodel(indices(:,b));
            benchB=bench(indices(:,b));
            tststatB(b,i)=mean(bsdata-benchB)/whichstd(bsdata,isStudentized);
        end
    else
        tststat(i)=nan;
        tststatB(:,i)=nan;
    end
end
end