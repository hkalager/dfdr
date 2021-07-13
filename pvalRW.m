% Calculating the efficient p-values based on Bootstap based on the paper
% of Romano and Wolf (2016) paper link:
% http://www.sciencedirect.com/science/article/pii/S0167715216000389
%  Script developed by Arman Hassanniakalager for the 3rd chapter of PhD,
%  Created on 04 Jul 2017 17:48 BST,
%  Last modified 07 Jul 2017 14:59 BST.
% ****************************************
%% Input arguments:
%   -tststat: test statistic of the hyposethes
%   -tststatB: Bootstrap test statistic of the hyposethes
%   -mondj: adjustment for monotonicity
function pvaladj=pvalRW(tststat,tststatBUN,monadj)
if nargin<3
    monadj=1;
end
tststatBCNT=tststatBUN-tststat;
[toptststat,rperm]=sort(tststat,'descend');
topboots=tststatBCNT(:,rperm);
M=size(tststatBCNT,1);
S=size(tststatBCNT,2);
for m=1:M
    for j=1:S
        max_t(m,j)=max(topboots(m,j:S));  
    end
end
pvalrperm=nan(1,numel(tststat));
excesstimes1=sum(max_t(:,1)>=toptststat(1));
excesstimes2=sum(max_t(:,1)<=toptststat(1));
excesstimes=min(excesstimes1,excesstimes2);

pvalrperm(1)=(2*excesstimes+1)/(M+1);
for s=2:S
    excesstimes1=sum(max_t(:,s)>=toptststat(s));
    excesstimes2=sum(max_t(:,s)<=toptststat(s));
    excesstimes=min(excesstimes1,excesstimes2);
    pvalrperm(s)=(2*excesstimes+1)/(M+1);
    pvalrperm(s)=monadj*max(pvalrperm(s),pvalrperm(s-1))...
        +(1-monadj)*pvalrperm(s);
end
pvaladj(rperm)=pvalrperm;
pvaladj(isnan(tststat)) = nan;
end


% function pvaladj=pvalRW(tststat,tststatBUN,monadj)
% if nargin<3
%     monadj=1;
% end
% tststatB=tststatBUN-tststat;
% [toptststat,rperm]=sort(tststat,'descend');
% topboots=tststatB(:,rperm);
% M=size(tststatB,1);
% S=size(tststatB,2);
% for m=1:M
%     for j=1:S
%         max_t(m,j)=max(topboots(m,j:S));        
%     end
% end
% pvalrperm=nan(1,numel(tststat));
% excesstimes=sum(max_t(:,1)>=toptststat(1));
% if excesstimes==0
%     excesstimes=nan;
% end
% pvalrperm(1)=(excesstimes+1)/(M+1);
% for s=2:S
%     excesstimes=sum(max_t(:,s)>=toptststat(s));
%     if excesstimes==0
%         excesstimes=nan;
%     end
%     pvalrperm(s)=(excesstimes+1)/(M+1);
%     pvalrperm(s)=monadj*max(pvalrperm(s),pvalrperm(s-1))...
%         +(1-monadj)*pvalrperm(s);
% end
% pvaladj(rperm)=pvalrperm;
% 
% end




