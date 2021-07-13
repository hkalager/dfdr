function pval=pvalsimp(tststat,tststatB)
M=size(tststatB,1);
S=size(tststatB,2);
pval=nan(1,numel(tststat));
for s=1:S
    excesstimes=sum(tststatB(:,s)>=tststat(s));
    pval(s)=excesstimes/M;
end

end