function [PORT,FDR]=portfolio_FDR(FDRtarget,Perfs,pvalues,pi_0hat)

pvalues(Perfs<=0)=2;

nbstrats=length(Perfs);
idx=(1:nbstrats)';

A=[Perfs pvalues idx];

A=sortrows(A,-1);
A=sortrows(A,2);

PORT=zeros(nbstrats,1);
s=1;
gamma=A(s,2);
oldFDR=-1;
FDR=compute_FDR(gamma,Perfs,pvalues,pi_0hat);
while(FDR<FDRtarget && A(s,2)<2 && s<=nbstrats)
    PORT(A(s,3))=1;
    s=s+1;
    gamma=A(s,2);
    oldFDR=FDR;
    FDR=compute_FDR(gamma,Perfs,pvalues,pi_0hat);
end

FDR=oldFDR;
















