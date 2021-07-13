function [PORT,FDR]=my_portfolio_FDR(FDRtarget,Perfs,pvalues,pi_0hat)
pvalues(Perfs<=0)=2;
nbins=20;
fdrbins=1/nbins;
nbstrats=length(Perfs);
idx=(1:nbstrats)';

A=[Perfs pvalues idx];

A=sortrows(A,-1);
A=sortrows(A,2);

PORT=zeros(nbstrats,1);
s=1;
%rejectionlevel=ceil(A(s,2)/fdrbins)*fdrbins;
rejectionlevel=max(A(s,2),fdrbins);
%rejectionlevel=A(s,2);

oldFDR=-1;
FDR=compute_FDR(rejectionlevel,Perfs,pvalues,pi_0hat);
if FDR<FDRtarget
    while(FDR<FDRtarget && A(s,2)<2 && s<=nbstrats)
        PORT(A(s,3))=1;
        s=s+1;
        %rejectionlevel=ceil(A(s,2)/fdrbins)*fdrbins;
        %rejectionlevel=A(s,2);
        rejectionlevel=max(A(s,2),fdrbins);
        oldFDR=FDR;
        FDR=compute_FDR(rejectionlevel,Perfs,pvalues,pi_0hat);
    end
    FDR=oldFDR;
else
    Rplus=find(A(:,2)<=rejectionlevel,1,'last');
    PORT(A(1:Rplus,3))=1;
end