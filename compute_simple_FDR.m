function FDR=compute_simple_FDR(gamma,Perfs,pvalues,pi_0)

nbstrats=length(Perfs);

Rejection_Set=sum(pvalues<=gamma);

if(Rejection_Set>0)
    Fplus=min(Rejection_Set,nbstrats*pi_0*gamma);
    FDR=Fplus/Rejection_Set;    
else
    FDR=2;    
end


