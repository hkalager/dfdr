function BRCpvalue=compute_BRCpvalue(Perfs,Perfs_B)

[nbstrats,B]=size(Perfs_B);

diffs=Perfs_B-kron(ones(1,B),Perfs);

maxperf=max(Perfs);

maxdiffs=max(diffs);

BRCpvalue=sum(maxdiffs>maxperf)/B;
