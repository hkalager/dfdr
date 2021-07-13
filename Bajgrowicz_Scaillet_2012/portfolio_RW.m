function PORT = portfolio_RW(alpha, Perfs, Perfs_B)

[nbstrats, ~] = size(Perfs_B);
PORT = zeros(nbstrats, 1);

Diffs = Perfs_B - Perfs;

idx = (1:nbstrats)';

A = [Perfs idx];
A = sortrows(A, -1);

olds = 0;
s = 1;
while(s > olds && s <= nbstrats-1)
    olds = s;
    MaxDiffs = max(Diffs(A(s:end, 2), :));
    c = quantile(MaxDiffs, 1 - alpha);
    while(A(s, 1) - c > 0 && s <= nbstrats-1)
        PORT(A(s, 2)) = 1;
        s = s + 1;
    end
end