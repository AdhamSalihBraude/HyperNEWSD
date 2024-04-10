function matingPool = selection(ind,popSize)
N = size(ind(1).SubProbelmsScore,2);

matingPool = [];
for i = 1 : popSize
    S = randi(popSize,1,2);
    S3 = randi(N,1,1);
    if ind(S(1)).SubProbelmsScore(S3) < ind(S(2)).SubProbelmsScore(S3)
        matingPool = [matingPool,ind(S(1))];
    else
        matingPool = [matingPool,ind(S(2))];
    end
end
end