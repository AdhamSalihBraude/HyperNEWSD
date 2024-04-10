function E = eliteSet(population,NSubProblems)
fitnessvec = nan(size(population,2),NSubProblems);
for i = 1:size(population,2)
    fitnessvec(i,:) = population(i).SubProbelmsScore;
end
[~ , I] = min(fitnessvec,[],1);
I = unique(I);
E = population(I);
end