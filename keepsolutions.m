function SolutionsToKeep = keepsolutions(population,Senario)
SolutionsToKeep = [];
SatFit = [Senario(:).SatFit]';
% N = length(population);
%     M = length(population(1).Fit);
%     F = nan(N,M);
%     for i = 1 :N
%         F(i,:) = population(i).Fit;
%     end
     F = [population(:).Fit]';
    [~,ia,~] = unique(F,'rows');
    population = population(ia);
for i = 1 : length(population)
    if all(population(i).Fit<1.1)
        SolutionsToKeep = [SolutionsToKeep,population(i)];
    end
end